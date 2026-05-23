#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import xgboost as xgb

from saafec_seq.utils.protseqfeature import (
    aa1_label,
    aa1_map,
    aa1_mttype_label,
    aa3to1,
    mutation_chemical,
    mutation_hbonds,
    mutation_hydrophobicity,
    mutation_polarity,
    mutation_size,
    mutation_type,
    net_flexibility,
    net_hydrophobicity,
    net_volume,
)

__version__ = "1.0"
PROGRAM = "SAAFEC-SEQ"

CONFIG_PATH = Path.home() / ".saafec-seqrc"
DEFAULT_BLAST_THREADS = 6
DEFAULT_DB_NAME = "uniref100"
SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_PATH = SCRIPT_DIR / "data" / "regression.model"

ALLOWED_AA = set(aa1_label.keys())


@dataclass(frozen=True)
class Config:
    psiblast_bin_dir: str
    uniref_db_dir: str
    uniref_db_name: str = DEFAULT_DB_NAME
    blast_num_threads: int = DEFAULT_BLAST_THREADS

    @property
    def psiblast_path(self) -> str:
        return str(Path(self.psiblast_bin_dir) / "psiblast")

    @property
    def psiblast_base(self) -> str:
        return str(Path(self.uniref_db_dir) / self.uniref_db_name)


def _parse_key_value_text(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip().upper()] = value.strip()
    return values


def read_config_file(path: Path) -> Config:
    values = _parse_key_value_text(path.read_text(encoding="utf-8"))
    return Config(
        psiblast_bin_dir=values.get("PSIBLAST_BIN_DIR", ""),
        uniref_db_dir=values.get("UNIREF_DB_DIR", ""),
        uniref_db_name=values.get("UNIREF_DB_NAME", DEFAULT_DB_NAME),
        blast_num_threads=int(values.get("BLAST_NUM_THREADS", str(DEFAULT_BLAST_THREADS))),
    )


def write_config_file(path: Path, cfg: Config) -> None:
    path.write_text(
        "\n".join(
            [
                f"PSIBLAST_BIN_DIR={cfg.psiblast_bin_dir}",
                f"UNIREF_DB_DIR={cfg.uniref_db_dir}",
                f"UNIREF_DB_NAME={cfg.uniref_db_name}",
                f"BLAST_NUM_THREADS={cfg.blast_num_threads}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def _blast_db_exists(base_path: str) -> bool:
    base = Path(base_path)
    candidates: Iterable[Path] = (
        base.with_suffix(ext)
        for ext in (".phr", ".pin", ".psq", ".pog", ".pto", ".pot", ".ptf")
    )
    return any(candidate.exists() for candidate in candidates)


def is_valid_runtime_config(cfg: Config) -> bool:
    return Path(cfg.psiblast_path).is_file() and _blast_db_exists(cfg.psiblast_base)


def prompt_for_config() -> Config:
    print("SAAFEC-SEQ configuration is required.")
    psiblast_bin_dir = input("Path to PSI-BLAST binary directory: ").strip()
    uniref_db_dir = input("Path to UniRef100 database directory: ").strip()
    uniref_db_name = input(f"UniRef100 database name [{DEFAULT_DB_NAME}]: ").strip() or DEFAULT_DB_NAME
    blast_num_threads_raw = input(f"Number of BLAST threads [{DEFAULT_BLAST_THREADS}]: ").strip() or str(DEFAULT_BLAST_THREADS)

    try:
        blast_num_threads = int(blast_num_threads_raw)
    except ValueError:
        blast_num_threads = DEFAULT_BLAST_THREADS

    return Config(
        psiblast_bin_dir=psiblast_bin_dir,
        uniref_db_dir=uniref_db_dir,
        uniref_db_name=uniref_db_name,
        blast_num_threads=blast_num_threads,
    )


def load_config(reconfigure: bool = False, config_path: Path | None = None) -> Config:
    """
    Load ~/.saafec-seqrc first. If it is missing or invalid, prompt the user
    and rewrite the config file for future runs.
    """
    path = config_path or CONFIG_PATH

    if not reconfigure and path.exists():
        cfg = read_config_file(path)
        if is_valid_runtime_config(cfg):
            return cfg

    if not sys.stdin.isatty():
        if path.exists():
            raise RuntimeError(
                f"Invalid configuration in {path}. Run interactively to reconfigure it."
            )
        raise RuntimeError(
            f"Missing configuration file {path}. Run interactively to create it."
        )

    cfg = prompt_for_config()
    write_config_file(path, cfg)

    if not is_valid_runtime_config(cfg):
        raise RuntimeError(
            f"Configuration written to {path}, but the runtime paths still do not look valid."
        )

    return cfg


def read_fasta(fastafile):
    with open(fastafile, "r", encoding="utf-8") as fh:
        sequence = "".join(
            [a.strip() for a in fh if not a.strip().startswith(">")]
        ).upper()
   
    for aa in sequence:
        if aa not in ALLOWED_AA:
            print(
                f"ERROR>> Invalid protein sequence, unknown amino acid '{aa}' in file {fastafile}"
            )
            sys.exit(1)
    return sequence


def sequences_features(target_seqn, mutation_resid, wild_aa):
    mutation_resid = int(mutation_resid) - 1
    features = []
    if wild_aa.upper() != target_seqn[mutation_resid]:
        print("Wild type is not same as input sequence!")
        sys.exit(1)
    for i in range(int(mutation_resid) - 5, int(mutation_resid) + 6):
        if i < 0 or i >= len(target_seqn):
            features.append("0")
            continue
        else:
            features.append(aa1_label.get(target_seqn[i], 0))
    return features


def run_psiblast(config, file_name, load_existing=False):
    if not (load_existing and os.path.isfile(file_name + ".pssm")):
        psiblast_found = os.path.isfile(config.psiblast_path)
        psiblastdb_found = os.path.isdir(config.uniref_db_dir)
        errors = []
        if not psiblast_found:
            errors.append(f"PSIBLAST exe '{config.psiblast_path}' is missing.")
        if not psiblastdb_found:
            try:
                subprocess.run(["ls", "-ltr", config.uniref_db_dir], check=False)
            except Exception:
                pass
            errors.append(f"PSIBLASTBASE dir '{config.uniref_db_dir}' is missing.")
        if psiblast_found and psiblastdb_found:
            cmd = [
                config.psiblast_path,
                "-query",
                file_name,
                "-num_threads",
                str(config.blast_num_threads),
                "-db",
                config.psiblast_base,
                "-num_iterations",
                "3",
                "-out",
                f"{file_name}.out",
                "-out_ascii_pssm",
                f"{file_name}.pssm",
            ]
            log_path = Path.cwd() / "blast.log"
            with open(log_path, "a", encoding="utf-8") as logfh:
                subprocess.run(cmd, stdout=logfh, stderr=logfh, check=False)
        else:
            print("ERROR>> " + "\nERROR>> ".join(errors))
            sys.exit(1)


def pssm_check(file_name):
    if os.path.exists(file_name + ".pssm"):
        pass
    else:
        print("Can't get the PSSM, please check the sequence")
        sys.exit(1)


def psepssm(file, lamda=7):
    pssm = []
    features = []
    length = 0
    for line1 in open(file + ".pssm"):
        info1 = line1.strip().split()
        if len(info1) > 43:
            if info1[1].isupper():
                length = int(info1[0])
                for i in range(2, 22):
                    pssm.append(float(1) / (1 + np.e ** (-int(info1[i]))))
    for i in range(20):
        sum_pssm = 0
        for j in range(int(length)):
            sum_pssm = sum_pssm + float(pssm[i + 20 * j])
        features.append("%.2f" % (float(sum_pssm) / length))
    for i in range(1, lamda + 1):
        for j in range(20):
            s_pssm = 0
            for k in range(int(length) - i):
                s_pssm = (
                    float(pssm[20 * k + j]) - float(pssm[20 * (k + i) + j])
                ) ** 2 + s_pssm
            features.append("%.2f" % (float(s_pssm) / (int(length) - i)))
    return features


def mpssmscores(file_name, position, windows=7):
    pssm = []
    line1 = []
    for line2 in open(file_name + ".pssm"):
        line1.append(line2.strip())
    for j in range(int(position) - (windows // 2), int(position) + windows // 2 + 1):
        if j < 0:
            for i in range(2, 22):
                pssm.append(0)
            continue
        index = 0
        for info1 in line1:
            info2 = info1.split()
            if len(info2) > 43:
                if int(info2[0]) == j:
                    for i in range(2, 22):
                        pssm.append(info2[i])
                    index = 1
                if info1 == line1[-7] and index == 0:
                    for i in range(2, 22):
                        pssm.append(0)
    return pssm


def file_loop(config, chainA, seqnA, mutation_list, verbose):
    results = []
    for mt_index, mt_info in enumerate(mutation_list):
        mutation_resid = mt_info[0]
        wild_aa, mutation_aa = mt_info[1], mt_info[2]
        if len(wild_aa) == 3:
            wild_aa = aa3to1[wild_aa]
        if len(mutation_aa) == 3:
            mutation_aa = aa3to1[mutation_aa]
        if len(mutation_list) == 1 and wild_aa == mutation_aa:
            return ["%.2f" % 0.0 + " Neutral"]
        label = []
        label.append(net_volume(aa1_map, wild_aa, mutation_aa))
        label.append(net_hydrophobicity(aa1_map, wild_aa, mutation_aa))
        label.append(net_flexibility(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_hydrophobicity(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_polarity(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_type(aa1_mttype_label, wild_aa, mutation_aa))
        label.append(mutation_size(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_hbonds(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_chemical(aa1_map, wild_aa, mutation_aa))
        label = label + sequences_features(seqnA, mutation_resid, wild_aa)
        if mt_index == 0:
            run_psiblast(config, chainA)
            pssm_check(chainA)
        label = label + mpssmscores(chainA, mutation_resid)
        label = label + psepssm(chainA)
        if verbose:
            print(
                "Mutation features: ",
                chainA,
                mutation_resid,
                wild_aa,
                mutation_aa,
                label,
            )
        results.append(pred_feature(label, wild_aa, mutation_aa))
    return results


def pred_feature(label, wild_aa, mutation_aa):
    if not MODEL_PATH.is_file():
        raise FileNotFoundError(f"Missing model file: {MODEL_PATH}")
    model = xgb.Booster(model_file=str(MODEL_PATH))
    if wild_aa == mutation_aa:
        return "%.2f" % 0.0 + " Neutral"
    x = np.array(label)
    x = x.reshape((1, len(label)))
    x = xgb.DMatrix(x)
    y_pred = model.predict(x)
    if y_pred[0] > 0:
        return "%.2f" % y_pred[0] + " Destabilizing"
    else:
        return "%.2f" % y_pred[0] + " Stabilizing"


def is_file(filepath):
    val = os.path.isfile(filepath)
    if not val:
        raise argparse.ArgumentTypeError("`%s` is not a valid filepath" % filepath)
    return filepath


def validate_aa_code(code):
    if code is None:
        return 0
    code = str(code).strip().upper()
    return 1 if code in ALLOWED_AA else 0


def argument_parser():
    parser = argparse.ArgumentParser(
        prog="saafec-seq",
        description=(
            f"SAAFEC-SEQ_v{__version__} : "
            "Predict the free energy change of folding due to "
            "point mutation for protein from sequence."
        ),
    )
    parser.add_argument(
        "-i",
        "--input-sequence",
        help="fasta filename of the protein sequence",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="output file name for the predictions (default: output.out)",
        default="output.out",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="print verbose messages",
        action="store_true",
    )
    parser.add_argument(
        "-c",
        "--config",
        help="path to configuration file (default: ~/.saafec-seqrc)",
        default=None,
    )
    parser.add_argument(
        "--reconfigure",
        action="store_true",
        help="prompt for configuration even if a config file already exists",
    )
    group1 = parser.add_argument_group(
        "single point-mutation", "single point-mutation information"
    )
    group1.add_argument(
        "-p",
        "--position",
        type=int,
        help="position of mutation (starting index is 1)",
    )
    group1.add_argument(
        "-w",
        "--wild-type",
        type=validate_aa_code,
        help="wild-type amino acid's one-letter code",
    )
    group1.add_argument(
        "-m",
        "--mutant",
        type=validate_aa_code,
        help="mutant amino acid's one-letter code",
    )
    group2 = parser.add_argument_group(
        "multiple point-mutations", "multiple point-mutations listing"
    )
    group2.add_argument(
        "-f",
        "--mutation-list-file",
        help="Mutation list file. Each line has: wild-type resid mutant",
    )
    args = parser.parse_args()
    return args


def validate_mutation_info(mt_seqn, mutation_list):
    errors, warns = [], []
    valided_mt_list = []
    for mt in mutation_list:
        position, wild_type, mutant = mt[0], mt[1], mt[2]
        is_valid_mt = True
        if position < 1 or position > len(mt_seqn):
            warns.append(
                (
                    f"position {position} is beyond "
                    f"mutation_sequence length: {len(mt_seqn)}"
                )
            )
            is_valid_mt = False
        elif validate_aa_code(wild_type) == 0:
            warns.append(f"wild-type amino acid {wild_type} is invalid")
            is_valid_mt = False
        elif mt_seqn[position - 1] != wild_type:
            warns.append(
                (
                    f"wild-type amino acid '{wild_type}' mis matches"
                    f"with mutation sequence '{mt_seqn[position - 1]}'"
                )
            )
            is_valid_mt = False
        if validate_aa_code(mutant) == 0:
            warns.append(f"mutant amino acid {mutant} is invalid")
            is_valid_mt = False
        if is_valid_mt:
            valided_mt_list.append((position, wild_type, mutant))
    if len(valided_mt_list) == 0:
        errors.extend(warns)
        warns = []
    return errors, warns, valided_mt_list


def validate_input(args):
    results = {}
    results["mt_seqn"] = ""
    mutations_list = []
    results["errors"], results["warns"] = [], []
    if is_file(args.input_sequence):
        results["mt_seqn"] = read_fasta(args.input_sequence)
    if (args.mutation_list_file is not None) and is_file(args.mutation_list_file):
        with open(args.mutation_list_file, "r", encoding="utf-8") as fh:
            for line in fh:
                if not line.strip() or line.strip().startswith("#"):
                    continue
                info = line.strip().split()
                if len(info) != 3:
                    results["errors"].append(
                        (
                            "Mutation list file format error. "
                            "Every line should have exactly three "
                            "information: 'position wild-type mutant'"
                        )
                    )
                    continue
                mutations_list.append((int(info[1]), info[0], info[2]))
    elif not (
        (args.position is None) or (args.wild_type is None) or (args.mutant is None)
    ):
        mutations_list.append((args.position, args.wild_type, args.mutant))
    if len(mutations_list) > 0:
        errors, warns, valided_mt_list = validate_mutation_info(
            results["mt_seqn"], mutations_list
        )
        results["errors"].extend(errors)
        results["warns"].extend(warns)
        results["mt_list"] = valided_mt_list
    else:
        results["errors"].append("Mutation information is require but missing.")
    if len(results["errors"]) > 0:
        print("ERROR>> " + "\nERROR>> ".join(results["errors"]))
        sys.exit(1)
    if len(results["warns"]) > 0:
        print("WARNING>> " + "\nWARNING>> ".join(results["warns"]), file=sys.stderr)

    return results


def main():
    args = argument_parser()
    config = load_config(reconfigure=args.reconfigure, config_path=Path(args.config).expanduser() if args.config else None)
    results = validate_input(args)
    with open(args.output, "w", encoding="utf-8") as f:
        print("Position Wild Mutant ddG Type", file=f)
        preds = file_loop(
            config,
            args.input_sequence,
            results["mt_seqn"],
            results["mt_list"],
            args.verbose,
        )
        for mt, pred in zip(results["mt_list"], preds):
            print(mt[0], mt[1], mt[2], pred, file=f)

    print()
    print("Please cite SAAFEC-SEQ when using its results in a publication:")
    print(
        "Li, G.; Panday, S.K.; Alexov, E. "
        "SAAFEC-SEQ: A Sequence-Based Method for Predicting the Effect of "
        "Single Point Mutations on Protein Thermodynamic Stability. "
        "Int. J. Mol. Sci. 2021, 22, 606."
    )


if __name__ == "__main__":
    main()
