#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import numpy as np
import xgboost as xgb

from utils.protseqfeature import *

__version__ = 1.0
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def configure():
    config = {}
    config["PSIBLASTDIR"] = "/zfs/compbio/DELPHI/Script_2020-08-05/SAAFEC-SEQ/blast/"
    config[
        "PSIBLASTBASE"
    ] = "/zfs/compbio/DELPHI/Script_2020-08-05/SAAFEC-SEQ/Uniref100/uniref100"
    config["BLAST_NUM_THREADS"] = 6
    return config


def read_fasta(fastafile):
    sequence = "".join(
        [a.strip() for a in open(fastafile) if not a.strip().startswith(">")]
    ).upper()
    for i in sequence:
        if i not in "ARNDCQEGHILKMFPSTWYV":
            print(
                (
                    f"ERROR>> Invalid protein sequence, unknown amino acid"
                    f" '{i}' in file {fastafile}"
                )
            )
            sys.exit()
    return sequence


def sequences_features(target_seqn, mutation_resid, wild_aa):
    mutation_resid = int(mutation_resid) - 1
    features = []
    if wild_aa.upper() != target_seqn[mutation_resid]:
        print("Wild type is not same as input sequence!")
        sys.exit()
    for i in range(int(mutation_resid) - 5, int(mutation_resid) + 6):
        if i < 0 or i >= len(target_seqn):
            features.append("0")
            continue
        else:
            features.append(aa1_label.get(target_seqn[i], 0))
    return features


def run_psiblast(config, file_name, load_existing=False):
    if not (load_existing and os.path.isfile(file_name + ".pssm")):
        psiblast_found = os.path.isfile(os.path.join(config["PSIBLASTDIR"], "psiblast")) 
        psibloatdb_found = os.path.isdir(config["PSIBLASTBASE"])
        errors = []
        if not psiblast_found:
            errors.append(f"PSIBLAST exe '{os.path.join(config['PSIBLASTDIR'], 'psiblast')}' is missing.")
        if not psibloatdb_found:
            errors.append(f"PSIBLASTBASE dir '{config['PSIBLASTBASE']}' is missing.")
        if psiblast_found and psibloatdb_found:
            os.system(
                config["PSIBLASTDIR"]
                + "/psiblast -query "
                + file_name
                + " -num_threads "
                + str(config["BLAST_NUM_THREADS"])
                + " -db "
                + config["PSIBLASTBASE"]
                + " -num_iterations 3 -out "
                + file_name
                + ".out -out_ascii_pssm "
                + file_name
                + ".pssm 2>/dev/null"
            )
        else:
            print("ERROR>> " + "\nERROR>> ".join(errors))
            sys.exit()


def pssm_check(file_name):
    if os.path.exists(file_name + ".pssm"):
        pass
    else:
        print("Can't get the PSSM, please check the sequence")
        sys.exit()


def psepssm(file, lamda=7):
    pssm = []
    features = []
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
    delete_index = 0
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
        label.append(net_volume2(aa1_map, wild_aa, mutation_aa))
        label.append(net_hydrophobicity2(aa1_map, wild_aa, mutation_aa))
        label.append(net_flexibility2(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_hydrophobicity2(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_polarity2(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_type2(aa1_mttype_label, wild_aa, mutation_aa))
        label.append(mutation_size2(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_hbonds2(aa1_map, wild_aa, mutation_aa))
        label.append(mutation_chemical2(aa1_map, wild_aa, mutation_aa))
        # label.append(mutation_ala(wild_aa,mutation_aa))
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
    model = xgb.Booster(model_file=f"{__location__}/regression.model")
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


def argument_parser():
    parser = argparse.ArgumentParser(
        prog=f"saafec-seq.py",
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
        help="".join(
            [
                "Mutation list file. The file should have one mutation per line ",
                "where every line has following three information seperated ",
                "by space. wild-type(1-letter-code) resid mutant(1-letter-code)",
            ]
        ),
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
    # If there is not a single valid mutation turn all warnings to errors
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
    if (not args.mutation_list_file is None) and is_file(args.mutation_list_file):
        for line in open(args.mutation_list_file):
            if not line.strip() or line.strip().startswith("#"):
                continue
            info = line.strip().split(" ")
            mutations_list.append((int(info[1]), info[0], info[2]))
            if len(info) != 3:
                errors.append(
                    (
                        "Mutation list file format error. "
                        "Every line should have exactly three "
                        "information: 'position wild-type mutant'"
                    )
                )
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
        errors.append("Mutation information is require but missing.")
    if len(results["errors"]) > 0:
        print("ERROR>> " + "\nERROR>> ".join(results["errors"]))
        sys.exit()
    if len(results["warns"]) > 0:
        print("WARNING>> " + "\nWARNING>> ".join(results["warns"]), file=sys.stderr)

    return results


def main():
    args = argument_parser()
    struct = None
    config = configure()
    results = validate_input(args)
    f = open(args.output, "w")
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
    f.close()


if __name__ == "__main__":
    main()
