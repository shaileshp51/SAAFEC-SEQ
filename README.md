# SAAFEC-SEQ
 Predicting change of folding free energy of a protein due to mutation using sequence only
1. blast software only need to copy to any directory
2. Please change the installation directory of blast
3. run the python script:
 1) single mutation
 python Mutation_pred.py -A protein_sequence -p position -w wild -m mutant -o outputfile(Optional)
 2) mutations list
 python Mutation_pred.py -A protein_sequence -f mutation_list -o outputfile

Example:
python Mutation_pred.py -A fastaA -p 182 -w C -m A