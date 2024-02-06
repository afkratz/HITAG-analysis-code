#This script runs 
python hitag.py --validate -g grnas/sg_grnas_validated.csv -l linkers/sg_linkers.fasta -i reads/HEK_SG/HEK_concat.fastq -o HEK_SG
python hitag.py --validate -g grnas/sg_grnas_validated.csv -l linkers/sg_linkers.fasta -i reads/HCT_SG/HCT_F2_junction.fastq -o HCT_SG
python hitag.py --validate -g grnas/tf_grnas_f1_validated.csv -l linkers/tf_linkers.fasta -i reads/HEK_TF/HEK293T_F1_set1.fastq -o HEK_TF_F1_S1
python hitag.py --validate -g grnas/tf_grnas_f1_validated.csv -l linkers/tf_linkers.fasta -i reads/HEK_TF/HEK293T_F1_set2.fastq -o HEK_TF_F1_S2
python hitag.py --validate -g grnas/tf_grnas_f2_validated.csv -l linkers/tf_linkers.fasta -i reads/HEK_TF/HEK293T_F2_set1.fastq -o HEK_TF_F2_S1
python hitag.py --validate -g grnas/tf_grnas_f2_validated.csv -l linkers/tf_linkers.fasta -i reads/HEK_TF/HEK293T_F2_set2.fastq -o HEK_TF_F2_S2

#This last call runs a script that executes a series of filtering steps, identifies off-targets, and reports summaries of results for each run above
python process_hitag_results.py