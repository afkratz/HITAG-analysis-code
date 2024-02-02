

import analyze_hitag
import os
def run_analysis_one():
    """
    HEK-SG dataset
    """
    input_files =['reads/HEK_SG/A252_R1.fastq','reads/HEK_SG/A258_R1.fastq','reads/HEK_SG/A272_R1.fastq','reads/HEK_SG/A279_R1.fastq']
    output_file = 'reads/HEK_SG/HEK_concat.fastq'
    if not os.path.exists('reads/HEK_SG/HEK_concat.fastq'):
        with open(output_file, 'w') as out_f:
            for file_name in input_files:
                with open(file_name, 'r') as in_f:
                    for line in in_f:
                        out_f.write(line)

    analyze_hitag.run_analysis(
        input_file_name='reads/HEK_SG/HEK_concat.fastq',
        linker_file_name='sg_linker_index',
        output='HEK_hitag',
        target_file_name='sg_grnas_validated.csv'
    )



def run_analysis_two():
    """
    HCT-SG dataset
    """
    analyze_hitag.run_analysis(
        input_file_name='reads/HCT_SG/HCT_F2_junction.fastq',
        linker_file_name='sg_linker_index',
        output='HCT',
        target_file_name='sg_grnas_validated.csv'  
        )

def run_analysis_three():
    """
    HEK-TF dataset
    """
    analyze_hitag.run_analysis(
        input_file_name='reads/HEK_TF/HEK293T_F1_set1.fastq',
        linker_file_name='tf_linker_index',
        output = 'HEK_TF_F1_S1',
        target_file_name='tf_grnas_f1_validated.csv',
    )
    
    analyze_hitag.run_analysis(
        input_file_name='reads/HEK_TF/HEK293T_F1_set2.fastq',
        linker_file_name='tf_linker_index',
        output = 'HEK_TF_F1_S2',
        target_file_name='tf_grnas_f1_validated.csv',
    )
    

    analyze_hitag.run_analysis(
        input_file_name='reads/HEK_TF/HEK293T_F2_set1.fastq',
        linker_file_name='tf_linker_index',
        output = 'HEK_TF_F2_S1',
        target_file_name='tf_grnas_f2_validated.csv',
    )
    
    analyze_hitag.run_analysis(
        input_file_name='reads/HEK_TF/HEK293T_F2_set2.fastq',
        linker_file_name='tf_linker_index',
        output = 'HEK_TF_F2_S2',
        target_file_name='tf_grnas_f2_validated.csv',
    )

if __name__=='__main__':
    run_analysis_one()
    #run_analysis_two()
    run_analysis_three()
    run_analysis_four()