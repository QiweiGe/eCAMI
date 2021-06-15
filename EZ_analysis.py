import argparse
import os

from Bio import SeqIO


def main():
    with open('path_summary.txt', 'r') as f:
        paths = f.readlines()
    f.close()
    for path in paths:
        path = path.split('\n')[0]
        target = path.split('/')[-2]
        file = path + target + '.faa'
        with open(file, 'r') as fasta_file:
            records = list(SeqIO.parse(file, "fasta"))
        fasta_file.close()
        with open("combined_tmhmm_result.faa", "a+") as combined_file:
            SeqIO.write(records, combined_file, "fasta")
        combined_file.close()


def EZ_number(folders):
    EZ_number_folder = "EZ_info/"
    try:
        os.mkdir(EZ_number_folder)
    except:
        pass
    path = 'examples/clustering/output/dbCAN3/'
    for folder in folders:
        folder_path = path + folder + '/kmer_for_each_cluster/*'
        command = "grep ^" + folder + " " + folder_path + " > " + EZ_number_folder + folder + '.txt'
        os.system(command)

def EZ_analysis(folders):
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-f', '--folders', nargs='+', required=True)
    # args = parser.parse_args()
    # folders = args.folders

    EZ_number_folder = "EZ_info/"
    for folder in folders:
        ez_file = EZ_number_folder + folder + '.txt'
        with open(ez_file, 'r') as f:
            refes = f.readlines()
        f.close()
        for refe in refes:
            refe = refe.split('\n')[0]
            refe_infos = refe.split()
            refe_path = refe_infos[0]
            refe_ez = refe_infos[1]
            #if '.' in refe_ez:
            processed_path = '/mnt/array2/smallproteins/eCAMI/hmm_analysis_combinded/hmm/' + folder + '/'
            file_info = refe_path.split('/')[-1].split(':')[0].split('.txt')[0]
            file_name = folder + '_cut_' + file_info + '.hmm'
            processed_file = processed_path + file_name

            try:
                with open(processed_file, 'r') as f:
                    hmm_info = f.read()
                f.close()
            except:
                file_info = 'cluster_0'
                file_name = folder + '_cut_' + file_info + '.hmm'
                processed_file = processed_path + file_name
                with open(processed_file, 'r') as f:
                    hmm_info = f.read()
                f.close()

            target_ezs = refe_ez.split('|')
            target = ''
            for target_ez in target_ezs:
                if '.' in target_ez:
                    target = target + '|' + target_ez

            name_info = hmm_info.split('NAME  ')[1].split('LENG')[0]
            print(name_info)
            if refe_ez in name_info:
                continue
            replace_name = folder + '_eCAMI_' + file_info.split('cluster_')[1] + '.hmm|' + refe_ez
            new_hmm_info = hmm_info.replace(file_name, replace_name)
            #new_hmm_info = hmm_info.replace(file_name + target, replace_name)
            with open(processed_file, 'w') as f:
                f.write(new_hmm_info)
            f.close()

def EZ_combine(folders):
    output_path = '/mnt/array2/smallproteins/eCAMI/hmm_refe_combine/'
    target_path = '/mnt/array2/smallproteins/eCAMI/hmm_analysis_combinded/hmm/'
    for folder in folders:
        combined_hmm = ''
        folder_path = target_path + folder + '/'
        files = os.listdir(folder_path)
        for file in files:
            if '.hmm' in file and '.aln' not in file and 'unclustered' not in file:
                file_path = folder_path + file
                with open(file_path, 'r') as f:
                    hmm = f.read()
                f.close()
                combined_hmm = combined_hmm + hmm + '\n'
        output_folder_path = output_path + folder + '/'
        rm_command = 'rm -rf ' + output_folder_path + '*'
        os.system(rm_command)

        output_file = folder + '.hmm'
        output_file_path = output_folder_path + output_file
        with open(output_file_path, 'w') as f:
            f.write(combined_hmm)
        f.close()

        hmmpress_command = 'hmmpress ' + output_file_path
        os.system(hmmpress_command)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folders', nargs='+', required=True)
    args = parser.parse_args()
    folders = args.folders
    EZ_number(folders)
    EZ_analysis(folders)
    EZ_combine(folders)
