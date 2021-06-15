import os
import pandas as pd
from Bio import SeqIO

def counter():
    source = '/mnt/array2/smallproteins/eCAMI/examples/clustering/input/dbCAN3_new/'
    files = os.listdir(source)
    cazy_family_EZ = []
    for file in files:
        cazy_family = file.split('.faa')[0]
        # if 'GH29' not in cazy_family:
        #     continue
        number_pro = len(list(SeqIO.parse(source+file, "fasta")))
        subfams_fasta_source = '/mnt/array2/smallproteins/eCAMI/examples/clustering/output/dbCAN3/' + cazy_family + '/fasta_for_each_cluster/'
        try:
            subfams_fasta_files = os.listdir(subfams_fasta_source)
            number_subfams = len(subfams_fasta_files)

            number_pro_in_subfams = 0
            for subfam_fasta in subfams_fasta_files:
                subfam_fasta_path = subfams_fasta_source + subfam_fasta
                number_pro_in_subfams += len(list(SeqIO.parse(subfam_fasta_path, "fasta")))

            number_EZ_subfams = 0
            number_EZ_pro = 0
            ez_file = "EZ_info/" + cazy_family + '.txt'
            with open(ez_file, 'r') as f:
                refes = f.readlines()
            f.close()
            for refe in refes:
                refe = refe.split('\n')[0]
                refe_infos = refe.split()
                refe_path = refe_infos[0]
                refe_ez = refe_infos[1]
                if '.' in refe_ez:
                    number_EZ_subfams += 1

                    target_ezs = refe_ez.split('|')
                    for target_ez in target_ezs:
                        if '.' in target_ez:
                            number_EZ_pro += int(target_ez.split(':')[1])
            cazyfamily_info = [cazy_family, number_pro, number_subfams, number_pro_in_subfams, number_EZ_subfams, number_EZ_pro]
            cazy_family_EZ.append(cazyfamily_info)

        except:
            continue
    df = pd.DataFrame(cazy_family_EZ)
    df.to_csv("cazy_family_EC.csv", index=None, header=['cazyfamily', '# of proteins', '# of eCAMI subfams (exclude unclassified)', '# of proteins in subfams', '# of subfams with EC', ' # of proteins in subfams with EC'])

def subfam_counter():
    source = '/mnt/array2/smallproteins/eCAMI/examples/clustering/input/dbCAN3/'
    files = os.listdir(source)
    cazysubfamily_EZ = []
    for file in files:
        cazy_family = file.split('.faa')[0]
        # if 'GH29' not in cazy_family:
        #     continue
        print(cazy_family)
        subfams_fasta_source = '/mnt/array2/smallproteins/eCAMI/examples/clustering/output/dbCAN3/' + cazy_family + '/fasta_for_each_cluster/'
        try:
            subfams_fasta_files = os.listdir(subfams_fasta_source)
            for subfam_fasta in subfams_fasta_files:
                subfam_name = subfam_fasta.split('.txt')[0]
                number_pro_in_subfam = 0
                number_EZ_pro = 0
                number_protein_hmmse = 0
                number_protein_domains_hmmse = 0
                number_protein_cdhit = 0
                print(subfam_name)

                cdhit_source = '/mnt/array2/smallproteins/eCAMI/hmm_analysis_combinded/hmm/' + cazy_family + '/'
                cdhit_files = os.listdir(cdhit_source)
                for cdhit_file in cdhit_files:
                    if 'nr' not in cdhit_file or 'clstr' in cdhit_file:
                        continue
                    print('cdhit_pre: ' + cdhit_file)
                    try:
                        flag = cdhit_file.split(subfam_name)[1]
                    except:
                        flag = 'None'
                    if flag != '.fa.nr':
                        continue
                    print("cdhit: " + cdhit_file)
                    cdhit_file_path = cdhit_source + cdhit_file
                    cdhit_content = list(SeqIO.parse(cdhit_file_path, "fasta"))
                    proteins = []
                    for seq_record in cdhit_content:
                        cazeFamilys = seq_record.id.split('|')[1:]
                        sequence_id = seq_record.id.split('|')[0]
                        if sequence_id not in proteins:
                            proteins.append(sequence_id)
                    number_protein_cdhit = len(proteins)
                print('cdhit done............')

                hmmse_source = '/mnt/array2/smallproteins/eCAMI/hmm_analysis_combinded/cut_domain_seq/' + cazy_family + '/'
                hmmse_files = os.listdir(hmmse_source)
                for hmmse_file in hmmse_files:
                    try:
                        flag = hmmse_file.split(subfam_name)[1]
                    except:
                        flag = 'None'
                    if flag != '.fa':
                        continue
                    print("hmmse: " + hmmse_file)
                    hmmse_file_path = hmmse_source + hmmse_file
                    print(hmmse_file_path)
                    hmmse_content = list(SeqIO.parse(hmmse_file_path, "fasta"))
                    number_protein_domains_hmmse = len(hmmse_content)
                    proteins = []
                    for seq_record in hmmse_content:
                        cazeFamilys = seq_record.id.split('|')[1:]
                        sequence_id = seq_record.id.split('|')[0]
                        if sequence_id not in proteins:
                            proteins.append(sequence_id)
                    number_protein_hmmse = len(proteins)
                print('hmmse done............')

                subfam_fasta_path = subfams_fasta_source + subfam_fasta
                number_pro_in_subfam = len(list(SeqIO.parse(subfam_fasta_path, "fasta")))
                ez_file = "EZ_info/" + cazy_family + '.txt'
                with open(ez_file, 'r') as f:
                    refes = f.readlines()
                f.close()

                for refe in refes:
                    if subfam_fasta in refe:
                        refe = refe.split('\n')[0]
                        refe_infos = refe.split()
                        refe_path = refe_infos[0]
                        refe_ez = refe_infos[1]
                        if '.' in refe_ez:
                            target_ezs = refe_ez.split('|')
                            for target_ez in target_ezs:
                                if '.' in target_ez:
                                    number_EZ_pro += int(target_ez.split(':')[1])
                print('subfam search done............')

                cazysubfamily_info = [cazy_family + '_' + subfam_name, number_pro_in_subfam, number_EZ_pro, number_protein_hmmse, number_protein_domains_hmmse, number_protein_cdhit]
                print(cazysubfamily_info)
                cazysubfamily_EZ.append(cazysubfamily_info)
        except:
            continue
    df = pd.DataFrame(cazysubfamily_EZ)
    df.to_csv("cazysubfamily_EC.csv", index=None, header=['CAZy subfam', '# of proteins', '# of proteins with EC', 'after hmmsearch # of remaining proteins', 'after hmmsearch # of remaining protein domains', 'after the usearch # of remaining proteins for mafft'])

if __name__ == '__main__':
    # create a excel with [cazyfamily,# of proteins,# of eCAMI subfams (exclude unclassified),# of proteins in subfams,# of subfams with EC, # of proteins in subfams with EC]
    counter()
    #create a excel with [ CAZy subfam,# of proteins,# of proteins with EC,after hmmsearch # of remaining proteins,after hmmsearch # of remaining protein domains,after the usearch # of remaining proteins for mafft]
    subfam_counter()
