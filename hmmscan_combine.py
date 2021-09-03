import argparse
import json
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folders', nargs='+', required=True)
    args = parser.parse_args()

    category = 'hmm_analysis_combinded/'
    try:
        os.mkdir(category)
    except:
        pass

    output_path = category + 'domain_result/'
    try:
        os.mkdir(output_path)
    except:
        pass

    cut_output_path = category + 'cut_domain_seq/'
    try:
        os.mkdir(cut_output_path)
    except:
        pass

    seq_summary = category + 'seq_summary'
    try:
        os.mkdir(seq_summary)
    except:
        pass

    hmm_path = 'ori_hmm_refe_combine/'

    path = 'examples/clustering/output/dbCAN3/'
    folders = args.folders

    try:
        with open(category + 'number_protein.json', 'r') as f:
            number_portein_all = json.load(f)
        f.close()
    except:
        number_portein_all = {}

    for folder in folders:
        number_portein = []
        output_path = category + 'domain_result/' + folder + '/'
        try:
            os.mkdir(output_path)
        except:
            pass


        cut_output_path = category + 'cut_domain_seq/' + folder + '/'
        try:
            os.mkdir(cut_output_path)
        except:
            pass

        folder_path = path + folder + '/fasta_for_each_cluster/'
        files = os.listdir(folder_path)
        for file in files:
            print(file)
            file_path = folder_path + file
            file_name = file.split('.txt')[0]
            output_file = output_path + folder + '_' + file_name + '.out.dm'
            hmm_file = hmm_path + folder + '/' + folder + '.hmm'
            command = 'hmmscan --domtblout ' + output_file + ' ' + hmm_file + ' ' + file_path
            os.system(command)

            parse_output_file = output_file + '.ps'
            command = 'sh hmmscan-parser.sh ' + output_file + ' > ' + parse_output_file
            os.system(command)

            # stringent_output_file = parse_output_file + '.stringent'
            # command = 'cat ' + parse_output_file + '| awk \'$5<1e-15&&$10>0.35\' > ' + stringent_output_file
            # os.system(command)

            print(file_path)
            record_dic = SeqIO.to_dict(SeqIO.parse(file_path, "fasta"))

            with open(parse_output_file, 'r') as file:
                results = file.readlines()
            file.close()

            cutted_fasta_seqs = []
            protein4count = []
            for result in results:
                result = result.split('\n')[0]
                attributes = result.split()
                query_id = attributes[2]
                start = int(attributes[5])
                end = int(attributes[6])
                original_seq = str(record_dic[query_id].seq)
                cutted_seq = original_seq[start:end+1]
                cutted_rec = SeqRecord(
                    Seq(cutted_seq),
                    id=record_dic[query_id].id,
                    name=record_dic[query_id].name,
                    description='start:' + str(start) + '|end:' + str(end)
                )
                cutted_fasta_seqs.append(cutted_rec)

                ###query_id
                protein4count.append(query_id.split('|')[0])
            number_portein += list(set(protein4count))
            SeqIO.write(cutted_fasta_seqs, cut_output_path + 'cut_' + file_name + '.fa', "fasta")
        number_portein_all[folder] = len(number_portein)

        with open(seq_summary + '/' + folder + '.json', 'w') as f:
            data = {}
            data['data'] = number_portein
            json.dump(data, f)
        f.close()
    with open(category + 'number_protein.json', 'w') as f:
        json.dump(number_portein_all, f)
    f.close()


if __name__ == '__main__':
    main()
