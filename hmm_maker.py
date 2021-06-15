import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folders', nargs='+', required=True)
    args = parser.parse_args()
    folders = args.folders

    output_path = "/mnt/array2/smallproteins/eCAMI/hmm_analysis_combinded/hmm/"
    try:
        os.mkdir(output_path)
    except:
        pass
    input_path = "/mnt/array2/smallproteins/eCAMI/hmm_analysis_combinded/cut_domain_seq/"
    for folder in folders:
        process_path = input_path + folder + '/'
        try:
            output_path = "/mnt/array2/smallproteins/eCAMI/hmm_analysis_combinded/hmm/" + folder + '/'
            os.mkdir(output_path)
        except:
            pass
       
        files = os.listdir(process_path)
        for file in files:
            brief_file = folder + '_' + file.split('.fa')[0]
            file_path = process_path + file
            
            cdhit_file = output_path + brief_file + '.fa.nr'
            command_0 = "cd-hit -i " + file_path + " -o " + cdhit_file + " -c 0.6 -n 2 -p 1 -d 200 -l 5 -s  0.95 -aL 0.95 -g 1"
            os.system(command_0)


            aln_file_1 = output_path + brief_file + '.fa.aln'
            command_1 = "mafft " + cdhit_file + ' > ' + aln_file_1
            os.system(command_1)

            aln_file_2 = output_path +  brief_file + '.hmm.aln'
            command_2 = "ln " + aln_file_1 + " -s " + aln_file_2
            os.system(command_2)

            hmm_file = output_path + brief_file + '.hmm'
            command_3 = "hmmbuild --informat afa " + hmm_file + " " + aln_file_2
            os.system(command_3)
            
if __name__ == '__main__':
    main()
