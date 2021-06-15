import json
import os

import pandas as pd

def diff_detail(cazyfamily):
    with open('hmm_analysis_no_combinded/seq_summary/' + cazyfamily + '.json', 'r') as f:
        data = json.load(f)
    f.close()
    nocom = data['data']
    with open('hmm_analysis_combinded/seq_summary/' + cazyfamily + '.json', 'r') as f:
        data_1 = json.load(f)
    f.close()
    com = data_1['data']
    diff = list(set(com).difference(set(nocom)))

    com_re = 'compare_result/'
    try:
        os.mkdir(com_re)
    except:
        pass

    re_file = open(com_re + cazyfamily + "_diff.csv", "w")
    re_file.write('sequece_id, hmmscan_result, hmmscan_path \n')
    path = '/mnt/array2/smallproteins/eCAMI/hmm_analysis_combinded/domain_result/' + cazyfamily + '/'
    files = os.listdir(path)
    diff_origion = []
    for file in files:
        if '.ps' in file:
            parse_output_file = path + file
            with open(parse_output_file, 'r') as file:
                results = file.readlines()
            file.close()

            for result in results:
                result = result.split('\n')[0]
                attributes = result.split()
                query_id = attributes[2].split('|')[0]
                if query_id in diff:
                    diff_origion.append(str(query_id) + ',' + result + ',' + parse_output_file + '\n')
    re_file.writelines(diff_origion)
    re_file.close()

def compare_result():
    with open('hmm_analysis_combinded/number_protein.json', 'r') as f:
        com_number = json.load(f)
    f.close()
    with open('hmm_analysis_no_combinded/number_protein.json', 'r') as f:
        nocom_number = json.load(f)
    f.close()

    with open('ori_seq.json', 'r') as f:
        ori_number = json.load(f)
    f.close()
    diff = []
    com = 0
    uncom = 0
    ori = 0
    for cazyfamily in com_number.keys():
        print(cazyfamily)
        ori += ori_number[cazyfamily]
        com += com_number[cazyfamily]
        uncom += nocom_number[cazyfamily]
        redu_value = ori_number[cazyfamily] - com_number[cazyfamily]
        diff_value = com_number[cazyfamily] - nocom_number[cazyfamily]
        if diff_value > 0:
            diff_detail(cazyfamily)
        diff.append([cazyfamily, ori_number[cazyfamily], com_number[cazyfamily], nocom_number[cazyfamily], redu_value, diff_value, round(redu_value/ori_number[cazyfamily], 2)])
    diff.append(['summary', ori, com, uncom, ori - com, com - uncom, round((ori - com) / ori, 2)])
    df = pd.DataFrame(diff)
    df.to_csv("diff.csv", index=None, header=['cazyfamily', 'original', 'combine', 'uncombine', 'reduce', 'diffe', 'reduce_rate'])

if __name__ == '__main__':
    compare_result()
