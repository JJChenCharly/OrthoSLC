from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

import sys#, getopt
import argparse
import os
from pathlib import Path

from itertools import combinations
from multiprocessing import Process, Manager

from OrthoSLC import __version__, mission_spliter

p = argparse.ArgumentParser(
    prog="python3 TK_SNPmat.py",
    description="Thanks for using OrthoSLC! (version: " + __version__ + ")\n"
)

p.add_argument("-i", "--input_path",
               metavar='',
                required=True,
                help="<dir> path/to/input/directory")
p.add_argument("-o", "--output_csv",
               metavar='',
                required=True,
                help="<dir> path/to/output/file.csv")
p.add_argument("-T", "--ID_TSV",
               metavar='',
                required=True,
                help="<txt> path/to/Step1.txt")
p.add_argument("-u", "--resource_total",
               metavar='',
                default="1",
                help="<int> total threads available")

args, extras = p.parse_known_args()

kaligned_path = args.input_path

def get_distance(str1, str2):
    if str1 == str2:
        return 0
    d = 0
    l = len(str1)
    
    for c in range(l):
        if str1[c] != str2[c]:
            d = d + 1
            
    return d

def distance_in_one_cluster(in_path_ls, 
                            share_ls, 
                            loop_ls, 
                            total_len):
    for in_path in in_path_ls:
        in_fasta = SeqIO.to_dict(SeqIO.parse(kaligned_path + '/' + in_path,
                                             'fasta'))


        spe_dict = {x[0:x.index('-')]:x for x in in_fasta.keys()} # !!!!!!!!!!!!!

        result_dict = {}

        for x in loop_ls:
            snp=get_distance(str(in_fasta[spe_dict[x[0]]].seq),
                             str(in_fasta[spe_dict[x[1]]].seq)
                            )
            result_dict[x] = snp

        cluster_len = len(str(in_fasta[spe_dict[x[1]]].seq))

        total_len.append(cluster_len)
        share_ls.append(result_dict)

if __name__ == '__main__':
    procc_num = int(args.resource_total)
    
    mission_ls = mission_spliter(os.listdir(kaligned_path), procc_num)
    
    manager = Manager()
    return_dict_ls = manager.list()
    len_ls = manager.list()
    
    a_fasta = SeqIO.parse(kaligned_path + '/' + os.listdir(kaligned_path)[0], 'fasta')
    a_ls = [x.id[0: x.id.index('-')] for x in a_fasta]
    a_ls.sort()
    
    par_ls = list(combinations(a_ls, 2))
    
    for submissions in mission_ls:
    
        jobs = []

        p = Process(target = distance_in_one_cluster,
                    args = (submissions, 
                            return_dict_ls,
                            par_ls,
                            len_ls
                           )

                   )
        p.start()
        jobs.append(p)


    for z in jobs:
        z.join()
    return_dict_ls = list(return_dict_ls)

#add up all the snp number
sum_snp = {}
for key in return_dict_ls[0]: #prepare a dict to save snp number
    sum_snp[key] = 0
    
for i in return_dict_ls:
    for key in sum_snp:
        sum_snp[key] = sum_snp[key] + i[key]

SLC_path_1 = args.ID_TSV

df_SLC_1 = pd.read_csv(SLC_path_1, sep = "\t"
                       , header=None
                       , index_col = 0)
dict_kalign_id = df_SLC_1[1].to_dict()

sum_snp_4_2 = {(dict_kalign_id[int(k[0])], dict_kalign_id[int(k[1])]): v for k, v in sum_snp.items()}
sum_snp_4_2 = dict(sorted(sum_snp_4_2.items()))

strain_naams = list(df_SLC_1[1])
df_mat = pd.DataFrame(0
                      , index=strain_naams
                      , columns=strain_naams
                     )
df_mat.index = strain_naams
df_mat.columns = strain_naams

cbn = list(combinations(strain_naams,
                        2
                   ))
for pairs in cbn:
    if (pairs[0], pairs[1]) in sum_snp_4_2.keys():
        df_mat.loc[pairs[1], pairs[0]] = sum_snp_4_2[(pairs[0], pairs[1])]
    elif (pairs[1], pairs[0]) in sum_snp_4_2.keys():
        df_mat.loc[(pairs[1], pairs[0])] = sum_snp_4_2[(pairs[1], pairs[0])]

df_mat.iloc[1:, 0: -1].to_csv(args.output_csv)