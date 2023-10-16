import sys, getopt
import os
from multiprocessing import Process
# import numpy as np

from Blast import BLAST

blast_entity = BLAST()
make_dbs = blast_entity.make_dbs

# parameter parsing
# directory path of output
dereped_dir_path = ''
# directory path of blastdb
blastdb_dir_path = ''
# how many jobs to parallel
process_number = 1
# set dbtype: nucl or prot
dbt = 'nucl'
# set path to makblastdb bin file
mbdb = 'makeblastdb'

argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv, 
                               "c:i:o:u:t:h", 
                               ["path_to_makeblastdb =",
                                "input_path =",
                                "output_path =",
                                "thread_number =", 
                                "dbtype =",
                                "help"]
                              )
    
except Exception as ex:
    print('Incorrect input command\n不正确指令')
    sys.exit()

for opt, arg in opts:
    if opt in ['-c', '--path_to_makeblastdb']:
        mbdb = arg
    elif opt in ['-i', '--input_path']:
        dereped_dir_path = arg
    elif opt in ['-o', '--output_path']:
        blastdb_dir_path = arg
    elif opt in ['-u', '--thread_number']:
        process_number = int(arg)
    elif opt in ['-t', '--dbtype']:
        dbt = arg
    elif opt in ['-h', '--help']:
        print("Thanks for using OrthoSLC! (version: " + blast_entity.version + ")\n")
        print("Usage: python Step4_makeblastdb.py -i input/ -o output/ [options...]\n")
        print("options:\n")
        print("  -i or --input_path -----------> <dir> path/to/input/directory of nr_genomes from Step 2")
        print("  -o or --output_path ----------> <dir> path/to/output/directory")
        print("  -c or --path_to_makeblastdb --> <cmd_path> path/to/makeblastdb, default: makeblastdb")
        print("  -u or --thread_number --------> <int> thread number, default: 1")
        print("  -t or --dbtype ---------------> <str> -dbtype <String, 'nucl', 'prot'>, default: nucl")
        print("  -h or --help -----------------> display this information")
        sys.exit()

if not os.path.exists(dereped_dir_path):
    print("Error: path provided to '-i or --input_path' do not exist. 路径不存在")
    sys.exit()
elif not os.path.exists(os.path.dirname(blastdb_dir_path)):
    print("Error: parent path provided to '-o or --output_path' do not exist. 路径不存在")
    sys.exit()

if __name__ == "__main__":
    
    # mkdir or not
    if os.path.exists(blastdb_dir_path):
        pass
    else:
        os.mkdir(blastdb_dir_path)
    
    mission_lst = os.listdir(dereped_dir_path)
    
    mission_lst = [os.path.join(dereped_dir_path, x) for x in mission_lst]
    
    mission_lst = BLAST.mission_spliter(mission_lst, 
                                        process_number
                                        )
    
    # mp
    jobs = []
    
    for sub_mission_lst in mission_lst:
    
        p = Process(target = make_dbs,
                    args = (sub_mission_lst, 
                            blastdb_dir_path, 
                            dbt,
                            mbdb
                           )
                   )
        p.start()
        jobs.append(p)


    for z in jobs:
        z.join()