import sys, getopt
import os
import subprocess
from multiprocessing import Process
import numpy as np

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
        print("Thanks for using OrthoSLC! (version: 0.1.1)\n")
        print("Usage: python Step5_makeblastdb.py -i input/ -o output/ [options...]\n")
        print("options:\n")
        print(" -i or --input_path ----------> path/to/input/directory")
        print(" -o or --output_path ---------> path/to/output/directory")
        print(" -c or --path_to_makeblastdb -> path/to/output/makeblastdb, default: makeblastdb")
        print(" -u or --thread_number -------> thread number, default: 1")
        print(" -t or --dbtype --------------> -dbtype <String, 'nucl', 'prot'>, default: nucl")
        print(" -h or --help ----------------> display this information")
        sys.exit()

if not os.path.exists(dereped_dir_path):
    print("Error: path provided to '-i or --input_path' do not exist. 路径不存在")
    sys.exit()
elif not os.path.exists(os.path.dirname(blastdb_dir_path)):
    print("Error: parent path provided to '-o or --output_path' do not exist. 路径不存在")
    sys.exit()

def make_a_db(in_abs_pth, 
              db_dir_abs_pth, 
              dbtype,
              makeblastdb_bin):
#   if db folder exist
    strain_naam = os.path.basename(in_abs_pth)
    strain_naam = os.path.splitext(strain_naam)[0]
    
    save_dest = os.path.join(db_dir_abs_pth, strain_naam)
    
    if os.path.exists(save_dest):
        pass
    else:
        os.mkdir(save_dest)
        
    save_naam = os.path.join(save_dest, strain_naam)
        
    subprocess.run([makeblastdb_bin, 
                    "-dbtype", dbtype, 
                    "-in", in_abs_pth,
                    "-out", save_naam,
                    "-parse_seqids"
                   ])
    
def make_dbs(in_abs_pth_lst,
             db_dir_abs_pth_,
             dbtype_,
             makeblastdb_bin_):
    for s in in_abs_pth_lst:
        
        make_a_db(s, 
                  db_dir_abs_pth_, 
                  dbtype_,
                  makeblastdb_bin_)

if __name__ == "__main__":
    
    # mkdir or not
    if os.path.exists(blastdb_dir_path):
        pass
    else:
        os.mkdir(blastdb_dir_path)
    
    mission_lst = os.listdir(dereped_dir_path)
    
    mission_lst = [os.path.join(dereped_dir_path, x) for x in mission_lst]
    
    mission_lst = list(np.array_split(mission_lst, 
                                      process_number
                                     )
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