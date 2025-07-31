import sys, getopt
import os
from multiprocessing import Process
# import numpy as np

from OrthoSLC import makeblastdb, __version__

# parameter parsing
# directory path of output
dereped_dir_path = None
# directory path of blastdb
blastdb_dir_path = None
# how many jobs to parallel
process_number = 1
# set dbtype: nucl or prot
dbt = 'nucl'
# set path to makblastdb bin file
mbdb_bin = 'makeblastdb'

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
        mbdb_bin = arg
    elif opt in ['-i', '--input_path']:
        dereped_dir_path = arg
    elif opt in ['-o', '--output_path']:
        blastdb_dir_path = arg
    elif opt in ['-u', '--thread_number']:
        process_number = int(arg)
    elif opt in ['-t', '--dbtype']:
        dbt = arg
    elif opt in ['-h', '--help']:
        print("Thanks for using OrthoSLC! (version: " + __version__ + ")\n")
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
    try:
        # mkdir or not
        if os.path.exists(blastdb_dir_path):
            pass
        else:
            os.mkdir(blastdb_dir_path)

        make_db = makeblastdb(op_path = blastdb_dir_path,
                            mbdb_bin_path = mbdb_bin,
                            dbt = dbt)
        
        # mp
        task_packs =[(os.path.join(dereped_dir_path, i), 1) for i in os.listdir(dereped_dir_path)]
        make_db.Pool_excute(task_packs, process_number)
    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        sys.exit(1)