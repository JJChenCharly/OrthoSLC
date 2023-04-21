import sys, getopt
import os
from multiprocessing import Process, Manager, Lock
# import numpy as np

from Blast import BLAST

blast_entity = BLAST()

# set your own
# path to concatenated fasta
cated_dereped_fasta_path = ''
# directory path of blastdb
blastdb_dir_path = ''
# how many jobs to parallel
process_number = 1
# set e_value
E_value = '1e-5'
# set path to blast bin file
blast_bin_ = 'blastn'
# set blast output directory
blast_op_dir = ''
# memory efficient mode 
mem_eff_mode = False
# blast output format
outfmt = str("6 qseqid sseqid score")
# blastp task, select from <'blastp' 'blastp-fast' 'blastp-short'>
blastp_task = ''

argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv, 
                               "c:i:d:o:e:u:m:f:t:h", 
                               ["path_to_blast =",
                                "query =",
                                "dir_to_dbs =",
                                "output_path =", 
                                "e_value =",
                                "blast_thread_num =",
                                "mem_eff_mode = ",
                                "outfmt = ",
                                "blastp_task = ",
                                "help"]
                              )
    
except Exception as ex:
    print('Incorrect input command\n不正确指令')
    sys.exit()

for opt, arg in opts:
    if opt in ['-c', '--path_to_blast']:
        blast_bin_ = arg
    elif opt in ['-i', '--query']:
        cated_dereped_fasta_path = arg
    elif opt in ['-d', '--dir_to_dbs']:
        blastdb_dir_path = arg
    elif opt in ['-o', '--output_path']:
        blast_op_dir = arg
    elif opt in ['-e', '--e_value']:
        E_value = arg
    elif opt in ['-u', '--blast_thread_num']:
        process_number = int(arg)
    elif opt in ['-m', '--mem_eff_mode']:
        mem_eff_mode_arg = arg

        if mem_eff_mode_arg == "on":
            mem_eff_mode = True
        elif mem_eff_mode_arg == "off":
            mem_eff_mode = False
        else:
            print("option ' -m or --mem_eff_mode' accept only <on> or <off>.")
            sys.exit()
    elif opt in ['-f', '--outfmt']:
        outfmt = arg
    elif opt in ['-t', '--blastp_task']:
        blastp_task = arg

    elif opt in ['-h', '--help']:
        print("Thanks for using OrthoSLC! (version: " + blast_entity.version + ")\n")
        print("Usage: python Step4_reciprocal_blast.py -i query.fasta -o output/ -d directory_of_dbs/ [options...]\n")
        print("options:\n")
        print("  -i or --query -------------> <fasta> path/to/dereped_cated.fasta from Step 3")
        print("  -d or --dir_to_dbs --------> <dir> path/to/directory/of/dbs by makeblastdb")
        print("  -o or --output_path -------> <dir> path/to/output/directory")
        print("  -c or --path_to_blast -----> <cmd_path> path/to/output/blastn or blastp, default: 'blastn'")
        print("  -e or --e_value -----------> <float> blast E value, default: 1e-5")
        print("  -u or --blast_thread_num --> <int> blast thread number, default: 1")
        print("  -m or --mem_eff_mode ------> <on/off> using memory efficient mode or not, select from <'on' or 'off'>, default: off")
        print("  -f or --outfmt ------------> <str> specify blast output format if needed, unspecified means `'6 qseqid sseqid score'` as default")
        print("  -t or --blastp_task  ------> <str> specify blastp_task, select from <'blastp' 'blastp-fast' 'blastp-short'>, unspecified means `'blastp'` as default")
        print("  -h or --help --------------> display this information")
        sys.exit()

if not os.path.exists(cated_dereped_fasta_path):
    print("Error: path provided to '-i or --input_path' do not exist. 路径不存在")
    sys.exit()
elif not os.path.exists(blastdb_dir_path):
    print("Error: path provided to '-d or --dir_to_dbs' do not exist. 路径不存在")
    sys.exit()
elif not os.path.exists(os.path.dirname(blast_op_dir)):
    print("Error: parent path provided to '-o or --output_path' do not exist. 路径不存在")
    sys.exit()

if 'blastn' in blast_bin_ and blastp_task != "":
    print("Error: set '-t or --blastp_task' for blastp mission only.")
    sys.exit()

# mkdir or not
if os.path.exists(blast_op_dir):
    pass
else:
    os.mkdir(blast_op_dir)
    
    
mission_lst = os.listdir(blastdb_dir_path)
mission_lst = [os.path.join(blastdb_dir_path, x) for x in mission_lst]

cmd_lst = [blast_bin_, 
            "-query", cated_dereped_fasta_path,
            # "-db", db_name,
            # "-out", save_dest,
            "-evalue", str(E_value),
            "-max_hsps", "1",
            # "-dust", "no",
            "-outfmt", outfmt,
            "-mt_mode", "1"
            # "-num_threads", str(t_num)
            ]

if 'blastn' in blast_bin_:
    cmd_lst.extend(["-dust", "no"])
elif 'blastp' in blast_bin_ and blastp_task == "":
    cmd_lst.extend(["-task", "blastp"])
else:
    cmd_lst.extend(["-task", blastp_task])

if mem_eff_mode:
    blast_low = blast_entity.RBB

    cmd_lst.extend(["-num_threads", str(process_number)])

    for db in mission_lst:
        
        blast_low(cmd_lst,
                  db, 
                  blast_op_dir
                  )
        
else:
    mission_lst = BLAST.mission_spliter(mission_lst, 
                                        process_number
                                        )
    
    blast_high = blast_entity.RBB_high

    if __name__ == "__main__":
        with Manager() as manager:
            smart_threads_d = manager.dict()
            
            st_lock = Lock()
            
            for x in range(len(mission_lst)):
                smart_threads_d[x] = 1
        
            # mp
            jobs = []
    
            pro_id = 0
            for db_lst in mission_lst:

                p = Process(target = blast_high,
                            args = (cmd_lst,
                                    db_lst,
                                    blast_op_dir,
                                    pro_id,
                                    smart_threads_d, 
                                    st_lock,
                                    process_number,
                                   )
                           )
                p.start()
                jobs.append(p)
                
                pro_id += 1


            for z in jobs:
                z.join()
