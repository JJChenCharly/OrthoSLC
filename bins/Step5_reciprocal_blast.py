import sys, getopt
import os
import subprocess
from multiprocessing import Process
import numpy as np

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

argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv, 
                               "c:i:d:o:e:u:m:h", 
                               ["path_to_blast =",
                                "query =",
                                "dir_to_dbs =",
                                "output_path =", 
                                "e_value =",
                                "blast_thread_num =",
                                "mem_eff_mode = "
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
            print("option ' -m or --mem_eff_mode' accept only <yes> or <no>.")
            sys.exit()

    elif opt in ['-h', '--help']:
        print("Thanks for using OrthoSLC! (version: 0.1)\n")
        print("Usage: python Step5_reciprocal_blast.py -i query.fasta -o output/ -d directory_of_dbs/ [options...]\n")
        print("options:\n")
        print(" -i or --query ---------------> path/to/concatenated.fasta")
        print(" -d or --dir_to_dbs ----------> path/to/directory/of/dbs")
        print(" -o or --output_path ---------> path/to/output/directory")
        print(" -c or --path_to_blast -------> path/to/output/blastn or blastp, default: 'blastn'")
        print(" -e or --e_value -------------> blast E value, default: 1e-5")
        print(" -u or --blast_thread_num ----> blast thread number, default: 1")
        print(" -m or --mem_eff_mode --------> using memory efficient mode or not, select from <'on' or 'off'>, default: off")
        print(" -h or --help ----------------> display this information")
        sys.exit()

if not os.path.exists(cated_dereped_fasta_path):
    print("Error: path provided to '-i or --input_path' do not exist. 路径不存在")
    sys.exit()
elif not os.path.exists(blastdb_dir_path):
    print("Error: path provided to '-d or --dir_to_dbs' do not exist. 路径不存在")
    sys.exit()

def sequential_blast_low(catted_fasta_pth,
                         db_pth, # alread abs
                         blast_res_dir, # already abs
                         blast_bin,
                         e_value, 
                         t_num):
#   get db name to blast
    strain_naam = db_pth.split('/')[-1]
    
    db_name = os.path.join(db_pth, strain_naam)
    save_dest = os.path.join(blast_res_dir, strain_naam + '.tab')
        
    subprocess.run([blast_bin, 
                    "-query", catted_fasta_pth,
                    "-db", db_name,
                    "-out", save_dest,
                    "-evalue", str(e_value),
                    "-max_hsps", "1",
                    "-dust", "no",
                    "-outfmt", str("6 qseqid sseqid score"),
                    "-mt_mode", "1",
                    "-num_threads", str(t_num)
                   ])
    
def sequential_blast_high(catted_fasta_pth,
                          db_pth, # alread abs
                          blast_res_dir, # already abs
                          blast_bin,
                          e_value):
#   get db name to blast
    strain_naam = db_pth.split('/')[-1]
    
    db_name = os.path.join(db_pth, strain_naam)
    save_dest = os.path.join(blast_res_dir, strain_naam + '.tab')
        
    subprocess.run([blast_bin, 
                    "-query", catted_fasta_pth,
                    "-db", db_name,
                    "-out", save_dest,
                    "-evalue", str(e_value),
                    "-max_hsps", "1",
                    "-dust", "no",
                    "-outfmt", str("6 qseqid sseqid score"),
                    "-mt_mode", "1",
                    "-num_threads", "1"
                   ])
def sequential_blast_high_s(catted_fasta_pth_,
                            db_pth_lst, # alread abs
                            blast_res_dir_, # already abs
                            blast_bin_,
                            e_value_):
    for db in db_pth_lst:
        sequential_blast_high(catted_fasta_pth_,
                              db,
                              blast_res_dir_, 
                              blast_bin_, 
                              e_value_)

# mkdir or not
if os.path.exists(blast_op_dir):
    pass
else:
    os.mkdir(blast_op_dir)
    
    
mission_lst = os.listdir(blastdb_dir_path)
mission_lst = [os.path.join(blastdb_dir_path, x) for x in mission_lst]

if mem_eff_mode:
    for db in mission_lst:
        sequential_blast_low(cated_dereped_fasta_path, 
                         db, 
                         blast_op_dir, 
                         blast_bin_, 
                         E_value, 
                         process_number)
        
else:
    mission_lst = list(np.array_split(mission_lst, 
                                      process_number
                                     )
                      )
    
    if __name__ == "__main__":
        # mp
        jobs = []

        for db_lst in mission_lst:

            p = Process(target = sequential_blast_high_s,
                        args = (cated_dereped_fasta_path, 
                                db_lst, 
                                blast_op_dir,
                                blast_bin_,
                                E_value
                               )
                       )
            p.start()
            jobs.append(p)


        for z in jobs:
            z.join()