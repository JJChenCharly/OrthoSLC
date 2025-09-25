import sys, getopt
import os
from multiprocessing import Process, Manager, Lock
# import numpy as np

from OrthoSLC import O_blast, mission_spliter, __version__

# set your own
# path to concatenated fasta
cated_dereped_fasta_path = ''
# directory path of blastdb
blastdb_dir_path = ''
# how many jobs to parallel
process_number = 1
# set e_value
E_value = '1e-5'
# set strand direction
strand = 'plus'
# set path to blast bin file
blast_bin_ = 'blastn'
# set blast output directory
blast_op_dir = ''
# memory efficient mode 
mem_eff_mode = False
# blast output format
outfmt = str("6 qseqid sseqid pident score evalue")
# blastp task, select from <'blastp' 'blastp-fast' 'blastp-short'>
blastp_task = ''
# blastp task <blastn' 'blastn-short' 'dc-megablast' 'megablast' 'rmblastn'>
blastn_task = ''

argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv, 
                               "c:i:d:o:e:s:u:f:t:T:h", 
                               ["path_to_blast =",
                                "query =",
                                "dir_to_dbs =",
                                "output_path =", 
                                "e_value =",
                                "strand =",
                                "blast_thread_num =",
                                "outfmt = ",
                                "blastp_task = ",
                                "blastn_task = ",
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
    elif opt in ['-s', '--strand']:
        strand = arg
    elif opt in ['-u', '--blast_thread_num']:
        process_number = int(arg)
    elif opt in ['-f', '--outfmt']:
        outfmt = arg
    elif opt in ['-t', '--blastp_task']:
        blastp_task = arg
    elif opt in ['-T', '--blastn_task']:
        blastn_task = arg

    elif opt in ['-h', '--help']:
        print("Thanks for using OrthoSLC! (version: " + __version__ + ")\n")
        print("Usage: python Step4_reciprocal_blast.py -i query.fasta -o output/ -d directory_of_dbs/ [options...]\n")
        print("options:\n")
        print("  -i or --query -------------> <fasta> path/to/dereped_cated.fasta from Step 3")
        print("  -d or --dir_to_dbs --------> <dir> path/to/directory/of/dbs by makeblastdb")
        print("  -o or --output_path -------> <dir> path/to/output/directory")
        print("  -c or --path_to_blast -----> <cmd_path> path/to/blastn or blastp, default: 'blastn'")
        print("  -e or --e_value -----------> <float> blast E value, default: 1e-5")
        print("  -s or --strand ------------> <str> select from <'both', 'minus', 'plus'>, default: plus")
        print("  -u or --blast_thread_num --> <int> blast thread number, default: 1")
        print("  -f or --outfmt ------------> <str> specify blast output format if needed, unspecified means `'6 qseqid sseqid pident score evalue'` as default")
        print("  -t or --blastp_task  ------> <str> select from <'blastp' 'blastp-fast' 'blastp-short'>, unspecified means `'blastp'` as default")
        print("  -T or --blastn_task  ------> <str> select from <'blastn' 'blastn-short' 'dc-megablast' 'megablast' 'rmblastn' >, unspecified means `'megablast'` as default")        
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

if 'blastp' in blast_bin_ and blastn_task != "":
    print("Error: set '-T or --blastn_task' for blastn mission only.")
    sys.exit()

# mkdir or not
if os.path.exists(blast_op_dir):
    pass
else:
    os.mkdir(blast_op_dir)
    
strain_naam_lst = os.listdir(blastdb_dir_path)
mission_lst = [os.path.join(blastdb_dir_path, x) for x in strain_naam_lst]

size = []
for strain_name, directory in zip(strain_naam_lst, mission_lst):
    representative = os.path.join(directory, f"{strain_name}.nsq")
    if not os.path.exists(representative):
        representative = os.path.join(directory, f"{strain_name}.psq")
    if os.path.exists(representative):
        size.append(os.path.getsize(representative))
    else:
        files = os.listdir(directory)
        size.append(sum(os.path.getsize(os.path.join(directory, filename)) for filename in files))
sorted_mission_lst = [(_[1], _[0]) for _ in sorted(zip(size, strain_naam_lst), reverse=True)]
total_size = sum([_[1] for _ in sorted_mission_lst])


mission_lst = mission_spliter(sorted_mission_lst, 
                                    process_number
                                    )

blast_ = O_blast(
    op_path = blast_op_dir,
    query_path = cated_dereped_fasta_path,
    blast_bin_path = blast_bin_,
    blastdb_dir_path = blastdb_dir_path,
    e_value = E_value,
    outfmt = outfmt,
    strand = strand,
    blastn_task = blastn_task,
    blastp_task = blastp_task
)

if __name__ == "__main__":
    try:
        blast_.Psudo_pool_excute(
            sorted_mission_lst = sorted_mission_lst,
            total_size = total_size,
            available_resource = process_number
        )  
    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        sys.exit(1)
        
