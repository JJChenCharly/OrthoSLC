# cat $wd/S7_op/* | mcl - --abc \
# -te $cpu -I 1.5 \
# -o - | \
# $bin_dir/cluster_fusion \
# $wd"/S3_op_pre_cluster.txt" \
# $wd/Final_cluster/

# -o $wd/Final_cluster/0.txt

import sys#, getopt
import argparse
import os
import subprocess

from OrthoSLC import __version__

argv = sys.argv[1:]

p = argparse.ArgumentParser(
    prog="python3 Step8_MCL.py",
    description="Thanks for using OrthoSLC! (version: " + __version__ + ")\n"
)
p.add_argument("-c", "--path_to_mcl",
               metavar='',
                default="mcl",
                help="<cmd_path> path/to/mcl")
p.add_argument("-i", "--input_path",
               metavar='',
                required=True,
                help="<dir> path/to/input/directory from Step 7")
p.add_argument("-o", "--output_path",
               metavar='',
                required=True,
                help="<dir> path/to/output/directory")
p.add_argument("-p", "--pre_cluster_path",
               metavar='',
                required=True,
                help="<txt> path/to/pre_cluster.txt from Step 3")
p.add_argument("-u", "--mcl_thread_num",
               metavar='',
                default="1",
                help="<int> number of threads for mcl")
# p.add_argument("-I", "--inflation",
#                 default="1.5",
#                 help="MCL inflation value")
# p.add_argument("--abc",
#                 action="store_true",
#                 help="pass the --abc flag to mcl")
# parse_known_args returns (Namespace, [extras])
args, extras = p.parse_known_args()

# sanity checks
if not os.path.isdir(args.input_path):
    p.error(f"input_path '{args.input_path}' does not exist or is not a directory")
if not os.path.isdir(os.path.dirname(args.output_path)):
    p.error(f"parent of output_path '{args.output_path}' does not exist")
if not os.path.isfile(args.pre_cluster_path):
    p.error(f"pre_cluster_path '{args.pre_cluster_path}' not found")


# elif opt in ['-h', '--help']:
#     print("Thanks for using OrthoSLC! (version: " + __version__ + ")\n")
#     print("Usage: python Step8_MCL.py -i input/ -o output/ [options...]\n")
#     print("options:\n")
#     print("  -i or --input_path --------> <dir> path/to/input/directory from Step 7")
#     print("  -o or --output_path -------> <dir> path/to/output/directory")
#     print("  -c or --path_to_mcl -------> <cmd_path> path/to/mcl, default: 'mcl'")
#     print("  -p or --pre_cluster_path --> <txt> path/to/pre_cluster.txt from Step 3")
#     print("  -u or --mcl_thread_num ----> <int> mcl thread number, default: 1")
#     print("  -h or --help --------------> display this information")
#     sys.exit()


# build the mcl command
# print(os.getcwd())
os.makedirs(args.output_path, exist_ok=True)

cmd = []
cmd.extend(["cat", os.path.join(args.input_path, "*"), '|'])
cmd.extend([args.path_to_mcl, '-', "-te", args.mcl_thread_num])

cmd.extend(extras)

script_dir = os.path.dirname(os.path.abspath(__file__))
# print(script_dir)
cmd.extend(['-o', '-', '|',
            os.path.join(script_dir, 'cluster_fusion'),  
            args.pre_cluster_path,
            args.output_path])
# print(cmd)

try:
    result = subprocess.run(' '.join(cmd), check=True, capture_output=True, text=True, shell=True)
except Exception as e:
    error_msg = f"Error: {cmd} -->\n Err CODE:\n {e.returncode}\n{e.stderr}"
    print(error_msg)