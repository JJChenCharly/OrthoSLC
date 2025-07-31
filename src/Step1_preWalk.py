import sys, getopt
import os

from OrthoSLC import __version__

def find_files(input_dir, ext):
    """Walk input_dir; yield (parent_basename, absolute_file_path) for each file ending with .ext"""
    input_dir = os.path.abspath(input_dir)
    for root, dirs, files in os.walk(input_dir):
        for fn in files:
            if fn.lower().endswith('.' + ext.lower()):
                full = os.path.abspath(os.path.join(root, fn))
                if root == input_dir:
                    # file lives in the top-level folder: use its own name (no ext) as “parent”
                    parent = os.path.splitext(fn)[0]
                else:
                    # file lives in a subdir: use the subdir’s basename
                    parent = os.path.basename(root)
                yield parent, full


def usage(version = __version__):
    print("Thanks for using OrthoSLC! (version: " + version + ")\n")
    print("Usage: python Step1_preWalk.py -i input/ -o output/ [options...]\n")
    print("options:\n")
    print("  -i or --input_path ---> <dir> path/to/input/directory")
    print("  -o or --output_path --> <txt> path/to/output_file")
    print("  -f or --format -------> <str> file extension to look for (e.g. 'sqn', 'ffn')")
    print("  -h or --help ---------> display this information")

    
input_path = None
output_path = None
fmt = None

argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(
        argv,
        "i:o:f:h",
        ["input_path=", 
         "output_path=", 
         "format=", "help"]
    )
except Exception as ex:
    print('Incorrect input command\n不正确指令')
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-i", "--input_path"):
        input_path = arg
    elif opt in ("-o", "--output_path"):
        output_path = arg
    elif opt in ("-f", "--format"):
        fmt = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(0)

if not input_path or not output_path or not fmt:
    usage()
    sys.exit(2)

if not os.path.isdir(input_path):
    print(f"Error: input path '{input_path}' does not exist or is not a directory.", file=sys.stderr)
    sys.exit(1)

out_dir = os.path.dirname(output_path) or "."
if not os.path.isdir(out_dir):
    print(f"Error: parent directory of output file '{out_dir}' does not exist.", file=sys.stderr)
    sys.exit(1)

# Collect matching files
entries = list(find_files(input_path, fmt))
if not entries:
    print(f"No *.{fmt} files found under {input_path}", file=sys.stderr)
    sys.exit(1)

# Write output
with open(output_path, 'w') as out:
    for idx, (parent, path) in enumerate(entries):
        num_str = str(idx)
        # if it's exactly 4 digits, pad to 5
        if len(num_str) == 4:
            num_str = num_str.zfill(5)
        out.write(f"{num_str}\t{parent}\t{path}\n")
