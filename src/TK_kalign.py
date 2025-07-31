# # run ----
# base='./'
# kalign_334='/data/share_data/Softwares/kalign/BIN/bin/kalign'

# # 提前make 输出的文件夹/ mkdir kalign
# time parallel -j 60 $kalign_334 \
# -i $base/orthoslc/S9_write_fasta/strict_core/{} \
# -o $base/kalign/{} --type dna -n 1 \
# ::: `ls $base/orthoslc/S9_write_fasta/strict_core` > $base/kaligned_log.txt 2>&1

import sys#, getopt
import argparse
import os
from typing import List
from pathlib import Path


from OrthoSLC import BASE_Tasker, __version__, mission_spliter

argv = sys.argv[1:]

p = argparse.ArgumentParser(
    prog="python3 TK_kalign.py",
    description="Thanks for using OrthoSLC! (version: " + __version__ + ")\n"
)
p.add_argument("-c", "--path_to_kalign",
               metavar='',
                default="kalign",
                help="<cmd_path> path/to/kalign")
p.add_argument("-i", "--input_path",
               metavar='',
                required=True,
                help="<dir> path/to/input/directory")
p.add_argument("-o", "--output_path",
               metavar='',
                required=True,
                help="<dir> path/to/output/directory")
p.add_argument("-u", "--resource_total",
               metavar='',
                default="1",
                help="<int> total threads available")

args, extras = p.parse_known_args()

# sanity checks
if not os.path.isdir(args.input_path):
    p.error(f"input_path '{args.input_path}' does not exist or is not a directory")
if not os.path.isdir(os.path.dirname(args.output_path)):
    p.error(f"parent of output_path '{args.output_path}' does not exist")
    sys.exit(2)

os.makedirs(args.output_path, exist_ok=True)

class Batch_kalign(BASE_Tasker):
    def __init__(self,
                 op_path: str,
                ):
        super().__init__(op_path)
        self.kalign_bin_path = args.path_to_kalign

        self.base_cmd_ = [
            self.kalign_bin_path
        ]

    def one_task_cmd_lst_modifier(self,
                                  task_pack: str, 
                                  assigned_resource: int,
                                  ) -> List[str]:
        strain_naam = Path(task_pack).stem
        strain_naam = Path(strain_naam)

        # print(task_pack + ' go: ---> ' + str(thread))

        # os.makedirs(self.op_path/strain_naam, exist_ok=True)
        
        cmd = self.base_cmd_.copy()
        cmd.extend(["-o", str(self.op_path/strain_naam.with_suffix('.fasta'))])
        cmd.extend(["-n", str(assigned_resource)])

        cmd.extend(["-i", task_pack])

        if extras != []:
            cmd.extend(extras)

        # if self.USE_SHELL:
        #     return ' '.join(cmd)
        # else:
        return cmd
        
strain_naam_lst = os.listdir(args.input_path)
mission_lst = [os.path.join(args.input_path, x) for x in strain_naam_lst]
size = [os.path.getsize(i) for i in mission_lst]
sorted_mission_lst = [(_[1], _[0]) for _ in sorted(zip(size, mission_lst), reverse=True)]
total_size = sum([_[1] for _ in sorted_mission_lst])

process_num = int(args.resource_total)

B_kalign = Batch_kalign(
    op_path = args.output_path,
)

if __name__ == "__main__":
    B_kalign.Psudo_pool_excute(
        sorted_mission_lst = sorted_mission_lst,
        total_size = total_size,
        available_resource = process_num
    )  
