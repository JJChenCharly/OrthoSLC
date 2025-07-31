import sys#, getopt
import subprocess
import argparse
import os
import random
from time import sleep
from typing import List, Dict
from pathlib import Path
from multiprocessing import Lock

from OrthoSLC import BASE_Tasker, __version__, mission_spliter

argv = sys.argv[1:]

p = argparse.ArgumentParser(
    prog="python3 TK_mafft.py",
    description="Thanks for using OrthoSLC! (version: " + __version__ + ")\n"
)
p.add_argument("-c", "--path_to_mafft",
               metavar='',
                default="mafft",
                help="<cmd_path> path/to/mafft")
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


class Batch_mafft(BASE_Tasker):
    def __init__(self,
                 op_path: str,
                 USE_SHELL: bool
                ):
        super().__init__(op_path)
        self.USE_SHELL = USE_SHELL
        self.mafft_bin_path = args.path_to_mafft

        self.base_cmd_ = [
            self.mafft_bin_path
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
        cmd.extend(["--thread", str(assigned_resource)])
        
        if extras != []:
            cmd.extend(extras)

        cmd.extend([task_pack])
        cmd.extend([">", str(self.op_path/strain_naam.with_suffix('.fasta'))])

        if self.USE_SHELL:
            return ' '.join(cmd)
        else:
            return cmd
        
    # override for use shell
    def Psudo_pool_one_sub_excute(
            self,
            task_queue, # each element is a task pack. task_pack[0] is the task string, task_pack[1] is the task size
            resource_pool: Dict[str, int],
            resource_lock_: Lock,
            total_size: int,
            available_resource: int
            ) -> None:
        """动态任务调度器 - 按任务大小比例分配线程"""
        while not task_queue.empty():
            try:
                # 动态获取任务
                with resource_lock_:
                    if task_queue.empty():
                        return
                    task_pack = task_queue.get_nowait()
            except:
                return
            
            # 计算任务大小比例
            size_ratio = task_pack[1] / total_size
            
            # 计算所需线程数 (至少1个)
            thread_should_be_given = max(1, int(size_ratio * available_resource))
            # threads_needed = max(1, int(size_ratio * available_resource))
            
            # 等待直到有足够资源
            while True:
                with resource_lock_:
                    if resource_pool['available'] >= thread_should_be_given:
                        if resource_pool['tasks_left']<available_resource:
                            thread_given = max(thread_should_be_given, resource_pool['available']//resource_pool['tasks_left'])
                            # thread_given = max(thread_should_be_given, resource_pool['available'])
                        else:
                            thread_given = thread_should_be_given
                        resource_pool['available'] -= thread_given
                        break
                # wait and check
                sleep(random.randint(1, 2))

            
            cmd = self.one_task_cmd_lst_modifier(task_pack[0], thread_given)
            # 执行任务
            try:
                result = subprocess.run(cmd, check=True, capture_output=True, text=True,
                                        shell=self.USE_SHELL)
            except subprocess.CalledProcessError as e:
                error_msg = f"Error: {cmd} -->\n Err CODE:\n {e.returncode}\n{e.stderr}"
                print(error_msg)
            except TypeError as e:
                error_msg = f"TypeError occurred while running command: {cmd}\nError message: {e}"
                print(error_msg)
            except Exception as e:
                error_msg = f"Unexpected error occurred while running command: {cmd}\nError message: {e}"
                print(error_msg)
            
            # 归还资源
            with resource_lock_:
                resource_pool['available'] += thread_given
                resource_pool['tasks_left'] -= 1
                
            # 标记任务完成
            task_queue.task_done()
        
strain_naam_lst = os.listdir(args.input_path)
mission_lst = [os.path.join(args.input_path, x) for x in strain_naam_lst]
size = [os.path.getsize(i) for i in mission_lst]
sorted_mission_lst = [(_[1], _[0]) for _ in sorted(zip(size, mission_lst), reverse=True)]
total_size = sum([_[1] for _ in sorted_mission_lst])

process_num = int(args.resource_total)

B_mafft = Batch_mafft(
    op_path = args.output_path,
    USE_SHELL=True
)

if __name__ == "__main__":
    B_mafft.Psudo_pool_excute(
        sorted_mission_lst = sorted_mission_lst,
        total_size = total_size,
        available_resource = process_num
    )  
