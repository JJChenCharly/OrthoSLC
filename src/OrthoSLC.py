import subprocess
import os
from pathlib import Path
import random
from time import sleep
import copy
from typing import List, Tuple, Dict, Any, Union, Optional
from multiprocessing import Pool, Process, Lock, Manager

__version__ = "1.0.0"

# psudo pool
# 任务分割函数 (保持原样，不在类内)
def mission_spliter(lst: List,  # 待分割列表(最困难任务在前)
                    num_splits: int  # 需要分割成的份数
                    ) -> List[List]:
    """
    任务分割函数 - 使用轮询方式分配任务
    
    参数:
        lst: 待分割的原始列表
        num_splits: 需要分割成的子列表数量
    
    返回:
        分割后的子列表集合
    
    说明:
        1. 如果请求的分割数 > 列表长度，则自动调整为列表长度
        2. 使用轮询分配算法确保任务均衡分配
        3. 自动移除产生的空子列表
    """
    n = len(lst)
    if num_splits > n:
        num_splits = n

    # 创建空结果列表
    result = [[] for _ in range(num_splits)]

    # 轮询方式分配任务
    for i, item in enumerate(lst):
        index = i % num_splits
        result[index].append(item)

    # 过滤空列表
    return [sublist for sublist in result if sublist]

class BASE_Tasker:
    def __init__(self, 
                 op_path: str,
                 V:str = __version__,
                 ):
        self.version = V
        self.op_path = Path(op_path)
        self.base_cmd_: List[str] = []  # 默认为空列表

    def one_task_cmd_lst_modifier(self,
                                  task_pack: Tuple[Any, ...], 
                                  assigned_resource: int
                                  ) -> List[str]:
        """
        修改基础命令以适配具体任务
        
        参数:
            task_pack: 任务参数元组
            assigned_resource: 分配的资源数
            
        返回:
            修改后的命令列表
        """
        # 在子类中应重写此方法
        return self.base_cmd_.copy() + list(task_pack) + [str(assigned_resource)]
    
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
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            except Exception as e:
                error_msg = f"Error: {cmd} -->\n Err CODE:\n {e.returncode}\n{e.stderr}"
                print(error_msg)
            
            # 归还资源
            with resource_lock_:
                resource_pool['available'] += thread_given
                resource_pool['tasks_left'] -= 1
                
            # 标记任务完成
            task_queue.task_done()

    def Psudo_pool_excute(self,
                        sorted_mission_lst: List[Tuple[str, int]],
                        total_size: int,
                        available_resource: int
                        ) -> None:
        with Manager() as manager:
            # 创建任务队列
            task_queue = manager.Queue()
            for task_pack in sorted_mission_lst:
                task_queue.put(task_pack)

           # 创建共享资源池
            resource_pool = manager.dict()
            resource_pool['available'] = available_resource  # 初始可用线程数
            resource_pool['tasks_left'] = len(sorted_mission_lst)
            st_lock = manager.Lock()
            # 创建进程 (进程数=min(总任务数, 进程数))
            jobs = []
            actual_procs = min(len(sorted_mission_lst), available_resource)
            
            for _ in range(actual_procs):
                p = Process(target=self.Psudo_pool_one_sub_excute,
                            args=(task_queue,
                                  resource_pool,
                                  st_lock,
                                  total_size,
                                  available_resource
                                  )
                                  )
                jobs.append(p)
                p.start()
            
            # 等待所有任务完成
            task_queue.join()
            
            # 终止所有进程
            for p in jobs:
                p.join(timeout=1.0)
                if p.is_alive():
                    p.terminate()

    def Pool_worker(self,
                    task_pack: Tuple[Any, ...],  # 任务参数包
                    thread: int = 1
                    ) -> str:
        """
        进程池工作函数
        
        参数:
            task: 任务参数元组
            thread: 分配的线程数
            
        返回:
            执行结果消息
        """
        cmd = self.one_task_cmd_lst_modifier(task_pack, thread)
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            error_msg = f"Error: {cmd} -->\n Err CODE:\n {e.returncode}\n{e.stderr}"
            print(error_msg)
    
    def Pool_excute(self,
                    tasks: List[Tuple[Any, ...]], 
                    processes: int) -> None:
        """
        并行执行任务（使用标准进程池）
        
        参数:
            tasks: 任务列表
            processes: 使用的进程数
        """
        # 使用进程池执行任务
        with Pool(processes=processes) as pool:
            # 使用starmap传递额外参数
            results = pool.starmap(self.Pool_worker, tasks)

class makeblastdb(BASE_Tasker):
    def __init__(self,
                 op_path: str,
                 mbdb_bin_path: str = 'makeblastdb',
                 dbt: str = 'nucl'
                 ):
        super().__init__(op_path)
        self.mbdb_bin_path = mbdb_bin_path

        self.base_cmd_ = [
            self.mbdb_bin_path,
            "-dbtype", dbt, 
            # "-in", ,
            # "-out", ,
            "-parse_seqids"
        ]
    def one_task_cmd_lst_modifier(self,
                                  task_pack: str, 
                                  thread: int
                                  ) -> List[str]:
        """
        修改基础命令以适配具体任务
        
        参数:
            task_pack: 任务参数元组
            assigned_resource: 分配的资源数
            
        返回:
            修改后的命令列表
        """
        strain_naam = Path(task_pack).stem

        os.makedirs(Path(self.op_path)/strain_naam, exist_ok=True)
        cmd = self.base_cmd_.copy()
        cmd.extend(["-in", task_pack])
        cmd.extend(["-out", self.op_path/strain_naam/strain_naam])
        return cmd
    
class O_blast(BASE_Tasker):
    def __init__(self,
                 op_path: str,
                 query_path:str,
                 blast_bin_path: str = 'blastn',
                 blastdb_dir_path:str = '',
                 e_value:str = '1e-5',
                 outfmt:str = str("6 qseqid sseqid pident score evalue"),
                 strand:str = 'plus',
                 blastn_task:str = '',
                 blastp_task:str = ''
                 ):
        super().__init__(op_path)
        self.blast_bin_path = blast_bin_path
        self.query_path = query_path
        self.blastdb_dir_path = Path(blastdb_dir_path)

        self.base_cmd_ = [
            self.blast_bin_path, 
            "-query", self.query_path,
            # "-db", db_name,
            # "-out", save_dest,
            "-evalue", e_value,
            "-max_hsps", "1",
            # "-dust", "no",
            "-outfmt", outfmt,
            "-mt_mode", "1"
            # "-num_threads", str(t_num)
            ]
        if 'blastn' in self.blast_bin_path:
            self.base_cmd_.extend(["-dust", "no", "-strand", strand])
            if blastn_task == '':
                pass
            else:
                self.base_cmd_.extend(["-task", blastn_task])
        elif 'blastp' in self.blast_bin_path:
            if blastp_task == '':
                pass
            else:
                self.base_cmd_.extend(["-task", blastp_task])

    def one_task_cmd_lst_modifier(self,
                                  task_pack: str, 
                                  thread: int
                                  ) -> List[str]:
        strain_naam = Path(task_pack)

        # print(task_pack + ' go: ---> ' + str(thread))
        
        cmd = self.base_cmd_.copy()
        cmd.extend(["-db", self.blastdb_dir_path/strain_naam/strain_naam])
        cmd.extend(["-out", self.op_path/strain_naam.with_suffix('.tab')])
        cmd.extend(["-num_threads", str(thread)])

        return cmd

            