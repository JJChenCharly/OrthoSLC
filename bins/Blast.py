import subprocess
import os
import random
from typing import List
import copy

__version__ = "0.2.2"

class BLAST:
    def __init__(self,
                 V = __version__):
        self.version = V

    @staticmethod
    def mission_spliter(lst: List, 
                        num_splits: int) -> List[List]:
        n = len(lst)
        random.shuffle(lst)
        
        split_size = n // num_splits
        remainder = n % num_splits
        start = 0
        result = []
        
        for i in range(num_splits):
            end = start + split_size
            if i < remainder:
                end += 1
            if lst[start:end] == []:
                break
            result.append(lst[start:end])
            start = end
        
        return result
    
    def make_a_db(self,
                  in_abs_pth, 
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
                        ], 
                        stdout = subprocess.DEVNULL
                        )
        
    def make_dbs(self,
                 in_abs_pth_lst,
                 db_dir_abs_pth_,
                 dbtype_,
                 makeblastdb_bin_):
        for s in in_abs_pth_lst:
            
            self.make_a_db(s, 
                           db_dir_abs_pth_, 
                           dbtype_,
                           makeblastdb_bin_)
            
    def RBB(self,
            cmd_lst_,
            db_pth, # alread abs
            blast_res_dir, # already abs
            ):
        strain_naam = db_pth.split('/')[-1]
    
        db_name = os.path.join(db_pth, strain_naam)
        save_dest = os.path.join(blast_res_dir, strain_naam + '.tab')
        
        cmd_lst_l = copy.deepcopy(cmd_lst_)

        cmd_lst_l.extend([
            "-db", db_name,
            "-out", save_dest
        ])

        subprocess.run(cmd_lst_l)
        
        
    def RBB_high(self,
                 cmd_lst_,
                 db_pth_lst, # alread abs
                 blast_res_dir_, # already abs
                 pro_id_, 
                 s_t_d,
                 st_lock_,
                 available_threads
                 ):
        for db in db_pth_lst:
            st_lock_.acquire()
            if sum(s_t_d.values()) >= available_threads:
                useable_th = s_t_d[pro_id_]
    #             print(str(pro_id_) + ' uses ' + str(useable_th))
            else:
                s_t_d[pro_id_] = s_t_d[pro_id_] + (available_threads - sum(s_t_d.values()
                                                                           )
                                                    )
                useable_th = s_t_d[pro_id_]
    #             print(str(pro_id_) + ' uses ' + str(useable_th))
            st_lock_.release()

            cmd_lst_l = copy.deepcopy(cmd_lst_)
                                                
            cmd_lst_l.extend(["-num_threads", str(useable_th)])

            self.RBB(
                cmd_lst_l,
                db, # alread abs
                blast_res_dir_, # already abs
            )
            
        st_lock_.acquire()
        s_t_d[pro_id_] = 0
    #     print(str(pro_id_) + ' finished')
        st_lock_.release()
            