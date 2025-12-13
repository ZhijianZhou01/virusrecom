# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Institution: Hunan University
Email: zjzhou@hnu.edu.cn
Time: 2022/5/20 20:18

"""
import sys
import os
import subprocess
import platform

from my_func import (make_dir, get_all_path, resolve_file_path)


class SeqAlign(object):

    def __init__(self,
                 query_lineage_path,
                 other_lineage_dir,
                 run_record,
                 run_id,
                 thread_num,
                 out_file):

        """
        Run the sequence alignment
        :param query_lineage_path: filepath of query sequence 
        :param other_lineage_dir:  dirpath of other lineages
        """

        super(SeqAlign, self).__init__()

        self.query_seq_path = query_lineage_path

        self.lineage_file_dir = other_lineage_dir

        self.run_record = run_record

        self.run_id = run_id

        self.thread_num = thread_num

        self.out_file = out_file

    def run(self):

        lineage_name_list = []

        query_seq_path = self.query_seq_path.replace("\\", "/")

        query_seq_dir, query_seq_prefix = resolve_file_path(query_seq_path)


        lineage_file_dir = self.lineage_file_dir.replace("\\", "/")

        lineage_file_list = get_all_path(lineage_file_dir)

        seq_for_mafft_path = (self.run_record + "/" + query_seq_prefix
                              + "_" + self.run_id
                              + "_merge.fasta")

        seq_for_mafft_file = open(seq_for_mafft_path,"w",encoding="utf-8")


        with open(query_seq_path) as query_seq:
            seq_for_mafft_file.write(query_seq.read().replace(">",">"
                                             + query_seq_prefix + "_")
                                     + "\n")


        for each_path in lineage_file_list:

            each_path = each_path.replace("\\","/")

            input_data_dir, out_prefix = resolve_file_path(each_path)
            lineage_name_list.append(out_prefix)

            with open(each_path,"r",encoding="utf-8") as lineage_file_input:

                seq_for_mafft_file.write(lineage_file_input.read().replace(">",">" + out_prefix + "_")
                                         + "\n")

        seq_for_mafft_file.close()



        print("Running MAFFT for sequence alignment..." + "\n")

        aligned_out_path = self.out_file 

        commd_list = [] 
        commd_list.append("mafft")

        commd_list.append("--inputorder")
        commd_list.append("--auto")
        commd_list.append("--thread " + str(self.thread_num))
        commd_list.append(seq_for_mafft_path + " > " + aligned_out_path)

        mafft_commd = " ".join(commd_list)

        print(mafft_commd)

        process = subprocess.Popen(mafft_commd,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   universal_newlines=True,
                                   shell=True)


        while True:

            output = process.stdout.readline()
            if process.poll() is not None:
                break

            elif output:
                print(output)

        process.terminate()

        print("Sequence alignment has been completed!" + "\n")


        return lineage_name_list
