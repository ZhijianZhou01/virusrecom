# -*- coding: utf-8 -*-

"""
Virus recombination(VirusRecom), detecting recombination of viral lineages using information theory
Author: Zhou Zhi-Jian
Institution: Hunan University
Email: zjzhou@hnu.edu.cn
Copyright：Copyright (c) Zhi-Jian Zhou
Time: 2022/4/17 17:29

"""

import sys
import os
import psutil

import numpy as np
import argparse

import pandas as pd

from datetime import datetime

from my_func import (resolve_file_path,make_dir,
                     mark_linegaes,export_seqname,
                     read_seq,calEnt, calEnt_gap)

from sequence_align import SeqAlign

from corecp import virus_infor_calculate


app_dir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))

app_dir = app_dir.replace("\\", "/")

sys.path.append(app_dir)



example_use = r'''
----------------☆ Example of use ☆-----------------

  (1) If the input sequence-data has been aligned:
      virusrecom -a alignment.fasta -q query_name -l lineage_name_list.txt -g n -m p -w 100 -s 20 -o outdir
      
  (2) If the input sequence-data was not aligned: 
      virusrecom -ua unalignment.fasta -at mafft -q query_name -l lineage_name_list.txt -g n -m p -w 100 -s 20 -t 2 -o outdir
      
  (3) If you have a *csv file with site’s WIC value:
      virusrecom -iwic *csv -g n (or -g y) -q query_name -w 100 -s 20 -o outdir
  
  Tip: the input-file and output-directory is recommended absolute path.
  
  Above is just a conceptual example, detailed usage in website: https://github.com/ZhijianZhou01/virusrecom
  
----------------------☆  End  ☆---------------------

'''


def calculate_memory():
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    memory = info.uss / 1024 / 1024
    return memory



if __name__ == "__main__":

    print("\n" + "-------------------------------------------------")

    print("  Name: Virus Recombination (VirusRecom)")

    print(
        "  Description: Detecting recombination of viral lineages (or subtypes) using information theory.")

    print("  Version: 1.1 (2023-02-28)")

    print("  Author: Zhi-Jian Zhou")

    print("  Citation: Brief Bioinform. 2023 Jan 19;24(1)")

    print("-------------------------------------------------" + "\n")


    def parameter():
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            prog="virusrecom",
            description="",
            epilog=example_use)


        parser.add_argument(
            "-a", dest="alignment",
            help="Aligned sequence file (*.fasta). Note, each sequence name requires containing lineage mark.",
            default="")

        parser.add_argument(
            "-ua", dest="unalignment",
            help="Unaligned (non-alignment) sequence file (*.fasta). Note, each sequence name requires containing lineage mark.",
            default = "")

        parser.add_argument(
            "-at",dest="align_tool",
            help="Program used for multiple sequence alignments (MSA). Supporting mafft, muscle and clustalo, such as ‘-at mafft’.",
            default="")

        parser.add_argument(
            "-iwic", dest="input_wic",
            help="Using the already obtained WIC values of reference lineages directly by a *.csv input-file.",
            default="")


        parser.add_argument(
            "-q", dest = "query",
            help="Name of query lineage (usually potential recombinant), such as ‘-q xxxx’. Besides, ‘-q auto’ can scan all lineages as potential recombinant in turn.",
            default="")


        parser.add_argument(
            "-l", dest="lineages",
            help="Path of a text-file containing multiple lineage marks.",
            default="")


        parser.add_argument("-g", dest="gap",
                            help="Reserve sites containing gap in subsequent analyses? ‘-g y’ means to reserve, and ‘-g n’ means to delete.",
                            type=str,
                            default="")


        parser.add_argument("-m", dest="method",
                            help="Method for scanning. ‘-m p’ means use polymorphic sites only, ‘-m a’ means use all the sites.",
                            type=str,
                            default="p")


        parser.add_argument("-w", dest="window",
                            help="Number of nucleotides sites per sliding window. Note: if the ‘-m p’ has been used, -w refers to the number of polymorphic sites per windows.",
                            type=int,
                            default=100)

        parser.add_argument("-s", dest="step",
                            help="Step size of the sliding window. Note: if the ‘-m p’ has been used, -s refers to the number of polymorphic sites per jump.",
                            type=int,
                            default=20)

        # max_recom_fragment
        parser.add_argument("-mr", dest="max_region",
                            help="The maximum allowed recombination region. Note: if the ‘-m p’ method has been used, it refers the maximum number of polymorphic sites contained in a recombinant region.",
                            type=int,
                            default=1000)

        parser.add_argument("-cp", dest="percentage",
                            help="The cutoff threshold of proportion (cp, default: 0.9) used for searching recombination regions when mWIC/EIC >= cp, the maximum value of cp is 1.",
                            type=float,
                            default=0.9)

        parser.add_argument("-cu", dest="cumulative",
                            help="Simply using the max cumulative WIC of all sites to identify the major parent. Off by default. If required, specify ‘-cu y’.",
                            type=str,
                            default="n")


        parser.add_argument("-b", dest="breakpoint",
                            help="Possible breakpoint scan of recombination. ‘-b y’ means yes, ‘-b n’ means no. Note: this option only takes effect when ‘-m p’ has been specified.",
                            type=str,
                            default="n")


        parser.add_argument("-bw", dest="breakwin",
                            help="The window size (default: 200) used for breakpoint scan, and step size is fixed at 1.",
                            type=int,
                            default=200)


        parser.add_argument(
            "-t", dest="thread",
            help="Number of threads (default: 1) used for MAS.",
            type=int,
            default=1)


        parser.add_argument(
            "-y", dest="y_start",
            help="Starting value (default: 0) of the Y-axis in plot diagram.",
            type=float,
            default=0.0)

        parser.add_argument(
            "-le", dest="legend",
            help="The location of the legend, the default is adaptive. '-le r' indicates placed on the right.",
            default="auto")


        # upper right


        parser.add_argument(
            "-owic", dest="only_wic",
            help="Only calculate site WIC value. Off by default. If required, please specify ‘-owic y’.",
            type=str,
            default="n")


        parser.add_argument(
            "-e",dest="engrave",
            help="Engraves file name to sequence names in batches. By specifying a directory containing one or multiple sequence files (*.fasta).",
            default="")


        parser.add_argument(
            "-en", dest="export_name",
            help="Export all sequence name of a *.fasta file.",
            default="")


        parser.add_argument(
            "-o", dest="outdir",
            help="Output directory to store all results.",
            type=str,
            default="")


        myargs = parser.parse_args(sys.argv[1:])

        return myargs


    myargs = parameter()

    start = datetime.today().now()

    start_memory = calculate_memory()


    #  parameter

    parameter_dic = {}

    parameter_dic["aligned_seq"] = myargs.alignment  # path of the aligned sequence file

    parameter_dic["unaligned_seq"] = myargs.unalignment   # path of the unaligned sequence file

    parameter_dic["align_tool"] = myargs.align_tool  #  alignment tool used for no-aligned sequence file

    parameter_dic["input_wic"] = myargs.input_wic            # input_wic, if you have

    parameter_dic["query_lineage_name"] = myargs.query         #  the name of query lineage

    parameter_dic["lineage_file"] = myargs.lineages          # other lineages

    parameter_dic["gaps_use"] = myargs.gap    #  strategy of handling gap

    parameter_dic["method"] = myargs.method   #  Scan method

    parameter_dic["windows_size"] = myargs.window    #  window size

    parameter_dic["step_size"] = myargs.step     #  step size

    parameter_dic["max_recom_fragment"] = myargs.max_region  #  maximum recombination region

    parameter_dic["recom_percentage"] = myargs.percentage  #  the specific cutoff threshold of proportion

    parameter_dic["cumulative_major"] = myargs.cumulative

    parameter_dic["breakpoints"] =  myargs.breakpoint  #  whether to search for recombination breakpoints

    parameter_dic["breakwins"] = myargs.breakwin  #  Window size used to search for recombination breakpoints

    parameter_dic["thread_num"] = myargs.thread  #  thread of MAS

    parameter_dic["y_start"] = myargs.y_start  #  Y-axis starting point when plotting

    parameter_dic["legend"] = myargs.legend

    parameter_dic["only_wic"] = myargs.only_wic

    parameter_dic["engrave_name"] = myargs.engrave

    parameter_dic["export_name"] = myargs.export_name

    parameter_dic["out_dir"] = myargs.outdir  # output directory



    # 处理不正确的输入

    # 检查是否指定了输出文件夹

    out_dir = parameter_dic["out_dir"].replace("\\", "/")   # 直接在这里就获取，并进行转义，不要改了

    if out_dir == "":
        print(
            ">>> Error！Please specify an output directory by '-o' options." + "\n")

        exit(0)


    make_dir(out_dir)  # 创建输出文件夹


    # 先判断用户是否只是简单标记谱系，而不是做分析
    if parameter_dic["engrave_name"] != "":   # 这里不要使用os.path.isdir()，否则不能输入相对路径

        mark_linegaes(parameter_dic["engrave_name"],
                      parameter_dic["out_dir"])

        exit(0)

    elif parameter_dic["export_name"] !="":  # 用户只是想提取序列的名称
        export_seqname(parameter_dic["export_name"],
                       parameter_dic["out_dir"])

        exit(0)


    print(myargs)

    print("\n")


    # 检查是否有输入的序列文件
    if (parameter_dic["unaligned_seq"] == ""
          and parameter_dic["aligned_seq"] == ""
          and parameter_dic["input_wic"] == ""):

        print(">>> Error！There is no input file." + "\n")

        exit(0)

    # 检查是否指定了查询谱系
    elif parameter_dic["query_lineage_name"] == "":
        print(">>> Error！Please specify the query lineage by '-q' options." + "\n")

        exit(0)


    elif parameter_dic["gaps_use"] == "":
        print(
            ">>> Error！Please specify the value of '-g' options." + "\n")

        exit(0)


    # query_seq_dir, query_prefix = resolve_file_path(query_seq_path)

    run_record = out_dir + "/" + "run_record"

    make_dir(run_record)

    seq_alignment_file = ""     # 对齐的序列文件



    ###  1.处理输入数据类型

    if parameter_dic["aligned_seq"] != "":        #  如果输入的序列是对齐的

        seq_alignment_file = parameter_dic["aligned_seq"].replace("\\", "/")


    elif parameter_dic["unaligned_seq"] != "":   # 如果序列事先没对齐

        if parameter_dic["align_tool"] !="":

            input_seq_dir, input_prefix = resolve_file_path(parameter_dic["unaligned_seq"])

            # 调用序列比对
            seq_alignment_file = (run_record + "/" + input_prefix
                                  + "_" + parameter_dic["align_tool"]
                                  + ".fasta")


            SeqAlign(parameter_dic["align_tool"],
                     parameter_dic["thread_num"],
                     parameter_dic["unaligned_seq"],
                     seq_alignment_file)

            if os.path.exists(seq_alignment_file):
                print(">>> Sequence alignment has been completed!" + "\n")

            else:

                print(">>> Error！Multiple sequence alignments did not succeed." + "\n")

                exit(0)


        else:
            print(">>> Error！Please specify a tool used for sequnece alignment." + "\n")

            exit(0)


    elif parameter_dic["input_wic"] != "":          # 如果是直接输入一个计算好位点的wic值文件

        out_dir = parameter_dic["out_dir"]

        query_prefix = parameter_dic["query_lineage_name"]

        virus_infor_calculate(parameter_dic,
                              pd.DataFrame(),
                              query_prefix,
                              [],
                              out_dir)

        duration = datetime.today().now() - start  # 打印总运行时间

        end_memory = calculate_memory()

        used_memory = end_memory - start_memory

        print("\n")

        print(">>> " + f"Occupied {used_memory}MB memory in total." + "\n")

        print(">>> " + "Take " + str(duration) + " seconds in total." + "\n")

        exit(0)


    else:   #  如果什么也没输入
        exit(0)


    lineage_mark_list = []   # 参考谱系名称列表

    query_prefix = parameter_dic["query_lineage_name"]    # 查询谱系的名称


    # 读取参考谱系名称
    with open(parameter_dic["lineage_file"]) as lineage_file:
        for line in lineage_file:
            line = line.strip()
            if line != "":
                lineage_mark_list.append(line)


    # 2.读取序列为数据框

    print(">>> " + "VirusRecom is importing sequence set, please wait patiently..." + "\n")

    seq_pd = read_seq(seq_alignment_file)   # 读取输入的对齐序列或者比对后的序列为数据框


    # 先从整体上考虑是否将含有gap的位点纳入重组分析

    seq_pd_clean = seq_pd     # 初始值考虑含有gap的位点

    parameter_dic["calEnt_use"] = calEnt_gap  # 初始值纳入gap进行分析


    if parameter_dic["gaps_use"].upper() == "N":  # 不保留gap

        print(">>> " + "VirusRecom is removing sites (columns) containing gap (-)..." + "\n")

        all_site_label = list(seq_pd.columns)

        seq_pd_clean = seq_pd[~seq_pd.isin(["-"])].dropna(axis=1)  # 去除含有gap的列（位点）


        parameter_dic["calEnt_use"] = calEnt  # 不考虑gap


        # 记录含有gap的位点（列）
        record_gap_path = (run_record + "/"
                           + "Record_of_deleted_gap_sites"
                           + ".txt")

        seq_pd_clean_site = list(seq_pd_clean.columns)

        with open(record_gap_path, "w",encoding="utf-8") as record_gap_file:

            record_gap_file.write("These sites with gap (-) in the file of "
                                  + seq_alignment_file + "\n")

            for each_site in all_site_label:
                if each_site not in seq_pd_clean_site:
                    record_gap_file.write("Site " + each_site + "\n")


        print(">>> " + "All sites containing gap have been cleared." + "\n")

    else:
        pass



    if parameter_dic["method"].upper() == "P":   # 只扫描多态位点(polymorphic)

        print(">>> " + "VirusRecom is extracting polymorphic sites..." + "\n")

        seq_pd_clean_site_old = list(seq_pd_clean)
        # print(seq_pd_clean_site_old)

        seq_pd_clean = seq_pd_clean.loc[:, (seq_pd_clean != seq_pd_clean.iloc[0]).any()]    # 去除没有突变的位点

        # 记录对齐序列中没有突变的位点
        seq_pd_clean_site_new = list(seq_pd_clean)


        record_same_sites_path = (run_record + "/"
                                  + "Record_of_same_sites_in_aligned_sequence"
                                  + ".txt")


        with open(record_same_sites_path, "w", encoding="utf8") as same_sites_file:
            same_sites_file.write("These same sites(no variation) in the file of "
                                  + seq_alignment_file + "\n")


            for each_site in seq_pd_clean_site_old:
                if each_site not in seq_pd_clean_site_new:
                    same_sites_file.write("Site " + each_site + "\n")

    else:
        pass



    analysis_site = np.size(seq_pd_clean,1)

    print(">>> " + "A total of " + str(analysis_site)

          + " sites will be included in the next analysis." + "\n")

    print(">>> " + "VirusRecom starts calculating weighted information content from each lineage..."
        + "\n")


    if query_prefix.upper() == "AUTO":    # 自动执行多个潜在重组体分析

        for each_lineage in lineage_mark_list:
            query_lineage_run = each_lineage

            new_lineage_name_list = []

            for lineage in lineage_mark_list:
                if lineage != query_lineage_run:
                    new_lineage_name_list.append(lineage)

            sub_outdir = out_dir + "/Results/" + query_lineage_run


            virus_infor_calculate(parameter_dic,
                                  seq_pd_clean,
                                  query_lineage_run,
                                  new_lineage_name_list,
                                  sub_outdir)


    else:

        sub_outdir = out_dir

        virus_infor_calculate(parameter_dic,
                              seq_pd_clean,
                              query_prefix,
                              lineage_mark_list,
                              sub_outdir)



    duration = datetime.today().now() - start  # 打印总运行时间

    end_memory = calculate_memory()

    used_memory = end_memory - start_memory

    print(">>> " + f"Occupied {used_memory}MB memory in total" + "\n")

    print(">>> " + "Take " + str(duration) + " seconds in total." + "\n")

    exit(0)