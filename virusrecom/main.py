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

app_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

app_dir = app_dir.replace("\\", "/")

sys.path.append(app_dir)

import psutil

import numpy as np
import argparse

from datetime import datetime


from my_func import (handle_input_file,make_dir,
                     calEnt, calEnt_gap,
                     handle_exceptions,
                     others_analysis)


from corecp import virus_infor_calculate


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
        default="")

    parser.add_argument(
        "-at", dest="align_tool",
        help="Program used for multiple sequence alignments (MSA). Supporting mafft, muscle and clustalo, such as ‘-at mafft’.",
        default="")

    parser.add_argument(
        "-iwic", dest="input_wic",
        help="Using the already obtained WIC values of reference lineages directly by a *.csv input-file.",
        default="")

    parser.add_argument(
        "-q", dest="query",
        help="Name of query lineage (usually potential recombinant), such as ‘-q xxxx’. Besides, ‘-q auto’ can scan all lineages as potential recombinant in turn.",
        default="")

    parser.add_argument(
        "-l", dest="lineages",
        help="Path of a text-file containing multiple lineage marks.",
        default="")

    parser.add_argument("-g", dest="gap",
                        help="Reserve sites containing gaps(-) in analyses? ‘-g y’ means to reserve, and ‘-g n’ means to delete.",
                        type=str,
                        default="")

    parser.add_argument("-m", dest="method",
                        help="Method for site scanning. ‘-m p’ uses polymorphic sites only, ‘-m a’ uses all the sites.",
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
        help="Number of threads (or cores) for calculations, default: 4.",
        type=int,
        default=4)

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
        "-e", dest="engrave",
        help="Write the file name to sequence names in batches. By specifying a directory containing one or multiple sequence files (*.fasta).",
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


    parser.add_argument(
        "--block",
        dest="block_size",
        help="Specifies the maximum number of sites per sub-block, different sub-blocks in sequence file will be sequentially loaded to calculate WIC. Default: 40000.",
        type=int,
        default=40000)


    parser.add_argument(
        "--no_wic_fig",
        dest="no_wic_figure",
        action="store_true",
        help="Do not draw the image of WICs.")


    parser.add_argument(
        "--no_mwic_fig",
        dest="no_mwic_figure",
        action="store_true",
        help="Do not draw the image of mWICs.")


    myargs = parser.parse_args(sys.argv[1:])

    return myargs




def starts():

    print("\n" + "-------------------------------------------------")

    print("  Name: Virus Recombination (VirusRecom)")

    print(
        "  Description: Detecting recombination of viral lineages (or subtypes) using information theory.")

    print("  Version: 1.3.2 (2024-07-23)")

    print("  Author: Zhi-Jian Zhou")

    print("  Citation: Brief Bioinform. 2023 Jan 19;24(1)")

    print("-------------------------------------------------" + "\n")


    myargs = parameter()

    start = datetime.today().now()


    #  parameter

    parameter_dic = {}

    parameter_dic["aligned_seq"] = myargs.alignment  # path of the aligned sequence file

    parameter_dic["unaligned_seq"] = myargs.unalignment   # path of the unaligned sequence file

    parameter_dic["align_tool"] = myargs.align_tool  #  alignment tool used for no-aligned sequence file

    parameter_dic["input_wic"] = myargs.input_wic.replace("\\", "/")     # input_wic, if you have

    parameter_dic["query_lineage_name"] = myargs.query         #  the name of query lineage

    parameter_dic["lineage_file"] = myargs.lineages          # other lineages

    parameter_dic["gaps_use"] = myargs.gap.strip().upper()    #  strategy of handling gap

    parameter_dic["method"] = myargs.method.strip().upper()   #  Scan method

    parameter_dic["windows_size"] = myargs.window    #  window size

    parameter_dic["step_size"] = myargs.step     #  step size

    parameter_dic["max_recom_fragment"] = myargs.max_region  #  maximum recombination region

    parameter_dic["recom_percentage"] = myargs.percentage  #  the specific cutoff threshold of proportion

    parameter_dic["cumulative_major"] = myargs.cumulative

    parameter_dic["breakpoints"] =  myargs.breakpoint.strip()  #  whether to search for recombination breakpoints

    parameter_dic["breakwins"] = myargs.breakwin  #  Window size used to search for recombination breakpoints

    parameter_dic["thread_num"] = int(myargs.thread)  #  thread

    parameter_dic["y_start"] = myargs.y_start  #  Y-axis starting point when plotting

    parameter_dic["legend"] = myargs.legend

    parameter_dic["only_wic"] = myargs.only_wic.strip().upper()

    parameter_dic["engrave_name"] = myargs.engrave

    parameter_dic["export_name"] = myargs.export_name

    parameter_dic["out_dir"] = myargs.outdir.replace("\\", "/")  # output directory

    parameter_dic["no_wic_figure"] = myargs.no_wic_figure

    parameter_dic["no_mwic_figure"] = myargs.no_mwic_figure

    parameter_dic["block_size"] = int(myargs.block_size)


    print(myargs)

    print("\n")


    make_dir(parameter_dic["out_dir"])


    with open(parameter_dic["out_dir"] + "/" + "input-parameter.txt",
              "w", encoding="utf-8") as input_para:

        input_para.write(
            "the used parameter of virusrecom in command-line interface."
            + "\n" * 2)

        for key in list(parameter_dic.keys()):
            input_para.write(key + "\t" + str(parameter_dic[key]) + "\n")

    # update dic

    if parameter_dic["gaps_use"].upper() == "N":

        parameter_dic["calEnt_use"] = calEnt
        parameter_dic["max_mic"] = 2

    else:
        parameter_dic["calEnt_use"] = calEnt_gap

        parameter_dic["max_mic"] = np.log2(5)


    others = others_analysis(parameter_dic)

    if others:
        sys.exit()


    exceptions_label = handle_exceptions(parameter_dic)

    if exceptions_label:
        sys.exit()


    query_prefix = parameter_dic["query_lineage_name"]


    if parameter_dic["input_wic"] != "":

        virus_infor_calculate("",
                              query_prefix,
                              [],
                              parameter_dic,
                              parameter_dic["out_dir"])


    else:

        seq_alignment_file = handle_input_file(parameter_dic)

        lineage_mark_list = []

        with open(parameter_dic["lineage_file"]) as lineage_file:
            for line in lineage_file:
                line = line.strip()
                if line != "":
                    lineage_mark_list.append(line)


        if query_prefix.upper() != "AUTO":

            out_dir = parameter_dic["out_dir"]

            virus_infor_calculate(seq_alignment_file,
                                  query_prefix,
                                  lineage_mark_list,
                                  parameter_dic,
                                  out_dir)


        elif query_prefix.upper() == "AUTO":

            for each_lineage in lineage_mark_list:
                query_lineage_run = each_lineage

                new_lineage_name_list = []

                for lineage in lineage_mark_list:
                    if lineage != query_lineage_run:
                        new_lineage_name_list.append(lineage)

                sub_outdir = parameter_dic["out_dir"] + "/Results/" + query_lineage_run

                make_dir(sub_outdir)

                virus_infor_calculate(seq_alignment_file,
                                      query_lineage_run,
                                      new_lineage_name_list,
                                      parameter_dic,
                                      sub_outdir)



    duration = datetime.today().now() - start

    print(">>> " + "Take " + str(duration) + " seconds in total." + "\n")


    sys.exit()



if __name__ == "__main__":
    starts()
