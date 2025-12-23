# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Institution: Hunan University
Email: zjzhou@hnu.edu.cn
Copyright：Copyright (c) Zhou Zhi-Jian
Time: 2022/5/20 20:18

"""

import subprocess


def SeqAlign(align_tool,thread_num,seq_file,out_file):

    print(
        "VirusRecom is running " + align_tool
        + " for sequence alignment..." + "\n")

    command_list = []

    command_list.append(align_tool)

    if align_tool.upper() == "MAFFT":

        command_list.append("--inputorder")
        command_list.append("--auto")
        command_list.append("--thread " + str(thread_num) + " ")
        command_list.append('"' + seq_file + '"' + " > "
                          + '"' + out_file + '"')


    elif align_tool.upper() == "MUSCLE":

        command_list.append("-maxiters 16")

        command_list.append(r" -in "
               + '"' + seq_file + '"'
               + " -out "
               + '"' + out_file + '"')


    elif align_tool.upper() == "CLUSTALO":

        command_list.append("--outfmt a2m")  # fasta format

        command_list.append("--threads " + str(thread_num))

        command_list.append("--force")

        command_list.append("--verbose")

        command_list.append(r" --in "
               + '"' + seq_file + '"'
               + " --out "
               + '"' + out_file + '"')


    align_command = " ".join(command_list)


    print(align_command)


    process = subprocess.Popen(align_command,
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

    return




