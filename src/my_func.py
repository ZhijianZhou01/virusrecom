# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Institution: Hunan University
Email: zjzhou@hnu.edu.cn
Time: 2022/4/17 21:32

"""

import os

import numpy as np
import pandas as pd


def get_all_path(open_dir_path):

    rootdir = open_dir_path

    path_list = []

    list = os.listdir(rootdir)

    for i in range(0, len(list)):
        com_path = os.path.join(rootdir, list[i])
        com_path = com_path.replace('\\', '/')
        #print(com_path)

        if os.path.isfile(com_path):
            path_list.append(com_path)
        if os.path.isdir(com_path):
            path_list.extend(get_all_path(com_path))

    return path_list



def make_dir(input_dir):
    """
    :param input_dir: directory which need to create
    :return:
    """
    if os.path.isdir(input_dir) == False:
        os.makedirs(input_dir)


def resolve_file_path(file_path):
    """
    Parse input file path
    :param file_path:
    :return:
    """
    input_data_dir = os.path.dirname(file_path)
    inputfile_name = file_path.split("/")[-1]
    file_farmat = inputfile_name.split(".")[-1]
    out_prefix = inputfile_name.replace(file_farmat,"").rstrip(".")

    return (input_data_dir,out_prefix)



def read_seq(file_path):
    """
    read sequence as data frame
    :param file_path:
    :return:
    """
    seq_tab = []
    max = 0
    with open(file_path, "r", encoding="utf-8") as gen_file_input:

        gen_file_list = gen_file_input.read().split(">")

        for n in range(1, int(len(gen_file_list))):
            seq = gen_file_list[n].strip()
            st = []
            seq_name = seq.split("\n")[0]
            # print(seq_name)

            seq_contain = seq.replace("\n", "").replace(seq_name, "").upper()

            st.append(seq_name)

            if len(seq_contain) >= max:
                max = len(seq_contain)

            for n in range(len(seq_contain)):
                st.append(seq_contain[n: n + 1])

            seq_tab.append(st)


    seq_sites_list = []
    seq_sites_list.append("SeqName")

    for x in range(max):
        seq_sites_list.append(str(x + 1))

    seq_pands = pd.DataFrame(seq_tab,columns=seq_sites_list).set_index("SeqName")


    return seq_pands



def calEnt(siteData):
    """
    Calculate the information content, no gaps
    :param siteData:Pass in a dataframe with only one column
    :return:
    """
    # print(siteData)
    n = siteData.shape[0]
    nt_count = siteData.value_counts()
    # print(nt_count)

    p = nt_count / n
    # print(p)

    ent = (-p * np.log2(p)).sum()

    IC = 2 - ent

    return IC


def calEnt_gap(siteData):
    """
    Calculate the information content, including gaps
    :param siteData:Pass in a dataframe with only one column
    :return:
    """
    # print(siteData)
    n = siteData.shape[0]
    nt_count = siteData.value_counts()
    # print(nt_count)

    p = nt_count / n
    # print(p)

    ent = (-p * np.log2(p)).sum()

    IC = np.log2(5) - ent

    return IC


