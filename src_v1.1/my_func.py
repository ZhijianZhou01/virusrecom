# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Institution: Hunan University
Email: zjzhou@hnu.edu.cn
Copyright：Copyright (c) Zhou Zhi-Jian
Time: 2022/4/17 21:32

"""

import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib

matplotlib.use("agg")  # avoid the  ModuleNotFoundError: No module named '_tkinter'

import matplotlib.pyplot as plt

from plt_corlor_list import plt_corlor


def get_all_path(open_dir_path):
    """
    获取某一文件夹里所有文件的路径 (包括子文件夹)
    :param open_dir_path: 文件夹
    :return: 所有文件的路径列表
    """

    rootdir = open_dir_path

    path_list = []

    list = os.listdir(rootdir)  # 列出文件夹下所有的目录与文件

    for i in range(0, len(list)):
        com_path = os.path.join(rootdir, list[i])
        com_path = com_path.replace('\\', '/')  #把\变成/
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
    解析输入文件路径
    :param file_path:
    :return:
    """
    file_path = file_path.replace("\\", "/")  # 一定要写啊

    input_data_dir = os.path.dirname(file_path)

    inputfile_name = file_path.split("/")[-1]

    file_farmat = inputfile_name.split(".")[-1]

    out_prefix = inputfile_name.replace(file_farmat,"").rstrip(".")

    return (input_data_dir,out_prefix)




def mark_linegaes(lineage_file_dir, outdir):

    print(">>> VirusRecom is running..." + "\n")

    lineage_file_list = get_all_path(lineage_file_dir)

    out_file = outdir.replace("\\", "/") + "/merge_with_mark.fasta"

    with open(out_file,"w",encoding="utf-8") as out_put:

        for each_path in lineage_file_list:
            # print(each_path)

            each_path = each_path.replace("\\", "/")

            input_data_dir, out_prefix = resolve_file_path(each_path)

            with open(each_path, "r", encoding="utf-8") as lineage_file_input:
                out_put.write(
                    lineage_file_input.read().replace(">",
                                                ">" + out_prefix + "_")+ "\n")

        print(">>> Finished, the output file is "
              + "'" + out_file + "'"
              + "\n")



def export_seqname(seq_file,outdir):

    print(">>> VirusRecom is running..." + "\n")

    seq_file = seq_file.replace("\\", "/")

    outdir = outdir.replace("\\", "/")

    out_file = outdir + "/" + "sequence_name.txt"

    with open(out_file, "w", encoding="utf-8") as output_file:

        with open(seq_file, "r", encoding="utf-8") as input_file:

            for line in input_file:
                line = line.strip()
                if line.startswith(">"):
                    output_file.write(line.replace(">", "") + "\n")

        print(">>> Finished, the output file is "
              + "'" + out_file + "'"
              + "\n")




def read_seq(file_path):
    """
    读取序列为数据框
    :param file_path: 输入文件路径
    :return:
    """
    seq_tab = []

    # max_length = 0
    with open(file_path, "r", encoding="utf-8") as gen_file_input:

        gen_file_list = gen_file_input.read().split(">")

        for n in range(1, int(len(gen_file_list))):
            seq = gen_file_list[n].strip()
            st = []
            seq_name = seq.split("\n")[0]
            # print(seq_name)

            seq_contain = seq.replace("\n", "").replace(seq_name, "").upper()

            st.append(seq_name)

            # if len(seq_contain) >= max_length:
            #     max_length = len(seq_contain)

            # for k in range(len(seq_contain)):
            #     st.append(seq_contain[k: k + 1])

            st += list(seq_contain.strip())

            seq_tab.append(st)  # 每条序列是一个st列表，seq_tab的每行储存不同的st（序列）


        seq_sites_list = []
        seq_sites_list.append("SeqName")

        seq_length = len(seq_tab[0]) -1    # 因为每条序列第一列为序列名称

        for x in range(seq_length):
            seq_sites_list.append(str(x + 1))

        seq_pands = pd.DataFrame(seq_tab,columns=seq_sites_list).set_index("SeqName")


        return seq_pands



def calEnt(siteData,mode):
    """
    计算信息量
    :param siteData:传入只有一列的数据框
    :return:
    """
    # print(siteData)

    if mode == "single":
        return 2

    else:
        n = siteData.shape[0]               # 总序列数（该位点总碱基数量）
        nt_count = siteData.value_counts()  # 每种碱基计数
        # print(nt_count)

        p = nt_count / n                    # 每种碱基占比
        # print(p)

        ent = (-p * np.log2(p)).sum()       # 计算香农熵

        IC = 2 - ent                        # 计算信息量

        return IC



def calEnt_gap(siteData, mode):
    """
    计算信息量
    :param siteData:传入只有一列的数据框
    :return:
    """

    # print(siteData)

    if mode == "single":
        return np.log2(5)

    else:

        n = siteData.shape[0]                # 总序列数（该位点总碱基数量）
        nt_count = siteData.value_counts()   # 每种碱基计数
        # print(nt_count)

        p = nt_count / n                    # 每种碱基占比
        # print(p)

        ent = (-p * np.log2(p)).sum()       # 计算香农熵

        IC = np.log2(5) - ent               # 计算信息量

        return IC




def wic_calculation(seq_data,lineage_name_list,used_calEnt,
                    site_dir,query_seq_prefix):
    """

    :param seq_data: 序列矩阵数据
    :param site_dir: 储存的文件路径
    :param query_seq_prefix: 查询谱系名称
    :return: 一个文件路径，生成的.csv文件
    """

    site_ic_table = (site_dir + "/"
                     + query_seq_prefix
                     + "_site_WIC_from_lineages.xlsx")

    site_ic_csv = (site_dir + "/"
                   + query_seq_prefix
                   + "_site_WIC.csv")


    # 先把查询谱系的数据抽取出来，为数据框
    query_seq = seq_data[
        seq_data.index.str.contains(
            query_seq_prefix) == True]    # 用于查询的潜在重组体的谱系

    query_seq_count = query_seq.shape[0]  # “查询谱系的”序列条数

    if query_seq_count == 0:

        print(">>> Error! The query lineage '" + query_seq_prefix + "'"
                          + " is not in sequence data." + "\n")

        exit(0)

    original_site_list = [int(x) for x in
                          query_seq.columns.tolist()]  # 用来记录碱基位点在原来基因组中的index（"Original_Index"）

    sites_count = len(original_site_list)  # 进入后续分析的位点的数量

    current_site_list = [i + 1 for i in range(sites_count)]  # 用来记录经过处理后，碱基位点实际的index（"Current_Index"）

    columns_name_list = ["Original_Index", "Current_Index"] + lineage_name_list

    sites_wic_data = pd.DataFrame(columns=columns_name_list,
                                          index=current_site_list)

    sites_wic_data["Original_Index"] = original_site_list

    sites_wic_data["Current_Index"] = current_site_list


    """

    本质上Current_Index不需要额外记录，因为循环数据框时本身就能获取位点顺序，
    但为了代码清晰，变量能及时释放，数据便于追寻，在之后生成的数据框增加这一列

    """

    # 采用双层循环，先遍历小的（外层循环每个谱系），然后再第二层遍历大的
    for each_lineage in lineage_name_list:

        # 抽取出这个谱系的序列数据
        lineage_seq_df = seq_data[
            seq_data.index.str.contains(each_lineage) == True]

        ent_probability = []  # 储存这个谱系的值的列表

        #  内层循环每个位点

        for (columnName, columnData) in lineage_seq_df.iteritems():
            # print(columnName)   # 列名
            # print(columnData)   # 该列数据

            lineage_site = list(columnData)  # 该列碱基，转成list。

            # 该谱系在该位点的信息量
            if len(lineage_site) == 1:  # 如果该位点只有一条序列
                IC = used_calEnt(columnData,"single")

            else:

                IC = used_calEnt(columnData,"")

            if query_seq_count == 1:  # 如果只有一条“查询序列”
                query_nt = query_seq[columnName][0]  # “查询序列”在该位点的碱基
                query_nt_ratio = 1  # “查询序列”的该位点碱基在“查询序列”中的占比

                try:
                    # “查询序列”该位点碱基在谱系中的占比
                    query_nt_lineage_ratio = lineage_site.count(query_nt) / \
                                             columnData.shape[0]

                    p_ent = query_nt_lineage_ratio * IC * query_nt_ratio  # 计算该位点加权后的信息量

                    ent_probability.append(p_ent)


                except:

                    print(">>> Error! The lineage '" + each_lineage + "'"
                          + " is not in sequence data." + "\n")

                    exit(0)



            else:  # 如果是多条查询序列（即查询序列也是谱系），计算该位点含量最多的碱基，以及比例

                query_seq_site_list = list(query_seq[columnName])  # 该位点碱基对应的列表

                query_seq_site_count = len(query_seq_site_list)  # 该位点碱基总数量

                maxpro_nt = max(query_seq_site_list,
                                key=query_seq_site_list.count)  # “查询序列”中出现次数最多的碱基，标记为maxpro_nt

                query_nt_ratio = query_seq_site_list.count(
                    maxpro_nt) / query_seq_site_count  # 该次数最多的碱基在“查询序列”中的占比



                try:
                    # 该次数最多的碱基在谱系中的占比
                    query_nt_lineage_ratio = lineage_site.count(maxpro_nt) / \
                                             columnData.shape[0]

                    p_ent = query_nt_lineage_ratio * IC * query_nt_ratio  # 计算该位点加权后的信息量

                    ent_probability.append(p_ent)


                except:

                    print(">>> Error! The lineage '" + each_lineage + "'"
                          + " is not in sequence data." + "\n")

                    exit(0)


        sites_wic_data[each_lineage] = ent_probability

        print("    " +
            "The calculation of " + each_lineage + "'s recombination contribution to "
            + query_seq_prefix + " has been completed!" + "\n")

    sites_wic_data.to_excel(excel_writer=site_ic_table,
                                    index=False,
                                    encoding="utf-8")

    sites_wic_data.to_csv(site_ic_csv, index=False, sep=",")

    return site_ic_csv




def mwic_calculation(sites_probability_data,lineage_name_list,
                     sites_count,site_map_dic,
                     windows_size,step_size,
                     output_path):

    step_probability_data = pd.DataFrame()  # 数据框，储存所有谱系在每个窗口的平均加权信息量

    slid_number = int(sites_count / step_size)  # 先计算出大致最多循环（滑动）多少次

    window_center_index = []     # 一个列表，储存将窗口中心位点index
    window_center_original = []  # 一个列表，储存将窗口中心映射到基因组原始位点

    slither_window_list = []  # 储存窗口的起始位点，结束位点，中心点  （以0开始）

    # 3.1. 计算出滑动时每个窗口的起始位点

    for n in range(slid_number):

        start_row = step_size * n

        end_row = min(start_row + windows_size, sites_count - 1)

        center_index = int((start_row + end_row) / 2)  # 标记该窗口的中心位点

        slither_window_list.append([start_row, end_row, center_index])

        center_site = center_index + 1  # 因为循环是从0开始的，所以要加上1

        window_center_index.append(center_site)

        # 定位该窗口中心位点在原基因组中的位点（考虑到可能使用的是处理后的多态性位点）,需要和sites_probability_data进行关联

        # 获取label_site行，"Original_Index"列的信息。

        _original_site = int(site_map_dic[str(center_site)])

        window_center_original.append(_original_site)

        if end_row == sites_count - 1:  # 不再滑动，这里减1是因为数据框是从0开始索引的
            break


    step_probability_data["Central_position(original)"] = window_center_original

    step_probability_data["Central_position(current)"] = window_center_index

    # 3.2.开始滑动计算平均WIC

    for each_lineage in lineage_name_list:  # 将谱系放在外循环，遍历谱系

        each_lineage_setp_pro = []  # 列表，储存这个谱系在每个窗口的平均加权信息量值

        lineage_ic_df = sites_probability_data[
            each_lineage]  # 抽取出这个谱系的加权位点信息量数据

        # 滑动抽取数据框

        for each_window in slither_window_list:
            start_row = each_window[0]

            end_row = each_window[1]

            step_df = list(lineage_ic_df[start_row:end_row])  # 数据框左闭右开

            mean_ic_per_win = np.mean(step_df)  # 计算该滑动窗口的平均加权信息量

            each_lineage_setp_pro.append(mean_ic_per_win)

        step_probability_data[each_lineage] = each_lineage_setp_pro

    step_probability_data.to_excel(
        excel_writer=output_path,
        index=False,
        encoding="utf-8")

    return step_probability_data,slither_window_list




def recomsplicing(sites_probability_data,sites_count,lineage_name_list,
                  recombination_frag,
                  step_size,max_mic,max_recom_fragment,recom_percentage):

    recom_region_dic = {}  # 一个字典，用于储存拼接后的重组区域

    for each_lineage in recombination_frag:
        lineage_frag_list = recombination_frag[each_lineage]

        if lineage_frag_list != []:

            frag_count = len(lineage_frag_list)  # 片段总个数
            # print(frag_count)

            cursor_site = 1  # 先将lineage_frag_list里第一个片段作为 起始片段（游标）
            cursor_center = lineage_frag_list[cursor_site - 1]  # 起始片段（游标）的中心

            detected_area = []  # 该谱系检测到的重组区域（拼接连续和非连续片段结果）

            while cursor_site <= frag_count:

                flage_label = False

                breakpoint_judgment = []

                # 判断游标起始片段后面的 每个断点（片段）到起始片段（游标）构成的重组区域Ri是否成立（和主要亲本的mWIC比较）
                for i in range(cursor_site - 1, frag_count):

                    # 第i个片段 到 起始片段（游标）的构成的 重组区域Ri
                    region_left = max(0, int(
                        cursor_center - step_size / 2))  # 重组区域Ri 左侧位点
                    region_right = min(sites_count, int(
                        lineage_frag_list[i] + step_size / 2))  # 重组区域Ri 右侧位点

                    region_sites_count = region_right - region_left + 1  # 该区间位点数量（这里要加上1）

                    # 计算该谱系在重组区域Ri内的的mWIC，注意，这里使用数据不是滑动窗口里的,是位点数据sites_probability_data里的

                    lineage_Ri_wic = sum(list(
                        sites_probability_data[each_lineage][
                        region_left:region_right]))

                    max_lineage_ic = 0  # 其他谱系最大值

                    for lineage in lineage_name_list:
                        if lineage != each_lineage:
                            lineage_ic = sum(list(
                                sites_probability_data[lineage][
                                region_left:region_right]))
                            if lineage_ic >= max_lineage_ic:
                                max_lineage_ic = lineage_ic

                                # print(max_lineage_ic)

                    # 如果该谱系的mWIC在该区间最大（大于其他谱系最大值），且大于指定的阈值比例（p）
                    if (lineage_Ri_wic > max_lineage_ic
                            and lineage_Ri_wic / (
                                    region_sites_count * max_mic) >= recom_percentage):
                        breakpoint_judgment.append([i + 1, region_left,
                                                    region_right,
                                                    "True"])  # 将该重组区域Ri标记为True全记录下来

                        flage_label = True

                if cursor_site == frag_count:
                    break

                if flage_label == False:

                    cursor_site = cursor_site + 1  # 起始片段（游标）后移一个

                    cursor_center = lineage_frag_list[cursor_site - 1]

                else:

                    # 如果最小的一个都超过了最大重组区间阈值范围
                    first_each_region = breakpoint_judgment[0]
                    first_ri = first_each_region[2] - first_each_region[1] + 1

                    if first_ri > max_recom_fragment:
                        cursor_site = cursor_site + 1  # 起始片段（游标）后移一个
                        cursor_center = lineage_frag_list[cursor_site - 1]

                    else:  # 找最近的一个

                        # 寻找在ri < mr 的条件下寻找最大

                        max_Ri = breakpoint_judgment[-1]  # breakpoint_judgment里最大的一个
                        max_Ri_sites = max_Ri[2] - max_Ri[1] + 1  # 该区间位点数量（这里要加上1）

                        if max_Ri_sites <= max_recom_fragment:  # 如果最大的一个都没有超过最大的片段阈值（mr）
                            detected_area.append([max_Ri[1], max_Ri[2]])

                            break

                        else:

                            for n in range(len(breakpoint_judgment)):
                                each_region = breakpoint_judgment[n]
                                ri = each_region[2] - each_region[
                                    1] + 1  # 该区间位点数量（这里要加上1）

                                if ri > max_recom_fragment:
                                    local_max_Ri = breakpoint_judgment[n - 1]
                                    detected_area.append([local_max_Ri[1],
                                                          local_max_Ri[2]])

                                    last_breakpoint = \
                                        breakpoint_judgment[n - 1][
                                            0]  # 注意，这里一定要写成breakpoint_judgment[n-1]，而不是breakpoint_judgment[n]

                                    if last_breakpoint == cursor_site:
                                        cursor_site = cursor_site + 1  # 起始片段（游标）后移一个
                                        cursor_center = lineage_frag_list[
                                            cursor_site - 1]

                                    else:
                                        cursor_site = last_breakpoint
                                        cursor_center = lineage_frag_list[
                                            cursor_site - 1]

                                    break

            recom_region_dic[each_lineage] = detected_area


    # 删除recom_region_dic里为空的谱系

    for i in list(recom_region_dic.keys()):
        if recom_region_dic[i] == []:
            del recom_region_dic[i]


    return recom_region_dic




def wic_plot(lineage_name,site_list,input_data,query_lineage_name,output_path):
    """
    wic 绘图

    :param lineage_name:  谱系名称列表
    :param site_list:     用于绘制X轴的位点index列表
    :param input_data:   输出的数据框数据
    :param output_path:   输出文件路径
    :return:
    """

    lineage_count = len(lineage_name)   # 谱系数量

    fig_high = int(max(site_list) * 3 / 10000) * 2

    fig, ax = plt.subplots(len(lineage_name), 1,
                           figsize=(
                           lineage_count, fig_high))  # 创建画布（figure画布）和多个ax(区域)

    fig.suptitle("Query seq: " + query_lineage_name, family="Arial")
    fig.tight_layout()

    cm1 = plt.cm.get_cmap("Reds")

    plt.subplots_adjust(top=0.95)

    for n in range(lineage_count):
        each_lineage = lineage_name[n]
        ax_n = ax[n]
        y = list(input_data[each_lineage])
        xx = ax_n.scatter(site_list, y, c=y,
                          label=each_lineage, s=5,
                          cmap=cm1)

        fig.colorbar(xx, ax=ax_n)

        ax_n.set_ylabel("WIC", family="Arial")

        ax_n.legend(loc="best")

        # 坐标轴刻度字体设置
        x1_label = ax_n.get_xticklabels()
        [x1_label_temp.set_fontname("Arial") for x1_label_temp in
         x1_label]
        y1_label = ax_n.get_yticklabels()
        [y1_label_temp.set_fontname("Arial") for y1_label_temp in
         y1_label]


    plt.xlabel("Site in alignment", family="Arial")  # X轴标签
    plt.savefig(output_path)
    plt.clf()   # 清空画布。后面还需要画图，这里要清空画布。



def mwic_plot(gap_used, lineage_name_list,
              window_center_original,
              step_probability_data,
              query_seq_prefix,
              output_file,
              y_start,
              legend_location):
    """
    绘制mWIC图

    :param gap_used: 是否将gap纳入了分析
    :param y_start: Y轴起始位点
    :param lineage_name_list: 参考谱系名称列表
    :param window_center_original: 滑动窗口的中点（index已定位到基因组上）
    :param step_probability_data: mwic数据框数据
    :param query_seq_prefix: 查询谱系名称
    :param output_file: 绘图后的输出文件路径
    :return:
    """

    fig, ax = plt.subplots()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    if gap_used.upper() == "N":

        plt.ylim((y_start, 2))

    else:
        plt.ylim((y_start, 2.5))

    for n in range(len(lineage_name_list)):
        each_lineage = lineage_name_list[n]

        try:
            plt.plot(window_center_original,
                     list(step_probability_data[each_lineage]),
                     label=each_lineage, color=plt_corlor[n], linewidth=1)
        except:
            plt.plot(window_center_original,
                     list(step_probability_data[each_lineage]),
                     label=each_lineage, color="black", linewidth=1)
        finally:
            pass

    # 图例
    if legend_location.upper() == "R":

        num1 = 1.05  # 水平位置
        num2 = 1  # 垂直位置
        num3 = 2  # 位于右上
        num4 = 0
        plt.legend(bbox_to_anchor=(num1, num2), loc=num3, borderaxespad=num4)
        plt.subplots_adjust(bottom=0.10, right=0.7)


    else:
        plt.legend()

        plt.subplots_adjust(bottom=0.10)


    plt.margins(0)
    # plt.subplots_adjust(bottom=0.10)
    plt.xlabel("Site in alignment", family="Arial")  # X轴标签
    plt.ylabel("Mean of weighted information content", family="Arial")  # Y轴标签
    plt.title("Query sequence: " + query_seq_prefix, family="Arial")  # 标题




    # 坐标轴刻度字体设置
    x1_label = ax.get_xticklabels()
    [x1_label_temp.set_fontname("Arial") for x1_label_temp in
     x1_label]
    y1_label = ax.get_yticklabels()
    [y1_label_temp.set_fontname("Arial") for y1_label_temp in
     y1_label]

    # plt.show()
    plt.savefig(output_file)

    plt.clf()  # 清空画布




def recombreak_plot(sites_probability_data,lineage_name_list,
                    sites_count,breakwins,
                    site_map_dic,
                    site_dir,query_seq_prefix):

    # 精确的重组断点搜寻
    # 外层循环每个谱系

    # 搜寻重组断点使用的窗口大小

    # 使用site_probability_data这个数据

    # 步长值默认是 1

    break_p_map = (site_dir + "/"
                   + query_seq_prefix
                   + "_ -lg(p-value)_for_potential_breakpoint.pdf")

    break_p_data = (site_dir + "/"
                    + query_seq_prefix
                    + "_ -lg(p-value)_for_potential_breakpoint.xlsx")

    breakpoint_data = pd.DataFrame()

    # 迭代次数 （滑动）多少次，等于序列总的位点数减去窗口的长度

    run_number = sites_count - breakwins + 1
    # 这里加1，是因为breakwins为包含的位点数量，比如第1位点到第5位点，实际间隔只有4，但是有5个位点

    central_pos_list = []

    lineage_count = len(lineage_name_list)  # 参考谱系的数量

    for lineage in lineage_name_list:

        central_pos_list = []
        negative_lg_p_list = []

        for i in range(run_number):
            start_site = i

            end_site = i + breakwins  # breakwins为搜寻重组断点的窗口

            central_pos = int((start_site + end_site) / 2)

            left_region = list(sites_probability_data[lineage][
                               start_site: central_pos - 1])
            right_region = list(
                sites_probability_data[lineage][central_pos + 1: end_site])

            # 定位原始位点
            original_site = int(site_map_dic[str(central_pos + 1)])

            central_pos_list.append(original_site)

            try:

                zihe_test = stats.mannwhitneyu(left_region, right_region,
                                               alternative="two-sided")
                p_value = zihe_test[1]
                negative_lg_p = - np.log10(p_value)
                negative_lg_p_list.append(negative_lg_p)

            except:
                negative_lg_p = 0
                negative_lg_p_list.append(negative_lg_p)

            finally:
                pass

        breakpoint_data["Original_Index"] = central_pos_list
        breakpoint_data[lineage] = negative_lg_p_list

    # print(breakpoint_data)

    breakpoint_data.to_excel(
        excel_writer=break_p_data,
        index=False,
        encoding="utf-8")

    fig_high2 = int(max(central_pos_list) * 3 / 10000) * 2

    figs, axs = plt.subplots(lineage_count, 1, figsize=(lineage_count,
                                                        fig_high2))  # 创建画布（figure画布）和多个ax(区域)

    figs.suptitle("Query seq: " + query_seq_prefix, family="Arial")
    figs.tight_layout()

    plt.subplots_adjust(top=0.95)

    for n in range(lineage_count):
        each_lineage = lineage_name_list[n]
        ax_n = axs[n]
        y = list(breakpoint_data[each_lineage])

        try:

            ax_n.plot(central_pos_list, y,
                      label=query_seq_prefix + " : " + each_lineage,
                      color=plt_corlor[n], )

        except:
            ax_n.plot(central_pos_list, y,
                      label=query_seq_prefix + " : " + each_lineage,
                      color="black", )

        finally:
            pass

        ax_n.set_ylabel("-lg(P)", family="Arial")

        ax_n.legend(loc="best")

        x1_label = ax_n.get_xticklabels()
        [x1_label_temp.set_fontname("Arial") for x1_label_temp in
         x1_label]
        y1_label = ax_n.get_yticklabels()
        [y1_label_temp.set_fontname("Arial") for y1_label_temp in
         y1_label]

    plt.xlabel("Site in alignment", family="Arial")  # X轴标签

    plt.savefig(break_p_map)

    plt.clf()  # 清空画布





