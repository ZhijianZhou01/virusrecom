# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Institution: Hunan University
Email: zjzhou@hnu.edu.cn
Copyright：Copyright (c) Zhou Zhi-Jian
Time: 2022/4/17 21:32

"""
import sys
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib

matplotlib.use("agg")  # avoid the  ModuleNotFoundError: No module named '_tkinter'

import matplotlib.pyplot as plt

from plt_corlor_list import plt_corlor


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

    file_path = file_path.replace("\\", "/")

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

            seq_tab.append(st)


        seq_sites_list = []
        seq_sites_list.append("SeqName")

        seq_length = len(seq_tab[0]) -1

        for x in range(seq_length):
            seq_sites_list.append(str(x + 1))

        seq_pands = pd.DataFrame(seq_tab,columns=seq_sites_list).set_index("SeqName")


        return seq_pands



def calEnt(siteData,mode):

    # print(siteData)

    if mode == "single":
        return 2

    else:
        n = siteData.shape[0]
        nt_count = siteData.value_counts()
        # print(nt_count)

        p = nt_count / n
        # print(p)

        ent = (-p * np.log2(p)).sum()

        IC = 2 - ent

        return IC



def calEnt_gap(siteData, mode):

    # print(siteData)

    if mode == "single":
        return np.log2(5)

    else:

        n = siteData.shape[0]
        nt_count = siteData.value_counts()
        # print(nt_count)

        p = nt_count / n
        # print(p)

        ent = (-p * np.log2(p)).sum()

        IC = np.log2(5) - ent

        return IC




def wic_calculation(seq_data,lineage_name_list,used_calEnt,
                    site_dir,query_seq_prefix):

    site_ic_table = (site_dir + "/"
                     + query_seq_prefix
                     + "_site_WIC_from_lineages.xlsx")

    site_ic_csv = (site_dir + "/"
                   + query_seq_prefix
                   + "_site_WIC.csv")


    query_seq = seq_data[
        seq_data.index.str.contains(
            query_seq_prefix) == True]

    query_seq_count = query_seq.shape[0]

    if query_seq_count == 0:

        print(">>> Error! The query lineage '" + query_seq_prefix + "'"
                          + " is not in sequence data." + "\n")

        sys.exit()

    original_site_list = [int(x) for x in
                          query_seq.columns.tolist()]

    sites_count = len(original_site_list)

    current_site_list = [i + 1 for i in range(sites_count)]

    columns_name_list = ["Original_Index", "Current_Index"] + lineage_name_list

    sites_wic_data = pd.DataFrame(columns=columns_name_list,
                                          index=current_site_list)

    sites_wic_data["Original_Index"] = original_site_list

    sites_wic_data["Current_Index"] = current_site_list


    for each_lineage in lineage_name_list:


        lineage_seq_df = seq_data[
            seq_data.index.str.contains(each_lineage) == True]

        ent_probability = []


        for (columnName, columnData) in lineage_seq_df.iteritems():

            lineage_site = list(columnData)

            if len(lineage_site) == 1:
                IC = used_calEnt(columnData,"single")

            else:

                IC = used_calEnt(columnData,"")

            if query_seq_count == 1:
                query_nt = query_seq[columnName][0]
                query_nt_ratio = 1

                try:

                    query_nt_lineage_ratio = lineage_site.count(query_nt) / \
                                             columnData.shape[0]

                    p_ent = query_nt_lineage_ratio * IC * query_nt_ratio

                    ent_probability.append(p_ent)


                except:

                    print(">>> Error! The lineage '" + each_lineage + "'"
                          + " is not in sequence data." + "\n")

                    sys.exit()


            else:

                query_seq_site_list = list(query_seq[columnName])

                query_seq_site_count = len(query_seq_site_list)

                maxpro_nt = max(query_seq_site_list,
                                key=query_seq_site_list.count)

                query_nt_ratio = query_seq_site_list.count(
                    maxpro_nt) / query_seq_site_count



                try:
                    query_nt_lineage_ratio = lineage_site.count(maxpro_nt) / \
                                             columnData.shape[0]

                    p_ent = query_nt_lineage_ratio * IC * query_nt_ratio

                    ent_probability.append(p_ent)


                except:

                    print(">>> Error! The lineage '" + each_lineage + "'"
                          + " is not in sequence data." + "\n")

                    sys.exit()


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

    step_probability_data = pd.DataFrame()

    slid_number = int(sites_count / step_size)

    window_center_index = []
    window_center_original = []

    slither_window_list = []


    for n in range(slid_number):

        start_row = step_size * n

        end_row = min(start_row + windows_size, sites_count - 1)

        center_index = int((start_row + end_row) / 2)

        slither_window_list.append([start_row, end_row, center_index])

        center_site = center_index + 1

        window_center_index.append(center_site)


        _original_site = int(site_map_dic[str(center_site)])

        window_center_original.append(_original_site)

        if end_row == sites_count - 1:
            break


    step_probability_data["Central_position(original)"] = window_center_original

    step_probability_data["Central_position(current)"] = window_center_index


    for each_lineage in lineage_name_list:

        each_lineage_setp_pro = []

        lineage_ic_df = sites_probability_data[
            each_lineage]


        for each_window in slither_window_list:
            start_row = each_window[0]

            end_row = each_window[1]

            step_df = list(lineage_ic_df[start_row:end_row])

            mean_ic_per_win = np.mean(step_df)

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

    recom_region_dic = {}

    for each_lineage in recombination_frag:
        lineage_frag_list = recombination_frag[each_lineage]

        if lineage_frag_list != []:

            frag_count = len(lineage_frag_list)
            # print(frag_count)

            cursor_site = 1
            cursor_center = lineage_frag_list[cursor_site - 1]

            detected_area = []

            while cursor_site <= frag_count:

                flage_label = False

                breakpoint_judgment = []


                for i in range(cursor_site - 1, frag_count):

                    region_left = max(0, int(
                        cursor_center - step_size / 2))
                    region_right = min(sites_count, int(
                        lineage_frag_list[i] + step_size / 2))

                    region_sites_count = region_right - region_left + 1

                    lineage_Ri_wic = sum(list(
                        sites_probability_data[each_lineage][
                        region_left:region_right]))

                    max_lineage_ic = 0

                    for lineage in lineage_name_list:
                        if lineage != each_lineage:
                            lineage_ic = sum(list(
                                sites_probability_data[lineage][
                                region_left:region_right]))
                            if lineage_ic >= max_lineage_ic:
                                max_lineage_ic = lineage_ic

                                # print(max_lineage_ic)


                    if (lineage_Ri_wic > max_lineage_ic
                            and lineage_Ri_wic / (
                                    region_sites_count * max_mic) >= recom_percentage):
                        breakpoint_judgment.append([i + 1, region_left,
                                                    region_right,
                                                    "True"])

                        flage_label = True

                if cursor_site == frag_count:
                    break

                if flage_label == False:

                    cursor_site = cursor_site + 1

                    cursor_center = lineage_frag_list[cursor_site - 1]

                else:

                    first_each_region = breakpoint_judgment[0]
                    first_ri = first_each_region[2] - first_each_region[1] + 1

                    if first_ri > max_recom_fragment:
                        cursor_site = cursor_site + 1
                        cursor_center = lineage_frag_list[cursor_site - 1]

                    else:

                        max_Ri = breakpoint_judgment[-1]
                        max_Ri_sites = max_Ri[2] - max_Ri[1] + 1

                        if max_Ri_sites <= max_recom_fragment:
                            detected_area.append([max_Ri[1], max_Ri[2]])

                            break

                        else:

                            for n in range(len(breakpoint_judgment)):
                                each_region = breakpoint_judgment[n]
                                ri = each_region[2] - each_region[
                                    1] + 1

                                if ri > max_recom_fragment:
                                    local_max_Ri = breakpoint_judgment[n - 1]
                                    detected_area.append([local_max_Ri[1],
                                                          local_max_Ri[2]])

                                    last_breakpoint = \
                                        breakpoint_judgment[n - 1][
                                            0]

                                    if last_breakpoint == cursor_site:
                                        cursor_site = cursor_site + 1
                                        cursor_center = lineage_frag_list[
                                            cursor_site - 1]

                                    else:
                                        cursor_site = last_breakpoint
                                        cursor_center = lineage_frag_list[
                                            cursor_site - 1]

                                    break

            recom_region_dic[each_lineage] = detected_area



    for i in list(recom_region_dic.keys()):
        if recom_region_dic[i] == []:
            del recom_region_dic[i]


    return recom_region_dic




def wic_plot(lineage_name,site_list,input_data,query_lineage_name,output_path):

    lineage_count = len(lineage_name)

    fig_high = int(max(site_list) * 3 / 10000) * 2

    fig, ax = plt.subplots(len(lineage_name), 1,
                           figsize=(
                           lineage_count, fig_high))

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

        x1_label = ax_n.get_xticklabels()
        [x1_label_temp.set_fontname("Arial") for x1_label_temp in
         x1_label]
        y1_label = ax_n.get_yticklabels()
        [y1_label_temp.set_fontname("Arial") for y1_label_temp in
         y1_label]


    plt.xlabel("Site in alignment", family="Arial")
    plt.savefig(output_path)
    plt.clf()



def mwic_plot(gap_used, lineage_name_list,
              window_center_original,
              step_probability_data,
              query_seq_prefix,
              output_file,
              y_start,
              legend_location):


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


    if legend_location.upper() == "R":

        num1 = 1.05
        num2 = 1
        num3 = 2
        num4 = 0
        plt.legend(bbox_to_anchor=(num1, num2), loc=num3, borderaxespad=num4)
        plt.subplots_adjust(bottom=0.10, right=0.7)


    else:
        plt.legend()

        plt.subplots_adjust(bottom=0.10)


    plt.margins(0)
    # plt.subplots_adjust(bottom=0.10)
    plt.xlabel("Site in alignment", family="Arial")
    plt.ylabel("Mean of weighted information content", family="Arial")
    plt.title("Query sequence: " + query_seq_prefix, family="Arial")


    x1_label = ax.get_xticklabels()
    [x1_label_temp.set_fontname("Arial") for x1_label_temp in
     x1_label]
    y1_label = ax.get_yticklabels()
    [y1_label_temp.set_fontname("Arial") for y1_label_temp in
     y1_label]

    # plt.show()
    plt.savefig(output_file)

    plt.clf()




def recombreak_plot(sites_probability_data,lineage_name_list,
                    sites_count,breakwins,
                    site_map_dic,
                    site_dir,query_seq_prefix):


    break_p_map = (site_dir + "/"
                   + query_seq_prefix
                   + "_ -lg(p-value)_for_potential_breakpoint.pdf")

    break_p_data = (site_dir + "/"
                    + query_seq_prefix
                    + "_ -lg(p-value)_for_potential_breakpoint.xlsx")

    breakpoint_data = pd.DataFrame()

    run_number = sites_count - breakwins + 1

    central_pos_list = []

    lineage_count = len(lineage_name_list)

    for lineage in lineage_name_list:

        central_pos_list = []
        negative_lg_p_list = []

        for i in range(run_number):
            start_site = i

            end_site = i + breakwins

            central_pos = int((start_site + end_site) / 2)

            left_region = list(sites_probability_data[lineage][
                               start_site: central_pos - 1])
            right_region = list(
                sites_probability_data[lineage][central_pos + 1: end_site])

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
                                                        fig_high2))

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

    plt.xlabel("Site in alignment", family="Arial")

    plt.savefig(break_p_map)

    plt.clf()





