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
from multiprocessing import Pool
from sequence_align import SeqAlign

matplotlib.use("agg")  # avoid the  ModuleNotFoundError: No module named '_tkinter'

import matplotlib.pyplot as plt

from plt_corlor_list import plt_corlor



def handle_exceptions(parameter_dic):

    if parameter_dic["out_dir"].strip() == "":
        print(
            ">>> Error！Please specify an output directory by '-o' options." + "\n")

        return True

    elif (parameter_dic["unaligned_seq"] == ""
          and parameter_dic["aligned_seq"] == ""
          and parameter_dic["input_wic"] == ""):

        print(">>> Error！There is no input file." + "\n")

        return True


    elif parameter_dic["query_lineage_name"] == "":
        print(">>> Error！Please specify the query lineage by '-q' options." + "\n")

        return True


    elif parameter_dic["gaps_use"] == "":
        print(
            ">>> Error！Please specify the value of '-g' options." + "\n")

        return True




def handle_input_file(parameter_dic):


    if parameter_dic["aligned_seq"] != "": 

        seq_alignment_file = parameter_dic["aligned_seq"].replace("\\", "/")

        return seq_alignment_file


    elif parameter_dic["unaligned_seq"] != "":

        if parameter_dic["align_tool"] != "":

            align_record = parameter_dic["out_dir"] + "/" + "align_record"

            make_dir(align_record)

            input_seq_dir, input_prefix = resolve_file_path(parameter_dic["unaligned_seq"])

            seq_alignment_file = (align_record + "/" + input_prefix
                                  + "_" + parameter_dic["align_tool"]
                                  + ".fasta")

            SeqAlign(parameter_dic["align_tool"],
                     parameter_dic["thread_num"],
                     parameter_dic["unaligned_seq"],
                     seq_alignment_file)

            file_size = os.path.getsize(seq_alignment_file)

            if os.path.exists(seq_alignment_file) and file_size > 0:

                file_size = os.path.getsize(seq_alignment_file)

                if file_size > 0:

                    print(">>> Sequence alignment has been completed!" + "\n")

                    return seq_alignment_file


            else:

                print(">>> Error！Multiple sequence alignments did not succeed." + "\n")

                sys.exit()


        else:
            print(">>> Error！Please specify a tool used for sequence alignment, "
                  "such as '-at mafft', or '-at muscle', or '-at clustalo'." + "\n")

            sys.exit()



    else:
        sys.exit()



def get_all_path(open_dir_path):

    if not os.path.exists(open_dir_path):
        raise FileNotFoundError(
            f"The directory {open_dir_path} does not exist.")

    path_list = [
        os.path.join(open_dir_path, f).replace("\\", "/")
        for f in os.listdir(open_dir_path)
        if os.path.isfile(os.path.join(open_dir_path, f))
    ]

    return path_list



def make_dir(input_dir):
    """
    :param input_dir: directory which need to create
    :return:
    """
    if not os.path.isdir(input_dir):
        os.makedirs(input_dir)


def resolve_file_path(file_path):

    file_path = file_path.replace("\\", "/")

    input_data_dir = os.path.dirname(file_path)

    inputfile_name = os.path.basename(file_path)

    file_prefix, _ = os.path.splitext(inputfile_name)

    return input_data_dir, file_prefix



def mark_lineages(lineage_file_dir, outdir):

    print(">>> VirusRecom is running..." + "\n")

    lineage_file_list = get_all_path(lineage_file_dir)

    out_file = outdir.replace("\\", "/") + "/merge_with_mark.fasta"

    with open(out_file,"w",encoding="utf-8") as out_put:

        for each_path in lineage_file_list:

            each_path = each_path.replace("\\", "/")

            input_data_dir, out_prefix = resolve_file_path(each_path)

            with open(each_path, "r", encoding="utf-8") as lineage_file_input:
                out_put.write(
                    lineage_file_input.read().replace(">",
                                                ">" + out_prefix + "_")+ "\n")

        print(">>> Finished, the output file is "
              + "'" + out_file + "'"
              + "\n")



def export_seqname(seq_file, outdir):

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


def others_analysis(parameter_dic):

    if parameter_dic["engrave_name"] != "":

        mark_lineages(parameter_dic["engrave_name"],
                      parameter_dic["out_dir"])

        return True

    if parameter_dic["export_name"] != "":
        export_seqname(parameter_dic["export_name"],
                       parameter_dic["out_dir"])

        return True

    else:
        return False


def read_seq(input_file_path):

    with open(input_file_path, "r",encoding="utf-8") as gen_file:

        seq_tab = []

        seq_dic = {}

        seq_name_list = []

        seq_name = ""

        for line in gen_file:
            line = line.strip()

            if line.startswith(">"):
                seq_name = line.strip(">")
                seq_name_list.append(seq_name)
                seq_dic[seq_name] = []

            elif line != "":
                seq_dic[seq_name].append(line)


        for each_seqname in seq_name_list:

            seq = "".join(seq_dic[each_seqname]).upper()

            seq_tab.append([each_seqname] + list(seq))

        seq_sites_list = ["SeqName"]

        seq_length = len(seq_tab[0]) - 1

        for k in range(seq_length):
            seq_sites_list.append(str(k + 1))

        seq_pands = pd.DataFrame(seq_tab, columns=seq_sites_list).set_index("SeqName")

        return seq_pands




def calEnt(siteData, mode):
 
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




def lineage_mwic_run(lineage_name,
                     slither_window_list,
                     lineage_wic_df):

    each_lineage_setp_pro = []


    for each_window in slither_window_list:
        start_row = each_window[0]

        end_row = each_window[1]

        step_df = list(lineage_wic_df[start_row:end_row])

        mean_ic_per_win = np.mean(step_df)

        each_lineage_setp_pro.append(mean_ic_per_win)


    return [lineage_name,each_lineage_setp_pro]



def mwic_calculation(sites_probability_data,
                     lineage_name_list,
                     sites_count,site_map_dic,
                     windows_size,step_size,
                     output_path,thread_num):

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


    p = Pool(int(thread_num))

    lineage_mwic_list = []


    for each_lineage in lineage_name_list:

        lineage_wic_df = sites_probability_data[each_lineage]

        mwic_result = p.apply_async(lineage_mwic_run,
                                       args=(each_lineage,slither_window_list,lineage_wic_df))

        lineage_mwic_list.append(mwic_result)

    p.close()
    p.join()

    for each_calculation in lineage_mwic_list:
        lineage_mwic = each_calculation.get()

        step_probability_data[lineage_mwic[0]] = lineage_mwic[1]


    # step_probability_data.to_excel(
    #     excel_writer=output_path,
    #     index=False)

    step_probability_data.to_csv(output_path,index=False, sep=",")


    return step_probability_data,slither_window_list




def recomsplicing(sites_probability_data,
                  sites_count,
                  lineage_name_list,
                  recombination_frag,
                  step_size,max_mic,
                  max_recom_fragment,
                  recom_percentage):

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
                            and lineage_Ri_wic / (region_sites_count * max_mic) >= recom_percentage):
                        breakpoint_judgment.append([i + 1, region_left,
                                                    region_right,"True"])

                        flage_label = True

                # print(breakpoint_judgment)
                
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

                            if breakpoint_judgment[-1][0] >= frag_count:

                                break

                            else:
                                cursor_site = breakpoint_judgment[-1][0] + 1
                                cursor_center = lineage_frag_list[cursor_site - 1]

                        else:

                            for n in range(len(breakpoint_judgment)):
                                each_region = breakpoint_judgment[n]
                                ri = each_region[2] - each_region[1] + 1

                                if ri > max_recom_fragment:
                                    local_max_Ri = breakpoint_judgment[n - 1]
                                    detected_area.append([local_max_Ri[1],
                                                          local_max_Ri[2]])

                                    last_breakpoint = \
                                        breakpoint_judgment[n - 1][0]

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
        if not recom_region_dic[i]:
            del recom_region_dic[i]


    return recom_region_dic




def wic_plot(ref_lineage_name_list,
             site_list,
             input_data,
             query_lineage_name,
             output_path):


    lineage_count = len(ref_lineage_name_list)

    fig_high = int(lineage_count * 2)

    fig, ax = plt.subplots(len(ref_lineage_name_list), 1,
                           figsize=(
                           lineage_count, fig_high))

    fig.suptitle("Query seq: " + query_lineage_name, family="Arial")
    fig.tight_layout()

    cm1 = plt.cm.get_cmap("Reds")

    plt.subplots_adjust(top=0.95)

    for n in range(lineage_count):
        each_lineage = ref_lineage_name_list[n]
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
    plt.close()



def mwic_plot(gap_used,
              ref_lineage_name_list,
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

    for n in range(len(ref_lineage_name_list)):
        each_lineage = ref_lineage_name_list[n]

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

    legend_font_size = 8

    if len(ref_lineage_name_list) >= 15:
        legend_font_size = 6

    elif len(ref_lineage_name_list) >= 25:
        legend_font_size = 4


    legend_font = {"family": "Arial",
             "weight": "normal",
             "size": legend_font_size}


    if legend_location.upper() == "R":

        num1 = 1.05
        num2 = 1
        num3 = 2
        num4 = 0
        plt.legend(bbox_to_anchor=(num1, num2), loc=num3, borderaxespad=num4,prop=legend_font)
        plt.subplots_adjust(bottom=0.10, right=0.7)


    else:
        plt.legend(prop=legend_font)

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

    plt.savefig(output_file)

    plt.clf()

    plt.close()




def recombreak_run(lineage,lineage_wic,
                   run_number,breakwins):

    negative_lg_p_list = []

    for i in range(run_number):
        start_site = i

        end_site = i + breakwins

        central_pos = int((start_site + end_site) / 2)

        left_region = list(lineage_wic[start_site: central_pos - 1])
        right_region = list(lineage_wic[central_pos + 1: end_site])

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


    return [lineage,negative_lg_p_list]



def recombreak_plot(sites_probability_data,
                    lineage_name_list,
                    sites_count,
                    breakwins,
                    site_map_dic,
                    site_dir,
                    query_seq_prefix,
                    thread_num):


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

    for k in range(run_number):
        start_site = k

        end_site = k + breakwins

        central_pos = int((start_site + end_site) / 2)

        original_site = int(site_map_dic[str(central_pos + 1)])

        central_pos_list.append(original_site)

    breakpoint_data["Original_Index"] = central_pos_list

    p2 = Pool(int(thread_num))

    ca_result_list = []

    for lineage in lineage_name_list:

        lineage_wic = sites_probability_data[lineage]

        ca_result = p2.apply_async(recombreak_run,
                                    args=(lineage, lineage_wic, run_number,breakwins))

        ca_result_list.append(ca_result)

    p2.close()
    p2.join()


    for each_calculation in ca_result_list:
        lineage_ca = each_calculation.get()

        breakpoint_data[lineage_ca[0]] = lineage_ca[1]

    # print(breakpoint_data)

    breakpoint_data.to_excel(
        excel_writer=break_p_data,
        index=False)


    fig_high2 = int(lineage_count * 2)

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

    plt.close()



def merge_intervals(intervals):
 
    if not intervals or len(intervals) == 1:
        return intervals

    intervals.sort(key=lambda x: x[0])

    merged = [intervals[0]]

    for current in intervals[1:]:
        prev = merged[-1]

        if current[0] <= prev[1]:
        
            merged[-1] = [prev[0], max(prev[1], current[1])]
        else:
            merged.append(current)

    return merged

