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
import pandas as pd
import numpy as np
import platform
import argparse
import time
import matplotlib
import scipy.stats as stats

matplotlib.use("agg")  # avoid the  ModuleNotFoundError: No module named '_tkinter'

import matplotlib.pyplot as plt

from datetime import datetime

# plt.style.use("ggplot")


from my_func import (resolve_file_path,get_all_path,read_seq, calEnt, calEnt_gap,make_dir)

from plt_corlor_list import plt_corlor

from sequence_align import SeqAlign

app_dir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
if platform.system().lower() == "windows":
    app_dir = app_dir.replace("\\", "/")

sys.path.append(app_dir)


example_use = r'''
-----------------------------------------------------
☆ Example of use ☆
  (1) If the input-sequence data was not aligned: 
      virusrecom -q query.fasta -l lineage_dir -g n -m p -w 100 -s 20 -t 2 -o outdir

  (2) If the input-sequence has been aligned:
      virusrecom -a alignment.fasta -q query_name -l lineage_name_list.txt -g n -m p -w 100 -s 20 -o outdir

  Note: sequence file or folder need to enter their absolute path, the above is just a conceptual example. For detailed usage, please check the website https://github.com/ZhijianZhou01/virusrecom.

-----------------------------------------------------

'''


def virus_infor_calculate(seq_matrix,
                          query_seq_prefix,
                          lineage_name_list,
                          sub_outdir,
                          windows_size,
                          step_size,
                          max_recom_fragment,
                          recom_percentage,
                          calibrate_majorparent,
                          breakpoints,
                          breakwins,
                          calEnt_use,
                          max_mic,
                          y_start):

    print("\n" + ">>> Treat " + query_seq_prefix + " as a potential recombination lineage...")

    site_dir = sub_outdir + "/" + "WICs of sites"

    slide_window_dir = sub_outdir + "/" + "WICs of slide_window"

    make_dir(site_dir)
    make_dir(slide_window_dir)

    site_ic_table = (site_dir + "/" + run_id + "_"
                     + query_seq_prefix
                     + "_WIC contribution from lineage in sites" + ".xlsx")

    site_ic_fig = (site_dir + "/" + run_id + "_"
                   + query_seq_prefix
                   + "_WIC contribution from lineage in sites" + ".pdf")

    window_ic_table = (slide_window_dir + "/" + run_id + "_"
                       + query_seq_prefix
                       + "_WIC contribution from lineage in sliding window"
                       + ".xlsx")

    window_ic_fig = (slide_window_dir + "/" + run_id + "_"
                     + query_seq_prefix
                     + "_WIC contribution from lineage in sliding window"
                     + run_id + ".pdf")


    sites_probability_data = pd.DataFrame()

    query_seq = seq_matrix[
        seq_matrix.index.str.contains(query_seq_prefix) == True]

    query_seq_count = query_seq.shape[0]

    site_list = []

    for each_lineage in lineage_name_list:

        site_list = []

        lineage_seq_df = seq_matrix[
            seq_matrix.index.str.contains(each_lineage) == True]

        ent_probability = []

        p_ent = 0

        for (columnName, columnData) in lineage_seq_df.iteritems():
            # print(columnName)
            # print(columnData)

            site_list.append(int(columnName))

            lineage_site = list(columnData)

            IC = calEnt_use(columnData)

            if query_seq_count == 1:
                query_nt = query_seq[columnName][0]
                query_nt_ratio = 1

                query_nt_lineage_ratio = lineage_site.count(query_nt) / columnData.shape[0]

                p_ent = query_nt_lineage_ratio * IC * query_nt_ratio

                ent_probability.append(p_ent)


            else:

                query_seq_site_list = list(query_seq[columnName])

                query_seq_site_count = len(query_seq_site_list)

                maxpro_nt = max(query_seq_site_list,
                                key=query_seq_site_list.count)

                query_nt_ratio = query_seq_site_list.count(
                    maxpro_nt) / query_seq_site_count


                query_nt_lineage_ratio = lineage_site.count(maxpro_nt) / columnData.shape[0]

                p_ent = query_nt_lineage_ratio * IC * query_nt_ratio

                ent_probability.append(p_ent)

        sites_probability_data["Site"] = site_list
        sites_probability_data[each_lineage] = ent_probability

        print("The calculation of " + each_lineage + "'s recombination contribution to "
              + query_seq_prefix + " has been completed!" + "\n")

    sites_probability_data.to_excel(excel_writer=site_ic_table,
                                    index=False,
                                    encoding="utf-8")


    lineage_count = len(lineage_name_list)

    fig_high = int(max(site_list) * 3 / 10000) * 2

    fig, ax = plt.subplots(len(lineage_name_list), 1,
                           figsize=(
                           lineage_count, fig_high))

    fig.suptitle("Query seq: " + query_seq_prefix, family="Arial")
    fig.tight_layout()

    cm1 = plt.cm.get_cmap("Reds")
    plt.subplots_adjust(top=0.95)

    for n in range(lineage_count):
        each_lineage = lineage_name_list[n]
        ax_n = ax[n]
        y = list(sites_probability_data[each_lineage])
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


    plt.savefig(site_ic_fig)
    plt.clf()


    print("VirusRecom starts scanning using sliding window ..." + "\n")


    sites_count = sites_probability_data.shape[0]

    step_probability_data = pd.DataFrame()

    # print(sites_probability_data.loc[0, "Site"])


    slid_number = int(sites_count / step_size)

    original_site_list = []


    for each_lineage in lineage_name_list:

        each_lineage_setp_pro = []

        lineage_ic_df = sites_probability_data[each_lineage]
        # print(lineage_ic_df)

        original_site_list = []

        for n in range(slid_number):

            start_row = step_size * n

            end_row = min(start_row + windows_size, sites_count)

            # print(start_row,end_row)

            step_df = list(lineage_ic_df[start_row:end_row])

            # print(step_df)


            mean_ic_per_win = np.mean(step_df)

            # print(mean_ic_per_win)

            each_lineage_setp_pro.append(mean_ic_per_win)

            label_site = int((start_row + end_row) / 2)

            original_site = sites_probability_data.loc[
                label_site, "Site"]

            original_site_list.append(original_site)

            if end_row == sites_count:
                break

        step_probability_data["Central position"] = original_site_list
        step_probability_data[each_lineage] = each_lineage_setp_pro

        # print(each_lineage + "'s scan has been completed!" + "\n")

    step_probability_data.to_excel(
        excel_writer=window_ic_table,
        index=False,
        encoding="utf-8")


    fig, ax = plt.subplots()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    if gaps_use.upper() == "N":

        plt.ylim((y_start, 2))

    else:
        plt.ylim((y_start, 2.5))

    for n in range(len(lineage_name_list)):
        each_lineage = lineage_name_list[n]

        try:
            plt.plot(original_site_list,
                     list(step_probability_data[each_lineage]),
                     label=each_lineage, color=plt_corlor[n], linewidth=1)
        except:
            plt.plot(original_site_list,
                     list(step_probability_data[each_lineage]),
                     label=each_lineage, color="black", linewidth=1)
        finally:
            pass

    plt.legend()

    plt.margins(0)
    plt.subplots_adjust(bottom=0.10)
    plt.xlabel("Site in alignment", family="Arial")
    plt.ylabel("Mean of weighted information content", family="Arial")
    plt.title("Query seq: " + query_seq_prefix, family="Arial")

    x1_label = ax.get_xticklabels()
    [x1_label_temp.set_fontname("Arial") for x1_label_temp in
     x1_label]
    y1_label = ax.get_yticklabels()
    [y1_label_temp.set_fontname("Arial") for y1_label_temp in
     y1_label]

    # plt.show()
    plt.savefig(window_ic_fig)

    plt.clf()

    recombination_frag = {}


    for each_lineage in lineage_name_list:

        # if each_lineage != major_parent:

        if not recombination_frag.__contains__(each_lineage):
            recombination_frag[each_lineage] = []

        potential_frag_list = []

        linegae_data_list = list(step_probability_data[each_lineage])

        for n in range(len(linegae_data_list)):

            start_row = step_size * n

            end_row = min(start_row + windows_size, sites_count)

            line_ic_all = list(step_probability_data.iloc[n,1:])

            if (linegae_data_list[n] == max(line_ic_all)
                    and linegae_data_list[
                        n] / max_mic >= recom_percentage):


                windows_center = int((start_row + end_row) / 2)

                potential_frag_list.append(windows_center)

        recombination_frag[each_lineage] = potential_frag_list

    recom_region_dic = {}

    for each_lineage in recombination_frag:
        lineage_frag_list = recombination_frag[each_lineage]

        if lineage_frag_list != []:

            frag_count = len(lineage_frag_list)

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
                                ri = each_region[2] - each_region[1] + 1

                                if ri > max_recom_fragment:
                                    local_max_Ri = breakpoint_judgment[n - 1]
                                    detected_area.append([local_max_Ri[1],
                                                          local_max_Ri[2]])

                                    last_breakpoint = breakpoint_judgment[n - 1][0]

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

    parents_region = {}

    for each_lineage in recom_region_dic:
        region_list = recom_region_dic[each_lineage]

        parents_region[each_lineage] = []

        range_list = []

        for each_region in region_list:
            region_range = (sites_probability_data.loc[each_region[1], "Site"]
                            - sites_probability_data.loc[each_region[0], "Site"])


            range_list.append(region_range)

        parents_region[each_lineage].append(sum(range_list))


    # print(parents_region)

    major_parent = ""

    try:

        major_parent = max(parents_region, key=parents_region.get)

    except:
        print("No significant recombination events were found in "
              + query_seq_prefix + "\n")

        recom_report_path = (sub_outdir + "/" + run_id + "_"
                             + "Possible recombination event in "
                             + query_seq_prefix + "_detailed.txt")

        with open(recom_report_path, "w", encoding="utf-8") as recom_report_file:

            recom_report_file.write(
                "No significant recombination events were found in "

                + query_seq_prefix)

        recom_report_jc_path = (sub_outdir + "/" + run_id + "_"
                                + "Possible recombination event in "
                                + query_seq_prefix + "_conciseness.txt")

        with open(recom_report_jc_path, "w", encoding="utf-8") as recom_report_jc:

            recom_report_jc.write(
                "No significant recombination events were found in "
                + query_seq_prefix)

        exit()

    finally:
        pass


    major_parent_ic = sum(list(sites_probability_data[major_parent]))
    mean_major_parent = major_parent_ic / sites_count



    if calibrate_majorparent.upper() != "N":

        accumulation_ic_list = []

        for lineage in lineage_name_list:
            accumulation_ic_list.append([lineage, sum(list(sites_probability_data[lineage]))/sites_count])

        hit_list_sort = sorted(accumulation_ic_list,
                               key=lambda x: x[1])


        max_ic_lineage = hit_list_sort[-1][0]


        if max_ic_lineage != major_parent:
            major_parent = max_ic_lineage


        major_parent_ic = sum(list(sites_probability_data[major_parent]))
        mean_major_parent = major_parent_ic / sites_count


    other_parental_markers = False



    recombination_dic = {}
    significant_recombination = {}

    for each_lineage in recom_region_dic:

        if each_lineage != major_parent:

            recombination_dic[each_lineage] = []
            significant_recombination[each_lineage] = []

            recom_region_list = recom_region_dic[each_lineage]

            for each_region in recom_region_list:
                each_region_start = each_region[0]
                each_region_end = each_region[1]

                this_lineage_ic_count = list(
                    sites_probability_data[each_lineage][
                    each_region_start:each_region_end])

                major_parent_ic_count = list(
                    sites_probability_data[major_parent][
                    each_region_start:each_region_end])

                region_mwic = sum(this_lineage_ic_count) / (
                            each_region_end - each_region_start + 1)


                left_start_site_original = sites_probability_data.loc[
                    each_region_start, "Site"]

                right_end_site_original = sites_probability_data.loc[
                    min(each_region_end, sites_count - 1), "Site"]


                try:

                    zihe_test = stats.mannwhitneyu(this_lineage_ic_count,
                                                   major_parent_ic_count,
                                                   alternative="two-sided")

                    p_value = zihe_test[1]

                    recombination_dic[each_lineage].append([str(
                        left_start_site_original) + " to " + str(
                        right_end_site_original) + "(mWIC: " + str(
                        region_mwic) + ")", "p_value: " + str(p_value)])

                    if p_value < 0.05:
                        other_parental_markers = True

                        significant_recombination[each_lineage].append([str(
                            left_start_site_original) + " to " + str(
                            right_end_site_original) + "(mWIC: " + str(
                            region_mwic) + ")", "p_value: " + str(p_value)])


                except:

                    recombination_dic[each_lineage].append([str(
                        left_start_site_original) + " to " + str(
                        right_end_site_original), "p_value: 1"])

                finally:
                    pass

    for i in list(significant_recombination.keys()):
        if significant_recombination[i] == []:
            del significant_recombination[i]

    if other_parental_markers == True:

        print("\n" + "Possible major parent: " + major_parent
              + "(global mWIC: " + str(mean_major_parent) + ")" + "\n")

        print("Other possible parents: " + "\n")

        print("Possible recombination region map at aligned genomes: " + "\n")

        for each_lineage in significant_recombination:
            print(each_lineage, significant_recombination[each_lineage])

    else:
        print("No significant recombination events were found in "
              + query_seq_prefix + "\n")


    recom_report_path = (sub_outdir + "/" + run_id + "_"
                         + "Possible recombination event in "
                         + query_seq_prefix + "_detailed.txt")

    with open(recom_report_path, "w", encoding="utf-8") as recom_report_file:

        if other_parental_markers == False:
            recom_report_file.write(
                "No significant recombination events were found in "

                + query_seq_prefix + "\n" * 2
                + "The most similar lineage: "
                + major_parent
                + "(global mWIC: " + str(mean_major_parent) + ")"
                + "\n" * 2
                + "The similar other lineage and not significant recombination regions (p>0.05):"
                + "\n")


        else:

            recom_report_file.write("Possible major parent: "
                                    + major_parent
                                    + "(global mWIC: "
                                    + str(mean_major_parent)
                                    + ")"
                                    + "\n" * 2
                                    + "Possible other parents and recombination regions:"
                                    + "\n")


        for key in recombination_dic:
            event_list = recombination_dic[key]

            recom_report_file.write(key + "\t")

            for each_envent in event_list:
                recom_report_file.write(", ".join(each_envent) + "\t")

            recom_report_file.write("\n")

        recom_report_file.write(
            "\n" + "Significance test of recombinant regions using Mann-Whitney-U test with two-tailed probabilities, "
                   "p-value less than 0.05 indicates a significant difference.")

    recom_report_jc_path = (sub_outdir + "/" + run_id + "_"
                            + "Possible recombination event in "
                            + query_seq_prefix + "_conciseness.txt")

    with open(recom_report_jc_path, "w", encoding="utf-8") as recom_report_jc:

        if other_parental_markers == False:
            recom_report_jc.write(
                "No significant recombination events were found in "
                + query_seq_prefix)

            recom_report_jc.write(
                "\n" + "Significance test of recombinant regions using Mann-Whitney-U test with two-tailed probabilities, "
                       "p-value less than 0.05 indicates a significant difference.")


        else:
            recom_report_jc.write("Possible major parent: "
                                  + major_parent
                                  + "(global mWIC: " + str(
                mean_major_parent) + ")"
                                  + "\n" * 2
                                  + "Other possible parents and significant recombination regions (p<0.05):"
                                  + "\n")


            for key in significant_recombination:
                event_list = significant_recombination[key]

                recom_report_jc.write(key + "\t")

                for each_envent in event_list:
                    recom_report_jc.write(", ".join(each_envent) + "\t")

                recom_report_jc.write("\n")

            recom_report_jc.write(
                "\n" + "Significance test of recombinant regions using Mann-Whitney-U test with two-tailed probabilities, "
                       "p-value less than 0.05 indicates a significant difference.")

    if method.upper() == "P" and breakpoints.upper() == "Y":

        print("\n" + "VirusRecom is running the algorithm of search for "
                     "recombination breakpoint..."
              + "\n")



        break_p_map = (site_dir + "/" + run_id + "_"
                       + query_seq_prefix
                       + "_ -lg(p-value) for potential breakpoint.pdf")

        break_p_data = (site_dir + "/" + run_id + "_"
                        + query_seq_prefix
                        + "_ -lg(p-value) for potential breakpoint.xlsx")

        breakpoint_data = pd.DataFrame()


        run_number = sites_count - breakwins + 1

        central_pos_list = []
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

                original_site = sites_probability_data.loc[
                    central_pos, "Site"]
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

            breakpoint_data["Site"] = central_pos_list
            breakpoint_data[lineage] = negative_lg_p_list

        # print(breakpoint_data)

        breakpoint_data.to_excel(
            excel_writer=break_p_data,
            index=False,
            encoding="utf-8")


        fig_high2 = int(max(central_pos_list) * 3 / 10000) * 2

        figs, axs = plt.subplots(lineage_count, 1, figsize=(lineage_count,fig_high2))

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


def calculate_memory():
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    memory = info.uss / 1024 / 1024
    return memory



if __name__ == "__main__":

    print("\n" + "-------------------------------------------------")

    print("  Name: Virus Recombination(VirusRecom)")

    print(
        "  Description: Detecting recombination of viral lineages (or subtypes) using information theory.")

    print("  Version: 1.0 (2022-04-18)")

    print("  Author: Zhi-Jian Zhou")

    print("-------------------------------------------------" + "\n")


    def parameter():
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            prog="virusrecom",
            description="",
            epilog=example_use)


        parser.add_argument(
            "-a", dest="alignment",
            help="FilePath of an aligned sequence set (*.fasta format) containing all sequences used for analysis, then the sequence alignment will be skipped. Default value is null. If “-a” parameter was used, the name of each sequence in aligned sequence set requires containing the mark (a unique string) of the lineage.",
            default="")


        parser.add_argument(
            "-q", dest = "query",
            help="FilePath of query lineage (usually potential recombinant, *.fasta format). Note, if the ‘-a’ parameter has been used, please enter the mark (a unique string) of query lineage here, such as ‘-q xxxx’, not a FilePath. Using ‘-q auto’ and all lineages will be scanned as potential recombinants in turn.",
            default="")


        parser.add_argument(
            "-l", dest="lineage",
            help="DirPath of reference lineages. One sequence file (*.fasta format) per lineage, and each lineage could contain multiple sequences. Note, if the ‘-a’ parameter has been used, please enter a file path of a text file containing the mark (a unique string) of lineage here, not a DirPath.",
            default="/home")


        parser.add_argument("-g", dest="gap",
                            help="Gaps (-) in the alignment were used in analysis? ‘-g y’ means to reserve gaps, and ‘-g n’ means to delete gaps.",
                            type=str,
                            default="n")


        parser.add_argument("-m", dest="method",
                            help="Scanning method of recombination analysis. ‘-m p’ means use polymorphic sites only, ‘-m a’ means all the monomorphic sites and polymorphic sites.",
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
                            help="The cutoff threshold of proportion (cp, default value was 0.9) used for searching recombination regions when mWIC/EIC >= cp, the maximum value of cp is 1.",
                            type=float,
                            default=0.9)

        parser.add_argument("-cm", dest="calibrate",
                            help="Whether to simply use the max cumulative WIC of all sites to identified the major parent. The default value is ‘n’ and means ‘no’. If required, please specify ‘-cm y’.",
                            type=str,
                            default="n")


        parser.add_argument("-b", dest="breakpoint",
                            help="Whe ther to run the breakpoint scan of recombination. ‘-b y’ means yes, ‘-b n’ means no. Note: this option only takes effect when ‘-m p’ has been specified.",
                            type=str,
                            default="n")


        parser.add_argument("-bw", dest="breakwin",
                            help="The window size (polymorphic sites, default value is 200) used for breakpoint scan. The step size is fixed at 1. Note: this option only takes effect when ‘-m p -b y’ has been specified.",
                            type=int,
                            default=200)


        parser.add_argument(
            "-o", dest="outdir",
            help="The path of the outdir of results.",
            type=str,
            default="/home")

        parser.add_argument(
            "-t", dest="thread",
            help="Number of threads used for the multiple sequence alignments (MSA), the default value is 1.",
            type=int,
            default=1)


        parser.add_argument(
            "-y", dest="y_start",
            help="Specify the starting value of the Y-axis scale in the picture, the default value is 0.",
            type=float,
            default=0.0)

        myargs = parser.parse_args(sys.argv[1:])

        return myargs


    myargs = parameter()


    # run_id = str(time.time()).split(".")[0]
    run_id = ""

    start = datetime.today().now()

    start_memory = calculate_memory()


    #  parameter
    seq_aligned_path = myargs.alignment             # path of the alignment

    query_seq_path = myargs.query

    lineage_file_dir = myargs.lineage

    gaps_use = myargs.gap                     #  strategy of handling gap

    method = myargs.method                    #  Scan method

    windows_size = myargs.window              #  window size

    step_size = myargs.step                   #  step size

    max_recom_fragment = myargs.max_region    #  maximum recombination region

    recom_percentage = myargs.percentage      #  the specific cutoff threshold of proportion

    calibrate_majorparent = myargs.calibrate

    breakpoints = myargs.breakpoint           #  whether to search for recombination breakpoints

    breakwins = myargs.breakwin               #  Window size used to search for recombination breakpoints

    thread_num = myargs.thread                #  thread of MAS

    y_start = myargs.y_start                  #  Y-axis starting point when plotting

    out_dir = myargs.outdir


    if gaps_use.upper() not in ["N","Y"]:
        print("Error, the parameter after '-g' is incorrect!")
        exit()

    if method.upper() not in ["P","A"]:
        print("Error, the parameter after '-m' is incorrect!")
        exit()

    print(myargs)



    aligned_out_path = ""


    query_seq_prefix = ""

    run_record = ""

    lineage_name_list = []



    # If the sequences were not aligned beforehand
    if seq_aligned_path == "":

        query_seq_path = query_seq_path.replace("\\", "/")

        query_seq_dir, query_seq_prefix = resolve_file_path(query_seq_path)

        run_record = out_dir + "/" + "run_record"

        make_dir(out_dir)
        make_dir(run_record)

        aligned_out_path = (run_record + "/" + query_seq_prefix
                            + "_" + run_id + "_merge_mafft.fasta")


        seq_align_task = SeqAlign(query_seq_path,
                                  lineage_file_dir,
                                  run_record,
                                  run_id,
                                  thread_num,
                                  aligned_out_path)

        lineage_name_list = seq_align_task.run()



    # If the sequences were aligned beforehand
    else:

        aligned_out_path = seq_aligned_path.replace("\\", "/")

        run_record = out_dir + "/" + "run_record"

        make_dir(out_dir)
        make_dir(run_record)

        query_seq_prefix = query_seq_path

        with open(lineage_file_dir) as lineage_file:
            for line in lineage_file:
                line = line.strip()
                if line != "":
                    lineage_name_list.append(line)




    print("VirusRecom starts calculating weighted information content from each lineage..."
          + "\n")

    seq_pd = read_seq(aligned_out_path)


    seq_pd_clean = seq_pd
    calEnt_use = calEnt_gap
    max_mic = np.log2(5)


    if gaps_use.upper() == "N":

        all_site_label = list(seq_pd.columns)

        seq_pd_clean = seq_pd[~seq_pd.isin(["-"])].dropna(axis=1)
        calEnt_use = calEnt
        max_mic = 2

        record_gap_path = (run_record + "/"
                           + "Record of deleted gap sites_"
                           + run_id + ".txt")

        seq_pd_clean_site = list(seq_pd_clean.columns)

        with open(record_gap_path, "w",encoding="utf-8") as record_gap_file:

            record_gap_file.write("These sites with gap(-) in the file of "
                                  + aligned_out_path + "\n")

            for each_site in all_site_label:
                if each_site not in seq_pd_clean_site:
                    record_gap_file.write("Site " + each_site + "\n")


    if method.upper() == "P":

        seq_pd_clean_site_old = list(seq_pd_clean)

        seq_pd_clean = seq_pd_clean.loc[:, (seq_pd_clean != seq_pd_clean.iloc[0]).any()]

        seq_pd_clean_site_new = list(seq_pd_clean)


        record_same_sites_path = (run_record + "/"
                                  + "Record of same sites in aligned sequence_"
                                  + run_id + ".txt")


        with open(record_same_sites_path, "w", encoding="utf8") as same_sites_file:
            same_sites_file.write("These same sites(no variation) in the file of "
                                  + aligned_out_path + "\n")


            for each_site in seq_pd_clean_site_old:
                if each_site not in seq_pd_clean_site_new:
                    same_sites_file.write("Site " + each_site + "\n")


    if seq_aligned_path != "" and query_seq_prefix.upper() == "AUTO":

        for each_lineage in lineage_name_list:
            query_lineage_run = each_lineage

            new_lineage_name_list = []

            for lineage in lineage_name_list:
                if lineage != query_lineage_run:
                    new_lineage_name_list.append(lineage)

            sub_outdir = out_dir + "/Results/" + query_lineage_run

            virus_infor_calculate(seq_pd_clean,
                                  query_lineage_run,
                                  new_lineage_name_list,
                                  sub_outdir,
                                  windows_size,
                                  step_size,
                                  max_recom_fragment,
                                  recom_percentage,
                                  calibrate_majorparent,
                                  breakpoints,
                                  breakwins,
                                  calEnt_use,
                                  max_mic,
                                  y_start)


    else:
        virus_infor_calculate(seq_pd_clean,
                              query_seq_prefix,
                              lineage_name_list,
                              out_dir,
                              windows_size,
                              step_size,
                              max_recom_fragment,
                              recom_percentage,
                              calibrate_majorparent,
                              breakpoints,
                              breakwins,
                              calEnt_use,
                              max_mic,
                              y_start)


    duration = datetime.today().now() - start

    end_memory = calculate_memory()

    print(f"Occupied {end_memory - start_memory}MB memory in total")

    print("Take " + str(duration) + " seconds in total." + "\n")
