# -*- coding: utf-8 -*-

"""
Virus recombination(VirusRecom), detecting recombination of viral lineages using information theory
Author: Zhou Zhi-Jian
Institution: Hunan University
Email: zjzhou@hnu.edu.cn
Copyright：Copyright (c) Zhi-Jian Zhou
Time: 2022/4/17 17:29

"""
import numpy as np
import pandas as pd

import scipy.stats as stats

from my_func import (make_dir, wic_calculation, mwic_calculation,
                     recomsplicing, wic_plot, mwic_plot, recombreak_plot)




def virus_infor_calculate(parameter_dic,
                          seq_matrix,
                          query_seq_prefix,
                          lineage_name_list,
                          sub_outdir):


    site_dir = sub_outdir + "/" + "WICs_of_sites"

    slide_window_dir = sub_outdir + "/" + "WICs_of_slide_window"

    site_ic_fig = (site_dir + "/"
                   + query_seq_prefix
                   + "_site_WIC_from_lineages.pdf")

    mwic_out_table = (slide_window_dir + "/"
                      + query_seq_prefix
                      + "_mWIC_from_lineages.xlsx")

    window_ic_fig = (slide_window_dir + "/"
                     + query_seq_prefix
                     + "_mWIC_from_lineages.pdf")



    is_use_gap = parameter_dic["gaps_use"]

    max_mic = 2

    if is_use_gap.upper() == "Y":
        max_mic = np.log2(5)


    windows_size = parameter_dic["windows_size"]

    step_size = parameter_dic["step_size"]

    max_recom_fragment = parameter_dic["max_recom_fragment"]

    recom_percentage = parameter_dic["recom_percentage"]

    cumulative_major = parameter_dic["cumulative_major"]

    y_start = parameter_dic["y_start"]

    legend_location = parameter_dic["legend"]



    make_dir(site_dir)

    site_ic_csv_path = parameter_dic["input_wic"]

    if site_ic_csv_path == "":

        print(">>> Treat " + query_seq_prefix
            + " as a potential recombination lineage..." + "\n")

        site_ic_csv_path = wic_calculation(seq_matrix, lineage_name_list,
                                           parameter_dic["calEnt_use"],
                                           site_dir, query_seq_prefix)

        if parameter_dic["only_wic"].upper() == "Y":

            print(">>> Finished, the output file is "
                  + "'" + site_ic_csv_path + "'"
                  + "\n")

            return

        else:
            pass

    else:
        pass


    sites_probability_data = pd.read_csv(site_ic_csv_path, sep=",",
                                         header=0)

    lineage_name_list = list(sites_probability_data.columns[2:])

    # print(lineage_name_list)

    current_site_list = sites_probability_data["Current_Index"]

    original_site_list = sites_probability_data["Original_Index"]

    sites_count = len(current_site_list)

    site_map_dic = {}

    for n in range(len(current_site_list)):
        original_site_index = str(original_site_list[n])

        current_site_index = str(current_site_list[n])

        site_map_dic[current_site_index] = original_site_index

    if parameter_dic["no_wic_figure"] != True:
        wic_plot(lineage_name_list, original_site_list,
                 sites_probability_data, query_seq_prefix, site_ic_fig)


    print("    " + "VirusRecom starts scanning using sliding window ..." + "\n")

    make_dir(slide_window_dir)

    step_probability_data, slither_window_list = mwic_calculation(
        sites_probability_data,
        lineage_name_list,
        sites_count, site_map_dic,
        windows_size, step_size,
        mwic_out_table)


    window_center_original = step_probability_data["Central_position(original)"]

    if parameter_dic["no_mwic_figure"] != True:
        mwic_plot(is_use_gap,lineage_name_list, window_center_original,
                  step_probability_data, query_seq_prefix, window_ic_fig,
                  y_start, legend_location)



    recombination_frag = {}


    for each_lineage in lineage_name_list:


        if not recombination_frag.__contains__(each_lineage):
            recombination_frag[each_lineage] = []

        potential_frag_list = []

        linegae_data_list = list(step_probability_data[
                                     each_lineage])

        for n in range(len(linegae_data_list)):

            line_ic_all = list(step_probability_data.iloc[n,
                               2:])

            linegae_mwic = linegae_data_list[n]

            if (linegae_mwic == max(line_ic_all)
                    and linegae_mwic / max_mic >= recom_percentage):

                windows_center = slither_window_list[n][2]

                potential_frag_list.append(windows_center)

        recombination_frag[each_lineage] = potential_frag_list

    # print(recombination_frag)

    recom_region_dic = recomsplicing(sites_probability_data,
                                     sites_count, lineage_name_list,
                                     recombination_frag,
                                     step_size, max_mic,
                                     max_recom_fragment,
                                     recom_percentage)


    parents_region = {}

    for each_lineage in recom_region_dic:
        region_list = recom_region_dic[each_lineage]

        parents_region[each_lineage] = []

        range_list = []

        for each_region in region_list:
            each_region_s = each_region[0] + 1
            each_region_n = each_region[1] + 1

            region_range = int(site_map_dic[str(each_region_n)]) - int(
                site_map_dic[str(each_region_s)])

            range_list.append(region_range)

        parents_region[each_lineage].append(sum(range_list))

    # print(parents_region)

    major_parent = ""



    try:

        major_parent = max(parents_region, key=parents_region.get)

    except:
        print("    " + "No significant recombination events were found in "
              + query_seq_prefix + "\n")

        print("    " +
            "Note: Please check whether it is a false negative caused by a higher cp value given!"
              + "\n")

        recom_report_path = (sub_outdir + "/"
                             + "Possible_recombination_event"
                             + "_detailed.txt")

        with open(recom_report_path, "w",
                  encoding="utf-8") as recom_report_file:

            recom_report_file.write(
                "No significant recombination events were found in "

                + query_seq_prefix + "\n"
                                     "Please check whether it is a false negative caused by a higher cp value given!")

        recom_report_jc_path = (sub_outdir + "/"
                                + "Possible_recombination_event"
                                + "_conciseness.txt")

        with open(recom_report_jc_path, "w",
                  encoding="utf-8") as recom_report_jc:

            recom_report_jc.write(
                "No significant recombination event were found in "
                + query_seq_prefix + "\n"
                + "Please check whether it is a false negative caused by a higher cp value given!")

        return

    finally:
        pass


    major_parent_ic = sum(list(sites_probability_data[major_parent]))

    mean_major_parent = major_parent_ic / sites_count

    if cumulative_major.upper() == "Y":

        accumulation_ic_list = []

        for lineage in lineage_name_list:
            accumulation_ic_list.append([lineage, sum(
                list(sites_probability_data[lineage])) / sites_count])

        hit_list_sort = sorted(accumulation_ic_list,
                               key=lambda x: x[1])

        max_ic_lineage = hit_list_sort[-1][0]

        if max_ic_lineage != major_parent:
            major_parent = max_ic_lineage


        major_parent_ic = sum(list(sites_probability_data[major_parent]))
        mean_major_parent = major_parent_ic / sites_count

    else:
        pass

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

                left_start_site_original = site_map_dic[
                    str(each_region_start + 1)]

                right_end_site_original = site_map_dic[
                    str(min(each_region_end + 1, sites_count))]

                try:

                    zihe_test = stats.mannwhitneyu(this_lineage_ic_count,
                                                   major_parent_ic_count,
                                                   alternative="two-sided")

                    p_value = zihe_test[1]

                    recombination_dic[each_lineage].append([
                        left_start_site_original + " to " +
                        right_end_site_original + "(mWIC: " + str(
                            region_mwic) + ")", "p_value: " + str(p_value)])

                    if p_value < 0.05:
                        other_parental_markers = True

                        significant_recombination[each_lineage].append([
                            left_start_site_original + " to " +
                            right_end_site_original + "(mWIC: " + str(
                                region_mwic) + ")", "p_value: " + str(p_value)])


                except:

                    recombination_dic[each_lineage].append([
                        left_start_site_original + " to " +
                        right_end_site_original, "p_value: 1"])

                finally:
                    pass




    for i in list(significant_recombination.keys()):
        if significant_recombination[i] == []:
            del significant_recombination[i]

    if other_parental_markers == True:

        print("\n" + "    "
              + "Possible major parent: " + major_parent
              + " (global mWIC: " + str(mean_major_parent) + ")" + "\n")

        print("    " + "Other possible parents and recombination region (map at the alignment): " + "\n")

        for each_lineage in significant_recombination:
            print("    " + each_lineage, significant_recombination[each_lineage])
            print("\n")

    else:
        print("    " + "No significant recombination events were found in "
              + query_seq_prefix + "\n")


    # print(recombination_dic)

    recom_report_path = (sub_outdir + "/"
                         + "Possible_recombination_event"
                         + "_detailed.txt")

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


    recom_report_jc_path = (sub_outdir + "/"
                            + "Possible_recombination_event"
                            + "_conciseness.txt")

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



    if parameter_dic["method"].upper() == "P" and parameter_dic[
        "breakpoints"].upper() == "Y":
        print("\n" + "    "
              + "VirusRecom is running the algorithm of search for "
                     "recombination breakpoint..."
              + "\n")

        breakwins = parameter_dic["breakwins"]

        recombreak_plot(sites_probability_data, lineage_name_list,
                        sites_count, breakwins,
                        site_map_dic,
                        site_dir, query_seq_prefix)




