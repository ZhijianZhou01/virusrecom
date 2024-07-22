# -*- coding: utf-8 -*-

"""
Author: Zhou Zhi-Jian
Email: zjzhou@hnu.edu.cn
Time: 2024/7/22 0:40

"""


import sys
import os
import psutil
import pandas as pd

from multiprocessing import Pool



def calculate_memory():
    pid = os.getpid()
    p = psutil.Process(pid)
    info = p.memory_full_info()
    memory = info.uss / 1024 / 1024
    return memory


def gap_processing(seq_df,gaps_use):

    gap_site = []

    if gaps_use.upper() == "N":

        print("    " + "VirusRecom is removing sites (columns) containing gap (-)..." + "\n")

        old_site_label = list(seq_df.columns)

        seq_pd_clean = seq_df[~seq_df.isin(["-"])].dropna(axis=1)  # 去除含有gap的列（位点)

        new_site_label = list(seq_pd_clean.columns)

        for each_site in old_site_label:
            if each_site not in new_site_label:
                gap_site.append(each_site)


        return seq_pd_clean,gap_site

    else:
        return seq_df,gap_site


def site_model(seq_df, method):

    invariant_site = []

    if method.upper() == "P":

        print("    " + "VirusRecom is extracting polymorphic sites..." + "\n")

        old_site_label = list(seq_df.columns)

        seq_pd_clean = seq_df.loc[:, (seq_df != seq_df.iloc[0]).any()]

        new_site_label = list(seq_pd_clean.columns)

        for each_site in old_site_label:
            if each_site not in new_site_label:
                invariant_site.append(each_site)

        return seq_pd_clean, invariant_site

    else:
        return seq_df, invariant_site




def piecewise_read(seq_file_path,
                   block_size,
                   data_blocks):

    print("    " + "VirusRecom is importing data blocks {}".format(data_blocks) + "\n")

    with open(seq_file_path, "r",encoding="utf-8") as gen_file:

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

        site_start = block_size * (data_blocks - 1)  #

        site_end = block_size * data_blocks


        if site_end >= seq_length:

            print("    Load sites: " + "{}".format(str(site_start + 1))
                  + " - {}".format(seq_length) + "\n")

            if data_blocks == 1:

                return seq_pands, False

            else:

                sub_pands = seq_pands.iloc[:, site_start:seq_length]

                return sub_pands,False

        else:
            sub_pands = seq_pands.iloc[:, site_start:site_end]

            print("    Load sites: " + "{}".format(str(site_start + 1))
                                  + " - {}".format(site_end) + "\n")

            return sub_pands,True





def lineage_wic_run(query_seq, each_lineage,
                    lineage_seq_df,used_calEnt,
                    query_seq_prefix):

    ent_probability = []

    query_seq_count = query_seq.shape[0]

    for (columnName, columnData) in lineage_seq_df.items():

        lineage_site = list(columnData)

        if len(lineage_site) == 1:
            IC = used_calEnt(columnData, "single")

        else:

            IC = used_calEnt(columnData, "")

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


    # print("    " +
    #       "The calculation of " + each_lineage + "'s recombination contribution to "
    #       + query_seq_prefix + " has been completed!" + "\n")

    return [each_lineage, ent_probability]



def wic_calculation(seq_data,
                    lineage_name_list,
                    query_seq_prefix,
                    parameter_dic,
                    site_idenx):


    used_calEnt = parameter_dic["calEnt_use"]

    thread_num = parameter_dic["thread_num"]


    query_seq = seq_data[seq_data.index.str.contains(query_seq_prefix) == True]

    query_seq_count = query_seq.shape[0]

    if query_seq_count == 0:
        print(">>> Error! The query lineage '" + query_seq_prefix + "'"
              + " is not in sequence data." + "\n")

        sys.exit()

    original_site_list = [int(x) for x in query_seq.columns.tolist()]

    sites_count = len(original_site_list)

    current_site_list = [site_idenx + i + 1 for i in range(sites_count)]

    columns_name_list = ["Original_Index", "Current_Index"] + lineage_name_list

    sites_wic_data = pd.DataFrame(columns=columns_name_list,
                                  index=current_site_list)

    sites_wic_data["Original_Index"] = original_site_list

    sites_wic_data["Current_Index"] = current_site_list

    p = Pool(int(thread_num))

    lineage_wic_list = []

    for each_lineage in lineage_name_list:

        lineage_seq_df = seq_data[seq_data.index.str.contains(each_lineage) == True]

        lineage_result = p.apply_async(lineage_wic_run,
                                       args=(query_seq, each_lineage, lineage_seq_df,
                                             used_calEnt, query_seq_prefix))

        lineage_wic_list.append(lineage_result)


    p.close()
    p.join()


    for each_calculation in lineage_wic_list:
        lineage_wic = each_calculation.get()

        sites_wic_data[lineage_wic[0]] = lineage_wic[1]


    return sites_wic_data,sites_count




def wic_compute_engine(seq_file_path,
                       query_seq_prefix,
                       lineage_name_list,
                       parameter_dic,
                       site_dir):

    # start_memory = calculate_memory()

    block_size = parameter_dic["block_size"]

    data_blocks = 0

    read_mark = True

    sites_wic_pd_list = []

    gap_site_list = []

    invariant_site_list = []

    site_idenx = 0

    while read_mark:

        data_blocks += 1

        seq_pd,read_mark = piecewise_read(seq_file_path, block_size, data_blocks)

        seq_pd_clean,gap_record = gap_processing(seq_pd, parameter_dic["gaps_use"])

        gap_site_list.append(gap_record)


        seq_pd_use,invariant_record = site_model(seq_pd_clean, parameter_dic["method"])

        invariant_site_list.append(invariant_record)


        if seq_pd_use.empty:
            print("    WIC for data_blocks {}".format(data_blocks) + " have been completed." + "\n")
            continue


        sites_wic_pd,last_site_count = wic_calculation(seq_pd_use,
                                       lineage_name_list,
                                       query_seq_prefix,
                                       parameter_dic,
                                       site_idenx)

        site_idenx += last_site_count


        sites_wic_pd_list.append(sites_wic_pd)


        print("    WIC for data_blocks {}".format(data_blocks) + " have been completed." + "\n")

    # print(sites_wic_pd_list)

    print(">>> " + "The WIC calculations of {} sites have been completed.".format(site_idenx)  + "\n")



    if not sites_wic_pd_list:
        print(">>> There are no variant sites in the data!")

        sys.exit()

    sites_probability_data = sites_wic_pd_list[0]

    for cookies_num in range(1,len(sites_wic_pd_list)):
        sites_probability_data = pd.concat([sites_probability_data,sites_wic_pd_list[cookies_num]],
                                           ignore_index= True)


    site_ic_table = (site_dir + "/"
                     + query_seq_prefix
                     + "_site_WIC_from_lineages.xlsx")

    site_ic_csv = (site_dir + "/"+ query_seq_prefix + "_site_WIC.csv")

    sites_probability_data.to_excel(excel_writer=site_ic_table,
                            index=False)

    sites_probability_data.to_csv(site_ic_csv, index=False, sep=",")


    # end_memory = calculate_memory()
    #
    # used_memory = end_memory - start_memory
    #
    #
    # print(">>> " + f"Occupied {used_memory}MB memory in the calculations of sites' WIC." + "\n")

    return site_ic_csv







