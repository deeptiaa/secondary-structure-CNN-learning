import os
import numpy as np
import pandas as pd
import requests
from Bio import SeqIO
from io import StringIO
import Bio.PDB.Polypeptide
import random
import itertools
import more_itertools as mit
import re


# parameters:
#      "uniprot_id" - string representing uniprot id of desired protein.
# This method uses a given uniprot id to query the uniprot data and
# return a string respresention of the protein sequence.
# E.g. MADIT
def get_protein_seq(uniprot_id):
    # importing fasta file from uniprot.org and getting protein sequence
    # taken from StackOverflow:
    # https://stackoverflow.com/questions/52569622/protein-sequence-from-uniprot-protein-id-python
    url = "http://www.uniprot.org/uniprot/"
    complete_url = url + uniprot_id + ".fasta"
    response = requests.post(complete_url)
    data = ''.join(response.text)
    sequence = StringIO(data)
    protein_seq = list(SeqIO.parse(sequence, 'fasta'))

    # protein sequence as string (single-letter amino acids)
    string_seq = str(protein_seq[0].seq)

    # protein sequence w/ three-letter convention
    protein_seq = get_expanded_seq(string_seq)
    return protein_seq


# parameter:
#      "seq" - string representing protein sequence in 1-letter convention.
# This method takes protein sequence string with 1-letter convention and returns
# a protein sequence with 3-letter convention.
# E.g. ADE -> ALA ASP GLU
def get_expanded_seq(seq):
    expanded_list = []
    split_seq = list(seq)
    for letter in split_seq:
        three_letter_abbr = Bio.PDB.Polypeptide.one_to_three(letter)
        expanded_list.append(three_letter_abbr)
    exanded_string = " ".join(expanded_list)
    return(exanded_string)


# parameters:
#      "full_protein_split" - list of amino acids in full protein in 3 letter convention.
#                             E.g. ["ALA", "GLY", "TYR"]
#      "domain_split" - list of amino acids in protein domain in 3 letter convention.
#                       E.g. ["ALA", "GLY"]
# This method prints the index of the given domain within the given protein.
# Starting value is inclusive and the ending value is exclusive.
# E.g. [(0, 3)]
def get_index_range(full_protein_split, domain_split):
    indexes = []
    for i in range(len(full_protein_split)):
        if full_protein_split[i:i+len(domain_split)] == domain_split:
            indexes.append((i, i+len(domain_split)))
    print(indexes)
    indexes.clear()


# parameter:
#      "split_mutation_column" - list of mutations, split by comma if there are multiple.
# This method returns a list with wild-type residue (first letter) from variant.
def get_wild_type(split_mutation_column):
    wild_type_list = []
    w_letters = []
    for string in split_mutation_column:
        if "wild-type" in string[0]:
            wild_type = "wild_type"
        elif "-" in string[0] or len(string) == 0:
            wild_type = np.nan
        else:
            for val in string:
                mutation_name = val.strip(" ")
                w_letters.append(mutation_name[0])
                wild_type = ",".join(w_letters)
        wild_type_list.append(wild_type)
        w_letters.clear()
    return wild_type_list


# parameter:
#      "split_mutation_column" - list of mutations, split by comma if there are multiple.
# This method returns a list with mutation residue (last letter) from variant.
def get_mutation_type(split_mutation_column):
    mutation_list = []
    m_letters = []
    for string in split_mutation_column:
        if "wild-type" in string[0]:
            mutation = "wild-type"
        elif "-" in string[0] or len(string) == 0:
            mutation = np.nan
        else:
            for val in string:
                mutation_name = val.strip(" ")
                m_letters.append(mutation_name[-1])
                mutation = ",".join(m_letters)
        mutation_list.append(mutation)
        m_letters.clear()
    return mutation_list


# parameter:
#      "split_mutation_column" - list of mutations, split by comma if there are multiple.
# This method returns a list with the position of mutation (number) from variant.
def get_position(split_mutation_column):
    position_list = []
    p_letters = []
    for string in split_mutation_column:
        if "wild-type" in string[0]:
            position = "wild-type"
        elif "-" in string[0] or len(string) == 0:
            position = np.nan
        else:
            for val in string:
                mutation_name = val.strip(" ")
                p_letters.append(mutation_name[1:-1])
                position = ",".join(p_letters)
        position_list.append(position)
        p_letters.clear()
    return(position_list)


# parameter:
#      "df" - dataframe of protein data with "MUTATED_RES" and "POSITION" columns.
# This method returns a list with the correctly formatted variant (mutation-position form).
def get_mutations_names_list(df):
    formatted_list = []
    expanded_abbv = []
    for mutation, position in zip(df["MUTATED_RES"], df["POSITION"]):
        split_mutations = mutation.split(",")
        split_positions = position.split(",")
        if "wild-type" in split_mutations[0].lower() or "wild-type" in split_positions[0].lower():
            abbv_names = "WT"
        else:
            for mut, pos in zip(split_mutations, split_positions):
                three_letter_mut = Bio.PDB.Polypeptide.one_to_three(mut.upper())
                position = str(int(pos))
                combined_name = position + three_letter_mut
                expanded_abbv.append(combined_name)
                abbv_names = ", ".join(expanded_abbv)
        expanded_abbv.clear()
        formatted_list.append(abbv_names)
    return(formatted_list)


# Parameters:
#      "df" - protein data dataframe with "POSITION" column
# This method takes the position column in the dataframe and splits it in order
# to help remove or keep mutatations depending on their position.
def get_positions_split(df):
    position_list_split = []

    for item in df["POSITION"]:
        item = item.split(",")  # splits positions into list
        int_item = [int(i) for i in item]
        position_list_split.append(int_item)

    return position_list_split


# Parameters:
#      "stride file" - stride file of protein
#      "is_sec_struc" - list of boolean values for each secondary structure value
#                       if it is, true, else false
# returns list of boolean values indicating if position is secondary strcuture or not
def get_sec_struc_boolean(stride_file):
    is_sec_struc = []
    sec_struc_assign = []

    for line in stride_file:
        if line.startswith('ASG'):
            split_line = line.split()
            sec_struc_assign.append(split_line[5])

    for sec_struc in sec_struc_assign:
        if (sec_struc == 'C' or sec_struc == 'T'):
            is_sec_struc.append(False)
        else:
            is_sec_struc.append(True)

    return is_sec_struc


# Parameters:
#      "orig_df" -
#      "start" -
#      "end" -
#      "not_included_list"
# This method does even more helpful stuff
def get_domain_dataset(orig_df, start, end, not_included_list):
    in_domain_list = []

    for val in orig_df["positions_split"]:
        for position in val:
            if not_included_list.count(position - start) == 0:  # if value is not in the list of values to exclude
                if position >= start and position < end:
                    in_domain = True
                else:
                    in_domain = False
            else:
                in_domain = False
        in_domain_list.append(in_domain)

    orig_df['in_domain'] = in_domain_list
    # print(in_domain_list)
    condition = orig_df['in_domain'] == True
    rows = orig_df.loc[condition, :]

    in_domain_df = pd.DataFrame(columns=orig_df.columns)
    in_domain_df = in_domain_df.append(rows, ignore_index=True)
    in_domain_df = in_domain_df.drop(['in_domain'], axis=1)
    return in_domain_df


# Parameters:
# - orig_df: original dataframe with all mutations and "positions_split" column which has mutation positions in split list
#            as ints
# - sec_st_df: new dataframe with all rows that have mutations in the secondary structure of protein
# - mixed_df: new dataframe with all rows that have mutations in both in and out of the secondary stucture of the protein
# - start: (inclusive) index where the domain of the protein in PDB file starts
# - end: (inclusive) index where the domain of the protein in PDB file ends
def get_ss_dataset(orig_df, bool_ss_list, domain_start_index):
    has_sec_str = []

    for val in orig_df["positions_split"]:
        # list of boolean values that are true if all mutation positions in line are sec. strc.
        all_pos_sec_struc = []

        for position in val:
            # print(position - domain_start_index)
            # print(str(position) + " " + str(domain_start_index))
            if (bool_ss_list[position - domain_start_index] == False):  # line up ss_indexes w/ position
                all_pos_sec_struc.append(False)
            else:
                all_pos_sec_struc.append(True)

        # all pos sec struc should match val list
        # if there's a value in all_pos_sec_struc that's false, append false
        # otherwise, append true
        # print("val")
        # print(val)
        # print("bool")
        # print(all_pos_sec_struc)
        if (all_pos_sec_struc.count(False) == 0):
            has_only_sec_str = True
        else:
            has_only_sec_str = False

        # print(has_only_sec_str)
        has_sec_str.append(has_only_sec_str)
        all_pos_sec_struc.clear()

    # print(len(has_sec_str)) # should match dataframe length
    orig_df['has_sec_str'] = has_sec_str

    condition = orig_df['has_sec_str'] == True
    rows = orig_df.loc[condition, :]

    sec_str_df = pd.DataFrame(columns=orig_df.columns)
    sec_str_df = sec_str_df.append(rows, ignore_index=True)
    sec_str_df = sec_str_df.drop(['has_sec_str'], axis=1)
    orig_df = orig_df.drop(['has_sec_str'], axis=1)

    return sec_str_df


def get_not_ss_dataset(orig_df, bool_ss_list, domain_start_index):
    is_not_sec_str = []

    for val in orig_df["positions_split"]:

        all_pos_sec_struc = []

        for position in val:
            # print(position - domain_start_index)
            # print(str(position) + " " + str(domain_start_index))
            if (bool_ss_list[position - domain_start_index] == False):
                all_pos_sec_struc.append(False)
            else:
                all_pos_sec_struc.append(True)

        if (all_pos_sec_struc.count(True) == 0):
            has_no_sec_str = True
        else:
            has_no_sec_str = False

        is_not_sec_str.append(has_no_sec_str)
        all_pos_sec_struc.clear()

    orig_df['is_not_sec_str'] = is_not_sec_str

    condition = orig_df['is_not_sec_str'] == True
    rows = orig_df.loc[condition, :]

    not_sec_str_df = pd.DataFrame(columns=orig_df.columns)
    not_sec_str_df = not_sec_str_df.append(rows, ignore_index=True)
    not_sec_str_df = not_sec_str_df.drop(['is_not_sec_str'], axis=1)
    orig_df = orig_df.drop(['is_not_sec_str'], axis=1)

    return not_sec_str_df


# parameters:
#      "txt_name" - desired name of formatted txt file for network. E.g. "pab1"
#      "protein_seq" - string of protein sequence in 3 letter convention. E.g. ALA GLU TYR
#      "df" - dataframe with cleaned protein data. Must contain "variant" and "score"
#             columns.
# This method cleans the protein data and formats it into a txt that can be processed by the
# network. It also prints the name of the file out for reference.
def write_data_file(txt_name, protein_seq, df):
    file_name = txt_name + ".txt"
    path_name = "../ML Script Data Files/" + file_name
    print("Filename: " + file_name)

    datafile = open(path_name, "w+")
    datafile.write(protein_seq + "\n")
    for index in range(len(df) - 1):
        datafile.write(df["variant"].iloc[index] + ": " + str(df["score"].iloc[index]) + "\n")
    datafile.write(df["variant"].iloc[len(df) - 1] + ": " + str(df["score"].iloc[len(df) - 1]))
    datafile.close()


# Parameters:
#      "stride file" - stride file of protein
#      "is_sec_struc" - list of boolean values for each secondary structure value
#                       if it is, true, else false
# returns list of boolean values indicating if position is in an alpha helix or not
def get_alpha_boolean(stride_file):
    # print('hi')
    is_alpha = []
    alpha_assign = []

    for line in stride_file:
        # print(line)
        # print("why isn't this working")
        if line.startswith('ASG'):
            split_line = line.split();
            # print(split_line[5])
            alpha_assign.append(split_line[5])

    #     print(alpha_assign)

    alpha_letters = ['H', 'G', 'I']
    for alpha in alpha_assign:
        if (alpha_letters.count(alpha) != 0):
            is_alpha.append(True)
        else:
            is_alpha.append(False)

    #     print(alpha_assign)
    #     print(is_alpha)

    return is_alpha


def get_beta_boolean(stride_file):
    is_beta = []
    beta_assign = []

    for line in stride_file:
        if line.startswith('ASG'):
            split_line = line.split();
            beta_assign.append(split_line[5])

    beta_letters = ['E', 'B', 'b']
    for beta in beta_assign:
        if (beta_letters.count(beta) != 0):
            is_beta.append(True)
        else:
            is_beta.append(False)

    #     print(beta_assign)
    #     print(is_beta)
    #     print(len(is_beta))

    return is_beta


def get_turns_boolean(stride_file):
    is_turn = []
    turn_assign = []

    for line in stride_file:
        if line.startswith('ASG'):
            split_line = line.split();
            turn_assign.append(split_line[5])

    for turn in turn_assign:
        if (turn == "T"):
            is_turn.append(True)
        else:
            is_turn.append(False)

    print(turn_assign)
    print(is_turn)

    return is_turn


# Parameters:
#      "stride file" - stride file of protein
#      "is_sec_struc" - list of boolean values for each secondary structure value
#                       if it is, true, else false
# returns list of boolean values indicating if position is secondary strcuture or not
def get_all_sec_struc_boolean(stride_file):
    is_sec_struc = []
    sec_struc_assign = []

    for line in stride_file:
        if line.startswith('ASG'):
            split_line = line.split();
            sec_struc_assign.append(split_line[5])

    for sec_struc in sec_struc_assign:
        if (sec_struc == 'C'):
            is_sec_struc.append(False)
        else:
            is_sec_struc.append(True)

    return is_sec_struc


# limit number of mutations to some number
# **use after get_domain dataset

# Parameters:
#    "indexes" - a boolean list indicating positions with secondary structure (True - in ss, False - not in ss)
# This method returns a list of indexes to exclude in order to match the number of positions in secondary structure
# and out of secondary structure
def get_excluded_res(indexes):
    # find the groups of secondary structure
    ss_ind = [i for i, val in enumerate(indexes) if val == True]
    ss_ind_groups = list(find_index_range(ss_ind))

    # find the groups of non secondary structure
    not_ss_ind = [i for i, val in enumerate(indexes) if val == False]
    not_ss_ind_groups = list(find_index_range(not_ss_ind))

    ind_to_remove = []

    num_false = indexes.count(False)
    num_true = indexes.count(True)

    if (num_false < num_true):  # is mostly ss
        ind_to_remove = remove_indices_helper(not_ss_ind_groups, ss_ind_groups)  # chunk with not_ss groups
    elif (num_false > num_true):  # NOT mostly ss
        ind_to_remove = remove_indices_helper(ss_ind_groups, not_ss_ind_groups)

    print("Num True Indices: " + str(num_true))
    print("Num False Indices: " + str(num_false))
    print("Difference: " + str(abs(num_true - num_false)))
    print("Num Indices to Remove: " + str(len(ind_to_remove)))

    return ind_to_remove
    # return list of indices to NOT include


# Parameters:
#    "chunked_list" - list of ints/tuples representing either ss/not-ss regions that should be matched by corresponding
#                     ss/not-ss regions
#    "to_chunk_list" - list of regions representing regions with excess values that is matched to regions in chunked list
# This method is a helper method that returns a list of indices to remove in order to match the groups of secondary
# structure and non-secondary structure
def remove_indices_helper(chunked_list, to_chunk_list):
    remainder = []
    count_to_remove = 0

    #     print("chunked len: " + str(len(chunked_list)))

    #     print("to_chunk len: " + str(len(to_chunk_list)))

    for chunk, to_chunk in zip(chunked_list, to_chunk_list):  # zip goes through the smallest of the lists

        chunk_exp_list = expand_list(chunk)
        to_chunk_exp_list = expand_list(to_chunk)

        if (len(chunk_exp_list) < len(to_chunk_exp_list)):
            remainder.append(to_chunk_exp_list[len(chunk_exp_list):])  # will add indices to remove to remainder list
        elif (len(chunk_exp_list) > len(to_chunk_exp_list)):
            count_to_remove = count_to_remove + (len(chunk_exp_list) - len(to_chunk_exp_list))

    if (len(chunked_list) > len(to_chunk_list)):  # idk if this works
        count_to_remove = count_to_remove + len(expand_list(chunked_list[-1]))

    remainder = list(itertools.chain.from_iterable(remainder))
    if (len(to_chunk_list) > len(chunked_list)):
        remainder_copy = remainder.copy()
        print("remainder before: " + str(len(remainder_copy)))
        remainder.extend(expand_list(to_chunk_list[-1]))
        print("remainder after: " + str(len(remainder)))

    remainder = delete_random_elems(remainder, count_to_remove)

    return remainder
    # returns indices of values that are not to be included


# Parameters:
#    "val" - integer or tuple to be cast as a list
# This method is a helper method that either casts a single integer as a list or expands the range of a tuple
# (inclusive, inclusive)
def expand_list(val):
    val_list = []
    if isinstance(val, int):
        val_list.append(val)
    else:
        val_list = list(range(val[0], val[-1]))
        val_list.append(val[-1])

    return val_list


# https://www.codegrepper.com/code-examples/python/python+remove+n+random+elements+from+a+list
# Parameters:
#    "input_list" - list of values
#    "n" - number of random elements to delete from the list
# This method is a helper method that removes a given number of random elements from a list
def delete_random_elems(input_list, n):
    to_delete = set(random.sample(range(len(input_list)), n))
    return [x for i,x in enumerate(input_list) if not i in to_delete]


# determining the ranges of false values
# https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list

# Parameters:
#    "int_indexes" - list containing a location values for a protein
# This method is a helper method which determines consecutive values in list in order to group regions of
# secondary structure and non-secondary structure. It returns a list with integers and tuples (inclusive, inclusive)
# representing where a given type of region starts and stops in the protein.
def find_index_range(int_indexes):
    for segment in mit.consecutive_groups(int_indexes):
        segment = list(segment)
        if len(segment) == 1:
            yield segment[0] # yield is like return, except that it
                             # retains state to enable function to resume where
                             # it left off (sequenve of vals vs. 1)
        else:
            yield segment[0], segment[-1]

def format_mavedb_variant(df, variant_col_name, offset):
    new_var_col = []
    for variant in df[variant_col_name]:
        wild_type = Bio.PDB.Polypeptide.three_to_one(variant[2:5].upper())
        position = int(re.findall("[0-9]+", variant)[0]) + offset
        mut_type = Bio.PDB.Polypeptide.three_to_one(variant[-3:].upper())
        new_var_col.append(wild_type + str(position) + mut_type)
    return new_var_col
