from collections import OrderedDict
import argparse
from utils import get_aln, get_tree_data, remove_identical_adjacent_columns, save_data

def IndelHistory(a, T, O):
    """parameter: a, a path of multiple alignment file
    T, a p-ath of tree file
    return a dictionary of indel history"""
    # input alignment and trees
    aln_dict, alignment = get_aln(a)
    tree, ancestor_list, v_list, order, v_dict = get_tree_data(T)
    # remove identical adjacent columns
    removed_columns, new_alignment = remove_identical_adjacent_columns(alignment,aln_dict)
    print("identical adjacent columns removed")

    opt_list = []# a value list
    opt_seq_list = []# a seq list

    # Calculate opt0
    min_opt0 = float('inf')# inital min_value
    min_opt0_index = None #inital index

    for i in range(2 ** len(ancestor_list)):# find the opt from all combination
        time = bin(i)[2:].zfill(len(ancestor_list))# generate each combination
        sum_sign = 0

        for v in v_list:
            father = tree[v][0]# get father inforamtion

            if v in new_alignment:
                v_value = new_alignment[v][0]# get the v_value
            else:
                v_value = time[ancestor_list.index(v)]

            father_value = time[ancestor_list.index(father)]#get father's value
            sign = abs(int(v_value) - int(father_value))
            sum_sign += sign

        if sum_sign < min_opt0:# compare with the current value
            min_opt0 = sum_sign
            min_opt0_index = i

    opt_list.append(min_opt0)# save the first col
    opt0_seq = bin(min_opt0_index)[2:].zfill(len(ancestor_list))
    opt_seq_list.append(opt0_seq)
    print(f'columns 1/{len(list(new_alignment.values())[0])} is finished')

    # Calculate opt_i for i from 1 to n
    for i in range(1, len(list(new_alignment.values())[0])):# cal all the col
        min_opt_i = float('inf')# inital min_value
        min_opt_i_index = None# inital index

        for j in range(2 ** len(ancestor_list)):# find the opt from all combination
            time = bin(j)[2:].zfill(len(ancestor_list))
            dist = 0

            for v in v_list:
                father = tree[v][0]

                if v in new_alignment:# the v is in edge
                    v_value = new_alignment[v][i]
                    brother_v = abs(
                        int(new_alignment[v][i - 1]) - int(opt_seq_list[i - 1][ancestor_list.index(father)]))# get the brother's value
                else:# not in edge
                    v_value = time[ancestor_list.index(v)]
                    brother_v = abs(int(opt_seq_list[i - 1][ancestor_list.index(v)]) -
                                    int(opt_seq_list[i - 1][ancestor_list.index(father)]))# get the brother's value

                father_value = time[ancestor_list.index(father)]# get father's value
                sign = abs(int(v_value) - int(father_value))#the differencebetween v and its father
                different = abs(sign - brother_v)# different is the difference brtyween sign and its i-1
                dist += (sign * different) # the distance its the sum of sign * fdifferent

            sign_i = dist + opt_list[i-1]
            if sign_i < min_opt_i:# compared to find the opt
                min_opt_i = sign_i
                min_opt_i_index = j

        opt_list.append(min_opt_i)
        opt_i_seq = bin(min_opt_i_index)[2:].zfill(len(ancestor_list))
        opt_seq_list.append(opt_i_seq)
        print(f"columns {i + 1}/{len(list(new_alignment.values())[0])} is finished")

    # Reset the opt-seq to full result
    result_history = [''.join(nums) for nums in zip(*opt_seq_list)]  # join each positions seeq into one
    for position in removed_columns:  # reinput he removed columns
        for i in range(len(result_history)):
            result_history[i] = (result_history[i][:position] + result_history[i][position - 1]
                                 + result_history[i][position:])

    result_history = [''.join(nums) for nums in result_history]  # join the list into seq
    result_history.extend(alignment)  # put the edge into the result
    ancestor_list.extend(list(aln_dict.keys()))
    all_sequence = {k: v for k, v in zip(ancestor_list, result_history)}  # dict the name and the seq
    all_sequence = OrderedDict((key, all_sequence[key]) for key in order)  # re order the result
    print("whole process is done")
    return save_data(O, all_sequence)

def main(args):
    print("Received the following arguments:")
    print(f"-a option: {args.a}")
    print(f"-n option: {args.n}")
    print(f"-o option: {args.o}")
    IndelHistory(args.a, args.n, args.o)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="indel")

    parser.add_argument("-a", required=True, help="fasta format alignment file")
    parser.add_argument("-n", required=True, help="ancestors annotated phylogenetic tree in newick format")
    parser.add_argument("-o", required=True, help="folder location where the output files will be stored")
    args = parser.parse_args()

    main(args)
