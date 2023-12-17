from ete3 import Tree
from collections import defaultdict
from collections import OrderedDict
import argparse

tree ='./toy_data/toy_data1/input_tree.nwk'#'./toy_data/CYP2U_165.nwk'
seqs ='./toy_data/toy_data1/extants.aln'#"./toy_data/CYP2U_165.aln"#

def readFasta(string, ignore = False):
    """ Read the given string as FASTA formatted data and return the list of
        sequences contained within it."""
    seqlist = []    # list of sequences contained in the string
    seqname = None  # name of *current* sequence
    seqinfo = None
    seqdata = []    # sequence data for *current* sequence
    for line in string.splitlines():    # read every line
        if len(line) == 0:              # ignore empty lines
            continue
        if line[0] == '>':  # start of new sequence
            if seqname:     # check if we've got one current
                try:
                    current = (seqdata, seqname,seqinfo)
                    seqlist.append(current)
                except RuntimeError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            # now collect data about the new sequence
            seqinfo = line[1:].split() # skip first char
            if len(seqinfo) > 0:
                try:
                    seqname = seqinfo[0]
                    if len(seqinfo) > 0: # more than a name
                        edited_info = ''
                        for infopart in seqinfo[1:]:
                            edited_info += infopart + ' '
                        seqinfo = edited_info
                    else:
                        seqinfo = ''
                except IndexError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            else:
                seqname = ''
                seqinfo = ''
            seqdata = []
        else:               # we assume this is (more) data for current
            cleanline = line.split()
            for thisline in cleanline:
                seqdata.extend(tuple(thisline.strip('*')))
    # we're done reading the file, but the last sequence remains
    if seqname:
        try:
            lastseq = (seqdata,seqname,seqinfo)
            seqlist.append(lastseq)
        except RuntimeError as errmsg:
            if not ignore:
                raise RuntimeError(errmsg)
    return seqlist

def get_aln(filename):
    fh = open(filename)
    seqlist = []
    batch = ''  # a batch of rows including one or more complete FASTA entries
    rowcnt = 0
    for row in fh:
        row = row.strip()
        if len(row) > 0:
            if row.startswith('>') and rowcnt > 0:
                more = readFasta(batch)
                if len(more) > 0:
                    seqlist.extend(more)
                batch = ''
                rowcnt = 0
            batch += row + '\n'
            rowcnt += 1
    if len(batch) > 0:
        more = readFasta(batch)
        if len(more) > 0:
            seqlist.extend(more)
    fh.close()


    aln_dict = {}
    for index in range(len(seqlist)):
        seq = ''.join('0' if char == '-' else '1' for char in seqlist[index][0])
        aln_dict.update({seqlist[index][1]: seq})
    alignment = list(aln_dict.values())
    return aln_dict, alignment


def get_tree_data(tree):
    """create neighbor dict """
    tree_file = open(tree, "r")
    my_tree = tree_file.read() + ";"
    tree = Tree(my_tree, format=1)
    # create neighbourhood object
    ancestor_list = []
    v_list = []
    order = []
    v_dict = defaultdict(list)
    V_father_dict = defaultdict(list)
    number=0
    for n in tree.traverse():
        # if the v have no name, give it a name
        if n.name=="":
            n.name='A'+str(number)
            number +=1
        order.append(n.name)
        if n.is_leaf() == False:
            ancestor_list.append(n.name)
            # if the v have no name, give it a name
            for c in n.children:
                if c.name == "":
                    c.name = 'A' + str(number)
                    number += 1
                V_father_dict[c.name] += [n.name]
                v_dict[n.name] += [c.name]
        if n.is_root() == False:
            v_list.append(n.name)
    return V_father_dict, ancestor_list, v_list, order,v_dict


def remove_identical_adjacent_columns(alignment, aln):
    removed_columns = []
    alignment = [list(seq) for seq in alignment]
    current_index = 0
    original_index = 0

    while current_index < len(alignment[0]) - 1:
        current_column = [seq[current_index] for seq in alignment]
        next_column = [seq[current_index + 1] for seq in alignment]

        # check and remove the identical columns
        if current_column == next_column:
            removed_columns.append(original_index + 1)
            for seq in alignment:
                seq.pop(current_index + 1)
        else:
            current_index += 1
        original_index += 1

    result_alignment = {list(aln.keys())[index]: ''.join(alignment[index]) for index in range(len(alignment))}
    return removed_columns, result_alignment

def save_data(seq):
    with open("result.aln", 'w') as file:
        for k, v in seq.items():
            file.write(f">{k}\n{v}\n")

def IndelHistory(a, T):
    """parameter: a, a path of multiple alignment file
    T, a p-ath of tree file
    return a dictionary of indel history"""
    # input alignment and trees
    aln_dict, alignment = get_aln(a)
    tree, ancestor_list, v_list, order,v_dict = get_tree_data(T)
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
    print(f"columns 1/{len(list(new_alignment.values())[0])} is finished")

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

            sign_i = dist
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
    return save_data(all_sequence)

def main(args):
    print("Received the following arguments:")
    print(f"-s option: {args.s}")
    print(f"-t option: {args.t}")
    IndelHistory(args.s, args.t)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="indel")

    parser.add_argument("-s", required=True, help="Specify seqs input")
    parser.add_argument("-t", required=True, help="Specify tree input")

    args = parser.parse_args()

    main(args)

