from ete3 import Tree
from collections import defaultdict
from collections import OrderedDict
import argparse

#'./toy_data/CYP2U_165.nwk'#'./toy_data/toy_data1/input_tree.nwk'

#"./toy_data/CYP2U_165.aln" #'./toy_data/toy_data1/extants.aln'

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
    #change charcter to 1 and gap to 0
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
    order = []
    v_dict = defaultdict(list)
    number=0
    for n in tree.traverse():#traverse
        # if the v have no name, give it a name
        if n.name=="":
            n.name='N'+str(number)
            number +=1
        #save the order of notes
        order.append(n.name)
        if n.is_leaf() == False:
            ancestor_list.append(n.name)
            # if the v have no name, give it a name
            for c in n.children:
                if c.name == "":
                    c.name = 'N' + str(number)
                    number += 1
                v_dict[n.name] += [c.name]
    return ancestor_list, order, v_dict


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
    aln, alignment = get_aln(a)
    ancestor_list, order, v_dict = get_tree_data(T)
    ancestor_list.reverse()
    # remove identical adjacent columns
    removed_columns, new_alignment = remove_identical_adjacent_columns(alignment, aln)

    print("identical adjacent columns removed")

    value_dict = defaultdict(list)# save all poissble value for each point
    opt_seq_list = []# the opt seqence list
    opt_list = [] #the opt value list
    opt0_seq='' #the opt seq in position0
    opt0_value='' # the opt value in position 0
    for v in ancestor_list: #for all the v in the ancestor list find its opt
        sons = v_dict[v] # find the child of the v
        sons.reverse()# reverse it to ensure its order itws correct
        value_dict[v] = [0,0]# to set the beginning value as 0
        clo_dict = defaultdict(list)# a dict save the opt seq of each child
        son_values = defaultdict(list)# a dict save the opt value of each child
        for son in sons:# to get each child contribute to the v c0lour
            if son in new_alignment:# if its a edge, just get the inofrmation from new_alignment
                v_value = new_alignment[son][0]
                value_dict[v][0] += abs(int(0)-int(v_value))#if the v is 0, cal the sign
                value_dict[v][1] += abs(int(1) - int(v_value))# if the v is 1, cal the sign
            else:# if child is a v, add the child's value to the sign
                v0 = [0+value_dict[son][0],1+value_dict[son][1]]# if child is 0 sign is 0, if child is 1 sign is 1
                clo_dict[0].append({son:v0.index(min(v0))})# find the min child
                son_values[0].append({son:value_dict[son][v0.index(min(v0))]})# find the min child's value
                value_dict[v][0] += min(v0)# add the the v's value
                v1 = [1+value_dict[son][0],0+value_dict[son][1]]# if child is 0 sign is 1, if child is 1 sign is 0 for v is 1
                clo_dict[1].append({son:v1.index(min(v1))})
                son_values[1].append({son: value_dict[son][v1.index(min(v1))]})
                value_dict[v][1] += min(v1)
        col= clo_dict[value_dict[v].index(min(value_dict[v]))] # find the opt child by getting the min v
        son_value=son_values[value_dict[v].index(min(value_dict[v]))]# also find its value

        for i0 in col: # save the value into a string
            opt0_seq+= str(list(i0.values())[0])
        for i1 in son_value:# save the value into a string
            opt0_value += str(list(i1.values())[0])

    opt0_seq += str(value_dict[v].index(min(value_dict[v]))) #add the root seq
    opt0_value += str(min(value_dict[v]))
    opt_seq_list.append(opt0_seq)# add the seq of 0 into the list
    opt_list.append(opt0_value)
    print(f"columns 1/{len(list(new_alignment.values())[0])} is finished")

    for index in range(1, len(list(new_alignment.values())[0])):# to get the rest of colums
        value_dict = defaultdict(list)
        opt_seq=''
        opt_value=''
        for v in ancestor_list:
            sons = v_dict[v]
            value_dict[v] = [0, 0]
            clo_dict = defaultdict(list)
            son_values = defaultdict(list)
            for son in sons:
                if son in new_alignment:
                    v_value = new_alignment[son][index]
                    brother_value=abs(int(new_alignment[son][index-1])-int(opt_seq_list[index-1][ancestor_list.index(v)]))# the i-1 sign
                    value_dict[v][0] += abs(int(0) - int(v_value))*abs(abs(int(0) - int(v_value))-brother_value)#times the difference between i-1 sign and current sign
                    value_dict[v][1] += abs(int(1) - int(v_value))*abs(abs(int(1) - int(v_value))-brother_value)
                else:
                    brother_value = abs(int(opt_seq_list[index - 1][ancestor_list.index(son)]) - int(opt_seq_list[index - 1][ancestor_list.index(v)]))
                    v0 = [0*brother_value + value_dict[son][0], 1*(1-brother_value) + value_dict[son][1]]#times the difference between i-1 sign and current sign
                    clo_dict[0].append((v0.index(min(v0))))
                    son_values[0].append(value_dict[son][v0.index(min(v0))])
                    value_dict[v][0] += min(v0)
                    v1 = [1*(1-brother_value) + value_dict[son][0], 0*brother_value + value_dict[son][1]]
                    clo_dict[1].append(v1.index(min(v1)))
                    son_values[1].append(value_dict[son][v1.index(min(v1))])
                    value_dict[v][1] += min(v1)

            value_dict[v][0] += int(opt_list[index-1][ancestor_list.index(v)])# adding the pervious opt value
            value_dict[v][1] += int(opt_list[index - 1][ancestor_list.index(v)])
            col = clo_dict[value_dict[v].index(min(value_dict[v]))]
            son_value = son_values[value_dict[v].index(min(value_dict[v]))]
            for i0 in col:
                opt_seq += str(i0)
            for i1 in son_value:
                opt_value += str(i1)
        opt_seq += str(value_dict[v].index(min(value_dict[v])))
        opt_seq_list.append(opt_seq)# adding the result to the list
        opt_value += str(min(value_dict[v]))
        opt_list.append(opt_value)

        print(f"columns {index+1}/{len(list(new_alignment.values())[0])} is finished")

    # reset the opt-seq to full result
    result_history = [''.join(nums) for nums in zip(*opt_seq_list)]# join each positions seeq into one
    for position in removed_columns:# reinput he removed columns
        for i in range(len(result_history)):
            result_history[i] = (result_history[i][:position] + result_history[i][position - 1]
                                 + result_history[i][position:])

    result_history = [''.join(nums) for nums in result_history]# join the list into seq
    result_history.extend(alignment)# put the edge into the result
    ancestor_list.extend(list(aln.keys()))
    all_sequence = {k: v for k, v in zip(ancestor_list, result_history)}# dict the name and the seq
    all_sequence = OrderedDict((key, all_sequence[key]) for key in order)# re order the result
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
