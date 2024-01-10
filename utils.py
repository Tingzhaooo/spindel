from ete3 import Tree
from collections import defaultdict

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

def save_data(o,seq):
    o_path = f"{o}/result.aln"
    with open(o_path, 'w') as file:
        for k, v in seq.items():
            file.write(f">{k}\n{v}\n")
