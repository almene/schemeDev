from ete3 import Tree
from ete3 import parser
from Bio import SeqIO
import vcfpy
from argparse import (ArgumentParser, FileType)
import numpy as np
import random
import pandas as pd
import csv


def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser_inner = ArgumentParser(
        description='BIO_Hansel Scheme development')
    parser_inner.add_argument('--in_vcf', type=str, required=True, help='VCF File of SNPs')
    parser_inner.add_argument('--in_nwk', type=str, required=False, default="none", help='Newick Tree of strains')
    parser_inner.add_argument('--reference', type=str, required=True, help='Reference fasta sequence from VCF')
    parser_inner.add_argument('--outdir', type=str, required=False, help='Output Directory to put results')
    parser_inner.add_argument('--min_snps', type=int, required=False, default=2,
                              help="Number of cannonical SNPs required to define a group")
    parser_inner.add_argument('--min_members', type=int, required=False, default=5,
                              help="Minimum number of members to define a group")
    parser_inner.add_argument('--group_info', type=str, required=False, default="none",
                              help="A tsv file that contains the leaf ids and group information,"
                                   " for specific formatting information please see the README")
    return parser_inner.parse_args()


def get_subtrees(tree_mem):
    """
       break the tree membership information into which subtrees of a given size are present

       Parameters
       ----------
       tree_mem : dictionary
        Contains the group information with the leaves as keys and the list of groups as values

       Returns
       -------
       subtree_sizes : dictionary
            contains the subtree information, keys represent the size of the sub tree, values are a tuple of leaf names,
            rank_id, group_id

       """
    # obtain the maximum rank size
    keys = list(tree_mem.keys())
    max_rank = len(tree_mem[keys[0]])
    subtree_sizes = dict()
    # for each rank
    for i in range(0, max_rank):
        # blank the dictionary for the current subtree
        subtree = dict()
        # for each of the leaves in the tree membership dictionary rearange the information so that the key is the
        # group_id and the value is a list of all the leaves that are in that subtree
        for key in tree_mem.keys():
            tree_id = tree_mem[key][i]
            if tree_id not in subtree.keys():
                subtree[tree_id] = [key]
            else:
                s_tree = list(subtree[tree_id])
                s_tree.append(key)
                subtree[tree_id] = s_tree
        # add to the dictionary of subtrees with the size of the subtree as a key and the values are
        # tuples that contain a list of the leaves and the rank and group ids
        for key in subtree.keys():
            size = len(subtree[key])
            if size not in subtree_sizes.keys():
                subtree_sizes[size] = [(subtree[key], i, key)]
            else:
                temp = list(subtree_sizes[size])
                temp.append((subtree[key], i, key))
                subtree_sizes[size] = temp
    return subtree_sizes


def search_st(ingroup, potential_rank, potential_group, subtrees):
    """
       given the membership of the partition dictated by a snp, what subtree has that membership

       Parameters
       ----------
        ingroup : list
            Itemized list of leaf members that either all have the ref allele or the alt allele
        potential_rank : list
            Currently identified potential subtree ranks that might be supported by the partition
        potential_group : list
            Currently identified potential subtree group ids that might be supported by the partition
       subtrees : dictionary
            contains the subtree information, keys represent the size of the sub tree, values are a tuple of leaf names,
            rank_id, group_id
       Returns
       -------
        potential_rank : list
            Currently identified potential subtree ranks that might be supported by the partition
        potential_group : list
            Currently identified potential subtree group ids that might be supported by the partition

       """
    # extract the size of the snp sharing group
    part_size = len(ingroup)
    # get the dictionary section that corresponds to the trees of the same size as the snp sharing group if it exists
    if part_size not in subtrees.keys():
        return potential_rank, potential_group
    potential = subtrees[part_size]
    ingroup.sort()
    # for each subtree of the target size check to see if it has the same members as the snp sharing group
    for item in potential:
        from_subtrees = item[0]
        from_subtrees.sort()
        if from_subtrees == ingroup:
            potential_rank.append(item[1])
            potential_group.append(item[2])
    # sort through the potential subtrees and order based on rank id
    if len(potential_rank) > 1:
        target_index = 0
        for i in range(0, len(potential_rank)):
            if potential_rank[i] == min(potential_rank):
                target_index = i
        potential_rank = [potential_rank[target_index]]
        potential_group = [potential_group[target_index]]
    return potential_rank, potential_group


def tsv_to_membership(infile):
    groups = dict()
    with open(infile) as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        max_len = 0
        for row in reader:
            leaf = row[0]
            if len(row)>max_len:
                max_len = len(row)
            for i in range(1, len(row)):
                if row[i] == "":
                    if i == len(row)-1:
                        break
                    else:
                        print("There is an error in your tsv file")
                else:
                    if leaf not in groups.keys():
                        groups[leaf] = [row[i]]
                    else:
                        groups[leaf].append(row[i])
    # ensure that all values are the same length for later matrixing
    for key in groups.keys():
        while len(groups[key])<max_len:
            groups[key].append(0)

    return groups


def parse_tree(tree_file):
    # Load a tree structure from a newick file.
    t = Tree(tree_file)
    # Need this otherwise groups derived from the tree are inaccurate
    t.resolve_polytomy()
    # Need to force root for consistency but should modify this behavior to support defined root
    root = t.get_midpoint_outgroup()
    # t.set_outgroup(root)
    return t


def get_tree_groups(ete3_tree_obj):
    """
    :param
        ete3_tree_obj: ete3 tree object
        Takes an ete3 format tree for processing into groupings indicated by the tree structure
    :return:
        memberships: dictionary
        Contains the group information with the leaves as keys and the list of groups as values
    """
    # initialize variables
    level_rank = 0
    memberships = dict()
    group = 1
    # visit all nodes by depth from the root and generate a list of group memberships for each leaf
    # a leaf is a part of the same membership as another leaf iff they share a common interior parental node
    # the root of the tree is ignored in this calculation
    for node in ete3_tree_obj.iter_descendants("levelorder"):
        names = node.get_leaf_names()
        # print(names)
        if node.dist < 0:
            node.dist = 0
        for n in names:
            if n == "Reference":
                continue
            if n not in memberships:
                memberships[n] = list()
            memberships[n].append(group)
        level_rank += 1
        group += 1
    # obtain the number of subtrees that contain a leaf
    n_ranks = list()
    for n in memberships:
        n_ranks.append(len(memberships[n]))
    # get the maximum depth of the tree, the largest number of subtrees that contain the same leaf
    max_ranks = max(n_ranks)
    # converts the group identifiers so that the numbering restarts at each level
    for i in range(0, max_ranks):
        group = 1
        lookup = dict()
        for n in memberships:
            length = len(memberships[n])
            if i < length:
                group_id = memberships[n][i]
                if group_id not in lookup:
                    lookup[group_id] = group
                    group += 1
                memberships[n][i] = lookup[group_id]
    # if the sample is a member of a smaller set of subtrees than the sample with the largest set of member subtrees
    # append zeros to the end of the list to make it the same size as the largest set of member subtrees
    for n in memberships:
        count_levels = len(memberships[n])
        if count_levels < max_ranks:
            for i in range(count_levels, max_ranks):
                memberships[n].append(0)
    # output the membership information in a tab delineated manner
    """for n in memberships:
        print("{}\t{}".format(n, "\t".join(str(v) for v in memberships[n])))"""
    return memberships


def main():
    # collect relevant user input and parse it into the appropriate variables
    global ref_seq
    global memberships
    args = parse_args()
    nwk_treefile = args.in_nwk
    reference_fasta = args.reference
    vcf_file = args.in_vcf
    min_snps = args.min_snps
    min_group_size = args.min_members
    group_info = args.group_info
    path = ""
    # double check that the user has input only one way to get group info
    if group_info == "none" and nwk_treefile == "none":
        print("Either a Newick tree file or a tsv with group information is required to proceed.")
        raise SystemExit(0)
    if group_info != "none" and nwk_treefile != "none":
        print("Either a Newick tree file or a tsv with group information is required to proceed. "
              "Please only use one of the two options")
        raise SystemExit(0)
    if group_info == "none":
        path = "tree"
    if nwk_treefile == "none":
        path = "groups"
    # hard coded size of positive and negative tiles surrounding SNP
    tile_size = 32
    tile_size = tile_size - 1
    mid_point = int(tile_size / 2)
    leaves = []
    # required spacing between at least two snps in the group support
    # size_cutoff = 200000
    # test that files exist
    # obtain the dictionary of leaf subtree membership after pre-processing the raw Newick Tree to resolve polytomys
    if path == "tree":
        try:
            tree = parse_tree(nwk_treefile)
            memberships = get_tree_groups(tree)
            for leaf in tree.traverse():
                if leaf.is_leaf():
                    leaves.append(leaf.name)
        except parser.newick.NewickError:
            print(
                f"There was an error in reading {nwk_treefile}.  "
                f"Verify that the file exists and is in the correct format")
            raise SystemExit(0)
    # if group information is provided as a tsv extract membership information
    if path == "groups":
        memberships = tsv_to_membership(group_info)
        leaves = list(memberships.keys())

    # process the memberships data to extract the subtree information from the membership information
    subtrees = get_subtrees(memberships)
    # read the reference fasta file by record
    with open(reference_fasta, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if not fasta:
            print("That is not a .fasta file")
            raise SystemExit(0)
    try:
        for seq_record in SeqIO.parse(reference_fasta, "fasta"):
            ref_seq = "{}".format(seq_record.seq)
            ref_length = len(ref_seq)
            if ref_length == 0:
                print("error in .fasta file")
                raise SystemExit(0)
    except FileNotFoundError:
        print(
            f"There was an error in reading {reference_fasta}.  "
            f"Verify that the file exists and is in the correct format")
        raise SystemExit(0)

    # collect the vcf file indicated in the arguments passed to the script
    try:
        reader = vcfpy.Reader.from_path(vcf_file)
    except FileNotFoundError:
        print(
            f"There was an error in reading {vcf_file}.  "
            f"Verify that the file exists")
        raise SystemExit(0)
    except vcfpy.exceptions.IncorrectVCFFormat:
        print(
            f"There was an error in reading {vcf_file}.  "
            f"Verify that the file is in the correct format")
        raise SystemExit(0)

    if min_snps < 1:
        print("Minimum SNP support for a group needs to be an integer greater than zero")
        raise SystemExit(0)
    if min_group_size < 1:
        print("Minimum group size needs to be an integer greater than zero")
        raise SystemExit(0)

    # Hardcoded header
    header = ['#CHROM', 'POS', 'REF', 'ALT'] + reader.header.samples.names
    samples = reader.header.samples.names
    discr = np.setdiff1d(samples, leaves)
    if len(discr) != 0:
        print("Something went wrong.  The leaf names in the provided Newick file or tsv file "
              "are not the same as the ones in the VCF file")
        print(f"The samples found in one of the files but not the other are : {discr}")
        raise SystemExit(0)
    scheme = dict()
    number_snps = 0
    number_groups = 0
    # print(subtrees)
    all_variable = []
    for record in reader:
        # skip over the metadata lines at the top of the file
        if not record.is_snv():
            continue

        # pull the positional information and the reference value and all the alt values from the vcf file
        line = [record.CHROM, record.POS, record.REF]
        line += [alt.value for alt in record.ALT]

        # skip any vcf calls that have multiple alternative alleles
        if len(record.ALT) > 1:
            continue

        # initialize variables and pull the ref and alt SNPs
        alt_base = record.ALT[0].value
        ref_base = record.REF
        position = record.POS
        all_variable.append(position)
        count_ref = 0
        count_alt = 0
        ref_group = list()
        alt_group = list()
        tracker = 0

        # go through the samples and divide them by alt vs ref bases
        for call in record.calls:
            state = int(call.data.get('GT'))
            sample_name = samples[tracker]
            if state == 0:
                count_ref += 1
                ref_group.append(sample_name)
            else:
                count_alt += 1
                alt_group.append(sample_name)
            tracker += 1

        valid_ranks = list()
        valid_groups = list()
        conflict_positions = dict()
        ref_candidate_partitions = dict()
        alt_candidate_partitions = dict()
        if count_ref >= 1 and count_alt >= 1:
            if len(ref_group) in subtrees.keys():
                valid_ranks, valid_groups = search_st(ref_group, valid_ranks, valid_groups, subtrees)
            if len(valid_ranks) > 0:
                for i in range(0, len(valid_ranks)):
                    if valid_ranks[i] not in ref_candidate_partitions:
                        ref_candidate_partitions[valid_ranks[i]] = list()
                    if valid_groups[i] == 0:
                        continue
                    for item in ref_candidate_partitions[valid_ranks[i]]:
                        if item["position"] == position:
                            """print(f"Warning duplication of {position}")"""
                    ref_candidate_partitions[valid_ranks[i]].append({"position": position, "ref_base": ref_base,
                                                                     "alt_base": alt_base,
                                                                     "ref_count": count_ref,
                                                                     "alt_count": count_alt,
                                                                     "positive_group": "ref",
                                                                     "group_id": valid_groups[i],
                                                                     "rank_id": valid_ranks[i],
                                                                     "qc_warnings": []})
            valid_ranks, valid_groups = search_st(alt_group, valid_ranks, valid_groups, subtrees)
            if len(valid_ranks) > 0:
                for i in range(0, len(valid_ranks)):
                    if valid_ranks[i] not in alt_candidate_partitions:
                        alt_candidate_partitions[valid_ranks[i]] = list()
                    if valid_groups[i] == 0:
                        continue
                    # get the tile start from the position and tile information
                    start = position - 1 - mid_point
                    alt_candidate_partitions[valid_ranks[i]].append({"position": position,
                                                                     "ref_base": ref_base,
                                                                     "alt_base": alt_base,
                                                                     "ref_count": count_ref,
                                                                     "alt_count": count_alt,
                                                                     "positive_group": "alt",
                                                                     "group_id": valid_groups[i],
                                                                     "rank_id": valid_ranks[i],
                                                                     "qc_warnings": []})
            positions = dict()

            for rank in ref_candidate_partitions:
                if rank not in scheme:
                    scheme[rank] = dict()
                for item in ref_candidate_partitions[rank]:

                    # filter out groups which are too small
                    if item['ref_count'] < min_group_size:
                        continue
                    if item['alt_count'] < min_group_size:
                        continue
                    if not item["group_id"] in scheme[rank]:
                        scheme[rank][item["group_id"]] = []
                    scheme[rank][item["group_id"]].append(item)

            for rank in alt_candidate_partitions:
                if rank not in scheme:
                    scheme[rank] = dict()
                for item in alt_candidate_partitions[rank]:
                    # filter out groups which are too small
                    if item['alt_count'] < min_group_size:
                        continue
                    if not item["group_id"] in scheme[rank]:
                        scheme[rank][item["group_id"]] = []
                    scheme[rank][item["group_id"]].append(item)

            """ for rank in ref_candidate_partitions:
                for item in ref_candidate_partitions[rank]:
                    pos = item["position"]
                    # filter out groups which are too small
                    if item['ref_count'] < min_group_size:
                        continue
                    if positions[pos] == 1:
                        print("{}\t{}".format(item, "REF-Single"))
                    else:
                        print("{}\t{}".format(item, "REF-Multi"))

            for rank in alt_candidate_partitions:
                for item in alt_candidate_partitions[rank]:
                    pos = item["position"]
                    # filter out groups which are too small
                    if item['alt_count'] < min_group_size:
                        continue
                    if positions[pos] == 1:
                        # print("Success")
                        print("{}\t{}".format(item, "ALT-Single"))
                    else:
                        print("{}\t{}".format(item, "ALT-Multi"))"""

    # filter out groups with less than minimum number of supporting snps
    required_tiles = []
    for rank in scheme:
        valid = dict()
        prev = []
        for group_id in scheme[rank]:
            if len(scheme[rank][group_id]) >= min_snps:
                valid[group_id] = scheme[rank][group_id]
                number_groups += 1
                duplicates = []
                for item in range(0, len(scheme[rank][group_id])):
                    if not prev:
                        prev = scheme[rank][group_id][item]["position"]
                        required_tiles.append(scheme[rank][group_id][item]["position"])
                        number_snps += 1
                    else:
                        if prev == scheme[rank][group_id][item]["position"]:
                            duplicates.append(item)
                            continue
                        else:
                            prev = scheme[rank][group_id][item]["position"]
                            required_tiles.append(scheme[rank][group_id][item]["position"])
                            number_snps += 1
                duplicates.sort(reverse=True)
                for i in range(0, len(duplicates)):
                    scheme[rank][group_id].pop(i)
        scheme[rank] = valid

    """ print(f"\n When the minimum group size is {min_group_size}, there are {number_snps} SNPs"
          f" that support a total of {number_groups} groups")"""
    # sort the lists of both all the variable positions and the ones of interest for the tile selection
    all_variable.sort()
    required_tiles.sort()
    # to reduce calculations store the location of the last position checked for conflict
    stored_place = 0
    for item in required_tiles:
        if stored_place < mid_point + 1:
            start = 0
        elif item - all_variable[stored_place] > mid_point:
            start = stored_place
        else:
            start = stored_place - mid_point
        for i in range(start, len(all_variable)):
            difference = item - all_variable[i]
            if difference > mid_point:
                continue
            elif difference < (-1 * mid_point):
                stored_place = i
                break
            elif mid_point >= difference >= (-1 * mid_point):
                if difference == 0:
                    continue
                if item not in conflict_positions.keys():
                    conflict_positions[item] = [all_variable[i]]
                else:
                    conflict_positions[item].append(all_variable[i])
            else:
                print(f"logic error\t : {difference}")
    # using the conflict position information pull the appropriate tiles from the reference genome
    for rank in scheme:
        for group_id in scheme[rank]:
            # pull the tiles without degenerate bases
            for item in range(0, len(scheme[rank][group_id])):
                position = scheme[rank][group_id][item]["position"]
                start = position - 1 - mid_point
                # which is the positive tile is determined by if the ref or alt base defines the in-group
                if scheme[rank][group_id][item]["positive_group"] == "ref":
                    scheme[rank][group_id][item]["positive_tile"] = ref_seq[start:position - 1] + \
                                                                    scheme[rank][group_id][item]["ref_base"] + \
                                                                    ref_seq[position:position + mid_point]
                    scheme[rank][group_id][item]["negative_tile"] = ref_seq[start:position - 1] + \
                                                                    scheme[rank][group_id][item]["alt_base"] + \
                                                                    ref_seq[position:position + mid_point]
                elif scheme[rank][group_id][item]["positive_group"] == "alt":
                    scheme[rank][group_id][item]["positive_tile"] = ref_seq[start:position - 1] + \
                                                                    scheme[rank][group_id][item]["alt_base"] + \
                                                                    ref_seq[position:position + mid_point]
                    scheme[rank][group_id][item]["negative_tile"] = ref_seq[start:position - 1] + \
                                                                    scheme[rank][group_id][item]["ref_base"] + \
                                                                    ref_seq[position:position + mid_point]
                else:
                    print("something went wrong")
                # if position has conflicts add degenerate bases
                if position in conflict_positions.keys():
                    for conflict in conflict_positions[position]:
                        index = conflict - start - 1

                        p_substring = scheme[rank][group_id][item]["positive_tile"][:index + 1]
                        n_substring = scheme[rank][group_id][item]["negative_tile"][:index + 1]
                        p_replace = scheme[rank][group_id][item]["positive_tile"][:index] + "N"
                        n_replace = (scheme[rank][group_id][item]["negative_tile"][:index] + "N")

                        scheme[rank][group_id][item]["positive_tile"] = scheme[rank][group_id][item]["positive_tile"]. \
                            replace(p_substring, p_replace, 1)
                        scheme[rank][group_id][item]["negative_tile"] = scheme[rank][group_id][item]["negative_tile"]. \
                            replace(n_substring, n_replace, 1)
                    if len(conflict_positions[position]) >= 1:
                        scheme[rank][group_id][item]["qc_warnings"].append("One or more degenerate bases")
    # double check that degenerate bases haven't removed support for a group
    # and that all support snps are not from the same region
    for rank in scheme:
        for group_id in scheme[rank]:
            good_snps = len(scheme[rank][group_id])
            positions = []
            for item in range(0, len(scheme[rank][group_id])):
                positions.append(scheme[rank][group_id][item]["position"])
                if len(scheme[rank][group_id][item]["qc_warnings"]) > 0:
                    good_snps -= 1
            if good_snps < min_snps:
                for item in range(0, len(scheme[rank][group_id])):
                    scheme[rank][group_id][item]["qc_warnings"].append("Lacks good tile support")
            """if max(positions) - min(positions) < size_cutoff:
                for item in range(0, len(scheme[rank][group_id])):
                    scheme[rank][group_id][item]["qc_warnings"].append("All support snps in the same genetic region")"""
    stripped_pos = []
    for rank in scheme:
        for group_id in scheme[rank]:
            for item in range(0, len(scheme[rank][group_id])):
                if len(scheme[rank][group_id][item]["qc_warnings"]) == 0:

                     stripped_pos.append(scheme[rank][group_id][item]["position"])
    """# write out stripped vcf file
    reader = vcfpy.Reader.from_path(vcf_file)
    writer = vcfpy.Writer.from_path('stripped_core.vcf', reader.header)
    for record in reader:
        print("I'm doing something!")
        if not record.is_snv():
            continue
        line = [record.CHROM, record.POS, record.REF]
        line += [alt.value for alt in record.ALT]
        position = record.POS
        print(position)
        if position in stripped_pos:
            writer.write_record(record)"""

    # print out scheme information
    """for rank in scheme:
        for group_id in scheme[rank]:
            # pull the tiles without degenerate bases
            for item in range(0, len(scheme[rank][group_id])):
                print("{}\t{}".format(scheme[rank][group_id][item], f"with group restrictions of {min_group_size}"
                                                                    f" members and {min_snps} snps"))"""
    # generate matrix representation to store masked rank, ID info
    # length of row
    rama = len(memberships[random.choice(list(memberships))])
    # length of column
    grouper = len(memberships)
    codes_wip = np.full((grouper, rama), 0)

    # mask unsupported groups with zeros
    for n in memberships:
        hierarchy = memberships[n]
        for i in range(0, len(hierarchy)):
            if i not in scheme:
                hierarchy[i] = 0
            elif not hierarchy[i] in scheme[i]:
                hierarchy[i] = 0
        # print("{}\t{}".format(n, "\t".join(str(v) for v in memberships[n])))
    # mask unsupported groups with zeros
    test = pd.DataFrame(memberships, dtype=object)
    test = test.T

    n_unique = test.apply(pd.Series.nunique)
    cols_to_drop = n_unique[n_unique == 1].index
    test.drop(cols_to_drop, axis=1)
    dropped = test.drop(cols_to_drop, axis=1)

    # make a vector with all the number of groups in each column
    n_unique = dropped.apply(pd.Series.nunique)
    current_codes = dict()
    col_groups = list(n_unique)
    biohansel_codes = dropped
    new_code = str()
    for i in range(1, biohansel_codes.shape[1]):
        current_max = 1
        used = dict()
        current_ids = dict()
        for j in range(0,col_groups[i]):
            for k in range(0,biohansel_codes.shape[0]):
                if dropped.values[k,i] not in current_ids.keys():
                    current_ids[dropped.values[k, i]] = [k]
                else:
                    current_ids[dropped.values[k, i]].append(k)
                    # this stores both the current group ids under consideration but also a list of all the rows that
                    # have that id
        else:
            for k in range(0, biohansel_codes.shape[0]): # go to every row
                # cycle through ids, skipping zero an identify which is connected to current line
                working_group = -1
                new_code=str()
                for identifier in current_ids.keys():
                    if k in current_ids[identifier]:
                        working_group = identifier
                # if that row is associated with group zero at column i just copy from left
                if working_group == 0:
                    biohansel_codes.values[k, i] = biohansel_codes.values[k, i - 1]
                else:
                    # if the group has already been used just pull the information
                    if working_group in used.keys():
                        biohansel_codes.values[k,i] = used[working_group]
                    else:
                        # if the value in the previous row is a zero start a new clade from the root
                        if biohansel_codes.values[k, i - 1] == 0:
                            for attempt in range(1, biohansel_codes.shape[0]):
                                if attempt in current_codes.keys():
                                    continue
                                else:
                                    new_code = str(attempt)
                                    break
                        # if there was a non-zero value in the previous column
                        else:
                            if biohansel_codes.values[k, i-1] in current_codes.keys():
                                current_max = max(current_codes[biohansel_codes.values[k, i-1]])+1
                            if biohansel_codes.values[k, i-1] not in current_codes.keys():
                                current_codes[biohansel_codes.values[k, i - 1]] = []
                            current_codes[biohansel_codes.values[k, i - 1]].append(current_max)
                            new_code = str(biohansel_codes.values[k, i - 1])+"."+str(current_max)
                        biohansel_codes.values[k,i] = new_code
                        used[working_group] = new_code
    """with pd.option_context('display.max_rows', None, 'display.max_columns',
                           None):  # more options can be specified also
        print(test)
        print(dropped)
        print(biohansel_codes)"""
    # for translating old column numbers to the new ones in the new matrix with dropped columns get the old indexs
    # and the indexes of the indexes is the new column number
    kept = n_unique[n_unique > 1].index
    rank_id = int()
    row_id = int()
    # print(biohansel_codes.shape)
    log = open(f"{vcf_file}_biohansel.log", "w+")
    fasta_file = open(f"{vcf_file}_biohansel.fasta", "w+")
    for rank in scheme:
        # translate the rank into the position in the table that has all the undefined columns are dropped
        for i in range(0, dropped.shape[1]):
            if int(kept[i]) == int(rank):
                rank_id = i
                # print(f"found : {int(kept[i])} equals {int(rank)}")
                break
            # else:
                # print(f"error finding {rank}")
                # print(f"{int(kept[i])} does not equal {int(rank)}")

        # print(rank_id)
        for group_id in scheme[rank]:
            if path == "tree":
                # get biohansel code for that rank, group combo
                for row in range(0, dropped.shape[0]):
                    if str(test.values[row, rank]) == str(group_id):
                        row_id = row
                        break
                    """else:
                        print(f"{test.values[row, rank]} does not equal {int(group_id)}")"""
                code = biohansel_codes.values[row_id, rank_id]
            # if the information from the groups in already provided use the old group codes
            if path == "groups":
                code = group_id

            for item in range(0, len(scheme[rank][group_id])):
                # exclude the the entries that have qc problems
                if len(scheme[rank][group_id][item]["qc_warnings"]) != 0:
                    continue
                # extract the required information from the entries that have good qc
                position = scheme[rank][group_id][item]["position"]
                pos_tile = scheme[rank][group_id][item]["positive_tile"]
                neg_tile = scheme[rank][group_id][item]["negative_tile"]

                # write out the biohansel fasta files and the accompanying logs
                log.write(f"{position}\t{code}\t{pos_tile}\n")
                fasta_file.write(f">{position}|{code}\n{pos_tile}\n")
                log.write(f"negative{position}\t{code}\t{neg_tile}\n")
                fasta_file.write(f">negative{position}|{code}\n{neg_tile}\n")
    log.close()
    fasta_file.close()


# call main function
if __name__ == '__main__':
    main()
