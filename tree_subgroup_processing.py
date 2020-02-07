"""
Biohansel requres a strain specific fasta format scheme to be able to classify samples.
Schemes exist for some common clonal pathogens but there are many untapped applications.
This module includes functions to aid in the generation of k-mers linked to predefined groups
to serve as the basis for automatic scheme development

"""
from argparse import ArgumentParser
import csv
import re
from ete3 import Tree
from ete3 import parser
from Bio import SeqIO
import vcfpy
import numpy as np
import pandas as pd


def parse_args():
    """
    Parse the input arguments, use '-h' for help

        Returns
        -------
            parser_inner.parse_args()

    """
    parser_inner = ArgumentParser(
        description='BIO_Hansel Scheme development')
    parser_inner.add_argument('--in_vcf', type=str, required=True, help='VCF File of SNPs')
    parser_inner.add_argument('--in_nwk', type=str, required=False, default="none",
                              help='Newick Tree of strains')
    parser_inner.add_argument('--reference', type=str, required=True,
                              help='Reference fasta sequence from VCF')
    parser_inner.add_argument('--outdir', type=str, required=False,
                              help='Output Directory to put results')
    parser_inner.add_argument('--min_snps', type=int, required=False, default=2,
                              help="Number of cannonical SNPs required to define a group")
    parser_inner.add_argument('--min_members', type=int, required=False, default=5,
                              help="Minimum number of members to define a group")
    parser_inner.add_argument('--min_parent', type=int, required=False, default=2,
                              help="Minimum size difference between new "
                                   "group and parent group to be considered valid")
    parser_inner.add_argument('--group_info', type=str, required=False, default="none",
                              help="A tsv file that contains the leaf ids and group information,"
                                   " for specific formatting information please see the README")
    return parser_inner.parse_args()


def get_subtrees(tree_mem):
    """
       break the tree membership information into which subtrees of a given size are present

       preparatory step for saving computational steps when
       comparing only group memberships of the same size

       Parameters
       ----------
       tree_mem : dictionary
        Contains the group information with the leaves as keys and the list of groups as values

       Returns
       -------
       subtree_sizes : dictionary
            contains the subtree information, keys represent the size of the sub tree,
            values are a tuple of leaf names, rank_id, g_id

       """
    # obtain the maximum rank size
    keys = list(tree_mem.keys())
    max_rank = len(tree_mem[keys[0]])
    subtree_sizes = dict()
    # for each rank
    for i in range(0, max_rank):
        # blank the dictionary for the current subtree
        subtree = dict()
        # for each of the leaves in the tree membership dictionary rearrange the
        # information so that the key is the g_id and the value is a list of \
        # all the leaves that are in that subtree
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


def search_st(ingroup: list, potential_rank: list, potential_group: list, subtrees):
    """
       given the membership of the partition dictated by a snp, what subtree has that membership

       Parameters
       ----------
        ingroup : list
            Itemized list of leaf members that either all have the ref allele or the alt allele
        potential_rank : list
            Currently identified potential subtree ranks that might be supported by the partition
        potential_group : list
            Currently identified potential subtree group ids that might be supported
            by the partition
       subtrees : dictionary
            contains the subtree information, keys represent the size of the sub tree,
            values are a tuple of leaf names, rank_id, g_id
       Returns
       -------
        potential_rank : list
            Currently identified potential subtree ranks that might be supported by the partition
        potential_group : list
            Currently identified potential subtree group ids that might
            be supported by the partition

       """
    # extract the size of the snp sharing group
    part_size = len(ingroup)
    # get the dictionary section that corresponds to the trees of the same size
    # as the snp sharing group if it exists
    if part_size not in subtrees.keys():
        return potential_rank, potential_group
    potential = subtrees[part_size]
    ingroup.sort()
    # for each subtree of the target size check to see if it has
    # the same members as the snp sharing group
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
    """
        obtain the hierarchical divisions that are present in a given
        tsv file for later functions to assess
        Parameters
        ----------
        infile : string
            indicates the file name of the group membership information
        Returns
        -------
        groups : dictionary
            Contains the group information with the leaves as keys and the list of groups as values

    """
    groups = dict()
    used = dict()
    try:
        with open(infile) as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            max_len = 0
            for row in reader:
                leaf = row[0]
                if len(row) > max_len:
                    max_len = len(row)
                for i in range(1, len(row)):
                    if row[i] == "":
                        if i == len(row) - 1:
                            break
                        print("There is an error in your tsv file "
                              "as there are empty internal columns")
                        raise SystemExit(0)
                    if leaf not in groups.keys():
                        groups[leaf] = [row[i]]
                    else:
                        groups[leaf].append(row[i])
                    if i not in used.keys():
                        used[i] = [row[i]]
                    else:
                        used[i].append(row[i])
    except FileNotFoundError:
        print("There was an error reading the tsv file.  Check the file name and try again")
        raise SystemExit(0)
    # ensure that all values are the same length for later matrix creation
    for key in groups.keys():
        while len(groups[key]) < max_len:
            groups[key].append(0)
    # ensure that all entries are numbers and replace to unique numbers if not already
    non_valid = r'[^0-9/.]'
    for i in range(0, max_len):
        switch = dict()
        new = str()
        for key in groups.keys():
            checking = str(groups[key][i])
            if i > 0:
                checking = str(groups[key][i - 1]) + "." + checking
            # check to see if the value contains only the accepted characters
            if re.search(non_valid, checking):
                if checking not in switch.keys():
                    if i == 0:
                        for j in range(1, len(groups.keys())):
                            if str(j) in used[i + 1]:
                                continue
                            new = str(j)
                            break
                        switch[checking] = new
                    else:
                        for j in range(1, len(groups.keys())):
                            new = str(groups[key][i - 1]) + "." + str(j)
                            if new in used[i + 1]:
                                continue
                            break
                        switch[checking] = new
                    used[i + 1].append(switch[checking])
                groups[key][i] = switch[checking]
    return groups


def parse_tree(tree_file):
    """
        resolves potential polytomys in the supplied newick tree
        Parameters
        ------
        tree_file: str
            user supplied tree filename
        Returns
        ------
        tree: ete3 tree object
            contains modified user tree
    """

    # Load a tree structure from a newick file.
    tree = Tree(tree_file)
    # Need this otherwise groups derived from the tree are inaccurate
    tree.resolve_polytomy()
    return tree


def get_tree_groups(ete3_tree_obj):
    """
        obtain the hierarchical divisions that are present in a
        given tree for later functions to assess
        Parameters
        ----------
            ete3_tree_obj: ete3 tree object
                Takes an ete3 format tree for processing into groupings
                 indicated by the tree structure
        Returns
        -------
            memberships: dictionary
                Contains the group information with the leaves as keys
                 and the list of groups as values
    """
    # initialize variables
    memberships = dict()
    group = 1
    # visit all nodes by depth from the root and generate a list of group memberships for each leaf
    # a leaf is a part of the same membership as another leaf if they share a common
    # interior parental node the root of the tree is ignored in this calculation
    for node in ete3_tree_obj.iter_descendants("levelorder"):
        names = node.get_leaf_names()
        # print(names)
        if node.dist < 0:
            node.dist = 0
        for sample in names:
            if sample == "Reference":
                continue
            if sample not in memberships:
                memberships[sample] = list()
            memberships[sample].append(group)
        group += 1
    # obtain the number of subtrees that contain a leaf
    n_ranks = list()
    for sample in memberships:
        n_ranks.append(len(memberships[sample]))
    # get the maximum depth of the tree, the largest number of subtrees that contain the same leaf
    max_ranks = max(n_ranks)
    # converts the group identifiers so that the numbering restarts at each level
    for i in range(0, max_ranks):
        # group numbering starts at one as zero groups are
        # placeholders and do not signify a real group
        group = 1
        lookup = dict()
        for sample in memberships:
            length = len(memberships[sample])
            if i < length:
                # obtain the old group numbering from the memberships dictionary
                g_id = memberships[sample][i]
                if g_id not in lookup:
                    # generate an entry in the lookup dictionary between
                    # the old and the new group numbers
                    lookup[g_id] = group
                    group += 1
                # replace the old group name with the renumbered group name
                # in the memberships dictionary
                memberships[sample][i] = lookup[g_id]
    # if the sample is a member of a smaller set of subtrees than the sample with
    # the largest set of member subtrees append zeros to the end of the list to make
    # it the same size as the largest set of member subtrees
    for sample in memberships:
        count_levels = len(memberships[sample])
        if count_levels < max_ranks:
            for i in range(count_levels, max_ranks):
                memberships[sample].append(0)
    return memberships


def add_tiles(scheme, ref_seq, mid_point, conflict_positions):
    """
    using the information from the vcf file and the requerence sequence generate the positive and
    negative k-mers and add them to the scheme object
    Parameters
    ----------
    scheme: dict
        storage location for the information needed to generate the biohansel scheme output
    ref_seq: seq
        contains the user supplied reference sequence
    mid_point: int
        user defined number of bases flanking the SNP site
    conflict_positions: dict
        information regarding ambiguous sites in the genome

    Returns
    -------
    scheme: dict
        storage location for the information needed to generate the biohansel scheme output
    """
    for rank in scheme:
        for g_id in scheme[rank]:
            # pull the tiles without degenerate bases
            for item in range(0, len(scheme[rank][g_id])):
                pos_tile = ""
                neg_tile = ""
                position = scheme[rank][g_id][item]["position"]
                start = position - 1 - mid_point
                # which is the positive tile is determined by if the ref or alt base defines
                # the in-group
                if scheme[rank][g_id][item]["positive_group"] == "ref":
                    pos_tile = ref_seq[start:position - 1] + scheme[rank][g_id][item]["ref_base"] +\
                        ref_seq[position:position + mid_point]
                    neg_tile = ref_seq[start:position - 1] + scheme[rank][g_id][item]["alt_base"] +\
                        ref_seq[position:position + mid_point]
                elif scheme[rank][g_id][item]["positive_group"] == "alt":
                    pos_tile = ref_seq[start:position - 1] + scheme[rank][g_id][item]["alt_base"] +\
                        ref_seq[position:position + mid_point]
                    neg_tile = ref_seq[start:position - 1] + scheme[rank][g_id][item]["ref_base"] +\
                        ref_seq[position:position + mid_point]
                else:
                    print("something went wrong")

                scheme[rank][g_id][item]["positive_tile"] = pos_tile
                scheme[rank][g_id][item]["negative_tile"] = neg_tile
                # if position has conflicts add degenerate bases
                if position in conflict_positions.keys():
                    for conflict in conflict_positions[position]:
                        index = conflict - start - 1

                        p_substring = scheme[rank][g_id][item]["positive_tile"][:index + 1]
                        n_substring = scheme[rank][g_id][item]["negative_tile"][:index + 1]
                        p_replace = scheme[rank][g_id][item]["positive_tile"][:index] + "N"
                        n_replace = (scheme[rank][g_id][item]["negative_tile"][:index] + "N")

                        scheme[rank][g_id][item]["positive_tile"] = \
                            scheme[rank][g_id][item]["positive_tile"]. \
                            replace(p_substring, p_replace, 1)
                        scheme[rank][g_id][item]["negative_tile"] = \
                            scheme[rank][g_id][item]["negative_tile"]. \
                            replace(n_substring, n_replace, 1)
                    if len(conflict_positions[position]) >= 1:
                        scheme[rank][g_id][item]["qc_warnings"].append(
                            "One or more degenerate bases")
    return scheme


def tile_generator(
        nwk_treefile, reference_fasta, vcf_file, min_snps: int, min_group_size: int,
        group_info, min_parent_size):
    """
        main function for generating potential biohansel tiles from group and
        varient call information

        Parameters
        ----------
            nwk_treefile, reference_fasta, vcf_file, group_info: str
                user supplied input file names obtained from the parse_args function
            min_snps, min_group_size, min_parent_size : int
                user supplied configuration settings from the parse_args function

        Raises
        ----------

    """
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
    else:
        path = "groups"
    # hard coded size of positive and negative tiles surrounding SNP
    tile_size = 32
    mid_point = int(tile_size / 2)
    leaves = []
    # test that files exist
    # obtain the dictionary of leaf subtree membership after pre-processing
    # the raw Newick Tree to resolve polytomys
    memberships = dict()
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
                f"Verify that the file exists and is in the correct"
                f" format for the ete3 python module")
            raise SystemExit(0)
    # if group information is provided as a tsv extract membership information
    if path == "groups":
        memberships = tsv_to_membership(group_info)
        leaves = list(memberships.keys())

    # process the memberships data to extract the subtree
    # information from the membership information
    subtrees = get_subtrees(memberships)
    # read the reference fasta file by record
    ref_seq = ""
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
            f"There was an error in reading {vcf_file} with vcfpy.  "
            f"Verify that the file exists")
        raise SystemExit(0)

    if min_snps < 1:
        print("Minimum SNP support for a group needs to be an integer greater than zero")
        raise SystemExit(0)
    if min_group_size < 1:
        print("Minimum group size needs to be an integer greater than zero")
        raise SystemExit(0)
    samples = []
    try:
        samples = reader.header.samples.names
    except vcfpy.exceptions.IncorrectVCFFormat:
        print(
            f"There was an error in reading {vcf_file}.  "
            f"Verify that the file format is correct")
        raise SystemExit(0)
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

        # pull the positional information and the reference
        # value and all the alt values from the vcf file
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
            elif state == 1:
                count_alt += 1
                alt_group.append(sample_name)
            else:
                print("there is a problem reading the state information")
            tracker += 1

        # generate places to put group and rank information
        valid_ranks = list()
        valid_ranks_alt = list()
        valid_groups = list()
        valid_groups_alt = list()
        ref_candidate_partitions = dict()
        alt_candidate_partitions = dict()

        # if there is variation at the snp location search to see
        # if the partitioning suggested by the snp is
        # one that is present in the provided tree
        if count_ref >= 1 and count_alt >= 1:
            if len(ref_group) in subtrees.keys():
                valid_ranks, valid_groups = search_st(ref_group, valid_ranks, valid_groups,
                                                      subtrees)
            if len(valid_ranks) < 1:
                continue
            for i in range(0, len(valid_ranks)):
                if valid_ranks[i] not in ref_candidate_partitions:
                    ref_candidate_partitions[valid_ranks[i]] = list()
                if valid_groups[i] == 0:
                    continue
                for item in ref_candidate_partitions[valid_ranks[i]]:
                    if item["position"] == position:
                        """print(f"Warning duplication of {position}")"""
                # if the partition is already supported add the information about that
                # snp if not created it
                ref_candidate_partitions[valid_ranks[i]].append({"position": position,
                                                                 "ref_base": ref_base,
                                                                 "alt_base": alt_base,
                                                                 "ref_count": count_ref,
                                                                 "alt_count": count_alt,
                                                                 "positive_group": "ref",
                                                                 "g_id": valid_groups[i],
                                                                 "rank_id": valid_ranks[i],
                                                                 "qc_warnings": []})
            if len(alt_group) in subtrees.keys():
                valid_ranks_alt, valid_groups_alt = search_st(alt_group, valid_ranks, valid_groups,
                                                              subtrees)
            if len(valid_ranks_alt) > 0:
                for i in range(0, len(valid_ranks_alt)):
                    if valid_ranks_alt[i] not in alt_candidate_partitions:
                        alt_candidate_partitions[valid_ranks_alt[i]] = list()
                    if valid_groups_alt[i] == 0:
                        continue
                    # if the partition is already supported add the information about that snp if
                    # not created it
                    alt_candidate_partitions[valid_ranks_alt[i]].append({"position": position,
                                                                         "ref_base": ref_base,
                                                                         "alt_base": alt_base,
                                                                         "ref_count": count_ref,
                                                                         "alt_count": count_alt,
                                                                         "positive_group": "alt",
                                                                         "g_id":
                                                                             valid_groups_alt[i],
                                                                         "rank_id": valid_ranks_alt[
                                                                             i],
                                                                         "qc_warnings": []})
            for rank in ref_candidate_partitions:
                if rank not in scheme:
                    scheme[rank] = dict()
                for item in ref_candidate_partitions[rank]:

                    # filter out groups which are too small
                    if item['ref_count'] < min_group_size:
                        continue
                    if item['alt_count'] < min_group_size:
                        continue
                    if not item["g_id"] in scheme[rank]:
                        scheme[rank][item["g_id"]] = []
                    scheme[rank][item["g_id"]].append(item)

            for rank in alt_candidate_partitions:
                if rank not in scheme:
                    scheme[rank] = dict()
                for item in alt_candidate_partitions[rank]:
                    # filter out groups which are too small
                    if item['alt_count'] < min_group_size:
                        continue
                    if not item["g_id"] in scheme[rank]:
                        scheme[rank][item["g_id"]] = []
                    scheme[rank][item["g_id"]].append(item)

    # filter out groups with less than minimum number of supporting snps
    required_tiles = []
    for rank in scheme:
        valid = dict()
        prev = []
        for g_id in scheme[rank]:
            if len(scheme[rank][g_id]) >= min_snps:
                valid[g_id] = scheme[rank][g_id]
                number_groups += 1
                duplicates = []
                for item in range(0, len(scheme[rank][g_id])):
                    if not prev:
                        prev = scheme[rank][g_id][item]["position"]
                        required_tiles.append(scheme[rank][g_id][item]["position"])
                        number_snps += 1
                    else:
                        if prev == scheme[rank][g_id][item]["position"]:
                            duplicates.append(item)
                            continue
                        prev = scheme[rank][g_id][item]["position"]
                        required_tiles.append(scheme[rank][g_id][item]["position"])
                        number_snps += 1
                duplicates.sort(reverse=True)
                for i in range(0, len(duplicates)):
                    scheme[rank][g_id].pop(i)
        scheme[rank] = valid

    # sort the lists of both all the variable positions and the ones of
    # interest for the tile selection
    all_variable.sort()
    required_tiles.sort()
    # to reduce calculations store the location of the last position checked for conflict
    stored_place = 0
    conflict_positions = dict()
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
            if difference < (-1 * mid_point):
                stored_place = i
                break
            if mid_point >= difference >= (-1 * mid_point):
                if difference == 0:
                    continue
                if item not in conflict_positions.keys():
                    conflict_positions[item] = [all_variable[i]]
                else:
                    conflict_positions[item].append(all_variable[i])
            else:
                print(f"logic error\t : {difference}")
    # using the conflict position information pull the appropriate tiles from the reference genome
    scheme = add_tiles(scheme, ref_seq, mid_point, conflict_positions)
    # double check that degenerate bases haven't removed support for a group
    # and that all support snps are not from the same region
    for rank in scheme:
        for g_id in scheme[rank]:
            good_snps = len(scheme[rank][g_id])
            positions = []
            for item in range(0, len(scheme[rank][g_id])):
                positions.append(scheme[rank][g_id][item]["position"])
                if len(scheme[rank][g_id][item]["qc_warnings"]) > 0:
                    good_snps -= 1
            if good_snps < min_snps:
                for item in range(0, len(scheme[rank][g_id])):
                    scheme[rank][g_id][item]["qc_warnings"].append("Lacks good tile support")

    stripped_pos = []
    for rank in scheme:
        for g_id in scheme[rank]:
            for item in range(0, len(scheme[rank][g_id])):
                if len(scheme[rank][g_id][item]["qc_warnings"]) == 0:
                    stripped_pos.append(scheme[rank][g_id][item]["position"])

    # mask unsupported groups with zeros
    for sample in memberships:
        hierarchy = memberships[sample]
        for i in range(0, len(hierarchy)):
            if i not in scheme:
                hierarchy[i] = 0
            elif not hierarchy[i] in scheme[i]:
                hierarchy[i] = 0
    # mask unsupported groups with zeros
    codes_start = pd.DataFrame(memberships, dtype=object)
    codes_start = codes_start.T

    n_unique = codes_start.apply(pd.Series.nunique)
    cols_to_drop = n_unique[n_unique == 1].index
    codes_start.drop(cols_to_drop, axis=1)
    dropped = codes_start.drop(cols_to_drop, axis=1)

    # make a vector with all the number of groups in each column
    n_unique = dropped.apply(pd.Series.nunique)
    current_codes = dict()
    biohansel_codes = dropped.copy()
    current_ids = dict()
    for i in range(0, biohansel_codes.shape[1]):
        current_max = 1
        used = dict()
        current_ids[i] = dict()
        for k in range(0, biohansel_codes.shape[0]):
            if dropped.values[k, i] not in current_ids[i].keys():
                current_ids[i][dropped.values[k, i]] = [k]
            else:
                current_ids[i][dropped.values[k, i]].append(k)
                # this stores both the current group ids under consideration but also a
                # list of all the rows that have that id
        if i == 0:
            continue
        for k in range(0, biohansel_codes.shape[0]):
            # if there was a non-zero value in the previous column build off of that node
            # go to every row cycle through ids,
            # skipping zero and identify which is connected to current line
            working_group = -1
            new_code = str()
            for thing in current_ids[i].keys():
                if k in current_ids[i][thing]:
                    working_group = thing
            # if the new group is only sample smaller than the previous group mask by copy from left
            this_group = 0
            last_group = min_parent_size
            if i > 0:
                last_group = len(current_ids[i - 1][dropped.values[k, (i - 1)]])
                this_group = len(current_ids[i][dropped.values[k, i]])
            if last_group - this_group < min_parent_size:
                biohansel_codes.values[k, i] = biohansel_codes.values[k, i - 1]
            # NOTE: NEED TO ADD QC FLAG TO SCHEME ENTRIES
            # if that row is associated with group zero at column i just copy from left
            elif working_group == 0:
                biohansel_codes.values[k, i] = biohansel_codes.values[k, i - 1]
            else:
                # if the group has already been used just pull the information
                if working_group in used.keys():
                    biohansel_codes.values[k, i] = used[working_group]
                else:
                    # if the value in the previous row is a zero start a new clade from the root
                    if biohansel_codes.values[k, i - 1] == 0:
                        for attempt in range(1, biohansel_codes.shape[0]):
                            if attempt in current_codes.keys():
                                continue
                            new_code = str(attempt)
                            break
                    # if the new group only shaves off one member ignore it and copy from the left
                    else:
                        if biohansel_codes.values[k, i - 1] in current_codes.keys():
                            current_max = max(current_codes[biohansel_codes.values[k, i - 1]]) + 1
                        else:
                            current_max = 1
                        if biohansel_codes.values[k, i - 1] not in current_codes.keys():
                            current_codes[biohansel_codes.values[k, i - 1]] = []
                        current_codes[biohansel_codes.values[k, i - 1]].append(current_max)
                        new_code = str(biohansel_codes.values[k, i - 1]) + "." + str(current_max)
                    biohansel_codes.values[k, i] = new_code
                    used[working_group] = new_code

    # for translating old column numbers to the new ones in the new matrix with dropped
    # columns get the old indexes and the indexes of the indexes is the new column number
    kept = n_unique[n_unique > 1].index
    rank_id = int()
    row_id = int()
    codes = open("codes.log", "w+")
    for i in range(0, len(leaves) - 1):
        codes.write(f"{biohansel_codes.index[i]} : {biohansel_codes.values[i, -1]}\n")
    codes.close()
    log = open(f"{vcf_file}S{min_snps}G{min_group_size}_biohansel.log", "w+")
    fasta_file = open(f"{vcf_file}S{min_snps}G{min_group_size}_biohansel.fasta", "w+")
    check = open("checking.txt", "w+")
    for i in range(0, biohansel_codes.shape[0]):
        check.write("{}\t{}".format(biohansel_codes.index.values[i],
                                    "\t".join(str(v) for v in biohansel_codes.values[i,])) + "\n")

    check.close()
    first_instance = dict()
    for rank in scheme:
        # translate the rank into the position in the table that has all
        # the undefined columns are dropped
        for i in range(0, dropped.shape[1]):
            if int(kept[i]) == int(rank):
                rank_id = i
                break
        for g_id in scheme[rank]:
            code = ""
            if path == "tree":
                # get biohansel code for that rank, group combo
                for row in range(0, dropped.shape[0]):
                    if str(codes_start.values[row, rank]) == str(g_id):
                        row_id = row
                        break
                code = biohansel_codes.values[row_id, rank_id]
                # mark the rank of first instance , if it shows up again it
                # is a masked internal node and tiles should not be output
                if code not in first_instance.keys():
                    first_instance[code] = rank_id
            # if the information from the groups in already provided use the old group codes
            if path == "groups":
                code = g_id

            for item in range(0, len(scheme[rank][g_id])):
                # exclude the the entries that have qc problems
                if len(scheme[rank][g_id][item]["qc_warnings"]) != 0:
                    continue
                if rank_id != first_instance[code]:
                    continue
                if rank_id < first_instance[code]:
                    print("There has been an error in logic")
                # extract the required information from the entries that have good qc
                position = scheme[rank][g_id][item]["position"]
                pos_tile = scheme[rank][g_id][item]["positive_tile"]
                neg_tile = scheme[rank][g_id][item]["negative_tile"]

                # write out the biohansel fasta files and the accompanying logs
                log.write(f"{position}\t{code}\t{pos_tile}\n")
                fasta_file.write(f">{position}-{code}\n{pos_tile}\n")
                log.write(f"negative{position}\t{code}\t{neg_tile}\n")
                fasta_file.write(f">negative{position}-{code}\n{neg_tile}\n")
    log.close()
    fasta_file.close()


def main():
    """
    The main body of the code that calls all the functions required to generate the biohansel tiles

    """

    # collect relevant user input and parse it into the appropriate variables
    args = parse_args()
    nwk_treefile = args.in_nwk
    reference_fasta = args.reference
    vcf_file = args.in_vcf
    min_snps = args.min_snps
    min_group_size = args.min_members
    group_info = args.group_info
    min_parent_size = args.min_parent
    tile_generator(nwk_treefile, reference_fasta, vcf_file, min_snps, min_group_size, group_info,
                   min_parent_size)


# call main function


if __name__ == '__main__':
    main()
