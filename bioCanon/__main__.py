"""
Biohansel requires a strain specific fasta format scheme to be able to classify samples.
Schemes exist for some common clonal pathogens but there are many untapped applications.
This module includes functions to aid in the generation of k-mers linked to predefined groups
to serve as the basis for automatic scheme development
"""
import csv
import re
import os
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import vcfpy
from Bio import SeqIO
from ete3 import Tree
from ete3 import parser


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
    parser_inner.add_argument('--group_info', type=str, required=False, default="none",
                              help="A tsv file that contains the leaf ids and group information,"
                                   " for specific formatting information please see the README")
    parser_inner.add_argument('--reference', type=str, required=True,
                              help='Reference fasta sequence from VCF')
    parser_inner.add_argument('--min_snps', type=int, required=False, default=2,
                              help="Number of cannonical SNPs required to define a group")
    parser_inner.add_argument('--min_members', type=int, required=False, default=5,
                              help="Minimum number of members to define a group")
    parser_inner.add_argument('--min_parent', type=int, required=False, default=2,
                              help="Minimum size difference between new "
                                   "group and parent group to be considered valid")
    parser_inner.add_argument('--flanking', type=int, required=False, default=15,
                              help="number of flanking positions on either side used in k-mer"
                                   "generation")
    parser_inner.add_argument('--outdir', type=str, required=False, default=os.getcwd(),
                              help='Output Directory to put results')
    return parser_inner.parse_args()


def get_subtrees(tree_mem):
    """
       break the tree membership information into which subtrees of a given size are present
       preparatory step for saving computational steps

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
        for key in subtree:
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


def generate_new_codes(checking, switch, groups, used, iterators: list):
    """
    if supplied group codes are not biohansel compatible generate new ones
    Parameters
    ----------
    checking: str
        old code that contains invalid characters
    switch: dict
        information on which new assignments have already been made
    groups: dict
        new group information after relabeling to biohansel compliant codes
    used: dict
        information on which codes have previously been used
    iterators: list
        contains the loop index and key under consideration
    """
    i = iterators[0]
    key = iterators[1]
    new = str()
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
    for key in groups:
        while len(groups[key]) < max_len:
            groups[key].append(0)
    # ensure that all entries are numbers and replace to unique numbers if not already
    non_valid = r'[^0-9/.]'
    warning = 0
    for i in range(0, max_len):
        switch = dict()
        for key in groups:
            checking = str(groups[key][i])
            if i > 0:
                checking = str(groups[key][i - 1]) + "." + checking
            # check to see if the value contains only the accepted characters,
            # if it contains illegal characters rename the groups
            if re.search(non_valid, checking):
                if warning == 0:
                    print("provided groups contain invalid characters for biohansel codes."
                          "New codes will be generated")
                    warning = 1
                iterators = [i, key]
                generate_new_codes(checking, switch, groups, used, iterators)
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
            Takes an ete3 format tree for processing into groupings indicated by the tree structure
        Returns
        -------
        memberships: dict
            Contains the group information with the leaves as keys and the list of groups as values
    """
    # initialize variables
    memberships = dict()
    group = 1
    # visit all nodes by depth from the root and generate a list of group memberships for each leaf
    # a leaf is a part of the same membership as another leaf if they share a common
    # interior parental node the root of the tree is ignored in this calculation
    for node in ete3_tree_obj.iter_descendants("levelorder"):
        names = node.get_leaf_names()
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


def mask_unsupported(memberships, scheme):
    """
    replace unsupported group identifiers with zeros
    Parameters
    ----------
    memberships: dict
        Contains the group information with the leaves as keys and the list of groups as values
    scheme: dict
        storage location for the information needed to generate the biohansel scheme output
    """
    for sample in memberships:
        hierarchy = memberships[sample]
        for i in range(0, len(hierarchy)):
            if i not in scheme:
                hierarchy[i] = 0
            elif not hierarchy[i] in scheme[i]:
                hierarchy[i] = 0


def replace_conflict(conflict, start, scheme, scheme_pieces: list):
    """
    change ambiguous positions to Ns in the k-mers
    Parameters
    ----------
    conflict, start: int
        genome positions relevant to ambiguous k-mer generation
    scheme: dict
        storage location for the information needed to generate the biohansel scheme output
    scheme_pieces: list
        contains the keys required to navigate to the appropriate location in the scheme

    """
    rank = scheme_pieces[0]
    g_id = scheme_pieces[1]
    item = scheme_pieces[2]
    index = conflict - start - 1

    p_substring = scheme[rank][g_id][item]["positive_tile"][:index + 1]
    n_substring = scheme[rank][g_id][item]["negative_tile"][:index + 1]
    p_replace = scheme[rank][g_id][item]["positive_tile"][:index] + "N"
    n_replace = (scheme[rank][g_id][item]["negative_tile"][:index] + "N")

    scheme[rank][g_id][item]["positive_tile"] = \
        scheme[rank][g_id][item]["positive_tile"].replace(p_substring, p_replace, 1)
    scheme[rank][g_id][item]["negative_tile"] = \
        scheme[rank][g_id][item]["negative_tile"].replace(n_substring, n_replace, 1)


def add_tiles(scheme, ref_seq, flanking, conflict_positions):
    """
    using the information from the vcf file and the reference sequence to generate the positive and
    negative k-mers and add them to the scheme object
    Parameters
    ----------
    scheme: dict
        storage location for the information needed to generate the biohansel scheme output
    ref_seq: seq
        contains the user supplied reference sequence
    flanking: int
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
                start = position - 1 - flanking
                # which is the positive tile is determined by if the ref or alt base defines
                # the in-group
                if scheme[rank][g_id][item]["positive_group"] == "ref":
                    pos_tile = ref_seq[start:position - 1] + scheme[rank][g_id][item]["ref_base"] +\
                               ref_seq[position:position + flanking]
                    neg_tile = ref_seq[start:position - 1] + scheme[rank][g_id][item]["alt_base"] +\
                               ref_seq[position:position + flanking]
                elif scheme[rank][g_id][item]["positive_group"] == "alt":
                    pos_tile = ref_seq[start:position - 1] + scheme[rank][g_id][item]["alt_base"] +\
                               ref_seq[position:position + flanking]
                    neg_tile = ref_seq[start:position - 1] + scheme[rank][g_id][item]["ref_base"] +\
                               ref_seq[position:position + flanking]
                else:
                    print("something went wrong")

                scheme[rank][g_id][item]["positive_tile"] = pos_tile
                scheme[rank][g_id][item]["negative_tile"] = neg_tile
                # if position has conflicts add degenerate bases
                if position in conflict_positions.keys():
                    for conflict in conflict_positions[position]:
                        scheme_pieces = [rank, g_id, item]
                        replace_conflict(conflict, start, scheme, scheme_pieces)
                    if len(conflict_positions[position]) >= 1:
                        scheme[rank][g_id][item]["qc_warnings"].append(
                            "One or more degenerate bases")
    return scheme


def filter_by_snps(scheme, min_snps, ignored_for_snps):
    """
    Takes the current scheme and removes any entries that do not meet the minimum snp support
    requirement.  Additionally it marks down the location of variable positions for
    tile generation.
    Parameters
    ----------
    scheme: dict
        contains the current working scheme information
    min_snps: int
        user specified desired number of snps required for a group to be validly supported
    ignored_for_snps: int
        number of snps that match a branch point but not the minimum snp support requirement
    Returns
    -------
    required_tiles: list
        List of positions for snps that support a valid grouping
    ignored_for_snps: int
        number of snps that match a branch point but not the minimum snp support requirement
    """
    # initialize required values
    required_tiles = list()
    number_snps = 0
    number_groups = 0

    # iterate through all of the pieces of the scheme
    for rank in scheme:
        valid = dict()
        prev = []
        for g_id in scheme[rank]:
            # the number of snps that support a group, rank combination can be obtained from the
            # length of the list that contains them.  It if is shorter than the user input cutoff
            # that group, rank combination is not valid
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
                        # there are occasional duplication of group rank support, this allows for
                        # the removal of the duplicate
                        if prev == scheme[rank][g_id][item]["position"]:
                            duplicates.append(item)
                            continue
                        prev = scheme[rank][g_id][item]["position"]
                        required_tiles.append(scheme[rank][g_id][item]["position"])
                        number_snps += 1
                duplicates.sort(reverse=True)
                for i in range(0, len(duplicates)):
                    scheme[rank][g_id].pop(i)
            else:
                ignored_for_snps += len(scheme[rank][g_id])
        scheme[rank] = valid
    return required_tiles, ignored_for_snps


def alt_or_ref(record, samples: list):
    """
    takes in a single record in a vcf file and returns the sample names divided into two lists:
    ones that have the reference snp state and ones that have the alternative snp state
    Parameters
    ----------
    record
        the record supplied by the vcf reader
    samples: list
        list of sample names
    Returns
    -------
    ref_group, alt_group : list
        lists of samples divided by ref or alt snp state
    """
    tracker = 0
    ref_group = []
    alt_group = []
    for call in record.calls:
        state = int(call.data.get('GT'))
        sample_name = samples[tracker]
        if state == 0:
            ref_group.append(sample_name)
        elif state == 1:
            alt_group.append(sample_name)
        else:
            print("there is a problem reading the state information")
            raise SystemExit(0)
        tracker += 1
    return ref_group, alt_group


def path_check(group_info, nwk_treefile):
    """
    checks that the sourcing on the group information is present and clear then generates the group
    information from the provided data
    Parameters
    ----------
    group_info, nwk_treefile: str
        one of the variables should include a file name from which the group information can be
        obtained

    Returns
    -------
    memberships: dict
        Contains the group information with the leaves as keys and the list of groups as values
    leaves: list
        list of sample names
    path: str
        identifier for which input was supplied
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
    return memberships, leaves, path


def id_conflict(required_tiles, flanking, all_variable: list):
    """
    search for the positions of snps that collide with other snps in the user defined k-mer window
    Parameters
    ----------
    required_tiles: list
        position information of all snps that support groups
    flanking: int
        size of flanking sequence on either side of the snp
    all_variable:list
        position information of all variable sites in the vfc file
    Returns
    -------
    conflict positions: dict
        data frame that contains the variable site positions around a target snp
    """
    stored_place = 0
    conflict_positions = dict()
    for item in required_tiles:
        if stored_place < flanking + 1:
            start = 0
        elif item - all_variable[stored_place] > flanking:
            start = stored_place
        else:
            start = stored_place - flanking
        for i in range(start, len(all_variable)):
            difference = item - all_variable[i]
            if difference > flanking:
                continue
            if difference < (-1 * flanking):
                stored_place = i
                break
            if flanking >= difference >= (-1 * flanking):
                if difference == 0:
                    continue
                if item not in conflict_positions.keys():
                    conflict_positions[item] = [all_variable[i]]
                else:
                    conflict_positions[item].append(all_variable[i])
            else:
                print(f"logic error\t : {difference}")
    return conflict_positions


def make_biohansel_codes(biohansel_codes, dropped, min_parent_size):
    """
    generate hierarchical codes for the groups identified
    Parameters
    ----------
    biohansel_codes : pandas data frame
        starts with the memberships and is then altered throughout the function into biohansel
        compatible group codes
    dropped : pandas data frame
        contains membership information excluding empty columns
    min_parent_size: int
        user defined requirement for differnce between parent node and child node size
    """
    current_codes = dict()
    current_ids = dict()
    for i in range(0, biohansel_codes.shape[1]):
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


def tile_generator(reference_fasta, vcf_file, numerical_parameters, groups, outdir):
    """
    main function for generating potential biohansel tiles from group and variant call information
    Parameters
    ----------
    reference_fasta, vcf_file: str
        user supplied input file names obtained from the parse_args function
    numerical_parameters: list
        user supplied configuration settings from the parse_args function
    groups: list
        contains information of memberships, leaf names and information source
    outdir: str
        user specified name of output directory
    """
    min_snps, min_group_size, min_parent_size, flanking = numerical_parameters
    memberships, leaves, path = groups
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
    all_variable = []
    total_snps = 0
    ignored_for_groupsize = 0
    ignored_for_multiple_states = 0
    ignored_for_degenerate_kmer = 0
    no_branchpoint = 0
    branchpoint_snps = 0
    ignored_for_snps = 0

    for record in reader:
        # skip over the metadata lines at the top of the file
        if not record.is_snv():
            continue
        total_snps += 1
        # pull the positional information and the reference
        # value and all the alt values from the vcf file
        line = [record.CHROM, record.POS, record.REF]
        line += [alt.value for alt in record.ALT]
        # skip any vcf calls that have multiple alternative alleles
        if len(record.ALT) > 1:
            ignored_for_multiple_states += 1
            continue
        # initialize variables and pull the ref and alt SNPs
        alt_base = record.ALT[0].value
        ref_base = record.REF
        position = record.POS
        all_variable.append(position)
        # go through the samples and divide them by alt vs ref bases
        ref_group, alt_group = alt_or_ref(record, samples)
        # generate places to put group and rank information
        valid_ranks = list()
        valid_ranks_alt = list()
        valid_groups = list()
        valid_groups_alt = list()
        ref_cand_part = dict()
        alt_cand_part = dict()
        # if there is variation at the snp location search to see
        # if the partitioning suggested by the snp is
        # one that is present in the provided tree
        if len(ref_group) >= 1 and len(alt_group) >= 1:
            if len(ref_group) in subtrees.keys():
                valid_ranks, valid_groups = search_st(ref_group, valid_ranks, valid_groups,
                                                      subtrees)
            if len(valid_ranks) >= 1:
                for i in range(0, len(valid_ranks)):
                    if valid_ranks[i] not in ref_cand_part:
                        ref_cand_part[valid_ranks[i]] = list()
                    if valid_groups[i] == 0:
                        continue
                    # if the partition is already supported add the information about that
                    # snp if not created it
                    ref_cand_part[valid_ranks[i]].append({"position": position,
                                                          "ref_base": ref_base,
                                                          "alt_base": alt_base,
                                                          "ref_count": len(ref_group),
                                                          "alt_count": len(alt_group),
                                                          "positive_group": "ref",
                                                          "g_id": valid_groups[i],
                                                          "rank_id": valid_ranks[i],
                                                          "qc_warnings": []})
            if len(alt_group) in subtrees.keys():
                valid_ranks_alt, valid_groups_alt = search_st(alt_group, valid_ranks, valid_groups,
                                                              subtrees)
            if len(valid_ranks_alt) > 0:
                for i in range(0, len(valid_ranks_alt)):
                    if valid_ranks_alt[i] not in alt_cand_part:
                        alt_cand_part[valid_ranks_alt[i]] = list()
                    if valid_groups_alt[i] == 0:
                        continue
                    # if the partition is already supported add the information about that snp if
                    # not created it
                    alt_cand_part[valid_ranks_alt[i]].append({"position": position,
                                                              "ref_base": ref_base,
                                                              "alt_base": alt_base,
                                                              "ref_count": len(ref_group),
                                                              "alt_count": len(alt_group),
                                                              "positive_group": "alt",
                                                              "g_id": valid_groups_alt[i],
                                                              "rank_id": valid_ranks_alt[i],
                                                              "qc_warnings": []})
            for rank in ref_cand_part:
                if rank not in scheme:
                    scheme[rank] = dict()
                for item in ref_cand_part[rank]:
                    # filter out groups which are too small
                    if item['ref_count'] < min_group_size:
                        ignored_for_groupsize +=1
                        continue
                    if item['alt_count'] < min_group_size:
                        ignored_for_groupsize += 1
                        continue
                    if not item["g_id"] in scheme[rank]:
                        scheme[rank][item["g_id"]] = []
                    scheme[rank][item["g_id"]].append(item)
            for rank in alt_cand_part:
                if rank not in scheme:
                    scheme[rank] = dict()
                for item in alt_cand_part[rank]:
                    # filter out groups which are too small
                    if item['alt_count'] < min_group_size:
                        continue
                    if not item["g_id"] in scheme[rank]:
                        scheme[rank][item["g_id"]] = []
                    scheme[rank][item["g_id"]].append(item)
            if not valid_ranks_alt and not valid_ranks:
                no_branchpoint += 1
        else:
            ignored_for_groupsize += 1
    # filter out groups with less than minimum number of supporting snps
    required_tiles, ignored_for_snps = filter_by_snps(scheme, min_snps, ignored_for_snps)
    # sort the lists of both all the variable positions and the ones of
    # interest for the tile selection
    all_variable.sort()
    required_tiles.sort()
    # to reduce calculations store the location of the last position checked for conflict
    conflict_positions = id_conflict(required_tiles, flanking, all_variable)
    # using the conflict position information pull the appropriate tiles from the reference genome
    scheme = add_tiles(scheme, ref_seq, flanking, conflict_positions)
    # double check that degenerate bases haven't removed support for a group
    # and that all support snps are not from the same region
    for rank in scheme:
        for g_id in scheme[rank]:
            good_snps = len(scheme[rank][g_id])
            positions = []
            for item in range(0, len(scheme[rank][g_id])):
                positions.append(scheme[rank][g_id][item]["position"])
                if len(scheme[rank][g_id][item]["qc_warnings"]) > 0:
                    ignored_for_degenerate_kmer +=1
                    good_snps -= 1
            if good_snps < min_snps:
                ignored_for_snps += good_snps
                for item in range(0, len(scheme[rank][g_id])):
                    scheme[rank][g_id][item]["qc_warnings"].append("Lacks good tile support")

            else:
                branchpoint_snps += good_snps

    stripped_pos = []
    for rank in scheme:
        for g_id in scheme[rank]:
            for item in range(0, len(scheme[rank][g_id])):
                if len(scheme[rank][g_id][item]["qc_warnings"]) == 0:
                    stripped_pos.append(scheme[rank][g_id][item]["position"])

    # mask unsupported groups with zeros
    mask_unsupported(memberships, scheme)
    codes_start = pd.DataFrame(memberships, dtype=object)
    codes_start = codes_start.T

    n_unique = codes_start.apply(pd.Series.nunique)
    cols_to_drop = n_unique[n_unique == 1].index
    codes_start.drop(cols_to_drop, axis=1)
    dropped = codes_start.drop(cols_to_drop, axis=1)

    # make a vector with all the number of groups in each column
    n_unique = dropped.apply(pd.Series.nunique)
    biohansel_codes = dropped.copy()
    make_biohansel_codes(biohansel_codes, dropped, min_parent_size)

    # for translating old column numbers to the new ones in the new matrix with dropped
    # columns get the old indexes and the indexes of the indexes is the new column number
    kept = n_unique[n_unique > 1].index
    rank_id = int()
    row_id = int()
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    out_path = os.path.join(os.getcwd(), outdir)
    codes = open(os.path.join(out_path, "codes.log"), "w+")
    for i in range(0, len(leaves) - 1):
        if len(biohansel_codes.values[i]) == 0:
            print("Something has gone wrong with the code assignments")
            raise SystemExit(0)
        codes.write(f"{biohansel_codes.index[i]} : {biohansel_codes.values[i, -1]}\n")
    codes.close()
    log = open(os.path.join(out_path, f"S{min_snps}G{min_group_size}_biohansel.log"),
               "w+")
    fasta_file = open(os.path.join(out_path, f"S{min_snps}G{min_group_size}"
                                             f"_biohansel.fasta"), "w+")
    snp_report = open(os.path.join(out_path, f"S{min_snps}G{min_group_size}"
                                             f"snp_report.txt"), "w+")

    not_mapped = np.setdiff1d(all_variable, required_tiles)
    print(not_mapped)
    print(len(not_mapped))
    for i in not_mapped:
        snp_report.write(f"The snp at {i} was not included in the scheme because it could not be"
                         f" matched to a branch point in the data.\n")
    first_instance = dict()
    if total_snps != (ignored_for_multiple_states + ignored_for_degenerate_kmer +
                      ignored_for_groupsize + branchpoint_snps + no_branchpoint + ignored_for_snps):
        print("Count logic off")
    print(f"There were {total_snps} in the vcf file.\n "
          f"{no_branchpoint} did not match a branch point\n "
          f"{ignored_for_groupsize} of these were ignored due to not matching the group "
          f"size requirement\n "
          f"{ignored_for_snps} of these because their group did not have "
          f"at least {min_snps} snps for support\n "
          f"{ignored_for_degenerate_kmer} of these were ignored because they produced "
          f"degenerate k-mer tiles\n {ignored_for_multiple_states} were ignored due to "
          f"having more than two alt states\n {branchpoint_snps} were found to match to a "
          f"branch point.")

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
                # mark the rank of first instance , if it shows up again it
                # is a masked internal node and tiles should not be output
                if code not in first_instance.keys():
                    first_instance[code] = rank_id
            for item in range(0, len(scheme[rank][g_id])):
                position = scheme[rank][g_id][item]["position"]
                # exclude the the entries that have qc problems
                if len(scheme[rank][g_id][item]["qc_warnings"]) != 0:
                    reasons = scheme[rank][g_id][item]["qc_warnings"]
                    snp_report.write(f"The snp at {position} was not included in the scheme due "
                                     f"to {reasons}\n")
                    continue
                if rank_id != first_instance[code]:
                    continue
                if rank_id < first_instance[code]:
                    print("There has been an error in logic")
                # extract the required information from the entries that have good qc
                pos_tile = scheme[rank][g_id][item]["positive_tile"]
                neg_tile = scheme[rank][g_id][item]["negative_tile"]

                # write out the biohansel fasta files and the accompanying logs
                snp_report.write(f"The snp at {position} was found to support the {code} group.\n")
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
    groups = path_check(args.group_info, args.in_nwk)
    numerical_parameters = [args.min_snps, args.min_members, args.min_parent, args.flanking]
    tile_generator(args.reference, args.in_vcf, numerical_parameters, groups, args.outdir)
# call main function


if __name__ == '__main__':
    main()
