from ete3 import Tree
from Bio import SeqIO
import vcfpy
from argparse import (ArgumentParser, FileType)


def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(
        description='BIO_Hansel Scheme development')
    parser.add_argument('--in_vcf', type=str, required=True, help='VCF File of SNPs')
    parser.add_argument('--in_nwk', type=str, required=True, help='Newick Tree of strains')
    parser.add_argument('--reference', type=str, required=True, help='Reference fasta sequence from VCF')
    parser.add_argument('--outdir', type=str, required=False, help='Output Directory to put results')
    parser.add_argument('--min_snps', type=int, required=False, default=2,
                        help="Number of cannonical SNPs required to define a group")
    parser.add_argument('--min_members', type=int, required=False, default=5,
                        help="Minimum number of members to define a group")
    return parser.parse_args()


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
        from_subtrees=item[0]
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
    for n in memberships:
        print("{}\t{}".format(n, "\t".join(str(v) for v in memberships[n])))
    return memberships


def main():

    """    # collect relevant user input and parse it into the appropriate variables
    args = parse_args()
    nwk_treefile = args.in_nwk
    reference_fasta = args.reference
    vcf_file = args.in_vcf
    min_snps = args.min_snps
    min_group_size = args.min_members"""

    nwk_treefile = "stripped.nwk"
    reference_fasta = "NZ_CP012921.fasta"
    vcf_file = "stripped_core.vcf"
    min_snps = 2
    min_group_size = 5

    # hard coded size of positive and negative tiles surrounding SNP
    tile_size = 32
    tile_size = tile_size - 1
    mid_point = int(tile_size / 2)
    # obtain the dictionary of leaf subtree membership after pre-processing the raw Newick Tree to resolve polytomys
    try:
        tree=parse_tree(nwk_treefile)
        memberships = get_tree_groups(tree)
        for leaf in tree.traverse():
            if leaf.is_leaf():
                leaves.append(leaf.name)
    except parser.newick.NewickError:
        print(
            f"There was an error in reading {nwk_treefile}.  "
            f"Verify that the file exists and is in the correct format")
        raise SystemExit(0)

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
            print(ref_length)
            if ref_length==0:
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
        print("Something went wrong.  The leaf names in the provided Newick file "
              "are not the same as the ones in the VCF file")
        print(f"The samples found in one of the files but not the other are : {discr}")
        raise SystemExit(0)
    scheme = dict()
    number_snps = 0
    number_groups = 0
    # print(subtrees)
    for record in reader:
        # skip over the metadata lines at the top of the file
        if not record.is_snv():
            continue

        # pull the positional information and the reference value and all the alt values from the vcf file
        # does not appear to be used later
        line = [record.CHROM, record.POS, record.REF]
        line += [alt.value for alt in record.ALT]

        # skip any vcf calls that have multiple alternative alleles
        if len(record.ALT) > 1:
            continue

        # initialize variables and pull the ref and alt SNPs
        alt_base = record.ALT[0].value
        ref_base = record.REF
        position = record.POS
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
                        if item["position"]==position:
                            print(f"Warning duplication of {position}")
                    ref_candidate_partitions[valid_ranks[i]].append({"position": position,
                                                                     "ref_base": ref_base,
                                                                     "alt_base": alt_base,
                                                                     "ref_count": count_ref,
                                                                     "alt_count": count_alt,
                                                                     "positive_tile": "positive_tile",
                                                                     "negative_tile": "negative_tile",
                                                                     "group_id": valid_groups[i],
                                                                     "rank_id": valid_ranks[i]})
            valid_ranks, valid_groups = search_st(alt_group, valid_ranks, valid_groups, subtrees)
            if len(valid_ranks) > 0:
                for i in range(0, len(valid_ranks)):
                    if valid_ranks[i] not in alt_candidate_partitions:
                        alt_candidate_partitions[valid_ranks[i]] = list()
                    if valid_groups[i] == 0:
                        continue
                    alt_candidate_partitions[valid_ranks[i]].append({"position": position,
                                                                     "ref_base": ref_base,
                                                                     "alt_base": alt_base,
                                                                     "ref_count": count_ref,
                                                                     "alt_count": count_alt,
                                                                     "positive_tile": "positive_tile",
                                                                     "negative_tile": "negative_tile",
                                                                     "group_id": valid_groups[i],
                                                                     "rank_id": valid_ranks[i]})
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
                    pos = item["position"]
                    if pos not in positions:
                        positions[pos] = 0
                    positions[pos] += 1

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
                    pos = item["position"]
                    if pos not in positions:
                        positions[pos] = 0
                    positions[pos] += 1
                    if positions[pos] > 1:
                        conflict_positions[pos] = item

            """            for rank in ref_candidate_partitions:
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
    for rank in scheme:
        valid = dict()
        for group_id in scheme[rank]:
            if len(scheme[rank][group_id]) >= min_snps:
                valid[group_id] = scheme[rank][group_id]
                number_snps += len(scheme[rank][group_id])
                number_groups += 1
                for item in range(0,len(scheme[rank][group_id])):
                    print("{}\t{}".format(scheme[rank][group_id][item], f"Min groups:{min_group_size},\t Min SNPs: {min_snps}"))
        scheme[rank] = valid

    print(f"\n When the minimum group size is {min_group_size}, there are {number_snps} SNPs"
          f" that support a total of {number_groups}\n groups")
    # mask unsupported groups with zeros
    for n in memberships:
        hierarchy = memberships[n]
        for i in range(0, len(hierarchy)):
            if i not in scheme:
                hierarchy[i] = 0
            elif not hierarchy[i] in scheme[i]:
                hierarchy[i] = 0
        print("{}\t{}".format(n, "\t".join(str(v) for v in memberships[n])))
    test = pd.DataFrame(memberships)
    test = test.T

    nunique = test.apply(pd.Series.nunique)
    cols_to_drop = nunique[nunique == 1].index
    test.drop(cols_to_drop, axis=1)
    dropped = test.drop(cols_to_drop, axis=1)

    # make a vector with all the number of groups in each column
    nunique = dropped.apply(pd.Series.nunique)
    col_groups = list(nunique)

    for i in range(0, len(col_groups))
# call main function
if __name__ == '__main__':
    main()
