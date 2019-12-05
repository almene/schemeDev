from ete3 import Tree
from ete3 import parser
from Bio import SeqIO
import vcfpy
from argparse import (ArgumentParser, FileType)
import numpy as np
import random
import pandas as pd

nwk_treefile = "stripped.nwk"
reference_fasta = "NZ_CP012921.fasta"
vcf_file = "stripped_core.vcf"
min_snps = 2
min_group_size = 5


tile_size = 32
tile_size = tile_size - 1
mid_point = int(tile_size / 2)
leaves = []
# required spacing between at least two snps in the group support
size_cutoff = 200000
# test that files exist
# obtain the dictionary of leaf subtree membership after pre-processing the raw Newick Tree to resolve polytomys

tree = parse_tree(nwk_treefile)
memberships = get_tree_groups(tree)
for leaf in tree.traverse():
	if leaf.is_leaf():
		leaves.append(leaf.name)


# process the memberships data to extract the subtree information from the membership information
subtrees = get_subtrees(memberships)
# read the reference fasta file by record


for seq_record in SeqIO.parse(reference_fasta, "fasta"):
	ref_seq = "{}".format(seq_record.seq)
	ref_length = len(ref_seq)
	if ref_length == 0:
		print("error in .fasta file")
		raise SystemExit(0)

# collect the vcf file indicated in the arguments passed to the script
reader = vcfpy.Reader.from_path(vcf_file)

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
				ref_candidate_partitions[valid_ranks[i]].append({"position": position, "ref_base": ref_base,\
																 "alt_base": alt_base,\
																 "ref_count": count_ref,\
																 "alt_count": count_alt,\
																 "positive_group": "ref",\
																 "group_id": valid_groups[i],\
																 "rank_id": valid_ranks[i],\
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
				alt_candidate_partitions[valid_ranks[i]].append({"position": position,\
																 "ref_base": ref_base,\
																 "alt_base": alt_base,\
																 "ref_count": count_ref,\
																 "alt_count": count_alt,\
																 "positive_group": "alt",\
																 "group_id": valid_groups[i],\
																 "rank_id": valid_ranks[i],\
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
