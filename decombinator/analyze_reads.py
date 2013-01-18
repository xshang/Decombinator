#! /usr/bin/env python
# encoding: utf-8
"""
Decombinator is a tool for the fast, efficient analysis of T cell receptor (TcR)
repertoire samples, designed to be accessible to those with no previous
programming experience. It is based on the Aho-Corasick algorithm which
uses a finite state automaton (FSA) to quickly assign a specific V and J
gene segment. From these assignments, it is then able to determine the number
of germline V and J deletions and the string of contiguous nucleotides which
lie between the 3' end of the V gene segment and the 5' end of the J gene
segment. These 5 variables form the identifier which uniquely categorises each
distinct TcR sequence. For more details, please see (Thomas et al.).
Decombinator assumes no prior programming experience.

original authors:
Niclas Thomas, James Heather, Wilfred Ndifon, John Shawe-Taylor, Benny Chain

>> Analyzing 1 file(s)
>> Importing known V, D and J gene segments and tags...
Importing sequences from testseqs.fastq  and assigning V and J regions...
17569 sequences were analysed
4066  sequences were successfully assigned
Time taken = 2.5801949501 seconds
"""
import sys
import string
import numpy as np
import decimal as dec
import Levenshtein as lev
from Bio import SeqIO
from time import time
from toolshed import nopen
from string import Template
from acora import AcoraBuilder

def analysis(fastqs, v_fasta, j_fasta, 
                with_reverse_complement_search=True, verbose=True):
    if verbose:
        sys.stderr.write('>> Analyzing %d file(s)\n' % len(fastqs))
        sys.stderr.write(">> Importing known V, and J gene segments and tags\n")

    # get the sequences per region
    v_genes = list(SeqIO.parse(nopen(v_fasta), "fasta"))
    j_genes = list(SeqIO.parse(nopen(j_fasta), "fasta"))
    v_regions = [v_genes[i].seq.upper() for i, v in enumerate(v_genes)]
    j_regions = [j_genes[i].seq.upper() for i, v in enumerate(j_genes)]

    ## Build keyword tries of V and J tags for fast assignment
    v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("tags_trbv.txt", "rU"), 10)
    j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("tags_trbj.txt", "rU"), 6)   

    v_builder = AcoraBuilder()
    # Add all V tags to keyword trie
    [v_builder.add(v_seqs[i]) for i, v in enumerate(v_seqs)]
    v_key = v_builder.build()
    
    j_builder = AcoraBuilder()
    # Add all J tags to keyword trie
    [j_builder.add(j_seqs[i]) for i, v in enumerate(j_seqs)]
    j_key = j_builder.build()

    # Build keyword tries for first and second halves of both V and J tags
    v_half1_builder = AcoraBuilder()
    for i in range(0,len(half1_v_seqs)):
        v_half1_builder.add(str(half1_v_seqs[i]))
    half1_v_key = v_half1_builder.build()

    v_half2_builder = AcoraBuilder()
    for i in range(0,len(half2_v_seqs)):
        v_half2_builder.add(str(half2_v_seqs[i]))
    half2_v_key = v_half2_builder.build()

    j_half1_builder = AcoraBuilder()
    for i in range(0,len(half1_j_seqs)):
        j_half1_builder.add(str(half1_j_seqs[i]))
    half1_j_key = j_half1_builder.build()

    j_half2_builder = AcoraBuilder()
    for i in range(0,len(half2_j_seqs)):
        j_half2_builder.add(str(half2_j_seqs[i]))
    half2_j_key = j_half2_builder.build()

    # correctly assign a seq read with all desired variables
    assigned_count = 0
    # this will simply track the number of sequences analysed in file
    seq_count = 0
    t0 = time() # Begin timer


    stemplate = Template('$v $j $del_v $del_j $nt_insert')

    for fastq in fastqs:
        #XXX
        print 'Importing sequences from', fastq, 'and assigning V and J regions...'
        
        for record in SeqIO.parse(nopen(fastq), "fastq"):
            
            found_seq_match = 0
            seq_count += 1
            
            hold_v = v_key.findall(str(record.seq))
            hold_j = j_key.findall(str(record.seq))

            if hold_v:
                # Assigns V
                v_match = v_seqs.index(hold_v[0][0])
                # Finds where the end of a full V would be
                temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1
                # If the number of deletions has been found
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ):
                    [ end_v, deletions_v] = get_v_deletions( record.seq, v_match, temp_end_v, v_regions )
            else:
                found_v_match = 0
                hold_v1 = half1_v_key.findall(str(record.seq))
                hold_v2 = half2_v_key.findall(str(record.seq))
                for i in range(len(hold_v1)):
                    indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
                    for k in indices:
                        if len(v_seqs[k]) == len(str(record.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                            if lev.hamming( v_seqs[k], str(record.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                                v_match = k
                                # Finds where the end of a full V would be
                                temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1
                                found_v_match += 1
                for i in range(len(hold_v2)):
                    indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                    for k in indices:
                        if len(v_seqs[k]) == len(str(record.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
                            if lev.hamming( v_seqs[k], str(record.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[k])] ) <= 1:
                                v_match = k
                                # Finds where the end of a full V would be
                                temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - 1
                                found_v_match += 1

            if hold_j:
                # Assigns J
                j_match = j_seqs.index(hold_j[0][0])
                # Finds where the start of a full J would be
                temp_start_j = hold_j[0][1] - jump_to_start_j[j_match]
                if get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
            else:
                found_j_match = 0
                hold_j1 = half1_j_key.findall(str(record.seq))
                hold_j2 = half2_j_key.findall(str(record.seq))
                for i in range(len(hold_j1)):
                    indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
                    for k in indices:
                        if len(j_seqs[k]) == len(str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                            if lev.hamming( j_seqs[k], str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                j_match = half1_j_seqs.index(hold_j1[i][0])
                                # Finds where the start of a full J would be
                                temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match]
                                found_j_match += 1
                for i in range(len(hold_j2)):
                    indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                    for k in indices:
                        if len(j_seqs[k]) == len(str(record.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
                            if lev.hamming( j_seqs[k], str(record.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[k])] ) <= 1:
                                j_match = half2_j_seqs.index(hold_j2[i][0])
                                # Finds where the start of a full J would be
                                temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - 6
                                found_j_match += 1

            if hold_v and hold_j:
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = record.seq[end_v+1:start_j])
                    # Write to analysis_file (text file) the classification of the sequence
                    print f_seq
                    assigned_count += 1
                    found_seq_match = 1
            elif hold_v and found_j_match == 1:
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                    print f_seq
                    assigned_count += 1
                    found_seq_match = 1
            elif found_v_match == 1 and hold_j:
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                    print f_seq
                    assigned_count += 1
                    found_seq_match = 1
            elif found_v_match == 1 and found_j_match == 1:
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                    print f_seq
                    assigned_count += 1
                    found_seq_match = 1

            if found_seq_match == 0 and with_reverse_complement_search:
                
                #####################
                # REVERSE COMPLEMENT
                #####################

                record_reverse = record.reverse_complement()
                hold_v = v_key.findall(str(record_reverse.seq))
                hold_j = j_key.findall(str(record_reverse.seq))

                if hold_v:
                    # Assigns V
                    v_match = v_seqs.index(hold_v[0][0])
                    # Finds where the end of a full V would be
                    temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1
                    # If the number of deletions has been found
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ):
                        [ end_v, deletions_v] = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                else:
                    found_v_match = 0
                    hold_v1 = half1_v_key.findall(str(record_reverse.seq))
                    hold_v2 = half2_v_key.findall(str(record_reverse.seq))
                    for i in range(len(hold_v1)):
                        indices = [y for y, x in enumerate(half1_v_seqs) if x == hold_v1[i][0] ]
                        for k in indices:
                            if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[half1_v_seqs.index(hold_v1[i][0])])]):
                                if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                                    v_match = k
                                    # Finds where the end of a full V would be
                                    temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1
                                    found_v_match += 1
                    for i in range(len(hold_v2)):
                        indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                        for k in indices:
                            if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
                                if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[k])] ) <= 1:
                                    v_match = k
                                    # Finds where the end of a full V would be
                                    temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - 1
                                    found_v_match += 1

                if hold_j:
                    # Assigns J
                    j_match = j_seqs.index(hold_j[0][0])
                    # Finds where the start of a full J would be
                    temp_start_j = hold_j[0][1] - jump_to_start_j[j_match]
                    if get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        [ start_j, deletions_j] = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                else:
                    found_j_match = 0
                    hold_j1 = half1_j_key.findall(str(record_reverse.seq))
                    hold_j2 = half2_j_key.findall(str(record_reverse.seq))
                    for i in range(len(hold_j1)):
                        indices = [y for y, x in enumerate(half1_j_seqs) if x == hold_j1[i][0] ]
                        for k in indices:
                            if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[half1_j_seqs.index(hold_j1[i][0])])]):
                                if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                    j_match = half1_j_seqs.index(hold_j1[i][0])
                                    # Finds where the start of a full J would be
                                    temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match]
                                    found_j_match += 1
                    for i in range(len(hold_j2)):
                        indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                        for k in indices:
                            if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
                                if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[k])] ) <= 1:
                                    j_match = half2_j_seqs.index(hold_j2[i][0])
                                    # Finds where the start of a full J would be
                                    temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - 6
                                    found_j_match += 1

                if hold_v and hold_j:
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record_reverse.seq[end_v+1:start_j]))
                        # Write to analysis_file (text file) the classification of the sequence
                        print f_seq
                        assigned_count += 1
                        found_seq_match = 1
                elif hold_v and found_j_match == 1:
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record_reverse.seq[end_v+1:start_j]))
                        print f_seq
                        assigned_count += 1
                        found_seq_match = 1
                elif found_v_match == 1 and hold_j:
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record_reverse.seq[end_v+1:start_j]))
                        print f_seq
                        assigned_count += 1
                        found_seq_match = 1
                elif found_v_match == 1 and found_j_match == 1:
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record_reverse.seq[end_v+1:start_j]))
                        print f_seq
                        assigned_count += 1
                        found_seq_match = 1
    if verbose:
        t = time() - t0
        sys.stderr.write('%d sequences were analysed\n' % seq_count)
        sys.stderr.write('%d sequences were successfully assigned\n' % assigned_count)
        sys.stderr.write('%s seconds elapsed\n' % t)

def get_v_deletions( rc, v_match, temp_end_v, v_regions_cut ):
    # This function determines the number of V deletions in sequence rc
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in rc.
    function_temp_end_v = temp_end_v
    pos = -1
    is_v_match = 0
    while is_v_match == 0 and 0 <= function_temp_end_v < len(rc):
        if str(v_regions_cut[v_match])[pos] == str(rc)[function_temp_end_v] \
                and str(v_regions_cut[v_match])[pos-1] == str(rc)[function_temp_end_v-1] \
                and str(v_regions_cut[v_match])[pos-2] == str(rc)[function_temp_end_v-2]:
            is_v_match = 1
            deletions_v = -pos - 1
            end_v = function_temp_end_v
        else:
            pos -= 1
            function_temp_end_v -= 1

    if is_v_match == 1:
        return [end_v, deletions_v]
    else:
        return []

def get_j_deletions( rc, j_match, temp_start_j, j_regions_cut ):
    # This function determines the number of J deletions in sequence rc
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in rc.
    function_temp_start_j = temp_start_j
    pos = 0
    is_j_match = 0
    while is_j_match == 0 and 0 <= function_temp_start_j+2 < len(str(rc)):
        if str(j_regions_cut[j_match])[pos] == str(rc)[function_temp_start_j] and \
                str(j_regions_cut[j_match])[pos+1] == str(rc)[function_temp_start_j+1] and \
                str(j_regions_cut[j_match])[pos+2] == str(rc)[function_temp_start_j+2]:
            is_j_match = 1
            deletions_j = pos
            start_j = function_temp_start_j
        else:
            pos += 1
            function_temp_start_j += 1
            
    if is_j_match == 1:
        return [start_j, deletions_j]
    else:
        return []

def get_v_tags(file_v, half_split):
    import string
    
    v_seqs = []
    jump_to_end_v = []
    for line in file_v:
        elements = line.rstrip("\n")
        v_seqs.append(string.split(elements)[0])
        jump_to_end_v.append(int(string.split(elements)[1]))

    half1_v_seqs = []
    half2_v_seqs = []

    for i in range(len(v_seqs)):
        half1_v_seqs.append(v_seqs[i][0:half_split])
        half2_v_seqs.append(v_seqs[i][half_split:])
    
    return [v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v]

def get_j_tags(file_j, half_split):
    import string
    
    j_seqs = []
    jump_to_start_j = []

    for line in file_j:
        elements = line.rstrip("\n")
        j_seqs.append(string.split(elements)[0])
        jump_to_start_j.append(int(string.split(elements)[1]))

    half1_j_seqs = []
    half2_j_seqs = []

    for j in range(len(j_seqs)):
        half1_j_seqs.append(j_seqs[j][0:half_split])
        half2_j_seqs.append(j_seqs[j][half_split:])

    return [j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j]

def main(args):
    analysis(args.fastqs, args.vfasta, args.jfasta, args.rev_comp, args.verbose)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastqs", metavar="FASTQ", nargs="+",
            help="single-end fastqs(s)")
    req = p.add_argument_group("required arguments")
    req.add_argument("-vf", "--v-fasta", dest="vfasta", required=True,
            help="v region fasta")
    req.add_argument("-jf", "--j-fasta", dest="jfasta", required=True,
            help="j region fasta")
    p.add_argument("--reverse-complement-search", dest="rev_comp",
            action="store_true", help="if not aligned, attempt to align the \
                                       reverse complement")
    p.add_argument("-v", "--verbose", action="store_true",
            help="prints alignment stats")
    main(p.parse_args())