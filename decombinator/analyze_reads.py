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
import Levenshtein as lev
from Bio import SeqIO
from time import time
from toolshed import nopen, reader
from string import Template
from acora import AcoraBuilder

def analysis(fastqs, vfasta, jfasta, vtags, jtags, rev_comp=False,
                verbose=False, sep=" "):
    if verbose:
        sys.stderr.write('>> Analyzing %d file(s)\n' % len(fastqs))
        sys.stderr.write(">> Importing known V, and J gene segments and tags\n")

    # get the sequences per region
    v_genes = list(SeqIO.parse(nopen(vfasta), "fasta"))
    j_genes = list(SeqIO.parse(nopen(jfasta), "fasta"))
    # XXX
    # classes to parse fasta, fastq, and method to reverse complement
    # get rid of biopython
    v_regions = [str(v_genes[i].seq.upper()) for i, v in enumerate(v_genes)]
    j_regions = [str(j_genes[i].seq.upper()) for i, v in enumerate(j_genes)]

    v_seqs, vleft_seqs, vright_seqs, v_ends = get_tags(vtags)
    j_seqs, jleft_seqs, jright_seqs, j_starts = get_tags(jtags)

    # full sequences
    builder = AcoraBuilder(v_seqs)
    v_key = builder.build()
    builder = AcoraBuilder(j_seqs)
    j_key = builder.build()
    
    # half sequences
    builder = AcoraBuilder(vleft_seqs)
    vleft_key = builder.build()
    builder = AcoraBuilder(vright_seqs)
    vright_key = builder.build()
    builder = AcoraBuilder(jleft_seqs)
    jleft_key = builder.build()
    builder = AcoraBuilder(jright_seqs)
    jright_key = builder.build()
    
    # correctly assigned sequences
    assigned_count = 0
    # number of sequences analysed
    seq_count = 0
    # begin clock
    t0 = time()
    
    # XXX
    stemplate = Template('$v $j $del_v $del_j $nt_insert')

    for fastq in fastqs:
        if verbose:
            sys.stderr.write(">> Starting %s...\n" % fastq)
        for i, record in enumerate(SeqIO.parse(nopen(fastq), "fastq")):
            # if i == 50:
            #     sys.exit()
            found_seq_match = 0
            seq_count += 1
            hold_v = v_key.findall(str(record.seq))
            hold_j = j_key.findall(str(record.seq))

            if hold_v:
                # the index position of the found sequence among known (v_seqs)
                v_match = v_seqs.index(hold_v[0][0])
                
                # new variable names
                # do not like lists for this task
                match_idx = v_seqs.index(hold_v[0][0])
                match_start_idx = hold_v[0][1]
                vseq_end = v_ends[match_idx] - 1
                end_of_v = match_start_idx + vseq_end
                
                # Finds where the end of a full V would be
                temp_end_v = hold_v[0][1] + v_ends[v_match] - 1
                
                # If the number of deletions has been found
                if get_v_deletions(record.seq, v_match, temp_end_v, v_regions):
                    end_v, deletions_v = get_v_deletions(record.seq, v_match, temp_end_v, v_regions)
            else:
                found_v_match = 0
                hold_v1 = vleft_key.findall(str(record.seq))
                hold_v2 = vright_key.findall(str(record.seq))
                for i in range(len(hold_v1)):
                    indices = [y for y, x in enumerate(vleft_seqs) if x == hold_v1[i][0] ]
                    for k in indices:
                        if len(v_seqs[k]) == len(str(record.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[vleft_seqs.index(hold_v1[i][0])])]):
                            if lev.hamming( v_seqs[k], str(record.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                                v_match = k
                                # Finds where the end of a full V would be
                                temp_end_v = hold_v1[i][1] + v_ends[v_match] - 1
                                found_v_match += 1
                for i in range(len(hold_v2)):
                    indices = [y for y, x in enumerate(vright_seqs) if x == hold_v2[i][0] ]
                    for k in indices:
                        if len(v_seqs[k]) == len(str(record.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[vright_seqs.index(hold_v2[i][0])])]):
                            if lev.hamming( v_seqs[k], str(record.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[k])] ) <= 1:
                                v_match = k
                                # Finds where the end of a full V would be
                                temp_end_v = hold_v2[i][1] + v_ends[v_match] - 1
                                found_v_match += 1

            if hold_j:
                # Assigns J
                j_match = j_seqs.index(hold_j[0][0])
                # Finds where the start of a full J would be
                temp_start_j = hold_j[0][1] - j_starts[j_match]
                if get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    [ start_j, deletions_j] = get_j_deletions( record.seq, j_match, temp_start_j, j_regions )
            else:
                found_j_match = 0
                hold_j1 = jleft_key.findall(str(record.seq))
                hold_j2 = jright_key.findall(str(record.seq))
                for i in range(len(hold_j1)):
                    indices = [y for y, x in enumerate(jleft_seqs) if x == hold_j1[i][0] ]
                    for k in indices:
                        if len(j_seqs[k]) == len(str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[jleft_seqs.index(hold_j1[i][0])])]):
                            if lev.hamming( j_seqs[k], str(record.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                j_match = jleft_seqs.index(hold_j1[i][0])
                                # Finds where the start of a full J would be
                                temp_start_j = hold_j1[i][1] - j_starts[j_match]
                                found_j_match += 1
                for i in range(len(hold_j2)):
                    indices = [y for y, x in enumerate(jright_seqs) if x == hold_j2[i][0] ]
                    for k in indices:
                        if len(j_seqs[k]) == len(str(record.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[jright_seqs.index(hold_j2[i][0])])]):
                            if lev.hamming( j_seqs[k], str(record.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[k])] ) <= 1:
                                j_match = jright_seqs.index(hold_j2[i][0])
                                # Finds where the start of a full J would be
                                temp_start_j = hold_j2[i][1] - j_starts[j_match] - 6
                                found_j_match += 1

            if hold_v and hold_j:
                f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = record.seq[end_v+1:start_j])
                # Write to analysis_file (text file) the classification of the sequence
                print f_seq
                assigned_count += 1
                found_seq_match = 1
            elif hold_v and found_j_match == 1:
                f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                print f_seq
                assigned_count += 1
                found_seq_match = 1
            elif found_v_match == 1 and hold_j:
                f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                print f_seq
                assigned_count += 1
                found_seq_match = 1
            elif found_v_match == 1 and found_j_match == 1:
                f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                print f_seq
                assigned_count += 1
                found_seq_match = 1
            
            #####################
            # REVERSE COMPLEMENT
            #####################
            if found_seq_match == 0 and rev_comp:

                record_reverse = record.reverse_complement()
                hold_v = v_key.findall(str(record_reverse.seq))
                hold_j = j_key.findall(str(record_reverse.seq))

                if hold_v:
                    # Assigns V
                    v_match = v_seqs.index(hold_v[0][0])
                    # Finds where the end of a full V would be
                    temp_end_v = hold_v[0][1] + v_ends[v_match] - 1
                    # If the number of deletions has been found
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ):
                        end_v, deletions_v = get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions )
                else:
                    found_v_match = 0
                    hold_v1 = vleft_key.findall(str(record_reverse.seq))
                    hold_v2 = vright_key.findall(str(record_reverse.seq))
                    for i in range(len(hold_v1)):
                        indices = [y for y, x in enumerate(vleft_seqs) if x == hold_v1[i][0] ]
                        for k in indices:
                            if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[vleft_seqs.index(hold_v1[i][0])])]):
                                if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v1[i][1]:hold_v1[i][1]+len(v_seqs[k])] ) <= 1:
                                    v_match = k
                                    # Finds where the end of a full V would be
                                    temp_end_v = hold_v1[i][1] + v_ends[v_match] - 1
                                    found_v_match += 1
                    for i in range(len(hold_v2)):
                        indices = [y for y, x in enumerate(vright_seqs) if x == hold_v2[i][0] ]
                        for k in indices:
                            if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[vright_seqs.index(hold_v2[i][0])])]):
                                if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[k])] ) <= 1:
                                    v_match = k
                                    # Finds where the end of a full V would be
                                    temp_end_v = hold_v2[i][1] + v_ends[v_match] - 1
                                    found_v_match += 1

                if hold_j:
                    # Assigns J
                    j_match = j_seqs.index(hold_j[0][0])
                    # Finds where the start of a full J would be
                    temp_start_j = hold_j[0][1] - j_starts[j_match]
                    if get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        start_j, deletions_j = get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions )
                else:
                    found_j_match = 0
                    hold_j1 = jleft_key.findall(str(record_reverse.seq))
                    hold_j2 = jright_key.findall(str(record_reverse.seq))
                    for i in range(len(hold_j1)):
                        indices = [y for y, x in enumerate(jleft_seqs) if x == hold_j1[i][0] ]
                        for k in indices:
                            if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[jleft_seqs.index(hold_j1[i][0])])]):
                                if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j1[i][1]:hold_j1[i][1]+len(j_seqs[k])] ) <= 1:
                                    j_match = jleft_seqs.index(hold_j1[i][0])
                                    # Finds where the start of a full J would be
                                    temp_start_j = hold_j1[i][1] - j_starts[j_match]
                                    found_j_match += 1
                    for i in range(len(hold_j2)):
                        indices = [y for y, x in enumerate(jright_seqs) if x == hold_j2[i][0] ]
                        for k in indices:
                            if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[jright_seqs.index(hold_j2[i][0])])]):
                                if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[k])] ) <= 1:
                                    j_match = jright_seqs.index(hold_j2[i][0])
                                    # Finds where the start of a full J would be
                                    temp_start_j = hold_j2[i][1] - j_starts[j_match] - 6
                                    found_j_match += 1

                if (hold_v and hold_j) or \
                        (hold_v and found_j_match == 1) or \
                        (found_v_match == 1 and hold_j) or \
                        (found_v_match == 1 and found_j_match == 1):
                    
                    f_seq = stemplate.substitute(v = v_match, j = j_match, 
                                del_v = deletions_v, del_j = deletions_j, 
                                nt_insert = str(record_reverse.seq[end_v + 1:start_j]))
                    fields = (v_match, j_match, deletions_v, deletions_j,
                                record_reverse.seq[end_v + 1:start_j])
                    assigned_count += 1
                    found_seq_match = 1
                    print sep.join(map(str, fields))
    if verbose:
        t = time() - t0
        sys.stderr.write('%d sequences were analysed\n' % seq_count)
        sys.stderr.write('%d sequences were successfully assigned\n' % assigned_count)
        sys.stderr.write('%s seconds elapsed\n' % t)

def get_v_deletions(seq, v_match, end_v, v_regions):
    # This function determines the number of V deletions in sequence seq
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in seq.
    seq = str(seq)
    pos = -1
    is_v_match = 0
    while is_v_match == 0 and 0 <= end_v < len(seq):
        if v_regions[v_match][pos] == seq[end_v] \
                and v_regions[v_match][pos - 1] == seq[end_v - 1] \
                and v_regions[v_match][pos - 2] == seq[end_v - 2]:
            is_v_match = 1
            deletions_v = -pos - 1
        else:
            pos -= 1
            end_v -= 1
    
    if is_v_match == 1:
        return end_v, deletions_v
    else:
        return

def get_j_deletions(seq, j_match, start_j, j_regions):
    # This function determines the number of J deletions in sequence seq
    # by comparing it to j_match, beginning by making comparisons at the
    # end of j_match and at position temp_end_j in seq.
    seq = str(seq)
    pos = 0
    is_j_match = 0
    while is_j_match == 0 and 0 <= start_j + 2 < len(seq):
        if j_regions[j_match][pos] == seq[start_j] and \
                j_regions[j_match][pos + 1] == seq[start_j + 1] and \
                j_regions[j_match][pos + 2] == seq[start_j + 2]:
            is_j_match = 1
            deletions_j = pos
        else:
            pos += 1
            start_j += 1
    
    if is_j_match == 1:
        return start_j, deletions_j
    else:
        return

def get_tags(tags):
    # tag file format is seq, length of something, name
    seqs = []
    ends = []
    # the halves of the sequence
    left = []
    right = []
    for l in reader(tags, header="seq end name".split(), sep=" "):
        seqs.append(l['seq'])
        ends.append(int(l['end']))
        left.append(l['seq'][:len(l['seq']) / 2])
        right.append(l['seq'][len(l['seq']) / 2:])
    return [seqs, left, right, ends]

def main(args):
    analysis(args.fastqs, args.vfasta, args.jfasta, args.vtags, args.jtags, 
                args.rev_comp, args.verbose, args.delimiter)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vfasta", help="V region fasta")
    p.add_argument("jfasta", help="J region fasta")
    p.add_argument("vtags", help="V tags")
    p.add_argument("jtags", help="J tags")
    p.add_argument("fastqs", metavar="FASTQ", nargs="+",
            help="single-end fastqs(s)")
    p.add_argument("-d", "--delimiter", default=" ", help="output delimiter")
    p.add_argument("--reverse-complement-search", dest="rev_comp",
            action="store_true", help="if not aligned, attempt to align the \
                                       reverse complement")
    p.add_argument("-v", "--verbose", action="store_true",
            help="prints alignment stats")
    main(p.parse_args())