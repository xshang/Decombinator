def analysis( Sequence_Reads, with_statistics=True, with_reverse_complement_search=True):
    import numpy as np
    import decimal as dec
    import string
    import operator as op
    import collections as coll
    from Bio import SeqIO
    from acora import AcoraBuilder
    from time import time, clock
    from string import Template
    from operator import itemgetter, attrgetter
    import Levenshtein as lev

    v_half_split, j_half_split = [10,6] # Do not change - V tags are split at position 10, J at position 6, to look for half tags if no full tag is found.

    ################

    print 'Commencing analysis on a total of', len(Sequence_Reads), 'file(s)'

    ## Create .txt file to store f=(v_index,j_index,v_deletions,j_deletions,nt_insert)
    analysis_file = open("DecombinatorResults.txt", "w")
    analysis_file.close()
    results = "DecombinatorResults.txt" # Name the .txt file to write to

    ################
    print ('Importing known V, D and J gene segments and tags...')

    handle = open("human_TRBV_region.fasta", "rU")
    v_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    handle = open("human_TRBJ_region.fasta", "rU")
    j_genes = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    v_regions = []
    for j in range(0, len(v_genes)):
        v_regions.append(string.upper(v_genes[j].seq))

    j_regions = []
    for j in range(0, len(j_genes)):
        j_regions.append(string.upper(j_genes[j].seq))

    ##############
    ## Build keyword tries of V and J tags for fast assignment
    v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v = get_v_tags(open("tags_trbv.txt", "rU"), v_half_split)
    j_seqs, half1_j_seqs, half2_j_seqs, jump_to_start_j = get_j_tags(open("tags_trbj.txt", "rU"), j_half_split)   

    v_builder = AcoraBuilder()
    for i in range(0,len(v_seqs)):
        v_builder.add(str(v_seqs[i])) # Add all V tags to keyword trie

    v_key = v_builder.build()

    j_builder = AcoraBuilder()
    for i in range(0,len(j_seqs)):
        j_builder.add(str(j_seqs[i])) # Add all J tags to keyword trie

    j_key = j_builder.build()

    ##############
    ## Build keyword tries for first and second halves of both V and J tags
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

    ###############
    ## Initialise variables
    assigned_count = 0 # this will just increase by one every time we correctly assign a seq read with all desired variables
    seq_count = 0 # this will simply track the number of sequences analysed in file
    t0 = time() # Begin timer

    ###############
    ## Open .txt file created at the start of analysis
    analysis_file = open(results, "a")
    stemplate = Template('$v $j $del_v $del_j $nt_insert') # Creates stemplate, a holder, for f. Each line will have the 5 variables separated by a space

    ###############
    ## Begin analysing sequences

    for i in range(len(Sequence_Reads)):
        
        print 'Importing sequences from', Sequence_Reads[i],' and assigning V and J regions...'
        handle = open(Sequence_Reads[i], "rU")
        
        for record in SeqIO.parse(handle, "fastq"):
            
            found_seq_match = 0
            seq_count += 1
            
            hold_v = v_key.findall(str(record.seq))
            hold_j = j_key.findall(str(record.seq))

            if hold_v:                
                v_match = v_seqs.index(hold_v[0][0]) # Assigns V
                temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ): # If the number of deletions has been found
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
                                temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                                found_v_match += 1
                for i in range(len(hold_v2)):
                    indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                    for k in indices:
                        if len(v_seqs[k]) == len(str(record.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
                            if lev.hamming( v_seqs[k], str(record.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[k])] ) <= 1:
                                v_match = k
                                temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                                found_v_match += 1

            if hold_j:
                j_match = j_seqs.index(hold_j[0][0]) # Assigns J
                temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
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
                                temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                                found_j_match += 1
                for i in range(len(hold_j2)):
                    indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                    for k in indices:
                        if len(j_seqs[k]) == len(str(record.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
                            if lev.hamming( j_seqs[k], str(record.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[k])] ) <= 1:
                                j_match = half2_j_seqs.index(hold_j2[i][0])
                                temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - 6 # Finds where the start of a full J would be
                                found_j_match += 1

            if hold_v and hold_j:
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                    print >> analysis_file, f_seq # Write to analysis_file (text file) the classification of the sequence
                    assigned_count += 1
                    found_seq_match = 1
            elif hold_v and found_j_match == 1:
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                    print >> analysis_file, f_seq
                    assigned_count += 1
                    found_seq_match = 1
            elif found_v_match == 1 and hold_j:
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                    print >> analysis_file, f_seq
                    assigned_count += 1
                    found_seq_match = 1
            elif found_v_match == 1 and found_j_match == 1:
                if get_v_deletions( record.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record.seq, j_match, temp_start_j, j_regions ):
                    f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record.seq[end_v+1:start_j]))
                    print >> analysis_file, f_seq
                    assigned_count += 1
                    found_seq_match = 1

            if found_seq_match == 0 and with_reverse_complement_search == True:
                
                #####################
                # REVERSE COMPLEMENT
                #####################

                record_reverse = record.reverse_complement()
                hold_v = v_key.findall(str(record_reverse.seq))
                hold_j = j_key.findall(str(record_reverse.seq))

                if hold_v:                
                    v_match = v_seqs.index(hold_v[0][0]) # Assigns V
                    temp_end_v = hold_v[0][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ): # If the number of deletions has been found
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
                                    temp_end_v = hold_v1[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                                    found_v_match += 1
                    for i in range(len(hold_v2)):
                        indices = [y for y, x in enumerate(half2_v_seqs) if x == hold_v2[i][0] ]
                        for k in indices:
                            if len(v_seqs[k]) == len(str(record_reverse.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[half2_v_seqs.index(hold_v2[i][0])])]):
                                if lev.hamming( v_seqs[k], str(record_reverse.seq)[hold_v2[i][1]:hold_v2[i][1]+len(v_seqs[k])] ) <= 1:
                                    v_match = k
                                    temp_end_v = hold_v2[i][1] + jump_to_end_v[v_match] - 1 # Finds where the end of a full V would be
                                    found_v_match += 1

                if hold_j:
                    j_match = j_seqs.index(hold_j[0][0]) # Assigns J
                    temp_start_j = hold_j[0][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
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
                                    temp_start_j = hold_j1[i][1] - jump_to_start_j[j_match] # Finds where the start of a full J would be
                                    found_j_match += 1
                    for i in range(len(hold_j2)):
                        indices = [y for y, x in enumerate(half2_j_seqs) if x == hold_j2[i][0] ]
                        for k in indices:
                            if len(j_seqs[k]) == len(str(record_reverse.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[half2_j_seqs.index(hold_j2[i][0])])]):
                                if lev.hamming( j_seqs[k], str(record_reverse.seq)[hold_j2[i][1]:hold_j2[i][1]+len(j_seqs[k])] ) <= 1:
                                    j_match = half2_j_seqs.index(hold_j2[i][0])
                                    temp_start_j = hold_j2[i][1] - jump_to_start_j[j_match] - 6 # Finds where the start of a full J would be
                                    found_j_match += 1

                if hold_v and hold_j:
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record_reverse.seq[end_v+1:start_j]))
                        print >> analysis_file, f_seq # Write to analysis_file (text file) the classification of the sequence
                        assigned_count += 1
                        found_seq_match = 1
                elif hold_v and found_j_match == 1:
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record_reverse.seq[end_v+1:start_j]))
                        print >> analysis_file, f_seq
                        assigned_count += 1
                        found_seq_match = 1
                elif found_v_match == 1 and hold_j:
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record_reverse.seq[end_v+1:start_j]))
                        print >> analysis_file, f_seq
                        assigned_count += 1
                        found_seq_match = 1
                elif found_v_match == 1 and found_j_match == 1:
                    if get_v_deletions( record_reverse.seq, v_match, temp_end_v, v_regions ) and get_j_deletions( record_reverse.seq, j_match, temp_start_j, j_regions ):
                        f_seq = stemplate.substitute( v = v_match, j = j_match, del_v = deletions_v, del_j = deletions_j, nt_insert = str(record_reverse.seq[end_v+1:start_j]))
                        print >> analysis_file, f_seq
                        assigned_count += 1
                        found_seq_match = 1
        handle.close()
    analysis_file.close()

    if with_statistics == True:
        timed = time() - t0
        print seq_count, 'sequences were analysed'
        print assigned_count, ' sequences were successfully assigned'
        print 'Time taken =', timed, 'seconds'

def get_v_deletions( rc, v_match, temp_end_v, v_regions_cut ):
    # This function determines the number of V deletions in sequence rc
    # by comparing it to v_match, beginning by making comparisons at the
    # end of v_match and at position temp_end_v in rc.
    function_temp_end_v = temp_end_v
    pos = -1
    is_v_match = 0
    while is_v_match == 0 and 0 <= function_temp_end_v < len(rc):
        if str(v_regions_cut[v_match])[pos] == str(rc)[function_temp_end_v] and str(v_regions_cut[v_match])[pos-1] == str(rc)[function_temp_end_v-1] and str(v_regions_cut[v_match])[pos-2] == str(rc)[function_temp_end_v-2]:
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
        if str(j_regions_cut[j_match])[pos] == str(rc)[function_temp_start_j] and str(j_regions_cut[j_match])[pos+1] == str(rc)[function_temp_start_j+1] and str(j_regions_cut[j_match])[pos+2] == str(rc)[function_temp_start_j+2]:
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
    
    v_seqs = [] # Holds all V tags
    jump_to_end_v = [] # Holds the number of jumps to make to look for deletions for each V region once the corresponding tag has been found
    for line in file_v:
        elements = line.rstrip("\n") # Turns every element in a text file line separated by a space into elements in a list
        v_seqs.append(string.split(elements)[0]) # Adds elements in first column iteratively
        jump_to_end_v.append(int(string.split(elements)[1])) # Adds elements in second column iteratively

    half1_v_seqs = []
    half2_v_seqs = []

    for i in range(len(v_seqs)):
        half1_v_seqs.append(v_seqs[i][0:half_split])
        half2_v_seqs.append(v_seqs[i][half_split:])
    
    return [v_seqs, half1_v_seqs, half2_v_seqs, jump_to_end_v]

def get_j_tags(file_j, half_split):
    import string
    
    j_seqs = [] # Holds all J tags
    jump_to_start_j = [] # Holds the number of jumps to make to look for deletions for each J region once the corresponding tag has been found

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

def get_distinct_clones( handle_results, with_count=False ):

    ## LOOKS THROUGH TEXT FILE OF CLASSIFIERS AND WRITES NEW FILE CONTAINING ALL DISTINCT CLASSIFIERS, OPTIONALLY WITH COUNT OF ALL DISTINCT CLASSIFIERS
    ## with_count=True writes file with counts of all distinct classifiers
    
    from string import Template
    import collections as coll
    from operator import itemgetter, attrgetter

    write_to = open("distinct_clones.txt", "w")

    if with_count == True:
        stemplate = Template('$count $element')
        d = coll.defaultdict(int)
        for line in handle_results:
            elements = line.rstrip("\n")
            d[elements] += 1
        d_sorted = sorted(d.items(), key=itemgetter(1), reverse=True)
        for k in d_sorted:
            f_seq = stemplate.substitute( count = k[1], element = k[0])
            print >> write_to, f_seq
    else:
        stemplate = Template('$element')
        d = coll.defaultdict(int)
        for line in handle_results:
            elements = line.rstrip("\n")
            if elements not in d:
                d[elements] = 1
                f_seq = stemplate.substitute(element = elements)
                print >> write_to, f_seq
    
    handle_results.close()
    write_to.close()

def get_translated_sequences( handle_results, chain="beta", with_outframe=False, fullaaseq=False, handle_vb=open("human_TRBV_region.fasta","rU"), handle_jb=open("human_TRBJ_region.fasta","rU"), handle_va=open("human_TRAV_region.fasta","rU"), handle_ja=open("human_TRAJ_region.fasta","rU") ):

    ## TRANSLATES CLASSIFIERS TO AA SEQUENCES VIA THEIR NT SEQUENCE
    ## Default settings are -
    ## chain = "beta" or chain = "alpha"
    ## with_outframe=True or False: writes all aa seqeunces to file, including those that are out-of-frame (with stop codon symbol *)
    ## fullaaseq=True or False: True writes the whole V(D)J aa sequence to file, False, writes only the CDR3 region.

    from Bio.Seq import Seq
    from Bio import SeqIO
    from Bio.Alphabet import generic_dna
    import string
    import re
    
    vb_raw = list(SeqIO.parse(handle_vb, "fasta"))
    handle_vb.close()
    jb_raw = list(SeqIO.parse(handle_jb, "fasta"))
    handle_jb.close()
    va_raw = list(SeqIO.parse(handle_va, "fasta"))
    handle_va.close()
    ja_raw = list(SeqIO.parse(handle_ja, "fasta"))
    handle_ja.close()

    vb_regions = []
    for i in range(0,len(vb_raw)):
        vb_regions.append(string.upper(vb_raw[i].seq))

    jb_regions = []
    for i in range(0,len(jb_raw)):
        jb_regions.append(string.upper(jb_raw[i].seq))

    va_regions = []
    for i in range(0,len(va_raw)):
        va_regions.append(string.upper(va_raw[i].seq))

    ja_regions = []
    for i in range(0,len(ja_raw)):
        ja_regions.append(string.upper(ja_raw[i].seq))

    write_to = open("translated_sequences.txt", "w")

    if chain == "beta":
        for line in handle_results:
            elements = line.rstrip("\n")
            classifier = string.split(elements)

            v = int(classifier[0])
            j = int(classifier[1])
            delv = int(classifier[2])
            delj = int(classifier[3])
            if len(classifier) == 5:
                ins = str(classifier[4])
            elif len(classifier) == 4:
                ins = ''

            if delv != 0:
                used_v = vb_regions[v][:-delv]
            elif delv == 0:
                used_v = vb_regions[v]

            if delj != 0:
                used_j = jb_regions[j][delj:]
            elif delj == 0:
                used_j = jb_regions[j]

            seq = str(used_v + ins + used_j)
            start = len(seq)%3
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, str(aaseq)
            else:     
                if re.findall('FG.G',str(aaseq)) and re.findall('C',str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G',str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[lower:upper+4]
                    if with_outframe == True:
                        print >> write_to, cdr3
                    elif '*' not in cdr3:
                        print >> write_to, cdr3

    if chain == "alpha":
        for line in handle_results:
            elements = line.rstrip("\n")
            classifier = string.split(elements)

            v = int(classifier[0])
            j = int(classifier[1])
            delv = int(classifier[2])
            delj = int(classifier[3])
            if len(classifier) == 5:
                ins = str(classifier[4])
            elif len(classifier) == 4:
                ins = ''

            if delv != 0:
                used_v = va_regions[v][:-delv]
            elif delv == 0:
                used_v = va_regions[v]

            if delj != 0:
                used_j = ja_regions[j][delj:]
            elif delj == 0:
                used_j = ja_regions[j]

            seq = str(used_v + ins + used_j)
            start = (len(seq)-1)%3
            aaseq = Seq(str(seq[start:]), generic_dna).translate()

            if fullaaseq == True:
                if with_outframe == True:
                    print >> write_to, str(aaseq)
                elif '*' not in aaseq:
                    print >> write_to, str(aaseq)
            else:     
                if re.findall('FG.G',str(aaseq))and re.findall('C',str(aaseq)):
                    indices = [i for i, x in enumerate(aaseq) if x == 'C']
                    upper = str(aaseq).find(re.findall('FG.G',str(aaseq))[0])
                    for i in indices:
                        if i < upper:
                            lower = i
                    cdr3 = aaseq[indices[lower]:upper+4]
                    if with_outframe == True:
                        print >> write_to, cdr3
                    elif '*' not in cdr3:
                        print >> write_to, cdr3
            
    handle_results.close()
    write_to.close()

def plot_v_usage( handle, savefilename="Vusage", tags=open("tags_trbv.txt", "rU"), order="frequency"):

    ## PLOTS V GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    num_genes = 0
    for line in tags:
        num_genes += 1

    freq_vector_v = [0]*num_genes
    for line in handle:
        elements = line.rstrip("\n")
        freq_vector_v[int(string.split(elements)[0])] += 1

    if order=="frequency":        
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_v)
        percent_usage_v = [0]*num_genes
        for i in range(num_genes):
            percent_usage_v[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
        gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/V12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27-1', 'V28-1', 'V29-1', 'V3-1', 'V30-1', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9')
        v_linked = [0]*len(percent_usage_v)
        for i in range(len(percent_usage_v)):
            v_linked[i] = (gene_list_v[i], percent_usage_v[i])
        sorted_v = sorted(v_linked, key=itemgetter(1))
        v_labels = [0]*len(sorted_v)
        v_percents = [0]*len(sorted_v)
        for j in range(len(sorted_v)):
            v_labels[j] = sorted_v[j][0]
            v_percents[j] = sorted_v[j][1]
        pos_v = np.arange(num_genes)+ 1
        plt.figure()
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.yticks( pos_v, v_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_v, v_percents, align = 'center', color = 'yellow', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)
        
    elif order=="chromosome":
        total = sum(freq_vector_v)
        fv = [0]*num_genes
        for i in range(num_genes):
            fv[i] = dec.Decimal(freq_vector_v[i]) / dec.Decimal(total)
        gene_list_v = ('10-1', '10-2', '10-3', '11-1', '11-2', '11-3', '12-3/4', '12-5', '13', '14', '15', '16', '18', '19', '2', '20-1', '24-1', '25-1', '27-1', '28-1', '29-1', '3-1', '30-1', '4-1', '4-2', '4-3', '5-1', '5-4', '5-5', '5-6', '5-8', '6-1', '6-4', '6-5', '6-6', '6-8', '6-9', '7-2', '7-3', '7-4', '7-6', '7-7', '7-8', '7-9', '9')
        chromosome_order = [ 14, 21, 23, 26, 31, 24, 25, 37, 32, 38, 44, 0, 3, 1, 4, 33, 39, 27, 34, 28, 40, 29, 35, 41, 36, 42, 30, 43, 8, 2, 5, 6, 7, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 22 ]
        gene_list_v = [ gene_list_v[i] for i in chromosome_order]
        fv = [fv[i] for i in chromosome_order]

        ind = np.arange(num_genes)
        width = 0.25

        fig = plt.figure()
        ax = fig.add_subplot(111)
        rects1 = ax.bar(ind, fv, width, color='yellow')

        ax.set_ylabel('Frequency', fontsize = 10)
        ax.set_xticks(ind+1*width)
        ax.set_xticklabels(gene_list_v)
        plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 10)
        plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 10)
        plt.grid(True)

        plt.savefig(str(savefilename)+'.png', dpi=300)

def plot_j_usage( handle, savefilename="Jusage", tags = open("tags_trbj.txt", "rU"), order="frequency"):

    ## PLOTS J GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec
    from operator import itemgetter, attrgetter

    num_genes = 0
    for line in tags:
        num_genes += 1

    freq_vector_j = [0]*num_genes
    for line in handle:
        elements = line.rstrip("\n")
        freq_vector_j[int(string.split(elements)[1])] += 1

    if order=="frequency":    
        plt.rcParams['figure.figsize'] = 10,10
        total = sum(freq_vector_j)
        percent_usage_j = [0]*num_genes
        for i in range(num_genes):
            percent_usage_j[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
        gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')
        j_linked = [0]*len(percent_usage_j)
        for i in range(len(percent_usage_j)):
            j_linked[i] = (gene_list_j[i], percent_usage_j[i])
        sorted_j = sorted(j_linked, key=itemgetter(1))
        j_labels = [0]*len(sorted_j)
        j_percents = [0]*len(sorted_j)
        for j in range(len(sorted_j)):
            j_labels[j] = sorted_j[j][0]
            j_percents[j] = sorted_j[j][1]
        pos_j = np.arange(num_genes)+ 1
        plt.figure()
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.yticks( pos_j, j_labels)
        plt.xlabel('Frequency Usage')
        plt.barh( pos_j, j_percents, align = 'center', color = 'red', height=0.2)
        plt.grid(True)
        plt.savefig(str(savefilename)+'.png', dpi=300)

    elif order=="chromosome":
        total = sum(freq_vector_j)
        fj = [0]*num_genes
        for i in range(num_genes):
            fj[i] = dec.Decimal(freq_vector_j[i]) / dec.Decimal(total)
        gene_list_j = ('1-1', '1-2', '1-3', '1-4', '1-5', '1-6', '2-1', '2-2', '2-3', '2-4', '2-5', '2-6', '2-7')

        ind = np.arange(num_genes)
        width = 0.25

        fig = plt.figure()
        ax = fig.add_subplot(111)
        rects1 = ax.bar(ind, fj, width, color='red')

        ax.set_ylabel('Frequency', fontsize = 10)
        ax.set_xticks(ind+width)
        ax.set_xticklabels(gene_list_j)
        plt.setp(ax.get_xticklabels(),rotation='vertical',fontsize = 10)
        plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 10)
        plt.grid(True)

        plt.savefig(str(savefilename)+'.png', dpi=300)

def plot_del_v( handle, savefilename="Vdels"):

    ## PLOTS V GERMLINE DELETIONS BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_v = [0]*50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_v[int(string.split(elements)[2])] += 1

    total = sum(deletions_v)
    for i in range(len(deletions_v)):
        deletions_v[i] = deletions_v[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_v, width, color='yellow')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of V germline deletions', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.ylim((0,0.2))
    plt.xlim((0,20))

    plt.savefig(str(savefilename)+'.png', dpi=300)


def plot_del_j( handle, savefilename="Jdels"):

    ## PLOTS J GERMLINE DELETIONS BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    deletions_j = [0]*50
    for line in handle:
        elements = line.rstrip("\n")
        deletions_j[int(string.split(elements)[3])] += 1

    total = sum(deletions_j)
    for i in range(len(deletions_j)):
        deletions_j[i] = deletions_j[i] / float(total)

    ind = np.arange(50)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, deletions_j, width, color='red')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of J germline deletions', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.ylim((0,0.2))
    plt.xlim((0,20))

    plt.savefig(str(savefilename)+'.png', dpi=300)

def plot_vj_joint_dist( handle, savefilename="VJusage", tags_v = open("tags_trbv.txt", "rU"), tags_j = open("tags_trbj.txt", "rU")):

    ## PLOTS VJ JOINT GENE USAGE BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string
    import decimal as dec

    num_v = 0
    for line in tags_v:
        num_v += 1

    num_j = 0
    for line in tags_j:
        num_j += 1
        
    joint_distribution = np.zeros((num_v,num_j))
    for line in handle:
        elements = line.rstrip("\n")

        v = int(string.split(elements)[0])
        j = int(string.split(elements)[1])

        joint_distribution[v,j] += 1

    joint_distribution = joint_distribution / sum(sum(joint_distribution))
    gene_list_v = ('V10-1', 'V10-2', 'V10-3', 'V11-1', 'V11-2', 'V11-3', 'V12-3/V12-4', 'V12-5', 'V13', 'V14', 'V15', 'V16', 'V18', 'V19', 'V2', 'V20-1', 'V24-1', 'V25-1', 'V27-1', 'V28-1', 'V29-1', 'V3-1', 'V30-1', 'V4-1', 'V4-2', 'V4-3', 'V5-1', 'V5-4', 'V5-5', 'V5-6', 'V5-8', 'V6-1', 'V6-4', 'V6-5', 'V6-6', 'V6-8', 'V6-9', 'V7-2', 'V7-3', 'V7-4', 'V7-6', 'V7-7', 'V7-8', 'V7-9', 'V9')
    gene_list_j = ('J1-1', 'J1-2', 'J1-3', 'J1-4', 'J1-5', 'J1-6', 'J2-1', 'J2-2', 'J2-3', 'J2-4', 'J2-5', 'J2-6', 'J2-7')

    pos_v = np.arange(num_v)+ 1
    pos_j = np.arange(num_j)+ 1
    
    plt.figure()
    plt.pcolor(joint_distribution)
    pos_ticks_v = pos_v-0.5
    pos_ticks_j = pos_j-0.5
    plt.yticks( pos_ticks_v, gene_list_v)
    plt.xticks( pos_ticks_j, gene_list_j)
    plt.colorbar()
    plt.pcolor(joint_distribution)
    yticklabels = plt.getp(plt.gca(), 'yticklabels')
    plt.setp(yticklabels, fontsize='8')
    xticklabels = plt.getp(plt.gca(), 'xticklabels')
    plt.setp(xticklabels, fontsize='8')
    plt.savefig(str(savefilename)+'.png', dpi=300)

def plot_insert_lengths( handle, savefilename="InsertLengths" ):

    ## PLOTS DISTRIBUTION OF INSERT LENGTH BASED ON A FILE OF CLASSIFIERS
    
    import numpy as np
    import matplotlib.pyplot as plt
    import string

    maxi = 100
    insert_lengths = [0]*maxi

    handle = open("male1_total_unique.txt", "rU")
    for line in handle:
        elements = line.rstrip("\n")

        classifier = string.split(elements)
        if len(classifier) == 5:
            insert_lengths[len(classifier[4])] += 1
        else:
            insert_lengths[0] += 1

    total = sum(insert_lengths)
    for i in range(len(insert_lengths)):
        insert_lengths[i] = insert_lengths[i] / float(total)

    ind = np.arange(maxi)
    width = 0.5

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rects1 = ax.bar(ind, insert_lengths, width, color='blue')

    ax.set_ylabel('Frequency', fontsize = 16)
    ax.set_xlabel('Number of nucleotides', fontsize = 16)
    ax.grid(True)

    plt.setp(ax.get_xticklabels(),rotation='horizontal',fontsize = 16)
    plt.setp(ax.get_yticklabels(),rotation='horizontal',fontsize = 16)
    plt.xlim((0,50))

    plt.savefig(str(savefilename)+'.png', dpi=300)
