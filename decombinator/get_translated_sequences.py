def get_translated_sequences( handle, chain="beta", with_outframe=False, fullaaseq=False ):

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

    handle_vb=open("human_TRBV_region.fasta","rU")
    handle_jb=open("human_TRBJ_region.fasta","rU")
    handle_va=open("human_TRAV_region.fasta","rU")
    handle_ja=open("human_TRAJ_region.fasta","rU")
    
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
        for line in handle:
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
                    elif '*' not in aaseq:
                        print >> write_to, cdr3

    if chain == "alpha":
        for line in handle:
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
                    elif '*' not in aaseq:
                        print >> write_to, cdr3
            
    handle.close()
    write_to.close()