#Decombinator

##Notes
This version is significantly different from the original while maintaining the
primary functionality.

####Original authors
Niclas Thomas, James Heather, Wilfred Ndifon, John Shawe-Taylor, Benny Chain.

####Their repo
https://github.com/uclinfectionimmunity/Decombinator

##Description
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

##Requires
NumPy
Biopython
matplotlib
acora
Levenshtein

##Usage
decombinator-analyze [OPTIONS] fastq
decombinator-distinct-clones
decombinator-translated-sequences

This stuff needs to be changed...
    Decombinator.plot_v_usage(handle=open("DecombinatorResults.txt","rU"),order="frequency")
    Decombinator.plot_j_usage(handle=open("DecombinatorResults.txt","rU"),order="frequency")
    Decombinator.plot_del_v(handle=open("DecombinatorResults.txt","rU"))
    Decombinator.plot_del_j(handle=open("DecombinatorResults.txt","rU"))
    Decombinator.plot_vj_joint_dist(handle=open("DecombinatorResults.txt","rU"))
    Decombinator.plot_insert_lengths(handle=open("DecombinatorResults.txt","rU"))