Decombinator: a tool for fast, efficient gene assignment in T cell receptor
sequences using a finite state machine

Niclas Thomas, James Heather, Wilfred Ndifon, John Shawe-Taylor, Benny Chain.

Introduction
############

Decombinator is a tool for the fast, efficient analysis of T cell receptor (TcR)
repertoire samples, designed to be accessible to those with no previous
programming experience. It is based on the Aho-Corasick algorithm which
uses a finite state automaton (FSA) to quickly assign a specific V and J
gene segment. From these assignments, it is then able to determine the number
of germline V and J deletions and the string of contiguous nucleotides which
lie between the 3' end of the V gene segment and the 5' end of the J gene
segment. These 5 variables form the identier which uniquely categorises each
distinct TcR sequence. For more details, please see (Thomas et al.).
Decombinator assumes no prior programming experience.

############
Decombinator has been tested on Python v2.6 and v2.7.
Decombinator requires the following Python modules:-

NumPy
Biopython
matplotlib
acora
Levenshtein

NumPy, BioPython and matplotlib install in a very straightforward manner - just follow the instructions,
accepting all defaults. To install acora, you will need GCC or something akin to it for your platform.
For Linux users, installation of GCC should be straightforward. For Windows users, something like MinGW
is needed - I recommend installing the GCC Win32 binaries from this website http://www.develer.com/oss/GccWinBinaries,
which is very straightforward to install - again, follow the instructions provided there. Installation on
a Mac should be possible with Xcode, though I haven't tested this yet. Installation of the Levenshtein
package will probably require Python setuptools to be downloaded and installed.

For anyone requiring further guidance for installation of any of the above packages, please see
www.ucl.ac.uk/innate2adaptive/software for a comprehensive set of instructions on installation and usage.

Once these modules have been successfully installed, Decombinator can be run
using Python from the command prompt (assuming Python has already been added
to your PATH, and that the downloaded folder Decombinator is located on your desktop) via, for example:-

cd C:\Users\...\Desktop\Decombinator
python
import Decombinator
sequencereads = ["mysequences.fastq"]
Decombinator.analysis( sequencereads )

#############
Further functionality is obtained by then using:-

Decombinator.plot_v_usage(handle=open("DecombinatorResults.txt","rU"),order="frequency")
Decombinator.plot_j_usage(handle=open("DecombinatorResults.txt","rU"),order="frequency")
Decombinator.plot_del_v(handle=open("DecombinatorResults.txt","rU"))
Decombinator.plot_del_j(handle=open("DecombinatorResults.txt","rU"))
Decombinator.plot_vj_joint_dist(handle=open("DecombinatorResults.txt","rU"))
Decombinator.plot_insert_lengths(handle=open("DecombinatorResults.txt","rU"))

*plot_v_usage and plot_j_usage have additional option:-
order="chromosome" # Plots V or J usage according to chromosome position

#############
Two further functions exist:-

Decombinator.get_distinct_clones(handle=open("DecombinatorResults.txt","rU"),with_count=False)
Decombinator.get_translated_sequences(handle=open("DecombinatorResults.txt","rU"),chain="beta",with_outframe=False,fullaaseq=False)

*option for get_distinct_clones() -
with_count=True # Includes the count of each distinct clone when writing classifiers to .txt file
*options for get_translated_seqeuences() -
chain="alpha" # when analysing alpha sequences
with_outframe=True # includes out of frame trascripts when writing to file
fullaaseq=True # includes the full translated sequence (CDR1, CDR2 and CDR3) when writing to file.