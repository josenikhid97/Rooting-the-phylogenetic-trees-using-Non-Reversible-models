import subprocess as sub
import os
import time
import sys
from Bio import SeqIO
import edlib
from edlib import align, getNiceAlignment
from multiprocessing.dummy import current_process



# 1. Align the sequences with edlib
# 3. Run the datasets with mst 
# 4. Run the datasets with mapple
# 5. convert datasets form fasta to mutdiff format
# 6. combine all the single files into a big file and visualise it 


#####################################
####### Bool Flags ##################
#####################################
remove_shortseq_with_80percent_gaps = False
align_sequences = True
convert_multiline_seq_to_singleline = False
replace_short_seq_with_gaps = False
trim_long_seq_to_ref_length = False
align_seq_with_mafft = False
run_with_iqtree = False




#####################################
###### Paths ########################
#####################################

####### input ##############

# reference seq fasta file
reference_file_path = "/project/exaptation/Jose_proj/data/covid_ref_seq/EPI_ISL_402124_uppercase.fasta"
# unaligned seqs fasta file (size n)
unaligned_seqs_file_path = "/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/unaligned_seq/full_variants_covid.fasta"
# aligned seq fasta file 
aligned_seq_file_path = "/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/Aligned_seq/fasta_format/aligned_100.fasta"
convert_mls_to_sls = "/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/Aligned_seq/mls_to_sls/all_variants_covid_sls.fasta"
replace_short_gaps = "/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/Aligned_seq/replaced_with_gaps/aligned_sars_all_10000.fasta"
removed_short_seq_path = "/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/Aligned_seq/removed_short_seq/removed_short__sars_aligned_10000.fasta"
length_seq = """ '/^>/{ seqlen=0; print; next; } seqlen < 29891 { if (seqlen + length($0) > 29891) $0 = substr($0, 1, 29891-seqlen); seqlen += length($0); print }' """
length_corrected_path = "/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/Aligned_seq/length_corrected/length_corrected_10000.fasta"
covid_proj_path = "/project/exaptation/Jose_proj/covid_results/"
parts_path = "/project/exaptation/Jose_proj/covid_results/parts/"
scripts_path = covid_proj_path + "scripts/"
mafft_result_path = "/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/Aligned_seq/mafft_result/aligned_covid_variants_full.fasta"

####### Tool paths #########
mafft_tool_path = "/project/exaptation/Tools/mafft --thread 5"
raxml_path = "/project/exaptation/Tools/raxml"
# iqtree_toolpath = 

if align_sequences:
    start = time.time()
    # unaligned_sequences_file_name = sys.argv[0]
    reference_sequence_file_name = '/project/exaptation/Jose_proj/data/covid_ref_seq/EPI_ISL_402124_lowercase.fasta'
    reference_sequence_file = open(reference_sequence_file_name,'r')
    ref_seq = ''
    for lines in reference_sequence_file:
        if lines.startswith('>'):
            ref_seq_name = lines.strip()[1:]
            # print(ref_seq_name)        
        else:
            ref_seq += lines.strip()
    reference_sequence_file.close()

    # print(ref_seq_name)
    # print(ref_seq)

    unaligned_sequences_file_name = '/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/unaligned_seq/unaligned_100.fasta'
    aligned_result_sequences_file_name = '/project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/Aligned_seq/edlib/jaligned_100.fasta'
    # sys.argv[1]
    unaligned_sequences_file = open(unaligned_sequences_file_name,'r')
    align_sequences_file = open(aligned_result_sequences_file_name,'w')
    ref_seq_length = len(ref_seq)
    unaligned_seq = ''

    for line in unaligned_sequences_file:
        if line.startswith('>'):
            if unaligned_seq != '':
                alignment_result = align(query=unaligned_seq, target=ref_seq, mode='NW', task='path')
                nice_alignment_result = getNiceAlignment(alignment_result, query=unaligned_seq, target=ref_seq)
                matched_aligned = nice_alignment_result['matched_aligned']
                query_aligned = nice_alignment_result['query_aligned']
                target_aligned = nice_alignment_result['target_aligned']
                cigar = alignment_result['cigar']
                edit_string = cigar.split('=')
                align_sequences_file.write('>' + seq_name +'\n'+ query_aligned +'\n')
            seq_name = line.strip()[1:]
            unaligned_seq = ''
        else:
            unaligned_seq += line.strip().lower()

    align_sequences_file.close()
    unaligned_sequences_file.close()
    time2 = time.time() - start
    print("Time to convert unaligned to aligned file: "+str(time2)) 

if convert_multiline_seq_to_singleline:
    # seqtk seq multi-line.fasta > single-line.fasta
    command_for_convert = sub.call("seqtk seq " + unaligned_seqs_file_path + " > " + convert_mls_to_sls , shell = True)
    print(command_for_convert, "The dataset has been converted form multiline sequences to singleline sequences")

if replace_short_seq_with_gaps:
    sequences = [s for s in SeqIO.parse(convert_mls_to_sls, 'fasta')]
    max_len = max([len(s.seq) for s in sequences])
    GAPs = "-"
    for seq in sequences:
        padding = GAPs*(max_len - len(seq.seq)) # creating the padding string
        seq.seq += padding

    SeqIO.write(sequences, replace_short_gaps, 'fasta')
    print("The short sequences have been replaced with the gaps to the max lenth of the sequence")

if remove_shortseq_with_80percent_gaps:
    # This will drop all the sequences that have gaps in >=80% of alignment positions.
    # (python fasta_drop.py file.fas trimmed.fas 0.8)
    FastaFile = open(replace_short_gaps, 'r')
    FastaDroppedFile = open(removed_short_seq_path, 'w')
    drop_cutoff = float(0.8)

    if (drop_cutoff > 1) or (drop_cutoff < 0):
        print('\n Sequence drop cutoff must be in 0-1 range !\n')
        sys.exit(1)

    for seqs in SeqIO.parse(FastaFile, 'fasta'):
        name = seqs.id
        seq = seqs.seq
        seqLen = len(seqs)
        gap_count = 0
        for z in range(seqLen):
            if seq[z]=='-' or seq[z]=='N':
                gap_count += 1
        if (gap_count/float(seqLen)) >= drop_cutoff:
            print(' %s was removed.' % name)
        else:
            SeqIO.write(seqs, FastaDroppedFile, 'fasta')

    FastaFile.close()
    FastaDroppedFile.close()
    print("Dropped all the sequences that have gaps in >=80percentage of alignment positions")

if trim_long_seq_to_ref_length:
    # awk '/^>/{ seqlen=0; print; next; } seqlen < 29891 { if (seqlen + length($0) > 29891) $0 = substr($0, 1, 29891-seqlen); seqlen += length($0); print }' aligned_1000.fasta > length_corrected_1000.fasta
    # sub.call("cat " + reference_file_path + " > " + raw_part_file_name, shell=True)
    # command_for_seq = sub.call("awk '/^>/{ seqlen=0; print; next; } seqlen " + " > " + " 29891 { if (seqlen + length($0) " + " > " + " 29891) $0 = substr($0, 1, 29891-seqlen); seqlen += length($0); print }' " + aligned_seq_file_path + " > " + length_corrected_path , shell=True)
    command_for_seq = sub.call("awk " + length_seq + removed_short_seq_path + " > " + length_corrected_path , shell=True)
    print(command_for_seq, "length of the sequences is corrected to the reference sequence length")


if align_seq_with_mafft:
    # ./mafft --thread 5 /project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/unaligned_seq/unaligned_5000.fasta > /project/exaptation/Jose_proj/data/gisaid/sequences_for_testing/unaligned_seq/aligned_5000.fasta
    # mafft_script_file = open(scripts_path + "mafft_batch.sh","w")
    command_for_mafft = sub.call(mafft_tool_path + " " + unaligned_seqs_file_path + " > " + mafft_result_path , shell = True)
    # mafft_script_file.write(command_for_mafft)
    # mafft_script_file.close()
    print(command_for_mafft, "The raw sequences has been aligned sucessfullly")

# if run_with_iqtree:
#     iqtree_command = tool_path + "iqtree2 -s " + test_data_path + test_prefix + ".fa -st DNA -te " + test_data_path + test_prefix + "_bi.fastree" + " -pre " + test_prefix + " -blfix -keep-ident -m GTR\{0.04,0.3,0.1,0.02,1.0,1.0\}+FQ -redo -nt 1  > " + result_path + test_prefix + "_iqtree.log"
#     # print(iqtree_command)
#     sub.call(iqtree_command, shell=True)       

## >hCoV-19/Wuhan/WIV04/2019