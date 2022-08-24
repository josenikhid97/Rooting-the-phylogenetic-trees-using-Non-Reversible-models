import subprocess as sub
import os
import sys
import time

# 1. Run the datasets with IQ-tree
# 2. Run the datasets with mst 
# 3. Run the datasets with mapple
# 4. Run the datasets with Rootdigger 


#################################
######## Bool Flags #############
#################################
convert_fasta_to_diff = False
Run_with_mstb = False
Run_with_iqt = False
Run_with_rtdg = True
Run_with_mpl = False

##################################
######### Input ##################
##################################

H3N2_fasta_input = "/project/exaptation/Jose_proj_bkp/data/H3N2/H3N2_156.fasta"
H3N2_newick_input = "/project/exaptation/Jose_proj_bkp/data/H3N2/H3N2_156.newick"
H3N2_diff_input = "/project/exaptation/Jose_proj_bkp/data/H3N2/H3N2_156.diff"
H3N2_ref_file = "/project/exaptation/Jose_proj_bkp/data/H3N2/H3N2_HA_ref_trimmed.fasta"

covid_ref_file = '/project/exaptation/covid_loop/EPI_ISL_402124_lowercase.fasta'
#no_dupli = "/project/exaptation/covid_loop/no_dupli_aligned_10covi_19-22_july.fasta"
covid_fasta_input = "/project/exaptation/covid_loop/aligned_10covi_19-22_july.fasta"
covid_newick_input = "/project/exaptation/covid_loop/aligned_10covi_19-22_july.newick"
covid_diff_input = "/project/exaptation/covid_loop/aligned_10covi_19-22_july.diff"

covid_fasta = "aligned_10covi_19-22_july.fasta"

covid_newick = "aligned_10covi_19-22_july.newick"

#################################
######## Tool Paths #############
#################################

mstb_tool = "/project/exaptation/Tools/mst-backbone --seq "
iqt_tool = "/project/exaptation/Tools/iqtree2 -s "
rtdg_tool = "/project/exaptation/Tools/root_digger/bin/rd --msa "
mpl_tool = "/project/exaptation/Tools/MAPLE/MAPLEv0.1.1.py --input "
pypy3_tool = "/project/exaptation/Tools/pypy3 "
diff_tool = "/project/exaptation/Tools/createDiffsFile.py --path /project/exaptation/covid_loop/ "

#################################
######## Result #################
#################################

rtdg_log = "/project/exaptation/Jose_proj_bkp/data/H3N2/H3N2_156_rtdg.log"
mpl_gtr_log = "/project/exaptation/Jose_proj_bkp/data/H3N2/H3N2_156_mpl_gtr.log"
mpl_jc_log = "/project/exaptation/Jose_proj_bkp/data/H3N2/H3N2_156_mpl_jc.log"
mpl_unrest_log = "/project/exaptation/Jose_proj_bkp/data/H3N2/H3N2_156_mpl_unrest.log"

cov_rtdg_log = "/project/exaptation/covid_loop/aligned_10covi_19-22_july_rtdg.log"
cov_mpl_gtr_log = "/project/exaptation/covid_loop/aligned_10covi_19-22_july_mpl_gtr.log"
cov_mpl_jc_log = "/project/exaptation/covid_loop/aligned_10covi_19-22_july_mpl_jc.log"
cov_mpl_unrest_log = "/project/exaptation/covid_loop/aligned_10covi_19-22_july_mpl_unrest.log"


# python3 createDiffsFile.py --path /pathToFolder/ --reference EPI_ISL_402124_lowercase.fasta --fasta 2021-03-31_unmasked.fa --output 2021-03-31_unmasked_differences.txt



if convert_fasta_to_diff:
    start = time.time()
    command_to_convert = sub.call("python3 " + diff_tool + "--reference EPI_ISL_402124_lowercase.fasta --fasta aligned_10covi_19-22_july.fasta --output aligned_10covi_19-22_july.diff", shell= True)
    print(command_to_convert, "The dataset ran on tool diff")
    time2 = time.time() - start
    print("Time to convert alignment file: "+str(time2))

start = time.time()

if Run_with_mstb:
    start = time.time()
    command_for_mstb = sub.call(mstb_tool + covid_fasta_input + " --out aligned_10covi_19-22_july", shell = True)
    print(command_for_mstb, "The dataset ran on tool mstb")
    time2 = time.time() - start
    print("Time to run mst-backbone: "+str(time2))

if Run_with_iqt:
    start = time.time()
    command_for_iqt = sub.call(iqt_tool + covid_fasta_input + " -m UNREST", shell = True)
    print(command_for_iqt, "The dataset ran on tool IQTREE")
    time2 = time.time() - start
    print("Time to run iqtree: "+str(time2))

if Run_with_rtdg:
    start = time.time()
    # command_for_rtdg = rtdg_tool + covid_fasta + " --tree " + covid_newick + " --prefix aligned_10covi_19-22_july" + " > " + cov_rtdg_log
    command_for_rtdg = sub.call(rtdg_tool + covid_fasta + " --tree " + covid_newick + " --prefix aligned_10covi_19-22_july" + " > " + cov_rtdg_log, shell = True)
    print(command_for_rtdg, "The dataset ran on tool RootDigger")
    time2 = time.time() - start
    print("Time to run rootdigger: "+str(time2))

if Run_with_mpl:
    start = time.time()
    # command_for_mpl_gtr1 = "python3 " + mpl_tool + covid_diff_input + " --reference " + covid_ref_file + " --output /project/exaptation/covid_loop/aligned_10covi_19-22_july_GTR --calculateLKfinalTree --model GTR" + " > " + mpl_gtr_log
    command_for_mpl_gtr = sub.call("python3 " + mpl_tool + covid_diff_input + " --reference " + covid_ref_file + " --output /project/exaptation/covid_loop/aligned_10covi_19-22_july_GTR --calculateLKfinalTree --model GTR" + " > " + cov_mpl_gtr_log, shell = True)
    print(command_for_mpl_gtr, "The dataset ran on tool MAPLE GTR")
    time2 = time.time() - start
    print("Time to run maple_gtr: "+str(time2))
    start1 = time.time()
    command_for_mpl_jc = sub.call("python3 " + mpl_tool + covid_diff_input + " --reference " + covid_ref_file + " --output /project/exaptation/covid_loop/aligned_10covi_19-22_july_JC --calculateLKfinalTree --model JC" + " > " + cov_mpl_jc_log, shell = True)
    print(command_for_mpl_jc, "The dataset ran on tool MAPLE JC")
    time3 = time.time() - start1
    print("Time to run maple_jc: "+str(time3))
    start2 = time.time()
    command_for_mpl_unrest = sub.call("python3 " + mpl_tool + covid_diff_input + " --reference " + covid_ref_file + " --output /project/exaptation/covid_loop/aligned_10covi_19-22_july_UNREST --calculateLKfinalTree --model UNREST" + " > " + cov_mpl_unrest_log, shell = True)
    print(command_for_mpl_unrest, "The dataset ran on tool MAPLE UNREST ")
    time4 = time.time() - start2
    print("Time to run maple_unrest: "+str(time4))


time2 = time.time() - start
print("Time to run all tools: "+str(time2))
