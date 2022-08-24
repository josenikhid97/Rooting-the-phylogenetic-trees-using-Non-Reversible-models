import subprocess as sub
import os
import sys
import time
import csv
from csv import reader
import pandas as pd
from ete3 import Tree
from scipy.stats import pearsonr
from sklearn import model_selection
from math import log



loglikelihood = True
pearson_correlation = False
Model_selection = False

if (os.path.exists("/project/exaptation/")):
    project_path = "/project/exaptation/"
elif (os.path.exists("/home/kalaghat/exaptation/")):
    project_path = "/home/kalaghat/exaptation/"
sys.path.insert(0,project_path+'Projects/MSTBasedForests/scripts/')
toolPath = project_path+"Tools/"

tool_list = ['IQ-tree_UNREST', 'GMM_MST-B', 'Maple_UNREST','Maple_GTR','Maple_JC'] 
prefix_list = ['Pearsons correlation', 'Maximum - likelihood', "AIC values ", 'AICc values ','BIC values']
test_prefix = 'aligned_10covi_19-22_july'
result_path = "/project/exaptation/covid_loop/"  

if loglikelihood:
    # def get_loglikelihood_for_tool_and_prefix(tool_name,test_prefix):
        T = {}
        T["Maple_GTR"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july_GTR_tree.tree")
        T["Maple_UNREST "] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july_UNREST_tree.tree")
        T["Maple_JC"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july_JC_tree.tree")
        T["GMM_MST-B"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july.fasta.newick_leafLabeledRooting")
        T["IQ-tree_UNREST"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july.fasta.bionj")
        # T["unrest_mst_rd"] = Tree(project_path + "/Jose_proj_bkp/ECCB/data/10000_aligned_sars_newick.rooted.tree")
        # T["unrest_mst_rd_fully"] = Tree(project_path + "/Jose_proj/ECCB/data/10000_aligned_sars_fully.rooted.tree")

        list = ['IQ-tree_UNREST','GMM_MST-B','Maple_UNREST ','Maple_GTR','Maple_JC']#,'unrest_mst_rd']#,'unrest_mst_rd_fully']#, 'unrest_iqtree','unrest_mst_rd']'unrest_maple','unrest_maple'GTR_maple',
        print("########### Pearson Correlation ##############")
        for tool_name in list:
            # print(file)
            root = T[tool_name].get_tree_root()
            # print(root)
            leaves = T[tool_name].get_leaves()
            dist_from_root = []
            sampling_time = []
            for l in leaves:
                dist_from_root.append(root.get_distance(l))
                sampling_time.append(float(l.name.split("/")[-1][0:4]))
            corr = pearsonr(dist_from_root,sampling_time)[0]
            
            print("Pearson's rho for SARS is " + tool_name,corr)
            
        for tool_name in tool_list:
            test_prefix = 'aligned_10covi_19-22_july'
        # def get_loglikelihood_for_tool_and_prefix(tool_name,test_prefix):
            result_path = "/project/exaptation/covid_loop/"  
            n = 10885
            k = 29891  
            if (tool_name == "IQ-tree_UNREST"):
                log_file_name =  result_path + test_prefix + ".fasta.log"
                # print(log_file_name)
                # iq_log_lik = ()
                result = sub.check_output('grep \'Optimal log-likelihood\' ' + log_file_name, shell=True).decode('utf8').strip()
                iq_log_lik = result.split(':')[1].strip()
                print('###############################################')
                print("IQ-Tree likelihood; ", iq_log_lik)

                n = 10885
                k = 29891
                def AIC(loglik, k):
                    return -2*loglik + 2*k
                # value = AIC(-138506.436, 1042)
                # value_1 = AIC(-112380, 1042)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                def BIC(loglik, n, k):
                    return -2*loglik + (4*(n-1))*log(k)
                
                AIC_LIK = AIC(-335454, k)
                BIC_LIK = BIC(-335454, n, k)
                AICc_LIK = AICc(AIC_LIK, n, k)

                print("Model:"+ tool_name +" for sars, AIC =", AIC_LIK)
                print("Model:"+ tool_name +" for sars, BIC  =", BIC_LIK )
                print(f"Model:"+ tool_name +" for sars, AICc  =", AICc_LIK)
                # corr = pearsonr(dist_from_root,sampling_time)[0]
                # print("Pearson's rho for SARS is " + tool_name,corr)
                print('###############################################')
            if (tool_name == "GMM_MST-B"):
                log_file_name =  result_path + test_prefix + ".mstbackbone_log"
                result = sub.check_output('grep \'Log likelihood of fitting a GM model to leaf-labeled T via restricted SEM is\' ' + log_file_name, shell=True).decode('utf8').strip()
                mst_log_lik = result.split('is')[1].strip()
                print("MST-Backbone likelihood: ", mst_log_lik)
                def AIC(loglik, k):
                    return -2*loglik + 2*k
                # value = AIC(-138506.436, 1042)
                # value_1 = AIC(-112380, 1042)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                def BIC(loglik, n, k):
                    return -2*loglik + (3+12*(n-1))*log(k)
                
                AIC_LIK = AIC(-1512300, k)
                BIC_LIK = BIC(-1512300, n, k)
                AICc_LIK = AICc(AIC_LIK, n, k)

                print("Model:"+ tool_name +" for sars, AIC =", AIC_LIK)
                print("Model:"+ tool_name +" for sars, BIC  =", BIC_LIK )
                print(f"Model:"+ tool_name +" for sars, AICc  =", AICc_LIK)
                # corr = pearsonr(dist_from_root,sampling_time)[0]
                # print("Pearson's rho for SARS is " + tool_name,corr)
                print('###############################################')

            # if (tool_name == "MST-B_RootDigger"):
            #     log_file_name =  result_path + test_prefix + "_rtdg.log"
            #     result = sub.check_output('grep \'Final LogLH\' ' + log_file_name, shell=True).decode('utf8').strip()
            #     log_lik = result.split(':')[1].strip()
            #     print(log_lik)
            #     def AIC(loglik, k):
            #         return -2*loglik + 2*k
            #     # value = AIC(-138506.436, 1042)
            #     # value_1 = AIC(-112380, 1042)
            #     def AICc(AIC, n, k):
            #         return AIC + 2*k*(k + 1)/(n-k-1)

            #     def BIC(loglik, n, k):
            #         return -2*loglik + (4*(n-1))*log(k)

            if (tool_name == "Maple_UNREST"):
                log_file_name =  result_path + test_prefix + "_mpl_unrest.log"
                result = sub.check_output('grep \'totalLK\' ' + log_file_name, shell=True).decode('utf8').strip()
                mplu_log_lik = result.split(':')[1].strip()
                print( "Maple UNREST likelihood:",mplu_log_lik)
                def AIC(loglik, k):
                    return -2*loglik + 2*k
                # value = AIC(-138506.436, 1042)
                # value_1 = AIC(-112380, 1042)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                def BIC(loglik, n, k):
                    return -2*loglik + (4*(n-1))*log(k)

                AIC_LIK = AIC(-332293.6110839789, k)
                BIC_LIK = BIC(-332293.6110839789, n, k)
                AICc_LIK = AICc(AIC_LIK, n, k)

                print("Model:"+ tool_name +" for sars, AIC =", AIC_LIK)
                print("Model:"+ tool_name +" for sars, BIC  =", BIC_LIK )
                print(f"Model:"+ tool_name +" for sars, AICc  =", AICc_LIK)
                # corr = pearsonr(dist_from_root,sampling_time)[0]
                # print("Pearson's rho for SARS is " + tool_name,corr)
                print('###############################################')

            if (tool_name == "Maple_GTR"):
                log_file_name =  result_path + test_prefix + "_mpl_gtr.log"
                result = sub.check_output('grep \'totalLK\' ' + log_file_name, shell=True).decode('utf8').strip()
                mplg_log_lik = result.split(':')[1].strip()
                print( "Maple GTR likelihood:",mplg_log_lik)
                def AIC(loglik, k):
                    return -2*loglik + 2*k
                # value = AIC(-138506.436, 1042)
                # value_1 = AIC(-112380, 1042)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                def BIC(loglik, n, k):
                    return -2*loglik + (3*(n-1))*log(k)

                AIC_LIK = AIC(-337818.80224763945, k)
                BIC_LIK = BIC(-337818.80224763945, n, k)
                AICc_LIK = AICc(AIC_LIK, n, k)

                print("Model:"+ tool_name +" for sars, AIC =", AIC_LIK)
                print("Model:"+ tool_name +" for sars, BIC  =", BIC_LIK )
                print(f"Model:"+ tool_name +" for sars, AICc  =", AICc_LIK)
                # corr = pearsonr(dist_from_root,sampling_time)[0]
                # print("Pearson's rho for SARS is " + tool_name,corr)
                print('###############################################')

            if (tool_name == "Maple_JC"):
                log_file_name =  result_path + test_prefix + "_mpl_jc.log"
                result = sub.check_output('grep \'totalLK\' ' + log_file_name, shell=True).decode('utf8').strip()
                mplj_log_lik = result.split(':')[1].strip()
                print( "Maple JC likelihood:",mplj_log_lik)
                def AIC(loglik, k):
                    return -2*loglik + 2*k
                # value = AIC(-138506.436, 1042)
                # value_1 = AIC(-112380, 1042)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                def BIC(loglik, n, k):
                    return -2*loglik + (n-1)*log(k)

                AIC_LIK = AIC(-357416.0980594669, k)
                BIC_LIK = BIC(-357416.0980594669, n, k)
                AICc_LIK = AICc(AIC_LIK, n, k)

                print("Model:"+ tool_name +" for sars, AIC =", AIC_LIK)
                print("Model:"+ tool_name +" for sars, BIC  =", BIC_LIK )
                print(f"Model:"+ tool_name +" for sars, AICc  =", AICc_LIK)
                # corr = pearsonr(dist_from_root,sampling_time)[0]
                # print("Pearson's rho for SARS is " + tool_name,corr)
                print('###############################################')        
            
            # return iq_log_lik,mst_log_lik,mplu_log_lik,mplg_log_lik,mplj_log_lik
    


if pearson_correlation:
    T = {}
    T["GTR_maple"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july_GTR_tree.tree")
    T["unrest_maple"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july_UNREST_tree.tree")
    T["JC_maple"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july_JC_tree.tree")
    T["gmm_mstb"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july.fasta.newick_leafLabeledRooting")
    T["unrest_iqtree"] = Tree(project_path + "covid_loop/aligned_10covi_19-22_july.fasta.bionj")
    # T["unrest_mst_rd"] = Tree(project_path + "/Jose_proj_bkp/ECCB/data/10000_aligned_sars_newick.rooted.tree")
    # T["unrest_mst_rd_fully"] = Tree(project_path + "/Jose_proj/ECCB/data/10000_aligned_sars_fully.rooted.tree")

    list = ['GTR_maple','unrest_maple','JC_maple','gmm_mstb','unrest_iqtree']#,'unrest_mst_rd']#,'unrest_mst_rd_fully']#, 'unrest_iqtree','unrest_mst_rd']'unrest_maple','unrest_maple'GTR_maple',

    for file in list:
        # print(file)
        root = T[file].get_tree_root()
        # print(root)
        leaves = T[file].get_leaves()
        dist_from_root = []
        sampling_time = []
        for l in leaves:
            dist_from_root.append(root.get_distance(l))
            sampling_time.append(float(l.name.split("/")[-1][0:4]))

        corr = pearsonr(dist_from_root,sampling_time)[0]

        print("########### Pearson Correlation##############")
        print("Pearson's rho for SARS is " + file,corr)
        # print(sampling_time)

# if Model_selection:
#     def AIC(loglik, k):
#         return -2*loglik + 2*k
#     # value = AIC(-138506.436, 1042)
#     # value_1 = AIC(-112380, 1042)
#     def AICc(AIC, n, k):
#         return AIC + 2*k*(k + 1)/(n-k-1)

#     def BIC(loglik, n, k):
#         return -2*loglik + (n-1)*log(k)

#     sars_Tool_names = ['Maple_UNREST','Maple_GTR','Maple_JC','MST-back_GMM','IQ-Tree_UNREST']
#     sars_tool_loglik = ['-332293.6110839789', '-337818.8022','-357416.0981','-1512300','-335453.75']

#     # unrest = 4
#     # gmm = 3 + 12
#     # jc = 0
#     # GTR = 3

#     #10_covi_sequences = 10855

#     # return -2*loglik + (3*(n-1))*log(k)
#     # for tool in Tool_names:
#     #     print(tool)
#     # for loglike in tool_loglik:
#     #     print(loglike)

#     # print(f"Model:"+ tool +" for sars, AIC =" AIC(loglike, 9499, 29891))
#     # print(f"Model:"+ tool +"  for sars, BIC  =" {BIC(loglike, 9499, 29891)})


#     # print(f"Model: MST-Backbone for SARS, BIC = {BIC(-121289.0, 9499, 29891)}")
#     # print(f"Model: IQTree for sars, BIC = {BIC(-417268.9775, 9499, 29891)}")
#     # print(f"Model: MST-Backbone for SARS, BIC = {BIC(-121289.0, 9499, 29891)}")
    # print(f"Model: MST-back_GMM for sars, AIC = {AIC(-357416.0981, 29891)}")
    # print(f"Model: MST-back_GMM for sars, BIC = {BIC(-357416.0981, 10855, 29891)}")
    # print(f"Model: MST-back_GMM for sars, AICc = {AICc(774614.1962, 10855, 29891)}")




if Model_selection:
    # def get_loglikelihood_for_tool_and_prefix(tool_name,test_prefix):
        for tool_name in tool_list:
            # test_prefix = 'aligned_10covi_19-22_july'
        # def get_loglikelihood_for_tool_and_prefix(tool_name,test_prefix):
            # result_path = "/project/exaptation/covid_loop/" 
            n = 10855
            k = 29891   
            if (tool_name == "IQ-tree_UNREST"):
                n = 10855
                k = 29891  
                # log_file_name =  result_path + test_prefix + ".fasta.log"
                # # print(log_file_name)
                # result = sub.check_output('grep \'Optimal log-likelihood\' ' + log_file_name, shell=True).decode('utf8').strip()
                # iq_log_lik = result.split(':')[1].strip()
                # print("IQ-Tree likelihood; ", iq_log_lik)
                def AIC(log_lik, k):
                    return -2*log_lik + 2*k
                def BIC(log_lik, n, k):
                    return -2*log_lik + (4*(n-1))*log(k)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                AIC_LIK = AIC(iq_log_lik, k)

                
                print(f"Model: for sars, AIC = ", AIC)
                print(f"Model: for sars, BIC  = " , AIC)
                # print(f"Model:"+ tool_name +" for sars, AICc  = " , AICc(AIC, n, k))
                
            

            if (tool_name == "GMM_MST-B"):
                # log_file_name =  result_path + test_prefix + ".mstbackbone_log"
                # result = sub.check_output('grep \'Log likelihood of fitting a GM model to leaf-labeled T via restricted SEM is\' ' + log_file_name, shell=True).decode('utf8').strip()
                # mst_log_lik = result.split('is')[1].strip()
                # print("MST-Backbone likelihood: ", mst_log_lik)
                def AIC(mst_log_lik, k):
                    return -2*mst_log_lik + 2*k
                def BIC(mst_log_lik, n, k):
                    return -2*mst_log_lik + (3+12*(n-1))*log(k)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                print(f"Model:"+ tool_name +" for sars, AIC =",AIC)
                print(f"Model:"+ tool_name +" for sars, BIC  =" ,BIC)
                print(f"Model:"+ tool_name +" for sars, AICc  =" ,AICc)

            # if (tool_name == "MST-B_RootDigger"):
            #     log_file_name =  result_path + test_prefix + "_rtdg.log"
            #     result = sub.check_output('grep \'Final LogLH\' ' + log_file_name, shell=True).decode('utf8').strip()
            #     log_lik = result.split(':')[1].strip()
            #     print(log_lik)
            #     def AIC(loglik, k):
            #         return -2*loglik + 2*k
            #     # value = AIC(-138506.436, 1042)
            #     # value_1 = AIC(-112380, 1042)
            #     def AICc(AIC, n, k):
            #         return AIC + 2*k*(k + 1)/(n-k-1)

            #     def BIC(loglik, n, k):
            #         return -2*loglik + (4*(n-1))*log(k)

            if (tool_name == "Maple_UNREST"):
                log_file_name =  result_path + test_prefix + "_mpl_unrest.log"
                result = sub.check_output('grep \'totalLK\' ' + log_file_name, shell=True).decode('utf8').strip()
                mplu_log_lik = result.split(':')[1].strip()
                print( "Maple UNREST likelihood:",mplu_log_lik)
                def AIC(mplu_log_lik, k):
                    return -2*mplu_log_lik + 2*k
                def BIC(mplu_log_lik, n, k):
                    return -2*mplu_log_lik + (4*(n-1))*log(k)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)
                
                print(f"Model:"+ tool_name +" for sars, AIC =" ,AIC)
                print(f"Model:"+ tool_name +" for sars, BIC  =", BIC)
                print(f"Model:"+ tool_name +" for sars, AICc  =" ,AICc)



            if (tool_name == "Maple_GTR"):
                # log_file_name =  result_path + test_prefix + "_mpl_gtr.log"
                # result = sub.check_output('grep \'totalLK\' ' + log_file_name, shell=True).decode('utf8').strip()
                # mplg_log_lik = result.split(':')[1].strip()
                # print( "Maple GTR likelihood:",mplg_log_lik)
                def AIC(mplg_log_lik, k):
                    return -2*mplg_log_lik + 2*k
                def BIC(mplg_log_lik, n, k):
                    return -2*mplg_log_lik + (3*(n-1))*log(k)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                print(f"Model:"+ tool_name +" for sars, AIC =" ,AIC)
                print(f"Model:"+ tool_name +" for sars, BIC  =", BIC)
                print(f"Model:"+ tool_name +" for sars, AICc  =" ,AICc)

            if (tool_name == "Maple_JC"):
                # log_file_name =  result_path + test_prefix + "_mpl_jc.log"
                # result = sub.check_output('grep \'totalLK\' ' + log_file_name, shell=True).decode('utf8').strip()
                # mplj_log_lik = result.split(':')[1].strip()
                # print( "Maple JC likelihood:",mplj_log_lik)
                def AIC(mplj_log_lik, k):
                    return -2*mplj_log_lik + 2*k
                def BIC(mplj_log_lik, n, k):
                    return -2*mplj_log_lik + (n-1)*log(k)
                def AICc(AIC, n, k):
                    return AIC + 2*k*(k + 1)/(n-k-1)

                print(f"Model:"+ tool_name +" for sars, AIC =" ,AIC)
                print(f"Model:"+ tool_name +" for sars, BIC  =" ,BIC)
                print(f"Model:"+ tool_name +" for sars, AICc  =" ,AICc)




# log_lik_dic = {}
# for tool in tool_list:
#     log_lik_dic[tool] = {}
#     for prefix in prefix_list:
#         log_lik_dic[tool][prefix] = 0
# for tool in tool_list:
#     for prefix in prefix_list:
#         print(tool,prefix)
#         log_lik = get_loglikelihood_for_tool_and_prefix(tool,prefix)
#         # log_lik = 42
#         log_lik_dic[tool][prefix] = log_lik
#         # dfObj[tool,prefix] = log_lik



# # writing csv
# csv_file = open(result_path + "results.csv","w")
# csv_file.write("Tool")
# for prefix in prefix_list:
#     csv_file.write( "," + prefix)
# csv_file.write("\n")
# for tool in tool_list:
#     csv_file.write(tool)
#     for prefix in prefix_list:
#         csv_file.write("," + str(log_lik_dic[tool][prefix]))
#     csv_file.write("\n")
# csv_file.close()