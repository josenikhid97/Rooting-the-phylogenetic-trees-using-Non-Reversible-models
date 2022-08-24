import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import csv
import networkx as nx

sars = '/project/exaptation/Jose_proj_bkp/ECCB/data/10000_aligned_sars_for_rd.csv'
H3n2 = '/project/exaptation/Jose_proj_bkp/ECCB/data/H3N2_HA_full.csv'

G = nx.Graph()
with open('/project/exaptation/Jose_proj_bkp/ECCB/data/H3N2_HA_full.csv', newline='') as csvfile:
    datareader = csv.reader(csvfile, delimiter='\t')
    for row in datareader:
        G.add_edge(row[0],row[1],d=float(row[2]))

fig = plt.figure(figsize=(125,123))
nx.draw_networkx(G,with_labels=True, 
                 font_size=30, 
                 cmap=plt.cm.coolwarm,
                 pos=nx.kamada_kawai_layout(G),vmin=0, vmax=1)
plt.savefig("H3N2_9504.pdf", format="pdf", bbox_inches="tight")
# plt.show()
print("H3N2 plot is done")

G1 = nx.Graph()
with open('/project/exaptation/Jose_proj_bkp/ECCB/data/10000_aligned_sars_for_rd.csv', newline='') as csvfile:
    datareader = csv.reader(csvfile, delimiter='\t')
    for row in datareader:
        G1.add_edge(row[0],row[1],d=float(row[2]))

fig = plt.figure(figsize=(225,223))
nx.draw_networkx(G1,with_labels=True, 
                 font_size=60, 
                 cmap=plt.cm.coolwarm,
                 pos=nx.kamada_kawai_layout(G1),vmin=0, vmax=1)
plt.savefig("sars_10000.pdf", format="pdf", bbox_inches="tight")
print("SARS plot is done")
# plt.show()
