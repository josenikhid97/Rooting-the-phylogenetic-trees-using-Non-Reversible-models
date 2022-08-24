import os

from ete3 import Tree
if (os.path.exists("/project/exaptation/")):
    project_path = "/project/exaptation/"
elif (os.path.exists("/home/kalaghat/exaptation/")):
    project_path = "/home/kalaghat/exaptation/"

# T = {}
# T["unrest_mst_rd"] = Tree(project_path + "/Jose_proj_bkp/ECCB/data/10000_aligned_sars_newick.rooted.tree")
# # T["unrest_mst_rd"].ladderize()
# T["unrest_mst_rd"].show()

# import matplotlib
# import matplotlib.pyplot as plt
# from Bio import Phylo
# ms_tree = Phylo.read("/project/exaptation/Jose_proj_bkp/ECCB/data/10000_aligned_sars.fasta.newick_leafLabeledRooting", "newick")
from ete3 import Tree, NodeStyle, TreeStyle, CircleFace, TextFace

T = Tree("/project/exaptation/Jose_proj_bkp/ECCB/data/10000_aligned_sars.fasta.newick_leafLabeledRooting")

print ("Hello Reading the input ..")
# Basic tree style
Ts = TreeStyle()

Ts.branch_vertical_margin = 10000

Ts.show_leaf_name = False
# Ts.margin_top = 50
# Ts.margin_bottom = 0
# Ts.margin_left = 0
# Ts.margin_right = 50
# Ts.legend.add_face(CircleFace(1,"red", label = "2019"), column=1)
# Ts.legend.add_face(CircleFace(1, "green", label = "2020"), column=1)
# Ts.legend.add_face(CircleFace(1, "black", label = "2021"), column=1)
# Ts.legend.add_face(CircleFace(1, "purple", label = "2022"), column=1)
# legend_margin_line = 5
# while legend_margin_line:
#     Ts.legend.add_face(TextFace(" "), column=0)
#     Ts.legend.add_face(TextFace(" "), column=1)
#     Ts.legend.add_face(TextFace(" "), column=2)
#     Ts.legend.add_face(TextFace(" "), column=3)
#     legend_margin_line -= 1

Ts.legend.add_face(CircleFace(100, "red"), column=0)
Ts.legend.add_face(TextFace("  2019 ", fsize=100, ftype="Arial"), column=1)

Ts.legend.add_face(CircleFace(100, "green"), column=0)
Ts.legend.add_face(TextFace("  2020 ", fsize=100, ftype="Arial"), column=1)

Ts.legend.add_face(CircleFace(100, "black"), column=0)
Ts.legend.add_face(TextFace("  2021 ", fsize=100, ftype="Arial"), column=1)

Ts.legend.add_face(CircleFace(100, "purple"), column=0)
Ts.legend.add_face(TextFace("  2022 ", fsize=100, ftype="Arial"), column=1)

# Position legend in lower right corner
Ts.legend_position=4

# Add title
Ts.title.add_face(TextFace("Sars 10000 Variants", fsize=300, ftype="Arial", fstyle="italic"), column=0)



#Ts.legend_position = 2

for n in T.get_leaves():
    for spf_name in n.name.split():
        new_name = spf_name.split("/")[3][:4]
        replaced_year = spf_name.replace(spf_name, new_name)
#         print(replaced_year)
        if replaced_year == "2019":
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "red"
            nstyle["size"] = 2500
            n.set_style(nstyle)
        elif replaced_year == "2020":
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "green"
            nstyle["size"] =2500
            n.set_style(nstyle)
        elif replaced_year == "2021":
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "black"
            nstyle["size"] = 2500
            n.set_style(nstyle)
        elif replaced_year == "2022":
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "purple"
            nstyle["size"] = 2500
            n.set_style(nstyle)
# T.render("sars_variants_10000_phylogenie_with_colors.png", w=600, tree_style = Ts)

#Ts.scale =  1200
Ts.tree_width = 50000
Ts.force_topology = True
Ts.allow_face_overlap = True
T.ladderize()
T.render("mytree_sars_19.svg", tree_style=Ts, h =80000, w=100000,dpi = 12000, units="px")
# T.show(tree_style=Ts)
print ("picture saved")


# fig = plt.figure(figsize=(550, 450), dpi=100) # create figure & set the size 
# matplotlib.rc('font', size=22)              # fontsize of the leaf and node labels 
# matplotlib.rc('xtick', labelsize=20)       # fontsize of the tick labels
# matplotlib.rc('ytick', labelsize=20)       # fontsize of the tick labels
# #turtle_tree.ladderize()
# axes = fig.add_subplot(1, 1, 1)
# # Phylo.draw(ms_tree, axes=axes)
# plt.savefig("sars_10000_phylogenie_with_colors.pdf", format="pdf", bbox_inches="tight")
# print("Good Bye")
# # fig.savefig("sars_cladogram")
#h =60000, w=80000, units="px"