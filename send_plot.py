from ete3 import Tree, NodeStyle, TreeStyle, CircleFace, TextFace

T = Tree("10000_aligned_sars.fasta.newick_leafLabeledRooting")

print ("Hello Reading the input ..")
# Basic tree style
Ts = TreeStyle()

Ts.branch_vertical_margin = 10000

Ts.show_leaf_name = False


Ts.legend.add_face(CircleFace(100, "red"), column=0)
Ts.legend.add_face(TextFace("  2019 ", fsize=100, ftype="Arial"), column=1)

Ts.legend.add_face(CircleFace(100, "green"), column=0)
Ts.legend.add_face(TextFace("  2020 ", fsize=100, ftype="Arial"), column=1)

Ts.legend.add_face(CircleFace(100, "black"), column=0)
Ts.legend.add_face(TextFace("  2021 ", fsize=100, ftype="Arial"), column=1)

Ts.legend.add_face(CircleFace(100, "purple"), column=0)
Ts.legend.add_face(TextFace("  2022 ", fsize=100, ftype="Arial"), column=1)

Ts.legend_position=4


Ts.title.add_face(TextFace("Variants", fsize=300, ftype="Arial", fstyle="italic"), column=0)

#Ts.legend_position = 2
for n in T.get_leaves():
    for spf_name in n.name.split():
        new_name = spf_name.split("/")[3][:4]
        replaced_year = spf_name.replace(spf_name, new_name)
        if replaced_year == "2019":
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "red"
            nstyle["size"] = 2000
            n.set_style(nstyle)
        elif replaced_year == "2020":
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "green"
            nstyle["size"] = 2000
            n.set_style(nstyle)
        elif replaced_year == "2021":
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "black"
            nstyle["size"] = 2000
            n.set_style(nstyle)
        elif replaced_year == "2022":
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere"
            nstyle["fgcolor"] = "purple"
            nstyle["size"] = 2000
            n.set_style(nstyle)

#Ts.scale =  1200
Ts.tree_width = 50000
Ts.force_topology = True
T.ladderize()
T.render("mytree_15.svg", tree_style=Ts, h =80000, w=100000, units="px")
# T.show(tree_style=Ts)