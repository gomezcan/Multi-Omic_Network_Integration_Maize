import matplotlib.cm as cm
import matplotlib.pyplot as plt
import igraph as ig
import pandas as pd
import numpy as np
import sys

## read network file
dfpdi = pd.read_csv(sys.argv[1], header = None, sep='\t') # network source target without headers 
dfpdi.columns = ["Source", "Target"]

# read name base
name = str(sys.argv[2])

'''
# Sample network to make printing test
Sample_dfpdi = dfpdi.sample(frac=0.01, replace=False, random_state=1)

# test with genes highly targeted
filtered = Sample_dfpdi.groupby('Source')['Target'].filter(lambda x: len(x) >= 200)

Sample_dfpdi = Sample_dfpdi[Sample_dfpdi['Target'].isin(filtered)]

print(" .. network and subnet size ..")
print(dfpdi.shape)
print(Sample_dfpdi.shape)

# make igraph from df
g = ig.Graph.DataFrame(Sample_dfpdi)
'''
# make ipgraph from df
g = ig.Graph.DataFrame(dfpdi)

### Identify communities ###

# Indirected version, requared to community_multilevel
g_un = g.as_undirected()

print(" .. Ready to calculate community ..")

# Method 1
community = g_un.community_multilevel()


print(" .. Community calculation done ..")

# get the membership vector
membersof = community.membership


# write df with vertex and membership index
nodes = [] # nodes

for i in range(len(g_un.vs)):
    nodes.append(str(g.vs[i]['name'])) 

out = pd.DataFrame(list(zip(nodes, membersof)), columns=['GeneID', "Community"]) 

out.to_csv("Community."+name+".txt",index=False, sep='\t')

'''
# plot net comunity
#layout = g.layout("lgl")
layout = g.layout("fr")

visual_style = {}
visual_style["edge_width"] = 0.05
visual_style["bbox"] = (300, 300)
visual_style["margin"] = 20
visual_style["layout"] = layout

fig, ax = plt.subplots(figsize=(20,14))
ax.axis("off")

ig.plot(community, target=ax, mark_groups = False, **visual_style, color='membership')



fig.savefig("Community."+name+".pdf", bbox_inches='tight')

#fig.savefig("PDI_comunities.pdf", bbox_inches='tight'i)
'''
