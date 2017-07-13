#!/usr/bin/env python
#
# 	This scripts plots the distribution of PFAM sizes wrt PDB structure annotation
#
###############################################################################

#===============================================================================
#== libraries
#===============================================================================

import pandas as pd
import numpy as np
import urllib2
from xml.dom import minidom


import plotly.graph_objs as go
from plotly.offline import plot as plotly_plot
from statsmodels.nonparametric.kde import KDEUnivariate
from sklearn.neighbors import KernelDensity
#===============================================================================
#== mapping from PFAM to PDB
#===============================================================================

#access to latest PFAM release
url = 'http://pfam.xfam.org/families?output=xml' #define XML location
dom = minidom.parse(urllib2.urlopen(url)) # parse the data

#prepare data frame for plotting
pfam_df = pd.DataFrame(columns=['accession', 'nr_sequences', 'nr_structures'])

#retrieve info from PFAM
entry_list = dom.getElementsByTagName("entry")
for id, entry in enumerate(entry_list):

    print("{0}/{1}".format(id, len(entry_list)))

    accession = str(entry.getAttribute('accession'))

    pfam_dom = minidom.parse(urllib2.urlopen("http://pfam.xfam.org/family/"+accession+"?output=xml"))
    nr_seq = int(pfam_dom.getElementsByTagName("full")[0]._get_firstChild().nodeValue)
    nr_struct = int(pfam_dom.getElementsByTagName("num_structures")[0]._get_firstChild().nodeValue)

    pfam_df.loc[id]= [accession, nr_seq, nr_struct]




#===============================================================================
#== plotting
#===============================================================================


plot_out = "/Users/Susann.Vorberg/pfam_pdb.html"


#remove erronous (?) entries with  family size = 0
pfam_df = pfam_df.drop(pfam_df.query('nr_sequences == 0').index)

#define counts for PFAmilies with and without annotated PDB structures
struct = np.log(pfam_df.query('nr_structures > 0')['nr_sequences'].values)
no_struct = np.log(pfam_df.query('nr_structures == 0')['nr_sequences'].values)




#prepare plot
fig={}
fig['data'] = []
fig['layout'] = {}
fig['layout']['xaxis'] = {}
fig['layout']['yaxis'] = {}

#define grid for kernel density estimation
x_grid = np.linspace(np.min(struct.tolist() +  no_struct.tolist()), np.max(struct.tolist() +  no_struct.tolist()), 500)
bandwidth=0.3
colors=['rgb(22, 96, 167)', 'rgb(205, 12, 24)']





#kernel density estimate for Pfamilies with annotated structure
kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(struct.reshape(-1,1))
struct_density = np.exp(kde.score_samples(x_grid.reshape(-1,1)))
struct_density_normalized_counts = len(struct)/np.sum(struct_density) * struct_density
fig['data'].append(
    go.Scatter(
        x = x_grid,
        y = struct_density_normalized_counts,
        mode='lines',
        line=dict(
            color=colors[0]
        ),
        name="PFAM families: "+str(len(struct))+" <br>with annotated structures"
    )
)


#kernel density estimate for Pfamilies without annotated structure
kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(no_struct.reshape(-1,1))
nostruct_density = np.exp(kde.score_samples(x_grid.reshape(-1,1)))
nostruct_density_normalized_counts = len(no_struct)/np.sum(nostruct_density) * nostruct_density
fig['data'].append(
    go.Scatter(
        x = x_grid,
        y = nostruct_density_normalized_counts,
        mode='lines',
        line=dict(
            color=colors[1]
        ),
        name="PFAM families: "+str(len(no_struct))+" <br>without annotated structures"
    )
)

#
# # #histogram for counts
# fig['data'].append(
#     go.Histogram(
#         x=struct,
#         histnorm="counts",
#         xbins=dict(start=np.min(struct), size=0.2, end=np.max(struct)),
#         opacity=0.4,
#         marker=dict(color=colors[0]),
#         name="Histogram of counts Pfamilies with structure",
#         showlegend=False
#     )
# )
# fig['data'].append(
#     go.Histogram(
#         x=no_struct,
#         histnorm="counts",
#         xbins=dict(start=np.min(no_struct), size=0.2, end=np.max(no_struct)),
#         opacity=0.4,
#         marker=dict(color=colors[1]),
#         name="Histogram of counts Pfamilies without structure",
#         showlegend=False
#     )
# )


#add vertical line for median of family size for families with structures
median_struct = np.median(struct)

fig['data'].append(
    go.Scatter(
        x=[median_struct, median_struct],
        y=[0, np.max([np.max(fig['data'][0]['y']), np.max(fig['data'][1]['y'])])],
        #y=[0, np.max([np.max(struct_density), np.max(nostruct_density)])],
        mode='lines+text',
        textfont=dict(
            family='sans serif',
            size=18,
            color=colors[0]
        ),
        text=["", "median family size: "+str(np.round(np.exp(median_struct), decimals=3))+" <br> (with annotated structure)"],
        textposition='right',
        line=dict(
            color=colors[0],
            width=4,
            dash='dash'),
        showlegend=False
    ))


#add vertical line for median of family size for families with NO structures
median_nostruct = np.median(no_struct)
fig['data'].append(
    go.Scatter(
        x=[median_nostruct, median_nostruct],
        y=[0, np.max([np.max(fig['data'][0]['y']), np.max(fig['data'][1]['y'])])],
        #y=[0, np.max([np.max(struct_density), np.max(nostruct_density)])],
        mode='lines+text',
        line=dict(
            color=colors[1],
            width=4,
            dash='dash'),
        textfont=dict(
            family='sans serif',
            size=18,
            color=colors[1]
        ),
        text=["", "median family size: "+str(np.round(np.exp(median_nostruct), decimals=3))+" <br> (without annotated structure)"],
        textposition="left",
        showlegend=False
    ))


#formatting of x-axis
fig['layout']['xaxis'].update(dict(
        title='family size [log number of sequences]',
        tickvals=np.log([10, 30, 100, 300, 1000, 10000, 100000, 500000]),
        ticktext=["10", "30", "100", "300", "1000", "10000", "100000", "500000"],
        exponentformat="e",
        showexponent='All',
    ))


#formatting of y-axis
fig['layout']['yaxis'].update(dict(
        title='Number of annotated structures in PDB',
        exponentformat="e",
        showexponent='All',
    ))

fig['layout']['title'] = "PFAM family sizes <br> Pfam 31.0 (March 2017, 16712 entries)"
fig['layout']['font'] = {}
fig['layout']['font']['size']= 18

#histogram bars will overlap
#fig['layout']['barmode'] = 'stack'

#finally plot that thing
plotly_plot(fig, filename=plot_out, auto_open=False)









from lxml import etree
import urllib


url="http://www.uniprot.org/statistics/TrEMBL%202017_05"


web = urllib.urlopen(url)
s = web.read()
html = etree.HTML(s)

## Get all 'tr'
tr_nodes = html.xpath('//table/tr')

## Get fourth element from table
print tr_nodes[4][0].text
print tr_nodes[4][1].xpath('a')[0].text

