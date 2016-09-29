import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
from plotly.tools import FigureFactory as FF


def plot_contact_map_someScore_plotly(plot_matrix, name, L, N, seqsep, plot_file):
    
    #sort matrix by confidence score
    plot_matrix.sort('confidence', ascending=False, inplace=True)
    
    #define triangle of plot that represents Predictions 
    heatmap_predicted = go.Heatmap(x=plot_matrix.residue_i.tolist(), 
                                   y=plot_matrix.residue_j.tolist(), 
                                   z=plot_matrix.confidence.tolist(), 
                                   name='predicted', 
                                   colorscale='Greys', reversescale=True, 
                                   colorbar=plotly.graph_objs.ColorBar(x=1.02, 
                                                                       y=0.4, 
                                                                       yanchor='bottom', 
                                                                       len=0.4,
                                                                       title="Score"))

    
    #if distances and class are available
    if 'class' in plot_matrix and 'distance' in plot_matrix:
        
        #define true and false positives among the L/5 highest scores
        coord_tp = [plot_matrix[:L/5].loc[plot_matrix['class'] == True, 'residue_i'].tolist() + plot_matrix[:L/5].loc[plot_matrix['class'] == True, 'residue_j'].tolist(), plot_matrix[:L/5].loc[plot_matrix['class'] == True, 'residue_j'].tolist()+ plot_matrix[:L/5].loc[plot_matrix['class'] == True, 'residue_i'].tolist()]
        size_tp  = plot_matrix[:L/5].loc[plot_matrix['class'] == True, 'confidence'].tolist() + plot_matrix[:L/5].loc[plot_matrix['class'] == True, 'confidence'].tolist()
        coord_fp = [plot_matrix[:L/5].loc[plot_matrix['class'] == False, 'residue_i'].tolist() + plot_matrix[:L/5].loc[plot_matrix['class'] == False, 'residue_j'].tolist(), plot_matrix[:L/5].loc[plot_matrix['class'] == False, 'residue_j'].tolist()+ plot_matrix[:L/5].loc[plot_matrix['class'] == False, 'residue_i'].tolist()]
        size_fp  = plot_matrix[:L/5].loc[plot_matrix['class'] == False, 'confidence'].tolist() + plot_matrix[:L/5].loc[plot_matrix['class'] == False, 'confidence'].tolist()

        #Mark TP and FP in the plot with little crosses
        tp = go.Scatter(x = coord_tp[0], 
                        y = coord_tp[1], 
                        mode = 'markers', 
                        marker = dict(symbol=134,  color="green", line=dict(width=2), size=12),#size_tp, sizeref=np.max([size_tp + size_fp])/15, sizemode = 'diameter'), 
                        name="TP (L/5)",
                        hoverinfo="none")
        fp = go.Scatter(x = coord_fp[0], 
                        y = coord_fp[1], 
                        mode = 'markers', 
                        marker = dict(symbol=134, color="red", line=dict(width=2), size=12),#size_fp, sizeref=np.max([size_tp + size_fp])/15, sizemode = 'diameter'), 
                        name="FP (L/5)",
                        hoverinfo="none")
        
        #define triangle on opposite site of Predictions 
        heatmap_observed = go.Heatmap(x=plot_matrix.residue_j.tolist(), 
                                      y=plot_matrix.residue_i.tolist(), 
                                      z=plot_matrix.distance.tolist(), 
                                      name='observed',  
                                      colorscale='RdBu' , reversescale=True,
                                      colorbar=plotly.graph_objs.ColorBar(x=1.02, 
                                                                          y=0, 
                                                                          yanchor='bottom', 
                                                                          len=0.4, 
                                                                          title="Distance [A]"))
        #put all plot elements in data list
        data = [heatmap_predicted,
                heatmap_observed,
                tp,
                fp,
                go.Scatter(x = [0,L], y = [0,L],                mode = 'lines', line = dict(color = ('rgb(0, 0, 0)'), width = 4),  showlegend=False),
                go.Scatter(x = [0,L-seqsep+1], y = [seqsep-1,L],mode = 'lines', line = dict(color = ('rgb(0, 0, 0)'), width = 2),  showlegend=False),
                go.Scatter(x = [seqsep-1,L], y = [0,L-seqsep+1],mode = 'lines', line = dict(color = ('rgb(0, 0, 0)'), width = 2),  showlegend=False)
                ]
    else:
        
        #if no pdb data is available, plot Predictions on both sies in the matrix
        heatmap_predicted_opposite = go.Heatmap(x=plot_matrix.residue_j.tolist(), 
                                                y=plot_matrix.residue_i.tolist(), 
                                                z=plot_matrix.confidence.tolist(), 
                                                name='predicted', 
                                                colorscale='Greys', reversescale=True, 
                                                showscale = False) #do not show the grey colorscale twice
        
        #put all plot elements in data list
        data = [heatmap_predicted,
                heatmap_predicted_opposite,
                go.Scatter(x = [0,L], y = [0,L],                mode = 'lines', line = dict(color = ('rgb(0, 0, 0)'), width = 4),  showlegend=False),
                go.Scatter(x = [0,L-seqsep+1], y = [seqsep-1,L],mode = 'lines', line = dict(color = ('rgb(0, 0, 0)'), width = 2),  showlegend=False),
                go.Scatter(x = [seqsep-1,L], y = [0,L-seqsep+1],mode = 'lines', line = dict(color = ('rgb(0, 0, 0)'), width = 2),  showlegend=False)
                ]

    

    #plot_contact_map_someScore
    plotly.offline.plot({
    "data": data,
    "layout": go.Layout(
        title = name + ", L: " + str(L) + ", N: " + str(N) + ", seq. separation: " + str(seqsep),
        width = 1000,
        height = 1000,
        xaxis1 = dict(range=[0.5,L+0.5], #removes white margin in plot
                   title='j'),  
        yaxis1 = dict(range=[0.5,L+0.5], #removes white margin in plot
                   title='i'),  
        legend = dict(x=1.02,y=1)       #places legend to the right of plot
    )
    }, filename=plot_file, auto_open=False)
    
    
