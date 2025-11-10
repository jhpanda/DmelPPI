
import csv,sys
import numpy as np
import pandas as pd
import networkx as nx
import markov_clustering as mc
import matplotlib.pyplot as plt
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout

import os,sys
fontfile = "/ru-auth/local/home/jpeng/.cache/matplotlib/fontlist-v330.json"
if os.path.isfile(fontfile):
    print(fontfile)
    os.remove("/ru-auth/local/home/jpeng/.cache/matplotlib/fontlist-v330.json")
import matplotlib.font_manager
matplotlib.get_cachedir()
matplotlib.rcParams['font.family'] = ['arial']

## important genes to label in networks ##

#important_genes = ["bys","pit","Sas10","CG9398","UQCR-C1","ND-75","Ance-2","Ance-3","Spn100A"]
important_genes = []
df = pd.read_csv("fbpp2fbgnlist.txt")
fbpp2symbol = dict(zip(df["FBpp"],df["Symbol"]))

def plot_graph(G,connectivity,figout,sizemultiper=1.5,width=6,height=5,node_labels=False,layout="spring"):
    colors_node = [connectivity[node]+3.5 for node in G.nodes]
    colors_edge = [G.edges[v]['weight'] for v in G.edges]
    #size_node = [sizemultiper+connectivity[node]*1.5 for node in G.nodes]
    #size_node = [sizemultiper*2 for node in G.nodes]
    if layout=="graphviz":
        size_node = [sizemultiper*connectivity[node] for node in G.nodes]
        #size_node = [sizemultiper+connectivity[node]*1.5 for node in G.nodes]
    else:
        size_node = [sizemultiper*connectivity[node] for node in G.nodes]

    if layout=="graphviz":
        #nodes_label = [fbpp2symbol[n] for n in G.nodes]
        nodes_label = {n:fbpp2symbol[n] for n in G.nodes if fbpp2symbol[n] in important_genes}
    else:
        nodes_label = {}

    if node_labels:
        nodes_label = {n:fbpp2symbol[n] for n in G.nodes}
        print(nodes_label)            

    nmax = max(connectivity.values())
    nmin = min(connectivity.values())

    if layout=="spring":
        pos = nx.spring_layout(G,k=0.03)
    elif layout=="graphviz":
        pos = graphviz_layout(G)


    plt.figure(figsize=(width,height))

    ec  = nx.draw_networkx_edges(G,pos,edgelist=G.edges,
                                    edge_color=colors_edge,
                                    #edge_color='gray',
                                    alpha=0.5,
                                    width=1.,
                                    edge_vmin=0,
                                    edge_vmax=1,
                                    edge_cmap=plt.cm.Greys)
                                    #)

    nc  = nx.draw_networkx_nodes(G,pos,nodelist=G.nodes,
                                    node_color=colors_node,
                                    #labels=nodes_label,
                                    #node_size=10,
                                    vmin=1,
                                    vmax=10,
                                    node_size=size_node,
                                    cmap=plt.cm.plasma)
    lc  = nx.draw_networkx_labels(G, pos, labels=nodes_label,font_size=6), 

    plt.axis('off')
    plt.tight_layout()
    #plt.savefig('environmental_adaptation_graph.pdf',dpi=300)
    plt.savefig(figout,dpi=300)
    plt.show()
    plt.close()

def network_mcl(network,ntop=10)
    # get the adjacency matrix (in sparse form)
    matrix = nx.to_scipy_sparse_matrix(network.G)
    result = mc.run_mcl(matrix,inflation=1.5) # run MCL with default parameters
    clusters = mc.get_clusters(result) # get clusters

    #print(result)
    sizes = []
    for cluster in clusters:
        #print(len(cluster))
        sizes.append(len(cluster))

    clusters = sorted(clusters,key=lambda x:len(x),reverse=True)
    nodes = [n for n in network.G.nodes]
    top_nodes = []
    top_clusters = []
    cluster_id = {}
    ntop = 10
    nclusters = len(clusters)
    #ntop = nclusters
    pref = f"{ntop}" if ntop<nclusters else "all"
    #rep_nodes = ["" for i in range(ntop)]
    #rep_nodes2 = ["" for i in range(ntop)]
    rep_nodes = []
    clusters_list = {}
    f = open(f"network.mcl.{pref}.txt","w")
    for i in range(ntop):
        cluster = [nodes[s] for s in clusters[i]]
        clusters_list[i] = cluster
        f.write("cluster %d %s\n"%(i+1,",".join(cluster)))
        print(i,len(cluster),", ".join([fbpp2symbol[c] for c in cluster]))
        top_clusters.append(cluster)

        cluster = sorted(cluster,key=lambda x: network.connectivity[x],reverse=True)
        #rep_nodes[i] = cluster[0]
        #rep_nodes2[i] = cluster[1]
        rep_nodes += cluster[0:2]
        for c in cluster:
            cluster_id[c] = i+1
            top_nodes.append(c)
        #print(n1,n2)
    f.close()

    subG = network.G.subgraph(top_nodes).copy()
    colors_node = [cluster_id[node] for node in subG.nodes]
    sizemultiper= 1.5
    size_node = [sizemultiper*network.connectivity[node] for node in subG.nodes]
    #size_node = [5 for node in subG.nodes]
    #nodes_label = {n:fbpp2symbol[n] for n in subG.nodes}
    nodes_label = dict([[n,""] if n not in rep_nodes else [n,fbpp2symbol[n]] for n in subG.nodes])

    plt.figure(figsize=(6,5))
    pos = graphviz_layout(subG)

    ec  = nx.draw_networkx_edges(subG,pos,edgelist=subG.edges,
                                    #edge_color=colors_edge,
                                    #edge_color='gray',
                                    alpha=0.66,
                                    width=.2,
                                    #edge_vmin=0,
                                    #edge_vmax=1,
                                    #edge_cmap=plt.cm.Greys)
                                    )

    nc  = nx.draw_networkx_nodes(subG,pos,nodelist=subG.nodes,
                                    node_color=colors_node,
                                    #labels=nodes_label,
                                    #node_size=10,
                                    vmin=1,
                                    vmax=ntop,
                                    node_size=size_node,
                                    cmap=plt.cm.cool)
    lc  = nx.draw_networkx_labels(subG,pos,labels=nodes_label,font_size=6),

    fig.tight_layout()
    plt.axis("off")
    plt.savefig(f"network.mcl.{pref}.pdf")

    return clusters


## the network class ##
class Network():
    def __init__(self):
        self.G = nx.Graph()
        self.connectivity = {}

    def load_alphafold(self,af_output,cut=0.23):
        df = pd.read_csv(af_output)

        for pair,pdockq in zip(df["pair"],df["pdockq"]):

            if pdockq>=cut:
                p1,p2 = pair.split("_")
                p1 = p1[5:]
                p2 = p2[5:]

                #w = 0.5 * np.exp(3*(pdockq-0.23))
                w = pdockq
                self.G.add_nodes_from([p1,p2])
                self.G.add_edge(p1, p2, weight=w )

    def report_connectivity(self,out):
        f = open(out,"w")
        f.write("## Connectivity of the predictions ##\n")
        mean_c = 0 ## mean connection
        for node in self.G.nodes:
            neighbors = [s for s in self.G.neighbors(node)]
            nc = len(neighbors)
            self.connectivity[node] = nc
            mean_c += nc
        mean_c = mean_c/len(self.G.nodes)
        self.connectivity = dict(sorted(self.connectivity.items(),
                                        key=lambda x:x[1],reverse=True))

        for node in self.connectivity:
            neighbors = [s for s in self.G.neighbors(node)]
            nc = self.connectivity[node]
            f.write("%s %d %s\n"%(node,nc,";".join(neighbors)))

        print("average number of partners: %.1f"%(mean_c))
        f.write("#average number of partners: %.1f\n"%(mean_c))
        f.close()

    def report_subnetworks(self,out,fmap=None):
        if fmap:
            dft = pd.read_csv("fbpp2fbgnlist.txt",sep=",")
            map_dict = dict(zip(dft["FBpp"],dft["FBgn"]))

        f = open(out,"w")
        subgraphs = nx.connected_components(self.G)
        i = 1
        
        #subgraphs = [tuple(s) for s in list(subgraphs)]
        #subnodes  = [len(s) for s in subgraphs]
        #g_dict = dict(zip(subgraphs,subnodes))
        subnodes = [tuple(s) for s in list(subgraphs)]
        subnodes_num = [len(s) for s in subnodes]
        subgraphs = [self.G.subgraph(s) for s in subnodes]
        #print("subnodes",subnodes)
        #print("subnum",subnodes_num)
        list2sort = []
        subgraphs_dict = {}
        for nodes,num,subG in zip(subnodes,subnodes_num,subgraphs):
            for n in nodes:
                subgraphs_dict[n] = subG
            list2sort += [[nodes,num,subG]]
        self.subgraphs_dict = subgraphs_dict ## easy to lookup subgraphs

        #print(list2sort)
        #print(list2sort[0])
        for nodes,num,subG in sorted(list2sort,key=lambda x:x[1],reverse=True):
            #print(nodes)
            #print(num)
            #print(subG)
            nodes = [s for s in nodes]
            nodes_connectivity = [self.connectivity[s] for s in nodes]

            list_connectivity = [[n,c] for n,c in zip(nodes,nodes_connectivity)]
            list_connectivity = dict(sorted(list_connectivity,key=lambda x:x[1],reverse=True))

            nodes = [s for s in list_connectivity]

            num_edges = subG.number_of_edges()
            #if fmap:
            #    nodes = [map_dict[s] for s in nodes]
            #    nodes = [s for s in nodes]

            print('cluster %4d %4d %4d'%(i,num,num_edges),','.join(nodes))
            f.write("cluster %4d %4d %4d%s\n"%(i,num,num_edges,','.join(nodes)))
            i += 1
        f.close()

    def plot(self,figout):
        G = self.G
        plot_graph(G,self.connectivity,figout,width=4,height=4)

    def plot_subnetwork(self,gene,pref,sizemultiper,node_labels=False):
        #subgraphs = nx.connected_components(self.G)
        #subgraphs = [tuple(s) for s in subgraphs]
        ##subnodes  = [len(s) for s in subgraphs]

        subG = None
        #for subg in subgraphs:
        #    subnodes = {}
        #    for p in subg:
        #        subnodes[p] = 1
        #    subnodes = [n for n in subnodes]
        #    if gene in subnodes:
        #        subG = self.G.subgraph(subnodes)
#
        #subgraphs = nx.connected_components(self.G)
        if gene in self.subgraphs_dict:
            subG = self.subgraphs_dict[gene]

        if subG==None:
            print("Error! %s not in any subclusters"%gene)
        else:
            figout = "subcluster%s_with_%s.png"%(pref,gene)
            G = subG
            connectivity = {}
            for node in G.nodes:
                neighbors = [s for s in G.neighbors(node)]
                nc = len(neighbors)
                connectivity[node] = nc
            plot_graph(G,connectivity,figout,sizemultiper=sizemultiper,
                        width=4,height=4,node_labels=node_labels,layout="graphviz")


if __name__ == "__main__":

    cut = 0.5
    af_output = "pae_summary.csv"
    subcluster_out = "result_subcluster_%.1f.txt"%cut
    connect_out = "result_connectivity_%.1f.txt"%cut
    figout = "result_network_%.1f.pdf"%cut
    network = Network()
    network.load_alphafold(af_output,cut=cut)
    print(network.G)
    network.report_connectivity(connect_out)
    network.report_subnetworks(subcluster_out,fmap="fbpp2fbgnlist.txt")
    network.plot(figout)

    # MCL clustering #
    network_mcl(network,ntop=10)