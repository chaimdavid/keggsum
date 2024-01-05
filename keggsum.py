#Import of necessary packages
from Bio.KEGG.KGML.KGML_parser import read
from Bio.KEGG.REST import *
import networkx as nx
from networkx.algorithms.approximation import *

#Definition of the KEGGSum algorithm, the user can select the percentage of the nodes the algorithm will keep and
#the centrality measure it will use. The defaults are 15% and PageRank respectively
def keggsum(kid,perc=15,cent="PageRank"):
    #Creation of an object that will store the kgml data downloaded from the KEGG database through their REST-like API
    pathway = read(kegg_get(kid, option="kgml"))
    #Creation of a list where all the nodes of the graph will be stored
    nodes_list = []
    #Storing of the nodes names avoiding to log the ones with "compound", "ortholog" type and "undefined" name
    for entry in pathway.entries.values():
        if entry.type != "compound" and entry.type != "ortholog" and entry.name != "undefined":
            nodes_list.append(entry.name)
    #Creation of lists where the nodes for each "relation"(edge) will be stored
    entries1 = []
    entries2 = []
    for relation in pathway.relations:
        entries1.append(relation.entry1)
        entries2.append(relation.entry2)
    #Creation of a list where the pairs of nodes for each edge will be stored
    tuple_lst = list(zip(entries1, entries2))
    #Exclusion of "compound" type nodes from the "relations" list
    tuple_lst_valid = []
    for (a, b) in tuple_lst:
        if a.type != "compound" and b.type != "compound":
            tuple_lst_valid.append((a, b))
    #Creation of the "graph" object and registration of the graph's edges
    G = nx.Graph()
    for (a, b) in tuple_lst_valid:
        G.add_edge(a.name, b.name)
    #Exclusion of the "undefined" nodes
    for i in list(G):
        if i == "undefined":
            G.remove_node(i)
    #Creation of a dictionary with acceptable centrality options
    cents_dict = {"Betweenness": nx.betweenness_centrality(G), "Degree": nx.degree_centrality(G),
                  "Closeness": nx.closeness_centrality(G),
                  "PageRank": nx.pagerank(G), "Katz": nx.katz_centrality(G),
                  "Eigenvector": nx.eigenvector_centrality(G, max_iter=600),
                  "Harmonic": nx.harmonic_centrality(G)}
    if cent in cents_dict:
        #Calculation of the graph's most important number of nodes according to the percentage of the
        # initial population of nodes the user selected.Else, the default percentage is set to 15%
        impor_perc = round(len(nodes_list) * perc / 100)
        print('\nThe node with the highest', cent, 'centrality is: ' + max(cents_dict[cent], key=cents_dict[cent].get))
        #Sorting the nodes by decreasing importance
        sorted_imps = sorted(cents_dict[cent].items(), key=lambda x: x[1], reverse=True)
        #Creation of the list where the most important nodes of the graph will be stored
        imp_nodes = []
        #Separation of the most important nodes of the graph
        for i in range(0, impor_perc):
            imp_nodes.append(sorted_imps[i][0])
        print('\nThe nodes of graph ' + kid + ' with the top', impor_perc,
              cent, 'centralities are: \n{}.'.format('\n'.join(imp_nodes)))
        #Creation of a subgraph that contains only the connected nodes
        largest_component = max(nx.connected_components(G), key=len)
        G2 = G.subgraph(largest_component)
        #Connection of the nodes with the Steiner Tree algorithm
        graph = steiner_tree(G2, imp_nodes)
        #Visualization of the summary graph
        nx.draw_networkx(graph)
    elif cent:
        raise ValueError("Invalid centrality option")
