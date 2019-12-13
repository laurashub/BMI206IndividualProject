import sys
import xml.etree.ElementTree as ET
import xlrd

class Reaction:
    def __init__(self, name, reactants, products):
        self.name = name
        self.reactants = reactants
        self.products = products
        self.out_edges = []
        self.in_edges = []
        self.betweenness_centrality = -1
        self.bridging_coefficient = -1
        self.clustering_coefficient = -1
        self.essential = None

def read_essential_rxns():
    ers = {}
    essential_rxns = xlrd.open_workbook('essential_rxns.xls')
    sheet = essential_rxns.sheet_by_index(0)
    for col_index in range(1, 8):
        model = sheet.cell(0, col_index).value
        rxns = [sheet.cell(row_index, col_index).value for row_index in range(1, sheet.nrows)]
        ers[model] = rxns
    return ers

def construct_graph(filename):
    rxns = []
    arcs = set()

    tree = ET.parse(filename)
    root = tree.getroot()

    for child in root[0][3]:
        reaction_id, reactants, products = '', [], []
        reaction_id = child.attrib['id'].strip('r_')
        for reactant in child[1]:
            reactants.append(reactant.attrib['species'].split('_')[1])
        for product in child[2]:
            products.append(product.attrib['species'].split('_')[1])
        if reactants != products:
            rxns.append(Reaction(reaction_id, reactants, products))

    for reaction1 in rxns:
        for reaction2 in rxns:
            if reaction1 != reaction2:
                if True in [product in reaction2.reactants for product in reaction1.products]:
                    reaction1.out_edges.append(reaction2)
                    reaction2.in_edges.append(reaction1)
                    arcs.add((reaction1, reaction2))
    return rxns, arcs

def betweenness_centrality(all_nodes):
    bcs = {} #empty dictionary, will hold betcens
    for s in all_nodes:
        bcs[s] = 0 #initalize all to 0
    for j in all_nodes:
        S = [] #empty stack
        P = {}
        sigma = {}
        d = {}
        delta = {}
        for k in all_nodes:
            P[k] = [] #empty list
            sigma[k] = 0
            d[k] = -1
            delta[k] = 0 #delta[v] <- 0 for v in V
        sigma[j] = 1 #sigma[s] <- 1
        d[j] = 0 #d[s] <- 0
        Q = [] # empty queue
        Q.append(j) #enqueue s -> Q
        while Q:
            v = Q[0] # dequeue v <- Q
            Q.remove(v)
            S.append(v) #push v -> S
            for w in v.out_edges:
                if d[w] < 0:
                    Q.append(w) #enqueue w -> Q
                    d[w] = d[v] + 1 #update distance
                if d[w] == (d[v] + 1):
                    sigma[w] += sigma[v]
                    P[w].append(v) #append v -> P[w]
        while S:
            w = S[-1] # pop w <- S
            S.remove(w)
            for v in P[w]:
                delta[v] += (sigma[v]/sigma[w])*(1+delta[w])
            if w != j:
                bcs[w] += delta[w]
    return bcs


#degree_total(i)^-1/ sum of total degrees of neighbors^-1
def bridging_coefficient(node):
    if len(node.in_edges) == 0 and len(node.out_edges) == 0:
        return 0
    num = total_degree(node)**-1
    denom = 0
    for neighbor in get_neighborhood(node):
        denom += total_degree(neighbor)**-1
    return num/denom
        
def total_degree(node):
    return len(node.in_edges) + len(node.out_edges)

def get_neighborhood(node):
    return set(node.out_edges + node.in_edges)

#C(i) = n_i/(k_i * (k_i - 1))
#n_i is the number of arcs between neighborhood of the node i
#k_i is the number of neighbors of node i
def clustering_coefficient(node, all_nodes):
    neighborhood = get_neighborhood(node)
    k_i = len(neighborhood)
    if k_i == 0 or k_i == 1:
        return None
    n_i = 0
    for j in neighborhood:
        for k in neighborhood:
            if j != k:
                if k in j.out_edges:
                    n_i +=1
    return n_i/(k_i * (k_i - 1))


#read in and populate graph
essential_rxns = read_essential_rxns()
print(essential_rxns.keys())

#calculate metrics
graphs = {}
for fname in essential_rxns.keys():
    print(fname)
    graphs[fname] = construct_graph(fname + '.xml')
    rxns, arcs = graphs[fname]
    print(len(rxns), len(arcs))
    betcens = betweenness_centrality(rxns)
    print("Bet cens done")
    for rxn in betcens.keys():
        rxn.betweenness_centrality = betcens[rxn]
        rxn.essential = rxn.name in essential_rxns[fname]
        rxn.bridging_coefficient = bridging_coefficient(rxn)
        rxn.bridging_centrality = rxn.betweenness_centrality * rxn.bridging_coefficient
        rxn.clustering_coefficient = clustering_coefficient(rxn, rxns)
    print(sum([rxn.essential for rxn in rxns]))

#write out to file
for fname in essential_rxns.keys():
    with open(fname + '_bridging_centrality.txt', 'w') as f:
        for rxn in sorted(graphs[fname][0], key=lambda x: x.bridging_centrality, reverse = True):
            f.write("{} {} {}\n".format(rxn.name, rxn.bridging_centrality, rxn.essential))
    with open(fname + '_betweenness_centrality.txt', 'w') as f:
        for rxn in sorted(graphs[fname][0], key=lambda x: x.betweenness_centrality, reverse = True):
            f.write("{} {} {}\n".format(rxn.name, rxn.betweenness_centrality, rxn.essential))
    with open(fname + '_clustering_coefficient.txt', 'w') as f:
        for rxn in sorted(list(filter(lambda x:x.clustering_coefficient is not None, graphs[fname][0])), key=lambda x: x.clustering_coefficient):
            f.write("{} {} {}\n".format(rxn.name, rxn.clustering_coefficient, rxn.essential))








