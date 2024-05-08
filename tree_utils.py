import numpy as np
import networkx as nx
from ete3 import Tree
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import random




def create_tree(n_leaves, scale=0.1):
    """
    :param n_leaves: Number of leaves
    :param branch: Branch information
    :param scale: used when branch is 'random'
    :return: returns the tree
    """
    n_nodes = 2 * n_leaves - 1
    tree = nx.Graph()
    node_names = [n for n in range(n_nodes)]
    leaves = node_names[:n_leaves]
    root = node_names[-1]
    next_child = root - 1
    opt = np.zeros((n_nodes, n_nodes))
    print("name:",node_names)
    for node in node_names[::-1]:  # iterate backwards
        if node in leaves:
            tree.add_node(node, type='leaf')
        else:
            # connect children with node
            child_1 = next_child
            t_1 = np.random.exponential(scale)

            opt[node][child_1]=t_1

            tree.add_edge(node, child_1, t=t_1)
            child_2 = child_1 - 1
            t_2 = np.random.exponential(scale)

            opt[node][child_2]=t_2

            tree.add_edge(node, child_2, t=t_2)
            if node == root:
                tree.add_node(node, type='root')
            else:
                tree.add_node(node, type='internal')
            # store ancestral data
            tree.add_node(node, children=[child_1, child_2])
            tree.add_node(child_1, parent=node, t=t_1)  # add branch length to parent
            tree.add_node(child_2, parent=node, t=t_2)
            # child for next node
            next_child -= 2
    for i in range(n_nodes):
        for j in range(n_nodes):
            if opt[i][j]==0:
                opt[i][j]= np.random.exponential(scale)#not sure if this should be completedly randomized
    #draw_dendrogram(tree, root)
    return tree, opt


def update_topology(tree, root):
    # add needed meta data to nodes: parent, children, type and branch length to parent
    n_nodes = len(tree)
    n_leaves = (n_nodes + 1) // 2
    postorder_nodes = list(nx.dfs_postorder_nodes(tree, root))
    edge_length=np.array([0.]*n_nodes)
    visited = []

    for node in postorder_nodes:
        if node < n_leaves:
            tree.add_node(node, type='leaf')
            visited.append(node)
            continue
        elif not node == root:
            children = [child for child in tree.adj[node] if child in visited]
            tree.add_node(node, type='internal', children=children)
        else:
            children = [child for child in tree.adj[node] if child in visited]
            tree.add_node(node, type='root', children=children)
        for child in children:
            t = tree.edges[(node, child)]['t']
            edge_length[child]=t
            tree.add_node(child, t=t, parent=node)
        visited.append(node)
    #draw_dendrogram(tree, root)
    return tree,edge_length


def draw_dendrogram(tree, root):
    print(nx2ete(tree, root))


def print_newick(tree, root):
    print(nx2ete(tree, root).write())


def nx2ete(graph, root):

    print("graph:",len(graph))
    print("root",root)
    nodes = len(graph)
    
    tree = Tree(f"{root};")
    topology = nodes * [None]
    has_parent = nodes * [False]
    has_parent[root] = True
    edges = list(graph.edges())
    for u, v in edges:        
        if has_parent[u]:
            if u == root:
                child = tree.add_child(name=v)
                child.dist = nx.shortest_path_length(graph, u, v, weight='t')
                topology[v] = child
            else:
                child = topology[u].add_child(name=v)
                child.dist = nx.shortest_path_length(graph, u, v, weight='t')
                topology[v] = child
            has_parent[v] = True
        elif has_parent[v]:
            if v == root:
                child = tree.add_child(name=u)
                child.dist = nx.shortest_path_length(graph, u, v, weight='t')
                topology[u] = child
            else:
                child = topology[v].add_child(name=u)
                child.dist = nx.shortest_path_length(graph, u, v, weight='t')
                topology[u] = child
            has_parent[u] = True
        else:
            edges.append((u, v))
    tree.children[0].get_tree_root().name = root
    # add branch lengths
    for node in tree.traverse('postorder'):
        if node.is_root():
            node.name = root
            node.dist = 0.
        else:
            node.dist = 1.
    return tree




