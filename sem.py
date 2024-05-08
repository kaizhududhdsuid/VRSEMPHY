import numpy as np
from tree_utils import create_tree, draw_dendrogram
from substitution_models import JukesCantor
from data_utils import get_char_data, simulate_seq
import computations as cmp
from mst_utils import get_mst, bifurcate_mst
import plot_utils as plot
import networkx as nx
from ete3 import Tree
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import random

def randomt_tree(n_leaves, scale=0.1):
    
    n_nodes = 2 * n_leaves - 1
    tree = nx.Graph()
    node_names = [n for n in range(n_nodes)]
    
    root = node_names[-1] 
    opt = np.zeros((n_nodes, n_nodes))
    for i in range(n_nodes):
        for j in range(n_nodes):
            if opt[i][j]==0:
                opt[i][j]= np.random.exponential(scale)#not sure if this should be completedly randomized
    
    weight_matrix = np.zeros((n_nodes, n_nodes))

    for node1 in range(n_nodes-1):     
        for node2 in range(n_nodes-1):
            weight_matrix[node1, node2] = random.uniform(-25, 0)

    weight_matrix[weight_matrix == 0] = -np.infty
    tree = get_mst(weight_matrix, opt)
    T_0,opt=create_tree(n_leaves)
    leaves = [n for n in T_0 if T_0.nodes[n]['type'] == 'leaf']
    T_l,edge = bifurcate_mst(tree, leaves, root)
    draw_dendrogram(tree, root)
    
    #print("root2",root)
    
    
    return T_l, opt

def exp_rate(opt, prob):
    n_nodes=len(opt[0])
    expected=np.zeros((n_nodes, n_nodes, len(prob[0][0])))
    for i in range(n_nodes):

        for j in range(n_nodes):
            t = opt[i][j]
            expected[i][j]=t*(0.1*prob[i][j]+10*1-(prob[i][j]))        
    return expected


def decode(opt, expect):
    n_nodes=len(opt[0])
    prob=np.zeros((n_nodes, n_nodes, len(expect[0][0])))
    for i in range(n_nodes):
        for j in range(n_nodes):
            if opt[i][j]==0:
                prob[i][j]=0.99
            else:
                prob[i][j]=((expect[i][j]/opt[i][j])-10)/(-9.9)
                # prob[i][j]=((expect[i][j]/opt[i][j])-10)/(-9.9)
    # for edge in opt:
    # 0.1x +10(1-x)=0.1x+10-10x=-9.9x+10
    #     prob[edge]=expect[edge]/ (tree.edges[edge]['t'])

    return prob

def sem(data, evol_mod, prob, perturb=True):
    n_leaves, n_sites = data.shape

    # create tree
    #T_0, opt = create_tree(n_leaves)  # init tree topology



    # T_0,opt=create_tree(n_leaves)


    T_0, opt=randomt_tree(n_leaves)
    leaves = [n for n in T_0 if T_0.nodes[n]['type'] == 'leaf']
    root = len(T_0) - 1

    # annealing variables
    sigma_l = 0.1  # annealing std
    rho = 0.95  # cooling factor

    # EM
    max_iter = 100
    T_l = T_0  # init topology
    ll_vec = []  # log-likelihood
    Tl_vec=[]
    for iter in range(max_iter):
        # E-step 
        Tl_vec.append(T_l)
        up_table = cmp.compute_up_messages(data, T_l, evol_mod,prob)

       
        down_table = cmp.compute_down_messages(data, T_l, evol_mod, up_table, prob)
        # compute log-likelihood
        ll_sites, ll = cmp.compute_loglikelihood(up_table, evol_mod.stat_prob)
        print(ll)
        if len(ll_vec)>0 and ll_vec[0]>ll:#figure out argmax maybe
            T_l=Tl_vec[0]
            ll_vec.append(ll_vec[0])
            continue
        ll_vec.append(ll)
        
        #obtain expected number of change trhough timing edge length and rate of change 
        expected_rate=exp_rate(opt, prob)
         
        # M-step
        W, t_opts = cmp.compute_weight_matrix(T_l, evol_mod, up_table, down_table)
        #work with t_opt, it has the length between each node
        

        old=prob
        prob=decode(t_opts, expected_rate)               
        if perturb:
            
            W_tilde = cmp.perturb_W(W, sigma_l)
        else:
            W_tilde = W
        mst = get_mst(W_tilde, t_opts)
        T_l,edge = bifurcate_mst(mst, leaves, root)
        t_opts=opt



       

        
        # print("edgeshape:",edge.shape)
        # print("probshp:",prob.shape)

        #Given a tree, I can find the edge and its corresponding loglikelihood, from there I can infer upward/downward messages->
        # update annealing temperature
        #Find Edge Length------Transition Matrix
        #obtain the edge length, calculate new expectd change at each site
        
        
        sigma_l *= rho
        if sigma_l <= 5e-3:
            break
        elif iter > 0 and np.abs(ll_vec[iter] - ll_vec[iter - 1]) < 1e-3:
            break

    return T_l, ll_vec


def main(load_data=False):
    #two bird with one stone
    # init substitution model
    jc = JukesCantor(alpha=0.1)
    true_ll = None
    
    
    
    n_leaves = 35
    tree, opt = create_tree(n_leaves)
    data,prob = simulate_seq(tree, jc, 1000)
        
        


        # u
        # p_table = cmp.compute_up_messages(data, tree, jc, prob)
        # ll_sites, true_ll = cmp.compute_loglikelihood(up_table, jc.stat_prob)
        #With the new topology u cab iterate through each nucleiotide  Psedoucount, a for loop for each
        #link length < 1 approximate # of change 
          

    T_l, ll_vec = sem(data, jc, prob, True)
    plot.plot_loglikelihood(ll_vec)


if __name__ == '__main__':
    main()
