# VRSEMPHY
Python implementation of the varied rate of evolution structural EM algorithm(VRSEMPHY), this implementation is thanks to Authors: Oskar Kviman & Niharika Gauraha KTH E-mail: {okviman, niharika}@kth.se. Through the python3 sem.py , the program can be run and the exact number of species and sites can be adjusted in the main method. 
computation.py stores methods for computing the marginal/joint probability(upward and downward messages), log-likelihood, etc... Intuitively, plot_util.py deals with plotting the log-likelihood graph, and data_util.py deals with simulating data, and mst_util deals with finding the maximum spanning tree and the conversion to the bifurcating graph. substituion_methods entail the assumed model of evolution and tree_util.py updates the tree topology

There is no specific requirement for the system, as most computers should suffice, as I am running on a particularly non-power device. Libraries needed for import would be networkx, random, and ete3.  
