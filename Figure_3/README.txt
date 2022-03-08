Code description

This folder contains the code to calculate the stress intensity factors K1 and K2 for vertical basal crevasses with specific crack lengths.

How to run it?

* Generate meshes with MeshGeneration.py in the folder "Notched_Mesh".

* Run the main program "Main_Notched.ipynb" in jupyter notebook of FEniCS in a docker container. This script computes the SIFs of vertical basal hydrofractures with varying cracklengths. The SIFs against cracklengths are saved in "SIF_deltatau02.txt" and "SIF_deltatau03.txt".

* The visualisation is done by the matlab code "SIF_Zc.m".
