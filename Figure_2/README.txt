Code Descriptions

This folder contains the code to compute the stresses in the cracked ice sheet to produce figure 2 in the manuscript.

How to run it

* Generate meshes with MeshGeneration.py in the folder "Mesh". Two meshes are already generated and saved here, with cracklengths = 0.3 and 0.5. In the main program the C=0.3 case is used.

* Run the main program "Main.ipynb" in jupyter notebook of FEniCS in a docker container. This script solves the displacement and stresses (saved in folder "Notched_Results").

* The visualisation is done in Make_Figures.ipynb in the folder "Notched_Results". Here the perturbation stress T is plotted with matplotlib.
