# Basal hydrofractures near sticky patches

This repository contains the code for the manuscript "Basal hydrofracture near sticky patches", which is planned to be submitted to [Journal of Galciology](https://www.cambridge.org/core/journals/journal-of-glaciology). The work aims to study the mode-I and mixed-mode hydrofracture propagation near the sticky patches that are regions with higher basal shear stress under glaciers or ice sheets.

## Description

Each folder named "Figure_n" (n=1,2,3...) contains the code to produce the corresponding figure in the manuscript. The computation relies on the open-source finite element library FEniCS (Logg and Wells, [2010](https://doi.org/10.1145/1731022.1731030)) and boundary element code developed by Davis, [2017](https://abdn.primo.exlibrisgroup.com/discovery/delivery/44ABE_INST:44ABE_VU1/12153058290005941?lang=en&viewerServiceCode=AlmaViewer). The mesh is generated by Gmsh (Geuzaine and Remacle, [2009](https://doi.org/10.1002/nme.2579)).

## Instructions to Run the Code

Each figure can be produced seperately by using the scripts and following the instructions attached in the corresponding folder.

### Dependencies

* [FEniCS](https://fenicsproject.org/) (Preferred in a Dokcer Container).
* [Gmsh](https://gmsh.info/)

## Authors

Hanwen Zhang
Timmothy Davis
Richard Katz
Laura Stevens
Dave May

## Acknowledgments

Thanks to
* [Brad Lipovsky](https://www.ess.washington.edu/people/profile.php?pid=lipovsky--brad)
* [Ching-Yao Lai](https://geosciences.princeton.edu/people/ching-yao-lai)

for their insightful suggestions for the LEFM model of the crevasses.

The basic idea that we explore arose during a conversation with
* [Poul Christoffersen](https://www.spri.cam.ac.uk/people/christoffersen/)
* [Robert Law](https://www.spri.cam.ac.uk/people/law/).

This research received funding from the European Research Council under Horizon 2020 research and innovation program grant agreement number 772255.
