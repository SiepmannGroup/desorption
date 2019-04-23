# Deep Neural Network Learning of Complex Binary Sorption Equilibria from Molecular Simulation Data

This repository contains the sample NN code used in the manuscript:

* Y. Sun, R. F. DeJaco, J. I. Siepmann, Deep Neural Network Learning of Complex Binary Sorption Equilibria from Molecular Simulation Data,  *Chem. Sci.* **2019**, 10, 4377–4388.  

The main 'executable' of this repository contain two IPython Notebook files, which train/load the neural networks used in the manuscript and validates the data presented in the figures and tables:
 * ```figures_prediction.ipynb```: validates Figure 2, Figure 3, and data in Table 1.
 * ```figures_transfer.ipynb```: validates Figure 4 and Table 2.
 
Other figures in the manuscript are not related to the neural network.

Dataset used in this manuscript can be found at <https://github.com/SiepmannGroup/MCCCS_DB/tree/master/diol-desorption>


**Please be advised that due to the randomness in neural network initialization and training, data obtained from the code in this repository may not be identical to those in the manuscript while general information will be consistent.**


The ```MCCCS-MN``` directory contains the Monte Carlo simulation program developed and used by the Siepmann research group. The source code, compilation and input files supplied correspond to the exact version and computer architecture on which the simulation results in the manuscript were obtained. Simulation details can be found in the Supplementry Information.
