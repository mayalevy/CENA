# CENA

CENA is a method for a joint identification of pairwise association together with the particular subset of cells in which the association is detected. In this implementation, CENA is limited to associations between the genes' expression levels (data from scRNA-sequencing) and an additional cellular meta-data of choice.


## Getting Started

CENA can be installed on Mac, Windows or Linux platforms by following instructions below

## Prerequisites

CENA is based on R but also requires python3 for running.
For Windows, we strongly recommend that python will be installed with conda (and not pure python) to ease the installation progress of the dependent libraries.
In addition to R, for Windows users, Rtools should be installed.

## Installing
Please make sure that the R and python that you install are compatible with your machine (32/64 bit).
### R installation
R can be downloaded and installed in the following link [download R](https://www.r-project.org/).
For windows users, Rtools should be installed [download R-tools](https://cran.r-project.org/bin/windows/Rtools).

### Python Installation
Python 3 should be downloaded and installed using miniconda [download miniconda](https://docs.conda.io/en/latest/miniconda.html). This is especially important for Windows users.
#### Packages of python
Dependent python libraries should be installed using the following commands:
```
conda install numpy
conda install -c vtraag python-igraph
conda install -c vtraag leidenalg
```
Please add python to the environment path variable, such that it will be recognized when you type python in the terminal.

### Environment validation
Before running CENA, please validate the installation process by checking the following sections:
#### Python is installed correctly
Please type the following commands in the terminal:
```
python # running python for making sure it works fine
(moved to the python console)
import igraph # making sure the pakcage igraph is installed correctly
import leidenalg # making sure the pakcage leidenalg is installed correctly
```
If python is not recognized, make sure you have added the python path to the evnironment path variable (PATH). If one of the modules is not recognize please repeat their installation and make sure to do so in the desired version of python.
#### R and python are working fine together
Please open R console and type the following R commands:
```
install.packages("reticulate") #reticulate library allows running python code
library(reticulate)
reticulate::import("igraph") #makes sure igraph is installed correctly
reticulate::import("leidenalg") # makes sure leidenalg is installed correctly
```
If one of the modules is not found, make sure the packages are installed on python, and that you are using the correct version of python. Reticulate uses the default python. In case you have many python versions on our computer, please install these packages on the default one (the one you get when you type which python), or change the reticulate to work with another version of python by typing.
```
py_config() # gives the direction of the current python
reticulate::use_python(**the python path**)# change the python path
```
Note: Reticulate has a bug. The function use_python does not work in Windows. Therefore, to use CENA on Windows please make sure to install all python packages on the default python version (you can change the default python by adding the python the at beginning of the environment PATH).

### CENA installation
For installing CENA from github, please open an R console, and type to following commands:
```
install.packages("devtools")
library("devtools")
install_github("mayalevy/CENA")
library("CENA")
```
## Running CENA

The required input data parameters for CENA are: 
A single cell gene expression matrix (genes X cells), a meta-data vector that we want to find association with, a 2-dim cell space (usually obtained using a dimension reduction), and the resolution parameter of the Leiden algorithm.
### The main function of the package:
```
CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter, no_cores = NULL, k1 = NULL, k2 = 10, Tw = 30, genesToRun = row.names(geneExpressionDataMatrix), python_path = NULL)
```
By changing the parameters of CENA function, the user may improve the results and the performance of the algorithm.
* **resolution_parameter** is the resolution parameter of the community detection algorithm Leiden which is responsible for more communities. Cutoff as desired; Higher values provide smaller clusters.
* **k1, k2** are responsible of the initial graph building, and hence for the connectivity of the graph and to the number of its edges. k1 should be kept relatively small. Default value of k1 is 1% of the cells. k2 was shown to have a relatively small effect on results. Therefore, to improve running times, the user may reduce the level of k2 for an initial run, and then, choose a higher k2 only for a genes with a good association with the meta-data. Increasing the k2 parameter should slightly improve the robustness of the results. Default value of k2 is 10.
* **Tw** is the subset size cutoff. Cutoff as desired; high values filter out small subsets. Default value is 30.
* **python_path** in case of not using the default python, you should specify the location of python by this parameter.

A full description of the parameters and return values can be found in the help page of the package.


#### Example:
Here is an example for the usage of CENA:
```
library("CENA")
data(cellSpace)
data(geneExpressionDataMatrix)
data(phenotypeData)
# running CENA on 5 genes
results = CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1)
```
Please notice that python path should be specified if not installed in the standard location
### Robustness Analysis:
After obtaining some interesting gene associations using CENA, it is advised to check their robustness. The robustness function runs CENA many times for a specific list of genes and returns their robustness scores across these runs. The robustness function applies the same parameters used in the original run (For more information, Levy et al., 2020).
#### Example:
```
library("CENA")
data(cellSpace)
data(geneExpressionDataMatrix)
data(phenotypeData)
results = CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1)
genesWithClusters = row.names(geneExpressionDataMatrix)[which(!is.na(rowSums(results$cluster_information)))]
robustnessResults = robustness(results, geneExpressionDataMatrix, phenotypeData, cellSpace, genesToRun = genesWithClusters, no_cores = 1)
```
### Specificity Analysis:
Another post analysis is applied by the specificity function. After obtaining some interesting gene associations using CENA, it is advised to check their specificity compared to permutations of the data. The specificity function applies the same parameters used in the original run (For more information, Levy et al., 2020).
#### Example:
```
data(cellSpace)
data(geneExpressionDataMatrix)
data(phenotypeData)
# running CENA on 5 genes
results = CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1)
genesWithClusters = row.names(geneExpressionDataMatrix)[which(!is.na(rowSums(results$cluster_information)))]
specificityResults = specificity(results, geneExpressionDataMatrix, phenotypeData, cellSpace, genesToRun = genesWithClusters, no_cores = NULL, numberOfRepeats =100)
```
## Authors

* **Maya Levy**
* **Amit Frishberg**
* **Irit Gat Viks**

## Citation
Maya Levy, Amit Frishberg, Irit Gat-Viks. ***Inferring cellular heterogeneity of associations from single cell genomics*** (submitted).