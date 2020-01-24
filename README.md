# CENA

CENA (CEllular Niche Association), is a method for a joint identification of pairwise association together with the particular subset of cells in which the association is detected. The algorithm relies on the input cell-state space to ensure a common cell state of all cells in the inferred cell subset. In this implementation, CENA tests association between multiple pairs of features: the expression of each gene against one additional meta-data.


## Getting Started

CENA can be installed on Mac, Windows or Linux platforms by the instructions below

## Prerequisites

CENA is based on R but also requires python3 for running.
For Windows, we strongly reccommend that python will be installed with conda (and not pure python) to ease the installation progress of the dependent libraries.
In addition to R, for Windows users, Rtools should bIn addition to R, for Windows users, Rtools should be installed.
e installed.



## Installing
Please make sure that the R and python that you install are compatible with your machine (32/64 bit).
### R installation
R can be downlowded and installed in the following link [download R](https://www.r-project.org/)
For windows users, Rtools should be installed [download R-tools](https://cran.r-project.org/bin/windows/Rtools)

### Python Installation
Python 3 should be downloaded and installed by miniconda (especially for Windows users) by the following link [download miniconda](https://docs.conda.io/en/latest/miniconda.html).
#### Packages of python
Dependent python libraries should be installed using the following commands:
```
conda install numpy
conda install -c vtraag python-igraph
conda install -c vtraag leidenalg
```
Please add the python command to the PATH, such that it will be recognized when you type python in the terminal.

### Environment validation
Before moving to running CENA, please validate by the following commands that everything is fine.
#### Python is installed correcly
Please type the following commands in the terminal:
```
python # running python for making sure it works fine
(moved to the python console)
import igraph # making sure the pakcage igraph is installed correctly
import leidenalg # making sure the pakcage leidenalg is installed correctly
```
If python is not recognized, make sure you have added the python path to the evnironment path variable (PATH). If one of the modules is not recognize please make sure that their installation succeed and was done on the current python.
#### R and python are working fine together
Please open R console and type the following R commands:
```
install.packages("reticulate") #reticulate library allows running python code
library(reticulate)
reticulate::import("igraph") #makes sure igraph is installed correctly
reticulate::import("leidenalg") # makes sure leidenalg is installed correctly
```
In case it returns that a module is not found, make sure the packages are installed on python, and that you are using the correct python.
Reticulate uses the default python, in case you have many python versions on our computer, please install these packages on the default one (the one you get when you type which python), or change the reticulate to work with another python by typing 
```
py_config() # gives the direction of the current python
reticulate::use_python(**the python path**)# change the python path
```
Note: Reticulate has a bug. use_python does not work in Windows, so Windows uses may make sure the installation of the python packages is done on the default python (you can change the default python by adding the python the the beginnig of the environment PATH)

### CENA installation
For installing CENA from github, please open an R console, and type to following commands:
```
install.packages("devtools")
library("devtools")
install_github("mayalevy/CENA")
library("CENA")
```
## Running CENA

The required input data for CENA should be the cell space (usually obtained using a dimension reduction), a gene expression matrix (genes X cells) and a phenotype vector that we want to find association with.
### The main function of the package:
```
CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter, no_cores = NULL, k1 = NULL, k2 = 10, Tw = 30, genesToRun = row.names(geneExpressionDataMatrix), python_path = NULL)
```
By changing the parameters of CENA function, the user may improve the results and the performance of the algorithm.
* **resolution_parameter** is the resolution parameter of the community detection algorithm Leiden which is responsible for more communities. Cutoff as desired; Higher values provide smaller clusters.
* **k1, k2** are responsible of the initial graph building, and hence for the connectivity of the graph and to the number of its edges. k1 should be kept relatively small. Default value of k1 is 1% of the cells. When running many genes, for leverage the running time, the user may take a rather small k2 for initial run which will take a sample of the whole graph, and then, for the more interesting genes choose a higher k2, which may take more time but will be more reliable. k2 have relatively small effect on results. Default value of k2 is 10.
* **Tw** is the cluster size that the user wants to discover. Cutoff as desired; high values filter out small subsets. Default value is 30
* **python_path** - in case of not using the default python, you should specify the location of python by this parameter.

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
After running CENA on many genes and choosing a few of them, one can run a robustness analysis for specific genes for checking the robustness of the gene results.
*robustness* function runs the analysis multiple times for a specific gene and returns the percentage of runs in which the cluster was found again in adittion to a score in the rang of 0-1 which represents the cluster score (low score means good score).
#### Example:
```
library("CENA")
data(cellSpace)
data(geneExpressionDataMatrix)
data(phenotypeData)
results = CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1)
robustnessResults = robustness(results, geneExpressionDataMatrix,phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1, genesToRun = row.names(geneExpressionDataMatrix)[4:5])
```
## Authors

* **Maya Levy**
* **Amit Frishberg**
* **Irit Gat Viks**

## Citation
Maya Levy, Amit Frishberg, Irit Gat-Viks. ***Inferring cellular heterogeneity of associations from single cell genomics*** (submitted).