# CENA

CENA (CEllular Niche Association), is a method for a joint identification of pairwise association together with the particular subset of cells in which the association is detected. The algorithm relies on the input cell-state space to ensure a common cell state of all cells in the inferred cell subset. In this implementation, CENA tests association between multiple pairs of features: the expression of each gene against one additional meta-data.


## Getting Started

CENA can be installed on Mac, Windows or Linux platforms by the instructions below

### Prerequisites

CENA is based on R but also requires python3 for running.
For Windows, we strongly reccommend that python will be installed with conda (and not pure python) to ease the installation progress of the dependent libraries.
In addition to R, for Windows users, Rtools should be installed.



### Installing
Please make sure that the R and python that you install are compatible with your machine (32/64 bit).
#### R installation
R can be downlowded and installed in the following link [download R](https://www.r-project.org/)
For windows users, Rtools should be installed [download R-tools](https://cran.r-project.org/bin/windows/Rtools)

#### Python Installation
Python 3 should be downloaded and installed by miniconda (especially for Windows users) by the following link [download miniconda](https://docs.conda.io/en/latest/miniconda.html).
##### Packages of python
Dependent python libraries should be installed using the following commands:
```
conda install numpy
conda install -c vtraag python-igraph
conda install -c vtraag leidenalg
```

## Running CENA

The required input data for CENA should be the cell space (usually obtained using a dimension reduction), a gene expression matrix (genes X cells) and a phenotype vector that we want to find association with.
#### The main functin of the package:
```
CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter, no_cores = NULL, k1 = NULL, k2 = 10, minClusterVolume = 30, genesToRun = row.names(geneExpressionDataMatrix), python_path = NULL)
```
By changing the parameters of the CENA function, the user may improve the results and the performance of the algorithm.
** resolution_parameter** is the resolution parameter of the community detection algorithm Leiden which is responsible for more communities.
** k1, k2 ** are responsible of the initial graph building, and hence for the connectivity of the graph and to the number of its edges. When running many genes, for leverage the running time, the user may take a rather small k2 for initial run which will take a sample of the whole graph, and then, for the more interesting genes choose a higher k2, which may take more time but will be more reliable.
In addition, the user may specify the cluster size he wants to discover by the parameter **minClusterVolume**.
A full description of the parameters and return values can be found in the help page of the package.


#### Example:
Here is an example for the usage of CENA:
```
data(cellSpace)
data(geneExpressionDataMatrix)
data(phenotypeData)
# running CENA on 5 genes
results = CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1)
```
Please notice that python path should be specified if not installed in the standard location

## Authors

* **Maya Levy**
* **Amit Frishberg**
* **Irit Gat Viks**

## Citation
Maya Levy, Amit Frishberg, Irit Gat-Viks. ***Inferring cellular heterogeneity of associations from single cell genomics*** (submitted).
