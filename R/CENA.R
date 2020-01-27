### code for package dependencies library(SDMTools) library(fields) library(igraph) library(parallel) library(foreach)
### library(doSNOW) library(reticulate)

# Sparsifies the matrix to have edges which are not too long (represented as distMatrix)
#' @keywords internal
sparsifingDistanceMatrix = function(distMatrix, cellSpace, alpha, beta, k1) {

  leafs = getBins(cellSpace, alpha, beta, k1)

  nonZeroEdges = do.call(rbind, lapply(leafs, function(leaf) {
    binNodes = leaf$points
    neighbors = getNeighborBins(leaf, leafs)
    binNames = row.names(binNodes)
    neigNames = c(row.names(binNodes), unlist(lapply(neighbors, function(currNeig) {
      row.names(currNeig$points)
    })))
    expand.grid(match(binNames, row.names(cellSpace)), match(neigNames, row.names(cellSpace)))
  }))
  nonDuplicatedEdges = as.matrix(unique(data.table::as.data.table(nonZeroEdges)))
  nonDuplicatedEdgesFinal = nonDuplicatedEdges[which((nonDuplicatedEdges[, 2] - nonDuplicatedEdges[, 1]) != 0), ]
  return(nonDuplicatedEdgesFinal)
}

# Keeps the order of the edges the same, no matter of the order - smaller the first
#' @keywords internal
get_edge_name = function(i, js, nodeNames, relevantEdges) {
  js = match(relevantEdges[js], nodeNames)
  smallers = js[which(js < i)]
  largers = js[which(js > i)]
  result = c()
  if (length(smallers) > 0) {
    result = c(result, paste(nodeNames[smallers], nodeNames[i], sep = "^"))
  }
  if (length(largers) > 0) {
    result = c(result, paste(nodeNames[i], nodeNames[largers], sep = "^"))
  }
  return(result)
}

#' @keywords internal
expand.grid.unique = function(x, y, include.equals = FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i - include.equals)])
    if (length(z))
      cbind(x[i], z, deparse.level = 0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

# Builds the edges network from the cells network
#' @keywords internal
getEdgesMatrix = function(ig, adjacencyNew, adjEdgeListRefined) {
  edgeMatrixForG = cbind(row.names(adjacencyNew)[adjEdgeListRefined[, 1]], row.names(adjacencyNew)[adjEdgeListRefined[, 2]])
  G = ig$Graph()
  vertex.attrs = list(name = row.names(adjacencyNew))
  G$add_vertices(reticulate::r_to_py(vertex.attrs$name))
  G$add_edges(reticulate::r_to_py(edgeMatrixForG))

  neig1List = G$neighborhood(G$vs$indices, as.integer(1))
  neig1Matrix = bigmemory::as.big.matrix(do.call(rbind, lapply(1:length(neig1List), function(i) {
    currList = neig1List[[i]] + 1
    currList = currList[currList != i]
    if (length(currList) > 1) {
      allCombinations = expand.grid.unique(currList, currList)
      combinationWithI = cbind(allCombinations, rep(i, dim(allCombinations)[1]), rep(i, dim(allCombinations)[1]))
      combinationWithI[which(combinationWithI[, 1] < combinationWithI[, 3]), c(1, 3)] = combinationWithI[which(combinationWithI[,
                                                                                                                                1] < combinationWithI[, 3]), c(3, 1)]
      combinationWithI[which(combinationWithI[, 2] < combinationWithI[, 4]), c(2, 4)] = combinationWithI[which(combinationWithI[,
                                                                                                                                2] < combinationWithI[, 4]), c(4, 2)]
      combinationWithI
    }
  })))
  # neig1Matrix = as.matrix(unique(data.table::as.data.table(neig1Matrix)))
  edge1Matrix = cbind(adjacencyNew[neig1Matrix[, c(1, 3)]], adjacencyNew[neig1Matrix[, c(2, 4)]])
  edge1Matrix = edge1Matrix/rowSums(edge1Matrix)
  # edgesScores = (1 / abs(1 - edge1Matrix[,1] / edge1Matrix[,2])) * (1 / abs(1 - edge1Matrix[,2] / edge1Matrix[,1]))
  edgesScores = 1/abs(edge1Matrix[, 2] - edge1Matrix[, 1])
  # edgesScores = exp(-abs(edge1Matrix[,2]-edge1Matrix[,1]))
  edgesScores[which(is.infinite(edgesScores))] = max(edgesScores[which(is.finite(edgesScores))])
  edgesMatrix = cbind(paste(row.names(adjacencyNew)[neig1Matrix[, 3]], row.names(adjacencyNew)[neig1Matrix[, 1]], sep = "^"),
                      paste(row.names(adjacencyNew)[neig1Matrix[, 4]], row.names(adjacencyNew)[neig1Matrix[, 2]], sep = "^"))

  # nodeNames = rownames(adjacency) edgeMatrix = do.call(rbind,sapply(seq(length(nodeNames)), function(currNodeName_ind) {
  # currNodeName = nodeNames[currNodeName_ind] relevantEdges = nodeNames[which(adjacency[currNodeName,] != 0)]
  # nodeValueVector = adjacency[currNodeName, which(adjacency[currNodeName,] != 0)] if(length(nodeValueVector)>1){
  # do.call(rbind,lapply(1:(length(nodeValueVector)-1),function(i){ neigValues =
  # nodeValueVector[(i+1):length(nodeValueVector)] scores = (1 / abs(1 - nodeValueVector[i] / neigValues)) * (1 / abs(1 -
  # neigValues / nodeValueVector[i])) edgeFirst = get_edge_name(i = currNodeName_ind, js = i,nodeNames,
  # relevantEdges)#paste(nonZeroEdges[1,relevantEdges[i]], nonZeroEdges[2,relevantEdges[i]], sep = '^') edgesSecond =
  # get_edge_name(i = currNodeName_ind, js = c((i+1):length(nodeValueVector)), nodeNames,
  # relevantEdges)#paste(nonZeroEdges[1,relevantEdges[(i+1):length(nodeValueVector)]],
  # nonZeroEdges[2,relevantEdges[(i+1):length(nodeValueVector)]], sep = '^') cbind(edgeFirst,edgesSecond,scores) })) } }))
  # scoreVector = as.numeric(edgeMatrix[,3]) scoreVector[which(is.infinite(scoreVector))] =
  # max(scoreVector[which(is.finite(scoreVector))]) return(list(edgeMatrix[,1:2],scoreVector))
  return(list(edgesMatrix = edgesMatrix, edgesScores = edgesScores))
}

# moving back to node network for deducing the cells which consist the subset
#' @keywords internal
getNodesForEachCluster = function(clusteringEdgesLabels, nonZeroEdges) {

  tab = table(clusteringEdgesLabels)
  # clusters = as.numeric(names(tab)[tab >= 30])
  clusters = as.numeric(names(tab))
  result = sapply(clusters[clusters != 0], function(cluster) {
    currLabels = which(clusteringEdgesLabels == cluster)
    currEdgesPairs = nonZeroEdges[currLabels]
    currNodeNames = unique(unlist(strsplit(currEdgesPairs, "\\^")))
    currNodeNames
  }, simplify = FALSE)

  # result = result[lapply(result, length) > 15]
  return(result)
}

#' @keywords internal
runJohnson = function(data) {
  selectedIndexes = which(data != "-")
  newData = rep("-", length(data))
  johnsonResults = try(Johnson::RE.Johnson(as.numeric(as.matrix(data[selectedIndexes])))$transformed, T)
  if (length(johnsonResults) == 1) {
    data
  } else {
    newData[selectedIndexes] = johnsonResults
    newData
  }
}
#' @keywords internal
normalizeData = function(data) {
  dataNew = apply(data, 2, runJohnson)
  row.names(dataNew) = row.names(data)
  apply(dataNew, 2, function(x) {
    selectedIndexes = which(x != "-")
    newVecor = rep("-", length(x))
    xNumbers = as.numeric(as.matrix(x[selectedIndexes]))
    newVecor[selectedIndexes] = (xNumbers - mean(xNumbers))/stats::sd(xNumbers)
    newVecor
  })
}

#' @keywords internal
computeNormalStats <- function(realDataScore, permDataScores, permsL) {
  Z = (realDataScore - mean(permDataScores))/(stats::sd(permDataScores))
  p = (-1) * stats::pnorm(realDataScore, mean = mean(permDataScores), sd = sd(permDataScores), lower.tail = FALSE, log.p = T)
  list(z = Z, p = p, ks = stats::ks.test(permDataScores, rnorm(permsL, mean = mean(permDataScores), sd = sd(permDataScores)))$p.value)
}

#' @keywords internal
getJohnsonStatsAlt = function(real, perms) {
  allScores = as.numeric(normalizeData(as.matrix(c(real, perms))))
  stats = computeNormalStats(allScores[1], allScores[2:length(allScores)], length(perms))
  p = stats$p
  ks = stats$ks
  z = stats$z
  list(z = z, p = p, ks = ks)
}

#' @keywords internal
calc_p_value_pnorm = function(perms, true_val) {
  result = tryCatch({
    getJohnsonStatsAlt(true_val, perms)
  }, warning = function(w) {

  }, error = function(e) {
    return(list())
  }, finally = {

  })
  return(result)
}

# calc the expliained percentage by he line.

# make phenotype permutations and check the p-value of the real subset vs. the permutated subset
#' @keywords internal
calc_r_and_p_value = function(currRegressionAnalysis, geneScoresPerGroup, phenScoresPerGroup, phenScoresPerGroupName) {

  if (is.null(currRegressionAnalysis) || is.null(stats::coefficients(currRegressionAnalysis)) || is.na(stats::coefficients(currRegressionAnalysis)) ||
      !(phenScoresPerGroupName %in% rownames(stats::coefficients(currRegressionAnalysis)))) {
    return(list(p_value_p_t = NULL, ks_p_t = NULL, perms_p_t = c(), r_value = NULL))
  }
  true_p_t = -log10(stats::coefficients(currRegressionAnalysis)[phenScoresPerGroupName, "Pr(>|t|)"])
  true_r_val = currRegressionAnalysis$r.squared

  # spearmanCor = cor.test(x=phenScoresPerGroup, y=geneScoresPerGroup, method = 'spearman')

  num_of_perms = 100  #TODO- move to another place

  # df = as.data.frame(cbind(phenScoresPerGroup,geneScoresPerGroup)) colnames(df) = c('phen','gene') true_MSE =
  # sqrt(mean((stats::lm(phen~gene,data = df)$residuals)^2)) rand_MSE = sapply(1:num_of_perms,function(i) { cellsForReg =
  # sample(1:length(phenScoresPerGroup),round(length(phenScoresPerGroup)*0.7)) cellsForValidate =
  # which(!(1:length(phenScoresPerGroup) %in% cellsForReg)) rand_regressionAnalysis = stats::lm(phen~gene,data =
  # df[cellsForReg,]) sqrt(mean((predict(rand_regressionAnalysis, df[cellsForValidate,]) - df$phen[cellsForValidate])^2)) })

  # phenotype permutations
  perms_p_t = sapply(1:num_of_perms, function(i) {
    rand_phenScoresPerGroup = phenScoresPerGroup[sample(1:length(phenScoresPerGroup))]
    rand_regressionAnalysis = stats::.lm.fit(x = cbind(1, as.matrix(rand_phenScoresPerGroup)), y = geneScoresPerGroup)

    rss <- sum(rand_regressionAnalysis$residuals^2)
    rdf <- length(geneScoresPerGroup) - 1
    resvar <- rss/rdf
    R <- chol2inv(rand_regressionAnalysis$qr)
    # R2 = 1 - sd(rand_regressionAnalysis$residuals)^2/sd(geneScoresPerGroup)^2
    se <- sqrt(diag(R) * resvar)

    # rand_intercept = rand_regressionAnalysis$coefficients[1] rand_slope = rand_regressionAnalysis$coefficients[2]

    rand_pvalue = -log10(2 * pt(abs(rand_regressionAnalysis$coef/se), rdf, lower.tail = FALSE)[2])

    # rand_regressionAnalysis = summary(stats::lm(geneScoresPerGroup~rand_phenScoresPerGroup)) rand_slope =
    # coefficients(rand_regressionAnalysis)['rand_phenScoresPerGroup','Estimate'] rand_intercept =
    # coefficients(rand_regressionAnalysis)['(Intercept)', 'Estimate']
    # -log10(coefficients(rand_regressionAnalysis)['rand_phenScoresPerGroup','Pr(>|t|)'])

    rand_pvalue
  })

  # try to use pnorm in order to calc p-value, and not pJohnson, because pJosnson returns zeros for extream real data In
  # addition, should change the groups file, in oder to sort by p-value which goes up + filter by p-value > 2 In group file :
  # order the results rev - if using non log, should rev the results. + change the p-value to be p-value >1.3 or p-value
  # <0.001

  stats_p_t = calc_p_value_pnorm(perms_p_t, true_p_t)
  # stats_p_t = pnorm(true_MSE,mean = mean(rand_MSE), sd = sd(rand_MSE), lower.tail = T, log.p = T)

  # return (list(p_value_p_t = stats_p_t$p, ks_p_t = stats_p_t$ks ,perms_p_t = perms_p_t, r_value = true_r_val))
  return(list(p_value_p_t = stats_p_t$p, r_value = true_r_val))
  # return (list(p_value_p_t = stats_p_t, r_value = true_r_val))
}

# Calculate all the subset characteries
#' @keywords internal
analyzeClustersFunction = function(adjacency, nodesPerClusters, geneVectorScores, phenVectorScores, cellSpace) {

  regressionAnalysis = lapply(nodesPerClusters, function(group) {


    chosenGroupIndexes = which(rownames(cellSpace) %in% group)
    geneGroup = cellSpace[chosenGroupIndexes, ]
    geneGroup_size = length(chosenGroupIndexes)
    ch_points_cords_geneGroup = grDevices::chull(geneGroup)
    non_ch_points = geneGroup[-ch_points_cords_geneGroup, ]
    ch_points = geneGroup[ch_points_cords_geneGroup, ]
    non_ch_points_cords = which(rownames(cellSpace) %in% rownames(non_ch_points))
    # ch_points_cords = which(rownames(xyLocations)%in% rownames(ch_points))
    ch_points_cords = rep(0, length(rownames(ch_points)))
    for (cord_ind in 1:length(rownames(ch_points))) {
      ch_points_cords[cord_ind] = which(rownames(cellSpace) == rownames(ch_points)[cord_ind])
    }
    ch_points_size = length(ch_points_cords)
    non_ch_points_size = length(non_ch_points_cords)

    out = SDMTools::pnt.in.poly(cellSpace[-c(ch_points_cords, non_ch_points_cords), ], cellSpace[ch_points_cords, ])

    invadors_points = out[which(out$pip == 1), 1:2]
    invadors_size = dim(invadors_points)[1]
    # percentage inside the convexhull
    percentageOfInvadors = invadors_size/(invadors_size + non_ch_points_size + ch_points_size)
    chosenGroupIndexes_with_invadors = which(rownames(adjacency) %in% c(group, c(rownames(invadors_points))))

    # non delta regression
    phenNonDelta_with_invadors = as.numeric(phenVectorScores[chosenGroupIndexes_with_invadors])
    geneNonDelta_with_invadors = as.numeric(geneVectorScores[chosenGroupIndexes_with_invadors])
    regressionAnalysisNonDelta = tryCatch({
      reg = stats::lm(geneNonDelta_with_invadors ~ phenNonDelta_with_invadors)
      list(r.squared = summary(reg)$r.squared, slope_nonDelta = as.numeric(reg$coefficients[2]), intercept_nonDelta = as.numeric(reg$coefficients[1]),
           p_t_general_nonDelta = coefficients(summary(reg))["phenNonDelta_with_invadors", "Pr(>|t|)"])
    }, error = function(e) {
      NULL
    })

    list(r_squared_nonDelta = regressionAnalysisNonDelta$r.squared, slope_nonDelta = regressionAnalysisNonDelta$slope_nonDelta,
         intercept_nonDelta = regressionAnalysisNonDelta$intercept_nonDelta, p_t_general_nonDelta = regressionAnalysisNonDelta$p_t_general_nonDelta,
         cluster = chosenGroupIndexes_with_invadors, invadorP = invadors_size/(invadors_size + geneGroup_size), p_value_p_t_general = regressionAnalysisNonDelta$r.squared)
  })

  return(regressionAnalysis)
}
## calc the numOfGenes genes which are significantlly differentially expressed inside and outside the cells subset using the
## t-test
#' @keywords internal
getClusterGenes = function(geneExpressionDataMatrix, cellClustersIndexes, numOfGenes) {
  # calc the t-test value for each gene
  clusterGenesResults = apply(geneExpressionDataMatrix, 1, function(geneRow) {
    insideCluster = geneRow[cellClustersIndexes]
    outsideCluster = geneRow[-cellClustersIndexes]
    stats::t.test(insideCluster, outsideCluster)$p.value
  })
  sortedclusterGenesResults = sort(clusterGenesResults)
  sortedclusterGenesResults[1:numOfGenes]
  # p.adjust(sortedclusterGenesResults,method = 'fdr')
}

# #DEMO:: geneName = 'GZMB' cellSpace =
# readRDS('/Users/mayalevy/Downloads/MarkerSelection/Runs/T-cell_1529229592/finalResults/temps/xyLocations')
# geneExpressionDataMatrix =
# readRDS('/Users/mayalevy/Downloads/MarkerSelection/Runs/T-cell_1529229592/finalResults/temps/geneDataForRegression')
# phenotypeData =
# readRDS('/Users/mayalevy/Downloads/MarkerSelection/Runs/T-cell_1529229592/finalResults/temps/phenVectorScoresOriginal')
# names(phenotypeData) = colnames(geneExpressionDataMatrix) resultsFolderPath =
# '/Users/mayalevy/Downloads/MarkerSelection/Runs/T-cell_1529229592/finalResults' discrete_phenotype = TRUE
# ratioToTiniesBin = 1 Tw = 15 geneValues = geneExpressionDataMatrix[geneName,] k1 = 5 k2 = 15
# CENA_single_gene(geneValues, phenotypeData, cellSpace, resultsFolderPath, discrete_phenotype, k1, k2, ratioToTiniesBin,
# Tw) # public

# CENA run for a single gene
#' @keywords internal
CENA_single_gene = function(geneValues, phenotypeData, cellSpace, k1, k2, Tw, nonZeroMatrix, resolution_parameter,
                            python_path) {
  # normalize geneExpression vector
  geneVectorScores = (geneValues - mean(geneValues))/sd(geneValues)

  # normalize phenotype vecotr - in case the phenotype is discretic, no need for normalization
  phenVectorScores = (phenotypeData - mean(phenotypeData))/sd(phenotypeData)

  # [cells X cells] - distance phenotype matrix of the cells

  # adjacency = bigmemory::big.matrix(nrow = dim(cellSpace)[1],ncol = dim(cellSpace)[1],init = 0)
  adjacency = Matrix::Matrix(0, nrow = dim(cellSpace)[1], ncol = dim(cellSpace)[1], sparse = TRUE)
  adjacency[nonZeroMatrix] = fields::rdist(geneVectorScores)[nonZeroMatrix]
  phenDistVector = fields::rdist(phenVectorScores)[nonZeroMatrix]
  phenDistVector[phenDistVector == 0] = 0.01
  adjacency[nonZeroMatrix] = adjacency[nonZeroMatrix]/phenDistVector
  # adjacency[is.na(adjacency)] = 0
  zeroGeneIndexes = which(geneValues == 0)
  if (length(zeroGeneIndexes) != 0) {
    adjacency[zeroGeneIndexes, zeroGeneIndexes] = 0
  }

  colnames(adjacency) = row.names(adjacency) = row.names(cellSpace)
  # colnames(adjacency) = row.names(adjacency) = 1:dim(cellSpace)[1]

  # avoid points in which both of the cells have geneExpression of 0.
  adjEdgeList = Matrix::which(adjacency > 0, arr.ind = T)
  if(!is.null(python_path)) {
    reticulate::use_python(python_path)
  }
  # initialGraph = igraph::graph(adjEdgeList, directed = FALSE) edgeAmmount = length(igraph::E(initialGraph)) bfsEdgeList =
  # do.call(rbind,lapply(sample(1:dim(adjacency)[2],k2),function(currRoot){
  # cbind(1:dim(adjacency)[2],igraph::dfs(initialGraph,currRoot,dist=TRUE)$dist)
  # #igraph::as_edgelist(igraph::mst(initialGraph, weights = sample(edgeAmmount,edgeAmmount))) })) bfsEdgeList =
  # bfsEdgeList[which(bfsEdgeList[,2]!=0),] bfsEdgeList = unique(data.table::as.data.table(bfsEdgeList))

  degrees = Matrix::colSums(sign(adjacency))
  probForNode = 1/degrees
  probForSampling = probForNode[adjEdgeList[, 1]] * probForNode[adjEdgeList[, 2]]
  probForSampling = probForSampling/sum(probForSampling)

  np <- reticulate::import("numpy", delay_load = TRUE)
  selectedEdges = np$random$choice(as.integer(dim(adjEdgeList)[1]), as.integer(k2 * length(degrees)), "False", probForSampling)
  adjEdgeListRefined = adjEdgeList[selectedEdges, ]
  if (dim(adjEdgeListRefined)[1] != 0) {
    adjEdgeListRefined = rbind(adjEdgeListRefined, adjEdgeListRefined[, 2:1])
    adjEdgeListRefined = as.matrix(unique(data.table::as.data.table(adjEdgeListRefined)))
  }

  ### support semi-linearity: we do no allow more than k1*k2 edges remove edges if there are more than k1*k2


  # for(gene_row_ind in seq(nrow(adjacency))) { gene_row = adjacency[gene_row_ind,] if(length(which(gene_row != 0)) > k2) {
  # to_remove = sample(which(gene_row != 0),length(which(gene_row != 0)) - k2) adjacency[gene_row_ind,to_remove] = 0
  # adjacency[to_remove,gene_row_ind] = 0 } } adjEdgeListRefined = Matrix::which(adjacency!=0,arr.ind = T)


  # move to adjancy matrix
  if (dim(adjEdgeListRefined)[1] != 0) {
    ig <- reticulate::import("igraph", delay_load = TRUE)
    adjacencyNew = Matrix::Matrix(0, nrow = dim(adjacency)[1], ncol = dim(adjacency)[1], sparse = TRUE)
    colnames(adjacencyNew) = row.names(adjacencyNew) = row.names(adjacency)
    adjacencyNew[adjEdgeListRefined] = adjacency[adjEdgeListRefined]
    ##### reticulate::source_python('/Users/mayalevy/Dropbox\ (Irit\ Gat\
    ##### Viks)/SC_ReferenceData/markerselection/AmitSimulations/python_try.py') aaa =
    ##### getEdgesMatrix_python(as.matrix(adjacencyNew),reticulate::r_to_py(as.matrix(adjEdgeListRefined)-1),
    ##### reticulate::r_to_py(row.names(adjacencyNew))) #aaa

    edgesMatrixResult = getEdgesMatrix(ig, adjacencyNew, adjEdgeListRefined)
    # edgesMatrixResult = getEdgesMatrix(adjacencyNew)
    gc(FALSE)

    leidenalg <- reticulate::import("leidenalg", delay_load = TRUE)

    G2 = ig$Graph()
    vertex.attrs = list(name = unique(c(edgesMatrixResult$edgesMatrix)))
    G2$add_vertices(reticulate::r_to_py(vertex.attrs$name))
    G2$add_edges(reticulate::r_to_py(edgesMatrixResult$edgesMatrix))

    clusteringEdges = leidenalg$find_partition(G2, leidenalg$RBConfigurationVertexPartition, weights = abs(edgesMatrixResult$edgesScores),
                                               resolution_parameter = resolution_parameter)

    # library(leiden) clusteringEdgesLabels = leiden::leiden(as_adjacency_matrix(G,attr = 'weight',sparse = F), partition_type
    # = 'ModularityVertexPartition', resolution_parameter = resolution_parameter, n_iterations = 4)

    optimiser = leidenalg$Optimiser()
    optimiser$merge_nodes(clusteringEdges)

    # clusteringEdges = igraph::cluster_louvain(G)
    clusteringEdgesLabels = clusteringEdges$membership

    # get the cluster nodes that had been found nodesPerClusters = getNodesForEachCluster(clusteringEdgesLabels, nonZeroEdges)
    nodesPerClusters = getNodesForEachCluster(clusteringEdgesLabels, vertex.attrs$name)
  } else {
    nodesPerClusters = list()
    clusteringEdgesLabels = c()
  }

  nodesPerClusters = nodesPerClusters[lapply(nodesPerClusters, length) > Tw]  ############################ CHANGE ACCORDING TO CLUSTERING RESULT
  # filter out groups with less that 3 phenotypes
  nodesPerClusters = lapply(nodesPerClusters, function(phens) {
    chosenGroupIndexes = which(rownames(adjacencyNew) %in% phens)
    phenNonDelta = as.numeric(phenVectorScores[chosenGroupIndexes])
    if (length(table(phenNonDelta)) > 2) {
      return(phens)
    } else {
      return(NULL)
    }
  })
  nodesPerClusters = Filter(Negate(is.null), nodesPerClusters)

  clustersAnalysis = analyzeClustersFunction(adjacencyNew, nodesPerClusters, geneVectorScores, phenVectorScores, cellSpace)
  # clustersAnalysis = clustersAnalysis[which(unlist(lapply(clustersAnalysis, function(x){x$invadorP<0.3})))]
  # ############################ CHANGE ACCORDING TO CLUSTERING RESULT

  relevantCommunetee = NA
  if (length(clustersAnalysis) != 0) {
    # take the best p_t_grneral community
    p_t_general = unlist(lapply(clustersAnalysis, function(x) {
      ifelse(!is.null(x$p_value_p_t_general), x$p_value_p_t_general, 0)
    }))

    # rRelevent = which(unlist(lapply(clustersAnalysis,function(x){x$r_squared_nonDelta}))>0.3) relevantCommunetee =
    # switch(all(is.na(p_t_general[rRelevent]))+1,
    # rRelevent[which.max(p_t_general[rRelevent])],NA)#ifelse(all(is.na(p_t_general)), NA, which.max(p_t_general) )#which.min

    relevantCommunetee = switch(all(is.na(p_t_general)) + 1, which.max(p_t_general), NA)  #ifelse(all(is.na(p_t_general)), NA, which.max(p_t_general) )#which.min
  }

  # relevantCellIndexes = which(colnames(geneExpressionDataMatrix) %in% nodesPerClusters[[relevantCommunetee]]) clusterGenes
  # = getClusterGenes(geneExpressionDataMatrix, relevantCellIndexes ,10)

  # geneName, slope, r2, list, geneScoresPerGroup, phenScoresPerGroup
  if (is.na(relevantCommunetee)) {
    result = NULL
  } else {
    result = clustersAnalysis[[relevantCommunetee]]
  }
  # result_p_t_general = c(result, list(geneName = geneName, cluster = result$cluster)) #cluster = relevantCellIndexes))
  # #clusterGenes = clusterGenes))

  result_gene = NULL
  if (!is.null(result) & length(result) != 0) {
    result_gene = result
  }
  return(list(res = result_gene, nodesPerClusters = nodesPerClusters))
}


# the run of the algorithm - run parrallely

# the run of CENA algorithm - run parrallely
#' @keywords internal
CENA_Main = function(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter, genesToRun = row.names(geneExpressionDataMatrix),
                     no_cores = NULL, k1 = NULL, k2 = 10, Tw = 30, python_path = NULL, clusterDetails = F) {

  if (is.null(k1)) {
    k1 = ceiling(0.01 * dim(cellSpace)[1])
  }
  colnames(geneExpressionDataMatrix) = rownames(cellSpace) = names(phenotypeData) = paste("cell", seq(ncol(geneExpressionDataMatrix)),
                                                                                          sep = "_")
  # if(!is.null(Tw) && Tw > 0 && Tw < 1) { Tw =
  # round(Tw * nrow(cellSpace)) } else { Tw = 0 }

  # sparsify the distMatrix by bins
  alpha = 1.2
  beta = 0.6
  nonZeroMatrix = sparsifingDistanceMatrix(bigmemory::as.big.matrix(fields::rdist(cellSpace, cellSpace)), cellSpace, alpha,
                                           beta, k1)
  gc(FALSE)

  ##### Main algorithm runs #####
  if (is.null(no_cores)) {
    no_cores = max(1, parallel::detectCores() - 1)
  }
  cl <- parallel::makeCluster(no_cores)

  parallel::clusterExport(cl = cl, varlist = c("cellSpace", "k1", "k2", "Tw", "CENA_single_gene", "resolution_parameter",
                                               "python_path", "expand.grid.unique", "getEdgesMatrix", "getNodesForEachCluster", "runJohnson", "normalizeData", "computeNormalStats",
                                               "getJohnsonStatsAlt", "get_edge_name", "calc_p_value_pnorm", "calc_r_and_p_value", "analyzeClustersFunction", "getClusterGenes",
                                               "nonZeroMatrix"), envir = environment())
  doSNOW::registerDoSNOW(cl)
  # progress bar

  pb <- utils::txtProgressBar(min = 1, max = length(genesToRun), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # t = proc.time()
  `%dopar2%` <- foreach::`%dopar%`
  `%do2%` <- foreach::`%do%`

  iterations = iterators::iter(geneExpressionDataMatrix[genesToRun, ], by = "row")
  # geneNameNum = NULL geneResults <- foreach::foreach(geneStat = iterations, .options.snow = opts) %do2% {
  geneResultsBasic <- foreach::foreach(geneStat = iterations, .options.snow = opts) %dopar2% {
    # geneResults <- for(i in 1:dim(geneExpressionDataMatrix)[1]) { geneStat = geneExpressionDataMatrix[i,] geneName =
    # genesToRun[i] print(geneName)
    gc(FALSE)
    CENA_single_gene(as.numeric(geneStat), phenotypeData, cellSpace, k1, k2, Tw, nonZeroMatrix, resolution_parameter,
                     python_path)
    # setTxtProgressBar(pb, geneNameNum) result_gene

  }
  # print(proc.time() - t)
  parallel::stopCluster(cl)
  close(pb)

  geneResults = lapply(geneResultsBasic, function(x) {
    x$res
  })
  geneAllNodes = lapply(geneResultsBasic, function(x) {
    x$nodesPerClusters
  })

  # arrange output
  geneNames = genesToRun
  # cluster information
  cluster_information = lapply(geneResults, function(x) {
    if (is.na(x) || is.null(x)) {
      return(c(NA, NA, NA, NA))
    } else {
      return(c(x$slope_nonDelta, x$intercept_nonDelta, x$r_squared_nonDelta, x$p_t_general_nonDelta))
    }
  })
  cluster_information = do.call(rbind, cluster_information)
  rownames(cluster_information) = geneNames
  colnames(cluster_information) = c("slope", "intercept", "r_squared", "p_value")
  # cluster
  cluster = lapply(geneResults, function(x) {
    if (is.na(x) || is.null(x)) {
      return(c(NA))
    } else {
      return(c(x$cluster))
    }

  })
  cluster[sapply(cluster, is.null)] <- NULL
  names(cluster) = geneNames
  # # cluster statistics
  # cluster_statistics = lapply(geneResults, function(x) {
  #     if (is.na(x) || is.null(x)) {
  #         return(c(NA))
  #     } else {
  #         return(c(x$p_value_p_t_general))
  #     }
  #
  # })
  # cluster_statistics[sapply(cluster_statistics, is.null)] <- NULL
  # names(cluster_statistics) = geneNames

  pars = list()
  pars[["resolution_parameter"]] = resolution_parameter
  pars[["k1"]] = k1
  pars[["k2"]] = k2
  pars[["Tw"]] = Tw
  pars[["python_path"]] = python_path

  if(clusterDetails){
    return(geneAllNodes)
  }
  final_geneResults = list(cluster_information = cluster_information, cluster = cluster, parameters = pars)
  return(final_geneResults)
}

#' Robustness check for CENA predictions
#'
#' Calculate the robustness scores for selected genes.
#'
#' @param prediction The result of the original CENA run. The robustness scores are calculated based on the parameters used in this run.
#' @param geneExpressionDataMatrix A matrix containing the single-cell RNA-seq data.
#' Each row corresponds to a certain gene and each column to a certain cell.
#' The algorithm assumes the order of the cells in this scRNA-seq data matrix is the same as the order in the meta-data ('phenotypeData') and the cell-state space ('cellSpace').
#' @param phenotypeData A vector containing the meta-data levels of each of the cells.
#' @param cellSpace Cell coordinates matrix for the single-cell data. Each row represent a cell and each column corresponds to a specific dimension. Only 2-dim spaces are allowed.
#' @param genesToRun A vector of genes for which robustness will be calculated. Due to relatively long running times, it is advised to choose only genes with informative clusters inferred by CENA.
#' @param no_cores The number of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.
#' @param numberOfRepeats Number of CENA runs used to calculate cluster robustness. The default value is 100.
#' @param minCoverage The minimum coverage required for matching between the predicted clusters and the clusters in each iteration. The default value is 0.5.
#' @return A matrix where rows are genes and columns are different robustness measures.
#' \describe{
#' \itemize{
#'  \item{Repetitivity: }{ A fraction representing how much the predicted gene clusters are shared across iterations. }
#'  \item{Normalized ranking mean: }{ Normalized ranking scores mean. Values range between 0 to 1, where a lower value means the perdicted gene cluster is more likely to be selected across iterations. }
#'  \item{Normalized ranking Sd: }{ Normalized ranking scores standard deviation across all robustness iterations. }
#'  }
#'}
#' @examples
#' data(cellSpace)
#' data(geneExpressionDataMatrix)
#' data(phenotypeData)
#' # running CENA on 5 genes
#' results = CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1)
#' genesWithClusters = row.names(geneExpressionDataMatrix)[which(!is.na(rowSums(results$cluster_information)))]
#' robustnessResults = robustness(results, geneExpressionDataMatrix, phenotypeData, cellSpace, genesToRun = genesWithClusters, no_cores = 1)
#' @export
robustness = function(prediction,geneExpressionDataMatrix, phenotypeData, cellSpace, genesToRun, no_cores = NULL,numberOfRepeats = 100, minCoverage = 0.5){
  resolution_parameter = prediction$parameters$resolution_parameter; k1 = prediction$parameters$k1; k2 = prediction$parameters$k2; Tw = prediction$parameters$Tw; python_path = prediction$parameters$python_path;
  indexesToRun = which(names(prediction$cluster) %in% genesToRun)
  prediction$cluster = prediction$cluster[indexesToRun]
  if(length(indexesToRun) != length(genesToRun) |
     length(prediction$cluster)==0|
     length(which(unlist(lapply(prediction$cluster,function(cluster){is.na(cluster)}))))>0) {
    print("All genes must have a predicted cluster!")
  }else{
    print("Running robustness check (This might take a while):")
    genesToRunRep = rep(genesToRun,numberOfRepeats)
    CENAForRobustness = CENA_Main(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter, genesToRunRep, no_cores = no_cores, k1 = k1, k2 = k2, Tw = Tw, python_path = python_path, clusterDetails = T)

    sortedAllNodes = lapply(1:length(genesToRunRep),function(currGeneIndex){
      currSetOfLists = CENAForRobustness[[currGeneIndex]]
      if(length(currSetOfLists)==0){
        currSetOfLists
      }else{
        currSetOfLists[order(unlist(lapply(currSetOfLists, function(currSelection){
          cellNumbersSelected = as.numeric(as.matrix(gsub("cell_","",currSelection)))
          currGeneExp = geneExpressionDataMatrix[genesToRunRep[currGeneIndex],cellNumbersSelected]
          currPhen = phenotypeData[cellNumbersSelected]
          summary(lm(currGeneExp~currPhen))$r.squared
        })),decreasing = T)]
      }
    })
    repScores = unlist(lapply(1:length(genesToRunRep), function(currGeneIndex){
      allNodesSelected = sortedAllNodes[[currGeneIndex]]
      if(length(allNodesSelected)==0){ return(NA) }
      originalCells = prediction$cluster[[genesToRunRep[currGeneIndex]]]
      scores = unlist(lapply(allNodesSelected,function(currSelection){
        cellNumbersSelected = as.numeric(as.matrix(gsub("cell_","",currSelection)))
        length(intersect(cellNumbersSelected,originalCells))/length(originalCells)
      }))
      if(max(scores)>minCoverage){
        which.max(scores)
      }else{
        NA
      }
    }))

    meanRatio = as.data.frame(t(do.call(cbind, lapply(unique(genesToRunRep), function(currGene){
      nonNAs = which(!is.na(repScores[currGene==genesToRunRep]))
      repetitivity = length(nonNAs)/numberOfRepeats
      ratioForClusters = repScores[which(genesToRunRep == currGene)[nonNAs]]/
        unlist(lapply(1:length(genesToRunRep), function(currGeneIndex){
          length(CENAForRobustness[[currGeneIndex]])
        }))[which(currGene==genesToRunRep)[nonNAs]]
      c(repetitivity, mean(ratioForClusters),sd(ratioForClusters))
    }))))
    colnames(meanRatio) = c("Repetitivity","Average rank mean","Average rank Sd")
    row.names(meanRatio) = genesToRun

    return(meanRatio)
  }
}


#' Specificity check for CENA predictions
#'
#' Calculate the specificity scores for selected genes.
#'
#' @param prediction The result of the CENA run.
#' @param geneExpressionDataMatrix A matrix containing the single-cell RNA-seq data.
#' Each row corresponds to a certain gene and each column to a certain cell.
#' The algorithm assumes the order of the cells in this scRNA-seq data matrix is the same as the order in the meta-data ('phenotypeData') and the cell-state space ('cellSpace').
#' @param phenotypeData A vector containing the meta-data levels of each of the cells.
#' @param cellSpace Cell coordinates matrix for the single-cell data. Each row represent a cell and each column corresponds to a specific dimension. Only 2-dim spaces are allowed.
#' @param genesToRun A vector of genes for which specificity will be calculated. Due to relatively long running times, it is advised to choose only genes with informative clusters inferred by CENA.
#' @param no_cores The number of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.
#' @param numberOfRepeats Number of CENA runs used to calculate cluster specificity. The default value is 100.
#' @return A vector of p-values (specificity scores) for each of the genes.
#' @examples
#' data(cellSpace)
#' data(geneExpressionDataMatrix)
#' data(phenotypeData)
#' # running CENA on 5 genes
#' results = CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1)
#' genesWithClusters = row.names(geneExpressionDataMatrix)[which(!is.na(rowSums(results$cluster_information)))]
#' specificityResults = specificity(results, geneExpressionDataMatrix, phenotypeData, cellSpace, genesToRun = genesWithClusters, no_cores = NULL, numberOfRepeats =100)
#' @export
specificity = function(prediction,geneExpressionDataMatrix, phenotypeData, cellSpace, genesToRun, no_cores = NULL,numberOfRepeats = 100) {
  resolution_parameter = prediction$parameters$resolution_parameter; k1 = prediction$parameters$k1; k2 = prediction$parameters$k2; Tw = prediction$parameters$Tw; python_path = prediction$parameters$python_path;
  indexesToRun = which(names(prediction$cluster) %in% genesToRun)

  if(any(is.na(indexesToRun)) | length(indexesToRun) != length(genesToRun) |
     length(prediction$cluster[indexesToRun])==0|
     length(which(unlist(lapply(prediction$cluster[indexesToRun],function(cluster){is.na(cluster)}))))>0) {
    print("All genes must have a predicted cluster!")
  } else{
    R2DiffReal = prediction$cluster_information[genesToRun,"r_squared"]-
      unlist(lapply(genesToRun,function(currGene){
        summary(lm(geneExpressionDataMatrix[currGene,]~phenotypeData))$r.squared
      }))
    genesToRunRep = rep(genesToRun,numberOfRepeats)
    permGEDataForRuns = t(apply(geneExpressionDataMatrix[genesToRunRep,],1,sample))
    print("Running the specificity check:")
    permCENARun = CENA_Main(permGEDataForRuns, phenotypeData, cellSpace, resolution_parameter, row.names(permGEDataForRuns), no_cores = no_cores, k1 = k1, k2 = k2, Tw = Tw, python_path = python_path, clusterDetails = F)
    R2DiffPerm = permCENARun$cluster_information[,"r_squared"]-unlist(lapply(1:dim(permGEDataForRuns)[1],function(i){
      summary(lm(permGEDataForRuns[i,]~phenotypeData))$r.squared
    }))
    p_values = unlist(lapply(genesToRun, function(currGene){
      currPermR2 = R2DiffPerm[genesToRunRep==currGene]
      currRealR2 = R2DiffReal[genesToRun == currGene]
      (length(which(currRealR2<=currPermR2))+1)/(length(currPermR2)+1)
    }))
    names(p_values) = genesToRun
    return(p_values)
  }
}

#' The Cell Niche Associations (CENA) algorithm
#'
#' CENA is a method for a joint identification of pairwise association together with
#' the particular subset of cells in which the association is detected.
#' In this implementation, CENA is limited to associations between the genes' expression levels
#' (data from scRNA-sequencing) and an additional cellular meta-data of choice.
#' Note that python3 should be installed on the machine, with the following libraries: numpy, igraph, leidenalg.
#'
#' @param geneExpressionDataMatrix A matrix containing the single-cell RNA-seq data.
#' Each row corresponds to a certain gene and each column to a certain cell.
#' The algorithm assumes the order of the cells in this scRNA-seq data matrix is the same as the order in the meta-data ('phenotypeData') and the cell-state space ('cellSpace').
#' @param phenotypeData A vector containing the meta-data levels of each of the cells.
#' @param cellSpace Cell coordinates matrix for the single-cell data. Each row represent a cell and each column corresponds to a specific dimension. Only 2-dim spaces are allowed.
#' @param resolution_parameter The resolution parameter of Leiden. Cutoff as desired; Higher values provide smaller clusters.
#' @param no_cores The number of cores which will be used for the analysis. The defalt (NULL) is total number of cores minus 1.
#' @param k1 The K1 parameter. The default value is 1 percentage of the cells. k1 should be kept relatively small.
#' @param k2 The K2 paremeter. The default value is 10. k2 have relatively small effect on results.
#' @param Tw The subset size cutoff (Tw). The default value is 30. Cutoff as desired; high values filter out small subsets.
#' @param genesToRun A vector of genes for which associations should be inferred. The default value is all genes in the input scRNA-seq matrix.
#' However, to reduce running times, it is recommended to use a small subset of the genes.
#' @param python_path The location of python on the machine. If not specified, the default python path will be selected.
#' Note that python3 and the following libraries should be installed: numpy, igraph, leidenalg.
#' @return
#' \describe{
#'  \item{cluster_information}{Association statistics for cell subsets obtained by CENA for each of the genes (rows).
#'  Columns corresponds to the r.squared, slope, intercept and p.value of these associations.}
#'  \item{cluster}{A list containing, for each gene, the identified cell subset in which the association occurs.}
#'  \item{parameters}{The initial parameters used by CENA.}
#'}
#' @examples
#' data(cellSpace)
#' data(geneExpressionDataMatrix)
#' data(phenotypeData)
#' # running CENA on 5 genes
#' results = CENA(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter = 8, no_cores = 1)
#' @export
CENA = function(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter, genesToRun = row.names(geneExpressionDataMatrix),
                no_cores = NULL, k1 = NULL, k2 = 10, Tw = 30, python_path = NULL) {
  print("Running CENA (This might take a while):")
  the_result = CENA_Main(geneExpressionDataMatrix, phenotypeData, cellSpace, resolution_parameter,genesToRun,
                         no_cores, k1, k2, Tw, python_path)
  ##### Combining cell predictions #####
  print("Finalizing...")
  return(the_result)
}
