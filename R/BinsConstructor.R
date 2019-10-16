

#' Reference Class to represent a bin for creating a tree of bins
#' In which the leaf bins are the final bins that will be used
#'
#' @field range A matrix which holds the start and the end values for each dimension
#' @field count The number of the points
#' @field parent The parent bin of the bin in the tree
#' @field isLeftSon is the son in the tree is left son
#' @importFrom methods setRefClass new
#' @keywords internal
BinNode <- setRefClass("BinNode",
                       fields = c("range", # matrix [dimsX2]
                                "count"#, # numeric
                                #"parent",#"BinNode"
                                #"isLeftSon"
                                )) # boolean

#' Reference Class to represent an internal bin in the bins tree
#'
#' @field left_son the left son of the internal bin
#' @field right_son the right son of the internal bin
#' @field midpoint the value that seperate the bin to left and right bins
#' @field discriminator_dim the dimension by which the bin was separated
#' @keywords internal
# InternalBinNode <- setRefClass("InternalBinNode",contains = "BinNode",
#                                fields=c(#"left_son", # BinNode
#                                         #"right_son", # BinNode
#                                         "midpoint", # numeric
#                                         "discriminator_dim")) # numeric
#' Reference Class to represent an leaf bin in the bins tree
#'
#' @field count_diff the different points number between the left and the right son
#' @field points a matrix of the points in this bin
#' @keywords internal
LeafBinNode <- setRefClass("LeafBinNode", contains = "BinNode",
                           fields = c("count_diff", # vector [dims]
                                    "points"), # matrix [pointsXdims]
                           methods = list(
                             split = function(dim) {
                               "split the bin into 2 sub-bins in the simension in which the count_diff is the largest"
                               midpoint = (max(points[,dim])+min(points[,dim]))/2 # the median point in the splitting dimension
                               ### left son
                               range_left = range
                               range_left[dim,2] = as.numeric(midpoint)

                               points_left = points[which(points[,dim]<midpoint),]
                               if(is.null(dim(points_left))) {
                                 points_left = t(as.matrix(points_left))
                               }

                               count_left = nrow(points_left)
                               count_diff_left = sapply(seq(1,ncol(points_left)), function(dim) {
                                 midpoint =  (max(points_left[,dim])+min(points_left[,dim]))/2
                                 return(abs(length(which(points_left[,dim] < midpoint)) -
                                              length(which(points_left[,dim] >= midpoint))))
                               })

                               left_son = LeafBinNode$new(range = range_left,
                                                          count = count_left,
                                                          points = points_left,
                                                          count_diff = count_diff_left#,
                                                          #isLeftSon = TRUE
                                                          )
                               #### right son
                               range_right = range
                               range_right[dim,1] = midpoint
                               points_right = points[which(points[,dim]>midpoint),]
                               if(is.null(dim(points_right))) {
                                 points_right = t(as.matrix(points_right))
                               }

                               count_right = nrow(points_right)

                               count_diff_right = sapply(seq(1,ncol(points_right)), function(dim) {
                                 midpoint =  (max(points_right[,dim])+min(points_right[,dim]))/2
                                 return(abs(length(which(points_right[,dim] < midpoint)) -
                                              length(which(points_right[,dim] >= midpoint))))
                               })

                               right_son = LeafBinNode$new(range = range_right,
                                                           count = count_right,
                                                           points = points_right,
                                                           count_diff = count_diff_right#,
                                                           #isLeftSon = FALSE
                                                           )




                               #### the substitution node
                               # substitute_node = InternalBinNode$new(#left_son=left_son,
                               #                                       #right_son=right_son,
                               #                                       #parent = parent,
                               #                                       count = count,
                               #                                       range = range,
                               #                                       discriminator_dim = dim,
                               #                                       midpoint = midpoint,
                               #                                       #isLeftSon = isLeftSon
                               #                                       )

                               #left_son$parent = substitute_node
                               #right_son$parent = substitute_node
                               # if(isLeftSon) {
                               #   parent$left_son <<- substitute_node
                               # } else {
                               #   parent$right_son <<- substitute_node
                               # }
                               #return(list(substitute_node, left_son, right_son))
                               return(list(left_son, right_son))
                             }))

# By a given cellSpace distribution, build bins which separate the points in such a awy that each bin contains approximetly the same number of points
#' @keywords internal
getBins = function(xyLocations, alpha = 1.2, beta = 0.6, P_target) {
  root = LeafBinNode$new(points = xyLocations,
                         range = t(sapply(seq(1, ncol(xyLocations)), function(dim) {
                           c(min(xyLocations[,dim]), max(xyLocations[,dim]))
                         })),
                         count = nrow(xyLocations),
                         #parent = NULL,
                         #isLeftSon = FALSE,
                         count_diff = sapply(seq(1,ncol(xyLocations)), function(dim) {
                           midpoint =  (max(xyLocations[,dim])+min(xyLocations[,dim]))/2
                           return(abs(length(which(xyLocations[,dim] < midpoint)) -
                                        length(which(xyLocations[,dim] >= midpoint))))
                         })
  )
  #firstTime = TRUE
  pool_to_check = c(root)

  leafs= c()
  while(length(pool_to_check) > 0) {
    curr_node = pool_to_check[[1]]
    pool_to_check = pool_to_check[-1]
    if(curr_node$count > alpha * P_target || any(curr_node$count_diff > beta * P_target)) {
      dim_to_split_by = which(curr_node$count_diff == max(curr_node$count_diff))
      if(length(dim_to_split_by)>1) {
        dim_to_split_by_range = sapply(dim_to_split_by, function(dim){
          curr_node$range[dim,2]-curr_node$range[dim,1]
        })
        max_range = which(dim_to_split_by_range == max(dim_to_split_by_range))
        dim_to_split_by = dim_to_split_by[max_range]
      }
      result = curr_node$split(dim_to_split_by)
      #left_son = result[[2]]
      #right_son = result[[3]]
      left_son = result[[1]]
      right_son = result[[2]]
      pool_to_check = c(pool_to_check, left_son)
      pool_to_check = c(pool_to_check, right_son)
      # if(firstTime) {
      #   root = result[[1]]
      #   firstTime = FALSE
      # }
    } else {
      leafs = c(leafs, curr_node)
    }
  }

  #return(list(root, leafs))
  return(leafs) ## NOTICE THAT THE RETURN IS DIFFERENT
}


# By a given bin, return its neighbor bins
#' @keywords internal
getNeighborBins = function(nodeLeaf, leafs) {
  bin_range = nodeLeaf$range
  neighbors = c()
  count = 0
  lapply(leafs, function(leaf) {
    count <<- count +1
    anydim = sapply(seq(1,nrow(bin_range)), function(dim_num) {
      if(bin_range[dim_num,1] == leaf$range[dim_num,2] ||
         bin_range[dim_num,2] == leaf$range[dim_num,1]) {
        # should be one btween the oder in the other dimensions
        the_oder_dimensions = seq(1,nrow(bin_range))[-dim_num]
        alldims = sapply(c(the_oder_dimensions), function(dim_num) {
          return(bin_range[dim_num,1] >= leaf$range[dim_num,1] && bin_range[dim_num,1] <= leaf$range[dim_num,2] ||
                   bin_range[dim_num,2] >= leaf$range[dim_num,1] && bin_range[dim_num,2] <= leaf$range[dim_num,2] ||
                   leaf$range[dim_num,1] >= bin_range[dim_num,1] && leaf$range[dim_num,1] <= bin_range[dim_num,2] ||
                   leaf$range[dim_num,2] >= bin_range[dim_num,1] && leaf$range[dim_num,2] <= bin_range[dim_num,2])
        })
        return(all(alldims==TRUE))
      } else {return(FALSE)}

    })
    if(any(anydim == TRUE)) {
      neighbors <<- c(neighbors, leaf)
    }
  })
  return(neighbors)
}

# plotTheBins = function(xyLocations, leafs) {
#   xmin = sapply(leafs, function(x) {
#     x$range[1,1]
#   })
#   xmax = sapply(leafs, function(x) {
#     x$range[1,2]
#   })
#   ymin = sapply(leafs, function(x) {
#     x$range[2,1]
#   })
#   ymax = sapply(leafs, function(x) {
#     x$range[2,2]
#   })
#
#   d = data.frame(x1 = xmin, x2 = xmax, y1 = ymin, y2 = ymax, r = seq(1,length(leafs)))
#   library(ggplot2)
#   ggplot() +
#     scale_x_continuous(name="x") +
#     scale_y_continuous(name="y") +
#     geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0.5) +
#     geom_point(data=as.data.frame(xyLocations), mapping=aes(x=V1,y=V2),size=1)+
#     geom_text(data=d, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=r), size=4)
# }

##main - melanoma
# runFolder = "/Users/mayalevy/Downloads/MarkerSelection/Runs/T.CD8_1543136145"
##main - influenza
# runFolder = "/Users/mayalevy/Downloads/MarkerSelection/Runs/T-cells_1543758932"
#
# finalResultsFolder = file.path(runFolder, "finalResults")
# finalResultsFolder_temps = file.path(finalResultsFolder, "temps")
# xyLocations = readRDS(file=file.path(finalResultsFolder_temps, "xyLocations"))
# result = getBins(xyLocations)
# plotTheBins(xyLocations, result)
# table(sapply(result, function(x){
#   x$count
# }))

# nodeLeaf = result[[2]][[71]]
# leafs = result[[2]]
# neighbors = getNeighborBins(nodeLeaf, leafs)

