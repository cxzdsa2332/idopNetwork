#' @title convert ODE results(ODE_solving2) to basic network plot table
#' @param result list result from qsODE_parallel
#' @return a list with basic information to plot network
#' @export
network_conversion <- function(result){
  n = ncol(result$fit)

  effect.mean = apply(result$fit,2,mean)[4:n]
  effect.predict.mean = apply(result$predict,2,mean)[4:n]
  effect.total = colSums(result$fit)[4:n]

  temp = matrix(NA,nrow = n-3, ncol=3)
  colnames(temp) = c("From", "To", "Effect")
  temp[,1] = colnames(result$fit)[4:n]
  temp[,2] = colnames(result$fit)[4]
  temp[,3] = effect.predict.mean
  if (nrow(temp)==2) {
    temp = t(data.frame(temp[-1,]))
  } else{
    temp = data.frame(temp[-1,])
  }

  temp2 = matrix(NA,nrow = n-3, ncol=3)
  colnames(temp2) = c("From", "To", "Effect")
  temp2[,1] = colnames(result$fit)[4:n]
  temp2[,2] = colnames(result$fit)[4]
  temp2[,3] = effect.total
  if (nrow(temp2)==2) {
    temp2 = t(data.frame(temp2[-1,]))
  } else{
    temp2 = data.frame(temp2[-1,])
  }
  output <- list(ind.name = colnames(result$fit)[4],
                 dep.name = colnames(result$fit)[5:n],
                 ODE.par = result$ODE.value,
                 ind.par = result$LOP_par[,3],
                 dep.par = result$LOP_par[,4:(n-1)],
                 effect.mean = effect.predict.mean,
                 effect.total = effect.total,
                 effect.all = result$fit,
                 edge = temp,
                 edge.total = temp2,
                 ind.effect = effect.predict.mean[1])
  return(output)
}

#' @title convert ODE results(ODE_solving2) to basic network plot table
#' @param result list result from qsODE_parallel
#' @return a list with basic information to plot network
#' @export
network_maxeffect <- function(result){
  module = result[[1]]
  after <- data.frame(do.call(rbind, lapply(module, "[[", "edge.total")))
  after$Effect = as.numeric(after$Effect)
  Module.maxeffect = max(abs(after$Effect))

  if (length(result) == 2) {
    after <- data.frame(do.call(rbind, lapply(result[[2]], "[[", "edge")))
    after$Effect = as.numeric(after$Effect)
    res.maxeffect = max(abs(after$Effect))
  } else{
    after = lapply(2:length(result),function(c) data.frame(do.call(rbind, lapply(result[[c]], "[[", "edge"))))
    after = do.call(rbind,after)
    after$Effect = as.numeric(after$Effect)
    res.maxeffect = max(abs(after$Effect))
  }
  maxeffect = max(c(Module.maxeffect,res.maxeffect))
  return(maxeffect)
}



#' @title generate network plot
#' @import igraph
#' @importFrom stats aggregate
#' @param result list result from network_conversion
#' @param title text for plot title
#' @param type select module effect or microbe effect
#' @param maxeffect control edge size when compare networks
#' @return network plot
#' @export
network_plot <- function(result, title = NULL, maxeffect = NULL, type = NULL){
  #ind effect control node size
  extra <- sapply(result,"[[", "ind.effect")

  if (is.null(type)) {
    after <- data.frame(do.call(rbind, lapply(result, "[[", "edge")))
  } else{
    after <- data.frame(do.call(rbind, lapply(result, "[[", "edge.total")))

  }
  after$Effect = as.numeric(after$Effect)
  rownames(after) = NULL


  after$edge.colour = NA
  for (i in 1:nrow(after)) {
    if(after$Effect[i]>=0){
      after$edge.colour[i] = "#FE433C"
    } else{
      after$edge.colour[i] = "#0095EF"
    }
  }

  #nodes
  nodes <- data.frame(unique(after[,2]),unique(after[,2]),extra)
  colnames(nodes) <- c("id","name","ind_effect")
  nodes$influence <- aggregate(Effect ~ To, data = after, sum)[,2]
  nodes$node.colour = NA
  for (i in 1:nrow(nodes)) {
    if(nodes$influence[i]>=0){
      nodes$node.colour[i] = "#FFC4C4"
    } else{
      nodes$node.colour[i] = "#89CFFD"
    }
  }

  #normalization
  if (is.null(maxeffect)) {
    after[,3] <- normalization(abs(after[,3]))
    nodes[,3:4] <- normalization(abs(nodes[,3:4]))
  } else{
    after[,3] <- (abs(after[,3]))/maxeffect*1.5+0.1
    nodes[,3:4] <- (abs(nodes[,3:4]))/maxeffect*1.5+0.1
  }

  net <- graph_from_data_frame( d=after,vertices = nodes,directed = T )

  #layout
  l <- layout_randomly(net)

  plot.igraph(net,
              vertex.label=V(net)$name,
              vertex.label.color="black",
              vertex.shape="circle",
              vertex.label.cex=V(net)$ind_effect*1.5,
              vertex.size=V(net)$ind_effect*20+5,
              edge.curved=0.05,
              edge.color=E(net)$edge.colour,
              edge.frame.color=E(net)$edge.colour,
              edge.width=E(net)$Effect*5,
              vertex.color=V(net)$node.colour,
              layout=l,
              main=title,
              margin=c(-.05,-.05,-.05,-.05)
  )
}

