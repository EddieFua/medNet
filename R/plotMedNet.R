library(grDevices)
library(igraph)

plotMedNet <- function(global.net, percentage = 1, spdf = TRUE) {
  # Visualize the global network generated from `findGlobalNet` given percentage of the predictive exposures
  #
  # Args:
  #   `global.net`: undirected global network generated from `findGlobalNet`
  #   `percentage`: percentage of predictive exposures to be shown based on their AUC
  #   `spdf`: save the figure in a PDF file or not
  #
  # Returns:
  #   `result`: a recordedplot object of the network plot

  if (percentage > 1 | percentage <= 0) {
    cat("The percentage should be between 0 and 1!\n")
  } else if(percentage != 1) {
    pname <- paste0("Top_", 100 * percentage, "%_predictive_exposures")

    # Indicate the top exposures
    rank <- unique(global.net$final_net[, 1])
    rank <- rank[1:floor(length(rank) * percentage)]

    # Generate the cutted network
    final_net <- global.net$final_net[global.net$final_net[, 1] %in% rank, ]
    index     <- igraph::V(global.net$net)$name %in% c(final_net[, 1], final_net[, 2])
    node      <- igraph::V(global.net$net)$name[index]
    vtype     <- igraph::V(global.net$net)$vtype[index]
    net       <- igraph::graph_from_data_frame(final_net, directed = F, vertices = data.frame(node, vtype))
  } else {

    pname <- "All_predictive_exposures"
    # Generate the full network
    #final_net <- global.net$final_net
    #index     <- igraph::V(global.net$net)$name %in% c(final_net[, 1], final_net[, 2])
    #node      <- igraph::V(global.net$net)$name[index]
    #vtype     <- igraph::V(global.net$net)$vtype[index]
    #net       <- igraph::graph_from_data_frame(final_net, directed = F, vertices = data.frame(node, vtype))
    final_net <- global.net$final_net
    index     <- igraph::V(global.net$net)$name %in% c(final_net[, 1], final_net[, 2])
    node      <- igraph::V(global.net$net)$name
    vtype     <- igraph::V(global.net$net)$vtype
    net       <- global.net$net
  }

  # Assign colors to vertices
  colrs <- c("Orange","gray70")
  igraph::V(net)$color <- colrs[igraph::V(net)$vtype]

  # Format the vertice label text size
  igraph::V(net)$label.cex <- ifelse(igraph::V(net)$vtype == 1, 6, 4) / (2 * max(nchar(node)))

  # Remove the frame of vertices
  igraph::V(net)$frame.color <- "white"

  # Format the node size
  igraph::V(net)$size <- ifelse(igraph::V(net)$vtype == 1, 10, 8)

  # Use Fruchterman and Reingold layout algorithm
  l <- igraph::layout.fruchterman.reingold(net)

  # Save plot object
  if(length(unique(final_net$type))==1){labels=c("Orange")}else{labels = c("Orange", "gray70")}
  plot(net, edge.color = as.character(factor(final_net$type, labels = labels)),
       edge.arrow.size = 0.2, edge.width = 2, layout = l, main = gsub("_", " ", pname))
  result <- grDevices::recordPlot()
  if(spdf) {
    pdf(paste0(gsub("%", "p", pname), ".pdf"), width = 10, height = 10, paper = 'special')
    plot(net, edge.color = as.character(factor(final_net$type, labels = labels)),
         edge.arrow.size = 0.2, edge.width = 2, layout = l, main = gsub("_", " ", pname))
    dev.off()
  }

  return(result)
}

