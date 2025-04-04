library(catmaid)
library(tidyverse)

source("~/R/conn.R")

ciliated_celltypes <- c("akrotroch", "prototroch", "metatroch", "paratroch", "telotroch", "nuchal_ciliated", "crescentcell")


get_bb_tagged_treenodes <- function(skid) {
  neuron <- read.neurons.catmaid(
    skid, 
    pid=11, 
    fetch.annotations = T)
  treenodes <- neuron[1][[1]]$tags$`basal body`
  return(treenodes)
}

get_bb_tagged_treenodes_for_annot <- function(annot) {
  annotation_query = paste("^", annot,"$", sep="")
  skids_with_annotation <- catmaid_skids(annotation_query, pid = 11)
  bb_taged_treenodes <- lapply(skids_with_annotation, get_bb_tagged_treenodes) %>%
    unlist()
  return(bb_taged_treenodes)
}

dist_to_colosest_other_bb_node <- function(query_node) {
  # this doesn't check if the other node is on same skeleton,
  # but in this particular case we can assume it is
  pos_node <-  catmaid_fetch(path = "11/nodes/location", 
                             body = list(node_ids=query_node))
  xyz <- c(pos_node[[1]][[2]],
           pos_node[[1]][[3]],
           pos_node[[1]][[4]])
  other_node <-  catmaid_fetch(path = "11/nodes/find-labels", 
                               body = list(
                                 x=xyz[1],
                                 y=xyz[2],
                                 z=xyz[3],
                                 label_regex="^basal body$"
                               ))
  # get second node, because the first node returned is the query node
  # technically this distance is to the query location, not to the query node itself,
  # but the difference is negligible
  distance <- other_node[[2]][[3]]
  return(distance)
}


bb_distances_by_celltype <- data.frame()
for (ciliated_celltype in ciliated_celltypes) {
  nodes_sampled <- get_bb_tagged_treenodes_for_annot(ciliated_celltype) %>%
    sample(100)
  distances <- vector()
  for (node in nodes_sampled) {
    distance <- dist_to_colosest_other_bb_node(node)
    distances <- c(distances, distance)
  }
  median_distance <- median(distances)
  df <- data.frame(celltype=ciliated_celltype, median_distance=median_distance)
  bb_distances_by_celltype <- rbind(bb_distances_by_celltype, df)
}





#############################################################################################
# method2: distances of all basal bodies of all skeletons

get_tag_treenodes <- function(skid) {
  neuron <- read.neurons.catmaid(
    167009, 
    pid=11, 
    fetch.annotations = T)
  treenodes <- neuron[1][[1]]$tags$`basal body`
  return(treenodes)
}

distance_between_2_nodes <- function(treenode1, treenode2) {
  pos1 <-  catmaid_fetch(path = "11/nodes/location", 
                         body = list(node_ids=treenode1))
  xyz1 <- c(pos1[[1]][[2]],
            pos1[[1]][[3]],
            pos1[[1]][[4]])
  pos2 <-  catmaid_fetch(path = "11/nodes/location", 
                         body = list(node_ids=treenode2))
  xyz2 <- c(pos2[[1]][[2]],
            pos2[[1]][[3]],
            pos2[[1]][[4]])
  x <- rbind(xyz1, xyz2)
  distance <- stats::dist(x, method = "euclidean") %>%
    as.vector()
  return(distance)
}

tagged_treenodes <- get_tag_treenodes(167009)
distance <- 999999999999
distances <- data.frame()
done_treenodes <- c()
for (treenode in tagged_treenodes) {
  done_treenodes <- c(treenode, done_treenodes)
  other_treenodes <- setdiff(tagged_treenodes, done_treenodes)
  for (other_treenode in other_treenodes) {
    dist_new <- distance_between_2_nodes(treenode, other_treenode)
    print(paste("node1 is ", treenode))
    print(paste("node2 is ", other_treenode))
    if (dist_new < distance) {
      distance <<- dist_new
      df <- data.frame(node1=treenode, node2=other_treenode, distance=distance)
    }
  }
  distances <- rbind(distances, df)
  # add other treenode to done treenodes
  # to not do the same pair twice
  done_treenodes <- c(treenode, df$node2)
}


