library(catmaid)
library(data.table)
library(tidyverse)

# get list of cell types in project --------------------------------------------
# assumes they have an annotation starting with "celltype:"
get_celltypes <- function(pid) {
  annotations <- catmaid_get_annotationlist(pid = 35)
  celltypes <- annotations$annotations |> filter(grepl("^celltype:", name)) |>
    select(name) |> pull() |> sub(pattern="^celltype:", replacement="")
  return(celltypes)
}

# get skids with multiple annotations ------------------------------------------
# note that second argument needs a vector if you are using more than 1 annotation!
get_skids_with_annot <- function(pid, annotations) {
  skids_list <- lapply(annotations, catmaid_skids, pid=35)
  intersection <- Reduce(intersect,skids_list)
  return(intersection)
}

# find skeletons without annotation --------------------------------------------
# if used without annotation name, returns skids without any annotation
get_skels_without_annot <- function(pid, annotationname = "") {
  skids <- unlist(catmaid_fetch(path = paste("/", pid, "/skeletons/", sep="")))
  # I'm using the above method of getting all skids instead of catmaid_skids(".*", pid)
  # because the above method is 100x faster (~0.06s vs ~6.05s for 988 skeletons)
  # it is not posible to omit the above line because 
  # catmaid_get_annotations_for_skeletons(".*", pid) doesn't return skeletons with no annotations at all
  skids_with_annot <- catmaid_get_annotations_for_skeletons(skids, pid = pid) |>
    filter(grepl(annotationname, annotation)) |> select(skid) |> 
    pull() |> unique()
  skids_without_annot <- setdiff(skids, skids_with_annot)
  return(skids_without_annot)
}


# crop substack ----------------------------------------------------------------
crop_substack <- function(tagname,
                         half_bb_size_x, half_bb_size_y, half_bb_size_z,
                         zoomlevel,
                         dest_dir,
                         pid, stackid) {
  # tag search only returns treenodes with tag, not connectors
  # however, we want to get either treenodes or connectors, so we have to use generic search
  tagname_edited <- gsub(" ", "%20", tagname)
  search_results <- catmaid_fetch(path = paste(pid, "/search?substring=", tagname_edited, sep = ""))
  # search doesn't accept regex, so we have to filter results
  ids=c() # keep track of the number of substacks we have to download from the server later
  # use ids rather than simple counter because it helps with debugging
  for (i in seq_along(search_results)) {
    if (search_results[[i]]$name == tagname && search_results[[i]]$class_name == "label") {
      for (j in seq_along(search_results[[i]]$nodes)) {
        ids <- c(ids, search_results[[i]]$nodes[[j]]$id)
        x <- search_results[[i]]$nodes[[j]]$x
        y <- search_results[[i]]$nodes[[j]]$y
        z <- search_results[[i]]$nodes[[j]]$z
        catmaid_fetch(
          path = paste(pid, "/crop/", sep = ""),
          body = list(
            stack_ids=stackid,
            min_x=x-half_bb_size_x,
            min_y=y-half_bb_size_y,
            min_z=z-half_bb_size_z,
            max_x=x+half_bb_size_x,
            max_y=y+half_bb_size_y,
            max_z=z+half_bb_size_z,
            zoom_level=zoomlevel,
            single_channel='true',
            rotationcw=0
          )
        )
      }
      for (j in seq_along(search_results[[i]]$connectors)) {
        ids <- c(ids, search_results[[i]]$connectors[[j]]$id)
        x <- search_results[[i]]$connectors[[j]]$x
        y <- search_results[[i]]$connectors[[j]]$y
        z <- search_results[[i]]$connectors[[j]]$z
        catmaid_fetch(
          path = paste(pid, "/crop/", sep = ""),
          body = list(
            stack_ids=21,
            min_x=x-half_bb_size_x,
            min_y=y-half_bb_size_y,
            min_z=z-half_bb_size_z,
            max_x=x+half_bb_size_x,
            max_y=y+half_bb_size_y,
            max_z=z+half_bb_size_z,
            zoom_level=zoomlevel,
            single_channel='true',
            rotationcw=0
          )
        )
      }
    }
  }
  
  if (length(ids) == 0) {
    print(paste("No treenodes or connectors tagged with", tagname, "were found"))
    return()
  }
  
  print(paste("Found", length(ids), "treenodes or connectors tagged with", tagname))
  print("Downloading the following ids:")
  print(ids)
  print("This will take some time")
  
  # server needs time to process this
  Sys.sleep(90)
  
  # download stacks ---------------------------------------------------------
  # the exact file names are needed for download,
  # to get them we need to check messages from the server
  
  # It would be cool if I can put the synapse ID as file name,
  # but it's not guaranteed that jobs will return in the same order
  
  cat_messages <- catmaid_fetch(
    path = "messages/list"
  )
  
  for (k in seq_along(ids)) {
    cat_message <- cat_messages[[k]]$text
    link <- regmatches(cat_message, 
                       regexpr("/crop/download/crop_.{6}.tiff", cat_message))
    filename <- regmatches(link, 
                           regexpr("crop_.{6}.tiff", link))
    full_link <- paste("https:/catmaid.ex.ac.uk", link, sep = "")
    destpath <- paste(dest_dir, "/", filename, sep = "")
    download.file(full_link, destfile = destpath)
  }
}


# number of postsynaptic sites of a connector ----------------------------------
n_pre_post <- function(connector_id, pid) {
  connector_info <- catmaid_fetch(path = paste(pid, "/connectors/", connector_id, sep = ""))
  # pre should always be 1, but maybe it's good to check if there is any weird stuff going on
  pre <- connector_info$partners %like% 'presynaptic_to'
  n_pre <- table(pre)["TRUE"][[1]]
  post <- connector_info$partners %like% 'postsynaptic_to'
  n_post <- table(post)["TRUE"][[1]]
  return(n_post)
}

# synapse polyadicity statistics for a project ---------------------------------
get_syn_polyadicity <- function(pid) {
  pre_connectors <- catmaid_fetch(
    path = paste(pid, "/connectors/", sep = ""),
    body = list(
      relation_type = "presynaptic_to",
      with_partners = "false"
    )
  )

  pre_connector_IDs <- lapply(pre_connectors$connectors,'[[',1)
  
  n_post <- lapply(pre_connector_IDs, n_pre_post, pid)
  
  n_post |> as.numeric() |> summary()
}


