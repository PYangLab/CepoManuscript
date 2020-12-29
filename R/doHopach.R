# code has been adapted from https://github.com/SydneyBioX/scClassify to include `labels` argument

runHOPACH <- function (data, plot = TRUE, kmax = 5, labels=labels) {
    
    dist <- hopach::distancematrix(t(data), d = "cor")
    clustresult <- hopach::hopach(t(data), dmat = dist, kmax = kmax)
    final_labels <- strsplit(as.character(format(clustresult$final$labels, 
                                                 scientific = FALSE)), "")
    level_list <- list()
    for (i in seq_len(nchar(clustresult$final$labels[1]))) {
        if (i != 1) {
            level_list[[i]] <- paste(level_list[[i - 1]], unlist(lapply(final_labels, 
                                                                        "[[", i)), sep = "")
        }
        else {
            level_list[[i]] <- unlist(lapply(final_labels, "[[", 
                                             i))
        }
    }
    cutree_list <- list()
    cutree_list[[1]] <- rep(1, length(colnames(data)))
    names(cutree_list[[1]]) <- colnames(data)
    cutree_list[2:(length(level_list) + 1)] <- lapply(level_list, 
                                                      function(x) {
                                                          x <- as.numeric(as.factor(x))
                                                          names(x) <- colnames(data)
                                                          x
                                                      })
    cutree_list[[length(cutree_list) + 1]] <- seq_len(length(clustresult$final$labels))
    names(cutree_list[[length(cutree_list)]]) <- colnames(data)
    if (plot) {
        treePlot <- plotCellTypeTree(cutree_list, labels=labels)
        return(list(cutree_list = cutree_list, plot = treePlot))
    }
    else {
        return(list(cutree_list = cutree_list))
    }
}

plotCellTypeTree <- function (cutree_list, group_level = NULL, labels = TRUE) {
    cutree_list <- cutree_list
    if (is.null(group_level)) {
        group_level <- 3
    }
    E_list <- list()
    for (i in seq_len(length(cutree_list))) {
        if (i != 1) {
            parent <- cutree_list[[i - 1]]
            E_list[[i]] <- paste(paste(i - 1, parent, sep = ""), 
                                 paste(i, cutree_list[[i]], sep = ""), sep = "_")
        }
    }
    V_list <- list()
    for (i in seq_len(length(cutree_list))) {
        V_list[[i]] <- paste(i, cutree_list[[i]], sep = "")
        names(V_list[[i]]) <- names(cutree_list[[i]])
    }
    g <- igraph::make_graph(unlist(strsplit(unique(unlist(lapply(E_list, 
                                                                 sort))), "_")))
    group_level <- group_level
    if (labels == FALSE) {
        group_col <- (V_list[[group_level]])
        names(group_col) <- V_list[[length(cutree_list)]]
        V_group <- rep(NA, length(igraph::V(g)))
        names(V_group) <- names(igraph::V(g))
        V_group[V_list[[length(V_list)]]] <- group_col[V_list[[length(V_list)]]]
        igraph::V(g)$group <- V_group
    } else {
        group_col <- color_annot$celltype_major[annot$celltype_major]
        names(group_col) <- V_list[[length(cutree_list)]]
        V_group <- rep(NA, length(igraph::V(g)))
        names(V_group) <- names(igraph::V(g))
        V_group[V_list[[length(V_list)]]] <- group_col[V_list[[length(V_list)]]]
        igraph::V(g)$group <- V_group
    }
    V_name <- rep(NA, length(igraph::V(g)))
    names(V_name) <- names(igraph::V(g))
    V_name[V_list[[length(V_list)]]] <- names(V_list[[length(V_list)]])
    igraph::V(g)$name <- V_name
    treePlot <- ggraph::ggraph(g, layout = "dendrogram") + ggraph::geom_edge_diagonal(colour = "black") + 
        ggraph::geom_node_text(aes(vjust = 1, hjust = -0.2, label = name), 
                               size = 2.7, alpha = 1) + ggraph::geom_node_point(aes(filter = leaf, 
                                                                                    alpha = 0.2, size = 3, colour = as.factor(igraph::V(g)$group))) + 
        ggplot2::theme_void() + ggplot2::scale_colour_hue() + 
        ggplot2::theme(legend.position = "none") + ggplot2::coord_flip() + 
        NULL
    treePlot
}