plotHeatmap <- function(sce, gene.list, celltype, slim=FALSE, verbose=F, scale=T) {
    
    #sce=sce
    #gene.list=gene.list.ds
    
    mat.ordered <- do.call(cbind, sapply(split(seq_len(ncol(sce)), celltype), function(i) logcounts(sce)[, i]))
    label.ordered <- do.call(c, sapply(split(seq_len(ncol(sce)), celltype), function(i) celltype[i]))
    
    exp.sel <- mat.ordered[unlist(gene.list),]
    rownames(exp.sel) <- paste(unlist(gene.list), 1:length(unlist(gene.list)), sep="_")
    annot <- data.frame(location = label.ordered)
    rownames(annot) <- colnames(mat.ordered)
    
    annot_row <- data.frame(genes = rep(names(gene.list), times=sapply(gene.list, length)))
    rownames(annot_row) <- paste(unlist(gene.list), 1:length(unlist(gene.list)), sep="_")
    
    require(pheatmap)
    require(RColorBrewer)
    breaksList = seq(0, 20, by = 1)
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(250)
    
    if (verbose==TRUE) {
            pheatmap(exp.sel,
                     cluster_cols = F, annotation_col = annot, annotation_row = annot_row,
                     cluster_rows = F, show_colnames = F, show_rownames = T, color = color)
    } else if (slim==FALSE) {
        pheatmap(exp.sel,
                 cluster_cols = F, annotation_col = annot, annotation_row = annot_row,
                 cluster_rows = F, show_colnames = F, show_rownames = F, color = color)
        } else {
        pheatmap(exp.sel,
                 cluster_cols = F, cluster_rows = F, show_colnames = F, show_rownames = F, 
                 legend=F, annotation_legend=F, annotation_names_col=F, annotation_names_row=F)
    }
}



make_dyno = function(sce, prop = 1L){
    ncells = ncol(sce)
    if(prop == 1){
        
    } else {
        sce = sce[,sample(ncells %>% seq_len, prop*ncells, replace = FALSE)]
    }
    wrap_expression(
        counts = counts(sce) %>% t,
        expression = logcounts(sce) %>% t,
        cell_info = colData(sce) %>% as.data.frame %>% tibble::rownames_to_column("cell_id"),
        group_info = colData(sce) %>% as.data.frame,
        feature_info = rowData(sce) %>% as.data.frame %>% tibble::rownames_to_column("feature_id")
    )
}

calcStability <- function(exprsMat, celltype, score=c("cv", "zprop")) {
    
    cty <- droplevels(factor(celltype))
    res <- lapply(levels(cty), function(type) {
        mat <- exprsMat[, celltype %in% type]
        cv <- apply(mat, 1, function(x) sd(x)/mean(x))
        zprop <- apply(mat, 1, function(x) length(x[x==0])/length(x))
        zprop[zprop == 0] <- 1/length(zprop) #add this line to replace any 0s
        return(list(cv=cv, zprop=zprop))
    })
    names(res) <- levels(cty)
    
    if (score=="cv") {
        cv_score <- lapply(levels(cty), function(anchor_type) {
            other_type <- levels(cty)[!levels(cty) %in% anchor_type]
            int_score <- lapply(other_type, function(other) {
                res[[anchor_type]]$cv/res[[other]]$cv
            })
            int_score <- do.call(cbind, int_score)
            stability_score <- apply(int_score, 1, function(x) 1 - log(mean(x)))
            stability_score <- sort(stability_score, decreasing=F)
            return(stability_score)
        })
        names(cv_score) <- levels(cty)
        score <- cv_score
    } else if (score=="zprop") {
        zprop_score <- lapply(levels(cty), function(anchor_type) {
            other_type <- levels(cty)[!levels(cty) %in% anchor_type]
            int_score <- lapply(other_type, function(other) {
                res[[anchor_type]]$zprop/res[[other]]$zprop
            })
            int_score <- do.call(cbind, int_score)
            stability_score <- apply(int_score, 1, function(x) 1 - log(mean(x)))
            stability_score <- sort(stability_score, decreasing=F)
            return(stability_score)
        })
        names(zprop_score) <- levels(cty)
        score <- zprop_score
    }
    
    return(score)
}

Cepo_components <- function(exprsMat, cellTypes, score=c("cv", "zprop")) {
    
    cts <- names(table(cellTypes))
    
    segIdx.list <- list()
    for(i in 1:length(cts)){
        idx <- which(cellTypes == cts[i])
        segIdx.list[[i]] <- segIndex(exprsMat[,idx], score=score)
    }
    names(segIdx.list) <- cts
    
    segMat <- segIdxList2Mat(segIdx.list)
    segGenes <- consensusSegIdx(segMat)
    names(segGenes) <- colnames(segMat)
    
    return(segGenes)
}

segIndex <- function(mat, score){
    # nz <- (rowSums(mat != 0) / ncol(mat))
    # ms <- apply(mat, 1, function(x){mean(x)})
    # cvs <- apply(mat, 1, function(x){
    #     return(sd(x)/mean(x))
    # })
    nz <- Matrix::rowMeans(mat != 0)
    ms <- Matrix::rowMeans(mat)
    sds <- matrixStats::rowSds(mat)
    cvs <- sds/ms
    names(nz) = names(sds) = names(cvs) = names(ms) = rownames(mat)
    
    x1 <- rank(nz)/(length(nz)+1)
    x2 <- 1 - rank(cvs)/(length(cvs)+1)
    
    if (score == "cv") {segIdx <- x2} 
    if (score == "zprop") {segIdx <- x1} 
    
    return(segIdx)
}

segIdxList2Mat <- function(segIdx.list) {
    allGenes <- unique(unlist(lapply(segIdx.list, names)))
    segMat <- matrix(0, nrow=length(allGenes), ncol=length(segIdx.list))
    rownames(segMat) <- allGenes
    colnames(segMat) <- names(segIdx.list)
    
    for(i in 1:length(segIdx.list)) {
        si <- segIdx.list[[i]]
        segMat[names(si),i] <- si
    }
    return(segMat)
}

consensusSegIdx <- function(mat) {
    tt <- mat
    
    CIGs <- list()
    for (i in 1:ncol(tt)) {
        avgRank <- c()
        for(j in 1:ncol(tt)) {
            if(i == j){next}
            avgRank <- cbind(avgRank, tt[,i] - tt[,j])
        }
        
        CIGs[[i]] <- sort(rowMeans(avgRank), decreasing = TRUE)
    }
    return(CIGs)
}

#adapted from http://genoweb.toulouse.inra.fr/~pmartin/pgpmartin/2018/11/14/nicer-scatterplot-in-gggally/
GGscatterPlot <- function(data, mapping, color_vec, ..., 
                          method = "spearman") {
    
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)
    cor <- cor(x, y, method = method)
    
    #Assemble data frame
    df <- data.frame(x = x, y = y)
    df$cols <- color_vec
    
    #Get 2D density for alpha
    dens2D <- MASS::kde2d(df$x, df$y)
    df$density <- fields::interp.surface(dens2D , 
                                         df[,c("x", "y")])
    
    if (any(df$density==0)) {
        mini2D = min(df$density[df$density!=0]) #smallest non zero value
        df$density[df$density==0] <- mini2D
    }
    
    #Prepare plot
    my_pal = colorRampPalette(viridis::magma(n=10)[9:3])(35)
    
    pp <- ggplot(df, aes(x=x, y=y, color = cols, alpha = 1/density)) +
        #ggrastr::geom_point_rast(size=5, alpha=0.6) + 
        ggplot2::geom_point(shape=16, size=0.8, show.legend = T) +
        scale_color_gradientn(colors=my_pal) +
        ggplot2::scale_alpha(range = c(.05, .6)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
        #ggplot2::geom_label(
           #data = data.frame(
               #xlabel = min(x, na.rm = TRUE),
               #ylabel = max(y, na.rm = TRUE),
               #lab = round(cor, digits = 3)),
           #mapping = ggplot2::aes(x = xlabel, 
                                  #y = ylabel, 
                                  #label = lab),
           #hjust = 0, vjust = 1,
           #size = 3, fontface = "bold",
           #inherit.aes = FALSE) +
    theme_classic()
    
    return(pp)
}


getStats <- function(res, method=c("Limma","Voom","MAST","ttest","EdgeR","Wilcoxon","DD","DP")) {
    
    if (method=="DP") {
        
        result <- lapply(res, function(x) {
            
            stats <- x$nonzero.pvalue.adj
            names(stats) <- rownames(x)
            x$dir[is.na(x$dir)] <- FALSE
            
            stats <- sapply(1:nrow(x), function(y) {
                if(x$dir[y]==T) {stats[y]} else {1-stats[y]}
            })
            
            stats <- -log10(stats)
            stats[is.na(stats)] <- 0
            stats[is.infinite(stats)] <- 0
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
            
        })
        return(result)
    }
    
    if (method=="DD") {
        
        result <- lapply(res, function(x) {
            #stats <- -log10(x$adj.pvalue)
            stats <- x$`stats.D^+`
            names(stats) <- rownames(x)
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method %in% c("Limma", "Voom")) {
        
        result <- lapply(res, function(x) {
            stats <- x$t
            names(stats) <- rownames(x)
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method=="EdgeR") {
        
        result <- lapply(res, function(x) {
            stats <- x$F
            names(stats) <- rownames(x)
            stats <- sapply(1:nrow(x), function(y) {
                if(x$logFC[y]>0) {stats[y]} else {-stats[y]}
            })
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method=="ttest") {
        
        result <- lapply(res, function(x) {
            stats <- x$stats.t
            names(stats) <- rownames(x)
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method=="Wilcoxon") {
        
        result <- lapply(res, function(x) {
            stats <- -x$stats.W
            names(stats) <- rownames(x)
            stats <- sort(stats, decreasing=T, na.last=T)
            return(stats)
        })
        return(result)
    }
    
    if (method=="MAST") {
        
        result <- lapply(res, function(x) {
            mast <- x$pval
            names(mast) <- rownames(x)
            x$logFC[is.na(x$logFC)] <- 0
            
            stats <- sapply(1:nrow(x), function(y) {
                if(x$logFC[y]>0) {mast[y]} else {1-mast[y]}
            })
            stats <- -log10(stats)
            stats[is.na(stats)] <- 0
            stats[is.infinite(stats)] <- 0
            
            #stats <- x$lambda
            #names(stats) <- rownames(x)
            #x$logFC[is.na(x$logFC)] <- 0
            #stats <- sapply(1:nrow(x), function(y) {
            #  if(x$logFC[y]>0) {stats[y]} else {-stats[y]}
            #})
            stats <- sort(stats, decreasing=T, na.last=T)
            
            return(stats)
        })
        return(result)
    }
    
}

FindDivergentGenes <- function(anchor.list, other.list, n_anchor=30, n_other=50) {
    
    # n_anchor: increase number to find more genes
    # n_other: decrease number to find more genes
    # anchor.list <- list(DS=DS.res)
    # other.list  <- list(DE=DE.res,
    #                    DP=DP.res,
    #                    BD=BD.res)
    #cv_score_rank <- lapply(cv_score, function(x) rank(x)) 
    
    other.list <- lapply(other.list, function(x) {
        res <- lapply(x, function(y) names(y)[1:n_other])
        return(res)
    })
    other.list <- lapply(1:length(other.list[[1]]), function(x) {
        unique(unlist(lapply(other.list, function(y) y[[x]])))
    })
    anchor.list <- lapply(anchor.list, function(x) {
        res <- lapply(x, function(y) names(y)[1:n_anchor])
        return(res)
    })
    anchor.unique.list <- lapply(1:length(other.list), function(x) {
        setdiff(anchor.list[["Cepo"]][[x]], other.list[[x]])
        
    })
    names(anchor.unique.list) <- names(anchor.list[["Cepo"]])
    return(anchor.unique.list)
    
}