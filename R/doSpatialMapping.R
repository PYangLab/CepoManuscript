df2numericMat <- function(df) {
    mat <- log2(apply(df, 2, as.numeric)+1)
    colnames(mat) <- colnames(df)
    rownames(mat) <- toupper(rownames(df))
    return(mat)
}


consensusRank <- function(mat, top=25) {
    zipCode <- list()
    for (i in 1:ncol(mat)) {
        avgRank <- c()
        for(j in 1:ncol(mat)) {
            if(i == j){next}
            avgRank <- cbind(avgRank, rank(mat[,i] - mat[,j]))
        }
        
        zipCode[[i]] <- names(sort(rowMeans(avgRank), decreasing = TRUE)[1:top])
    }
    return(zipCode)
}

doSpatialMap <- function(singleCells, spatial, spatial.zip, coordinate, k=3) {
    #singleCells <- pijuan_E6.5
    #spatial <- spatial_E6.5
    #spatial.zip <- intersect(DS, rownames(spatial_E6.5))
    #coordinate <- coordinate.E6.5
    
    map <- matrix(0,nrow=ncol(singleCells), ncol=4)
    for(i in 1:ncol(singleCells)) {
        #i=221
        o <- intersect(rownames(singleCells), spatial.zip)
        cs <- (cor(singleCells[o,i], spatial[o,], method = "spearman") + 1)/2
        tmp <- cbind(coordinate[colnames(cs),], cor=as.numeric(cs))
        # top k
        tops <- tmp[order(tmp[,ncol(tmp)], decreasing = TRUE)[1:k],]
        vote <- names(which.max(sapply(split(tops$cor, as.character(tops$labels)), function(x){sum(x)})))
        
        if(length(vote) == 0) {vote="unknown"}
        
        if (vote %in% c("E1", "E2", "E3")) {
            idx <- which(tops$labels == vote)
            sideA <- grepl("A", rownames(tops[idx,]))
            sideP <- grepl("P", rownames(tops[idx,]))
            if (sum(sideA) > sum(sideP)) { 
                tops <- tops[idx,]
                avgCs <- as.numeric(colMeans(matrix(apply(tops[which(sideA),1:3], 2, as.numeric), ncol=3)))
                map[i,]  <- c(avgCs, tops$cor[which(sideA)[1]])
            } else if (sum(sideA) < sum(sideP)) {
                tops <- tops[idx,]
                avgCs <- as.numeric(colMeans(matrix(apply(tops[which(sideP),1:3], 2, as.numeric), ncol=3)))
                map[i,]   <- c(avgCs, tops$cor[which(sideP)[1]])
            } else {
                tops <- tops[idx,]
                idx <- which.max(tops$cor)
                avgCs <- as.numeric(tops[idx,1:3])
                map[i,]   <- c(avgCs, tops$cor[idx[1]])
            } 
        } else {  
            idx <- which(tops$labels == vote)
            avgCs <- as.numeric(colMeans(matrix(apply(tops[idx,1:3], 2, as.numeric), ncol=3)))
            map[i,]  <- c(avgCs, tops$cor[idx[1]])
        }
    }
    # finding the average location
    ys <- apply(map[,1:3], 2, as.numeric)
    cellabs <- c()
    spalabs <- c()
    for (i in 1:nrow(ys)) {
        idx <- which.min(apply(apply(coordinate[,1:3], 2, as.numeric),1,function(x)sqrt(sum((x-ys[i,])^2))))
        cellabs <- c(cellabs, as.character(coordinate[idx,"labels"]))
        spalabs <- c(spalabs, rownames(coordinate[idx,]))
    }
    map.mat <- data.frame(map, cellabs, spalabs)
    colnames(map.mat) <- c("xcoord", "ycoord", "zcoord", "cor", "labels", "spatialLabel")
    return(map.mat)
}


purity <- function(map, j, pm, stage=c("E6.5", "E7.5")){
    
    if (stage == "E6.5") {
        for(i in 1:length(locs)) {
            if (locs[i] %in% names(map)) {
                p <- c(sum(map[[locs[i]]] == "Epiblast"), sum(map[[locs[i]]] == "Primitive Streak"), sum(map[[locs[i]]] == "Visceral endoderm")) / length(map[[locs[i]]])
                pm[i,j] <- paste(paste(p, collapse = " "), length(map[[locs[i]]]), sep=" ")
            }
        }
        return(pm)
    } else {
        for(i in 1:length(locs)) {
            if (locs[i] %in% names(map)) {
                p <- c(sum(map[[locs[i]]] == "Ectoderm"), sum(map[[locs[i]]] == "Endoderm"), sum(map[[locs[i]]] == "Mesoderm")) / length(map[[locs[i]]])
                pm[i,j] <- paste(paste(p, collapse = " "), length(map[[locs[i]]]), sep=" ")
            }
        }
        return(pm)
    }
    
}

