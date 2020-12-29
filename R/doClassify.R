calculateSimilarity <- function(exprsMat_train,
                                exprsMat_test,
                                similarity = c("pearson",  "spearman",
                                               "cosine", "jaccard", "kendall",
                                               "weighted_rank","manhattan")) {
    
    
    similarity <- match.arg(similarity, c("pearson",  "spearman",
                                          "cosine", "jaccard", "kendall",
                                          "weighted_rank","manhattan"))
    
    if ("dgCMatrix"  %in% is(exprsMat_test) &
        !"dgCMatrix" %in% is(exprsMat_train)) {
        exprsMat_train <- methods::as(exprsMat_train, "dgCMatrix")
    }
    
    if (!"dgCMatrix" %in% is(exprsMat_test) &
        "dgCMatrix" %in% is(exprsMat_train)) {
        exprsMat_test <- methods::as(exprsMat_test, "dgCMatrix")
    }
    
    
    if (similarity == "cosine") {
        
        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- proxyC::simil(Matrix::t(exprsMat_train),
                                    Matrix::t(exprsMat_test),
                                    method = "cosine")
        } else {
            corMat <- 1 - as.matrix(proxy::dist(t(as.matrix(exprsMat_train)),
                                                t(as.matrix(exprsMat_test)),
                                                method = "cosine"))
        }
        
        corMat[is.na(corMat) | is.infinite(corMat)] <- min(corMat)
        
    } else if (similarity == "kendall") {
        
        
        
        corMat <- stats::cor(as.matrix(exprsMat_train),
                             as.matrix(exprsMat_test), method = "kendall")
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1
        
    } else if (similarity == "jaccard") {
        
        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- proxyC::simil(Matrix::t(exprsMat_train),
                                    Matrix::t(exprsMat_test),
                                    method = "jaccard")
        } else {
            corMat <- 1 -
                as.matrix(proxy::dist(t(as.matrix(exprsMat_train > 0)),
                                      t(as.matrix(exprsMat_test > 0)),
                                      method = "Jaccard"))
            
        }
        corMat[is.na(corMat) | is.infinite(corMat)] <- min(corMat)
    }else if (similarity == "weighted_rank") {
        
        corMat <- wtd_rank2(as.matrix(exprsMat_train),
                            as.matrix(exprsMat_test), method = "pearson")
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1
        
    }else if (similarity == "manhattan") {
        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- 1 - as.matrix(proxy::dist(t(as.matrix(exprsMat_train)),
                                                t(as.matrix(exprsMat_test)),
                                                method = "Manhattan"))
        } else {
            corMat <- 1 - as.matrix(proxy::dist(t(as.matrix(exprsMat_train)),
                                                t(as.matrix(exprsMat_test)),
                                                method = "Manhattan"))
        }
        corMat[is.na(corMat) | is.infinite(corMat)] <- min(corMat)
        
    }else if (similarity == "spearman") {
        
        corMat <- stats::cor(as.matrix(exprsMat_train),
                             as.matrix(exprsMat_test), method = "spearman")
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1
        
    }else if (similarity == "pearson") {
        
        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- proxyC::simil(Matrix::t(exprsMat_train),
                                    Matrix::t(exprsMat_test),
                                    method = "correlation")
        } else {
            corMat <- stats::cor(as.matrix(exprsMat_train),
                                 as.matrix(exprsMat_test), method = "pearson")
        }
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1
        
    }else{
        
        if ("dgCMatrix" %in% is(exprsMat_test) &
            "dgCMatrix" %in% is(exprsMat_train)) {
            corMat <- proxyC::simil(Matrix::t(exprsMat_train),
                                    Matrix::t(exprsMat_test),
                                    method = "correlation")
        } else {
            corMat <- stats::cor(as.matrix(exprsMat_train),
                                 as.matrix(exprsMat_test), method = "pearson")
        }
        corMat[is.na(corMat) | is.infinite(corMat)] <- -1
    }
    
    return(corMat)
}


kNN <- function(exprsMat_train, exprsMat_test, trainLabels, k=10, similarity = "pearson") {
    corMat <- calculateSimilarity(exprsMat_train, exprsMat_test, similarity = similarity)
    corMat[is.na(corMat) | is.infinite(corMat)] <- -1
    topKNN <- apply(corMat, 2, function(x) trainLabels[order(x, decreasing = TRUE)][seq_len(k)])
    predRes <- apply(topKNN, 2, function(x) {
        tab <- table(x)/length(x)
        names(tab)[which(tab == max(tab, na.rm = TRUE))][1]})
    return(predRes)
}

mIntersect <- function (x, y, ...){
    if (missing(...)) 
        intersect(x, y)
    else intersect(x, mIntersect(y, ...))
}

avgAcc <- function(pred, testLabels) {
    cls <- unique(pred)
    accs <- c()
    for(i in 1:length(cls)) {
        idx <- which(testLabels == cls[i])
        accs <- c(accs, sum(testLabels[idx] == pred[idx]) / length(idx))
    }
    return(accs)
}
