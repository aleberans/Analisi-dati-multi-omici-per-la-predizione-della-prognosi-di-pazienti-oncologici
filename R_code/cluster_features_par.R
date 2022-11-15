library(doParallel)
library(genieclust)
library(PINSPlus)
cluster_features_par <- function(df, dim_split = 1000, 
                             maxD = 250,  M = 10, maxit = 3, method = "pearson", perturbation_based_clustering = FALSE){
  X = t(df)
  cat("CLUSTERING FEATURES\n")
  cat("maxD = ", maxD, "\n")
  
  nc_per_split = min(0.5, (maxD/nrow(X))^(1/maxit))
  
  print(nc_per_split)
  
  iter = 0
  while((nrow(X)>maxD) & (iter<maxit)){
    iter = iter+1
    cat("ITERATION=", iter, "\n")
    cat("dim_split =", dim_split, "nc_per_split=", nc_per_split , "\n")
    
    clustered_X = NULL
    #new_feats = NULL
    
    
    cl <- makePSOCKcluster(4)
    registerDoParallel(cl)
    clustered_X = NULL
    
    clustered_X = foreach(nR = seq(1, nrow(X), by=dim_split), .combine='rbind', .packages = c("genieclust", "PINSPlus")) %dopar% {
      subX = X[nR:min(nrow(X), (nR+dim_split-1)), ]
      if (nrow(subX)>M*2){
        #cat('dimensions of subX = ', dim(subX), '\n')
        k = min(max(round(nrow(subX)*nc_per_split),3),nrow(subX))
        if (k > 1){
          cc = cor(t(subX), use = "pairwise.complete.obs", method = method)
          dist_mat = 1- cor(t(subX), use = "pairwise.complete.obs")
          
          if (!perturbation_based_clustering){
            h <- gclust(dist_mat, M = min(M, round(nrow(subX)/2)))
            res = cutree(h, k)
          }else{
            
            res <- PerturbationClustering(data = data, kMax = k, kMin = k-10, clusteringFunction = function(dist_mat, k){
              # this function must return a vector of cluster
              h <- gclust(dist_mat, M = min(M, round(nrow(subX)/2)))
              res = cutree(h, k)
              return(res)            
            })
          
          }
          clusters_centroids = NULL
          rnames = NULL
          for (nc in min(res):max(res)){
            idx_cluster = which(res == nc)
  
  
            #new_feats[nc] = paste(features_in_cluster[idx_cluster], collapse = ' + ')
            if (length(idx_cluster)>1){
              
              centroidval = colMeans(subX[idx_cluster, ], na.rm = TRUE)
              xxx = rbind(centroidval,subX[idx_cluster, ])
              
              #prendi solo la prima riga per avere tutte e sole le correlazioni con il centroide
              xcor = cor(t(xxx), use = "pairwise.complete.obs")[1, ]
              
              #elimina la correlazione del centroide con se stesso
              idx_best_feat = which.max(abs(xcor[-1]))
              idx_cluster = idx_cluster[idx_best_feat]
              
              rnames = c(rnames, idx_cluster)
            }else{
              # TENGO O BUTTO?? per ora tengo
              #if (keep_isolated){
              rnames = c(rnames, idx_cluster)
              #}
            }
          }
  
          print('********************-------- *************************')
          print(rnames)
          clusters_centroids = subX[rnames, ]
          rownames(clusters_centroids) = names(rnames)
  
        }else{
          #cat('no clusters = ', k, 'too little; selecting all elements\n')
          clusters_centroids = subX
        }
        return(clusters_centroids)
      }else{
        return(subX)
      }
    }

    
    
    X = clustered_X[sample(nrow(clustered_X)), ]
    clustered_X = NULL
    
    cat("dim di X", dim(X), "\n")
    stopCluster(cl)
    
  }
  # cat('Too many features, (' , ncol(projected_data), ') I will cluster them!\n')
  # h <- gclust(emst_mlpack(X), distance = 'cosine')
  # #h <- gclust(X, distance = 'cosine')
  # 
  
  # 
  # cat('with clustering I have ', ncol(projected_data),  ' clustered features\n') 
   
   X = t(X)
   cat('final dimension = ', dim(X), '\n')
   
  
   return(X)
}