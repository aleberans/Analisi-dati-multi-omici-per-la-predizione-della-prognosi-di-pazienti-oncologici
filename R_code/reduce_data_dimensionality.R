reduce_data_dimensionality <- function(mat_data = NULL, str_desc = 'data', data_path = NULL, dim_red = dim_red, 
                         minD = 100, best_red = NULL, fagg_twonn = 'median',
                         dist_fun_twoNN = 'euclidean', ID_col = NULL, str_par = NULL,
                         perc_points = 0.95){
  library(rsvd)  
  
  #minD = min(minD, nrow(mat_data))
  if ((!is.null(dim_red)) & (ncol(mat_data)>minD)){
    
      ID_twonn = estimate_ID_twoNN(mat_data, fagg_twonn = fagg_twonn, dist_fun_twoNN = dist_fun_twoNN, 
                                   maxit = 51, perc_points = perc_points)
      cat('ID of input data = ', ID_twonn, '\n')
      
      if (is.null(best_red)){
        cat('dimensionality reduction - choosing between ', names(dim_red), '\n ')
      }else{
        cat('dimensionality reduction with ', best_red, '- evaluationg also ', 
            names(dim_red)[!(names(dim_red) %in% best_red)], '\n ')
        
      }
      red_data = list()
      IDs = NULL

      for (algo in names(dim_red)){
        cat('*****', algo, '***** \n')
        
        red_data[[algo]] = do.call(dim_red[[algo]], list(x= mat_data, k = minD))
      
        IDs = c(IDs, estimate_ID_twoNN(red_data[[algo]], fagg_twonn = fagg_twonn, dist_fun_twoNN = dist_fun_twoNN, 
                                       maxit = 51, perc_points = perc_points, verbose = 0))
      
        if (!is.null(data_path)){
          
          do_scatter(mat_data = red_data[[algo]], df_pt = df_pt, ID_col_in_pt = ID_col[['pt']], outcome_col = outcome_col, 
                     str_desc = str_desc, task ='red', save_path = save_path, str_par = str_par)
            
        }      
      }    
      names(IDs) = names(red_data)
      
      cat('estimated IDs = \n')
      print(IDs)
      
      if (is.null(best_red)){
        dists = abs(IDs-ID_twonn)
          
        gt_0 = which(dists>0)
        if (length(gt_0)>0){
          dists = dists[gt_0]
          IDs = IDs[names(dists)]
          best_red = names(dists)[which(dists == min(dists))]
          
          
        }else{
          dists = abs(dists)
          best_red = names(IDs)[which(IDs == max(IDs))]

        }
      }
      
      projected_data = red_data[[best_red]]
      ID_twonn = IDs[best_red]
      cat(' choosing', best_red, ' for last dimensionality reduction\n')
      cat('estimated ids on ', str_desc, ' after last dim red  = ', ID_twonn, '\n')
   }else{
     projected_data = mat_data
   }
   
 #  
   return(projected_data) 
}