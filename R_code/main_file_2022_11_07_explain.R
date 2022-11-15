library(dplyr)
library(stringr)
library(caret)
library(intrinsicDimension)
library(intRinsic)
library(impute)
library(tsne)
library(umap)
library(ggplot2)
library(genieclust)

#https://genieclust.gagolewski.com/#


setwd('C:/DATI/Anacleto/Progetti/MS_Finding/')
source('./utilities.R')


dati_rinominati <<- 'dati_elena_Discovery_rinominati'
dir_with_pvals <<- file.path(dati_rinominati, 'association_statistics')
save_path <<- 'results'
if (!dir.exists(save_path)) dir.create(save_path)
data_path <<- file.path(save_path,'data')
if (!dir.exists(data_path)) dir.create(data_path)



file_p <<- list(miRna =  file.path(dir_with_pvals, 'miRna_p.txt'),
                RnaSeq = file.path( dir_with_pvals, 'RnaSeq_p.tsv'),
                SNP =  file.path( dir_with_pvals,'Genotyping_p.txt'),
                methy = file.path( dir_with_pvals, 'methy_p.tsv'))

fn_orig_data <<- list(miRna = file.path( dati_rinominati, 'miRna.txt'),
                      RnaSeq = file.path( dati_rinominati, 'RnaSeq.tsv'),
                      SNP =  file.path( dati_rinominati,'Genotyping.traw'),
                      methy = file.path( dati_rinominati, 'methy.tsv'))


var_with_rownames <<- list(miRna = 'miRna',
                           RnaSeq = 'GeneID',
                           SNP =  'SNP',
                           methy = NULL)

drop_columns <<- list(miRna = NULL,
                      RnaSeq = NULL,
                      SNP =   c("CHR","X.C.M","POS","COUNTED","ALT"),
                      methy = NULL)


dim_red <<- list(
  "umap" = function(x=NULL, k=NULL, method = 'euclidean'){
    cat('umap with ', method, ' distance\n')
    if ((sum(is.na(x)))>0){ 
      x_dist  = compute_distance(x, method = method)
      proj_data = umap(x_dist, n_components =  k,  n_neighbors = 9, input="dist")
    }else{
      proj_data = umap(x, n_components =  k,  n_neighbors = 9)
    }
    return(as.matrix(proj_data$layout))},
  # "tsne" = function(x=NULL, k=NULL){
  #   return(tsne(imputa_con_knn(x), k = k))},
  "rpca" = function(x=NULL, k=NULL){
    s <- rpca(imputa_con_knn(x),k = k,
              center = TRUE,
              scale = TRUE,
              retx = TRUE,
              rand = TRUE)
    return(s$x)
  }
)

open_p_fun <<- list(
  miRna = function(fn){
    df_p = read.csv(file = fn, header = TRUE, sep = '\t')
    return(df_p)
  },
  RnaSeq = function(fn){
    # INNANZI TUTTO CONSIDERA SOLO I RNAseq con pvalue (di ferdinando) sotto una certa soglia
    df_p = read.csv(file = fn, 
                    header = TRUE, sep = '\t')
    row.names(df_p) = df_p$Symbol
    return(df_p)
  },
  SNP = function(fn){
    df_p = read.csv(file = fn, 
                    header = TRUE, sep = '\t')
    df_p[['padj']] = df_p$adj.P.Val
    row.names(df_p) = df_p$SNP
    return(df_p)
  },
  methy= function(fn){
    df_methy_p = read.csv(file = fn, 
                          header = TRUE, sep = '\t')
    df_methy_p[['padj']] = df_methy_p$adj.P.Val
    row.names(df_p) = df_p$ID
    return(df_p)
  }
)


#fracVar <<- 0.5 # se < 1, fracVar ? il è percent di variabili che voglio tenere
#fracVar > 1 è il numero di variabili che voglio tenere, se fracVar = Inf prendo tutto 
fracVar <<- Inf
activity_col <<- 'Activity'
center_col <<- 'CENTER'

name_pt_file <<- 'df_pt.R'

outcome_col <<- paste(center_col, activity_col, sep ='_')



create_pt_df <- function(fn){
  df_pt = read.csv(file= fn, 
                     header = TRUE, sep = '\t')
  
  
  colNames = names(df_pt)
  cols_ID = colNames[grepl('ID', colNames)]
  for (cc_ID in cols_ID){
    df_pt[[cc_ID]] =  gsub('X', '',gsub('-', '_', gsub('\\.','_',df_pt[[cc_ID]])))
  }
  df_pt = df_pt[order(df_pt$Patient_ID), ]
  ID_levels = df_pt$Patient_ID
  df_pt$Patient_ID = factor(df_pt$Patient_ID, levels = ID_levels)
  
  cat('creating new column = ', outcome_col, '\n')
  
  df_pt[[outcome_col]] = paste( df_pt[[center_col]], df_pt[[activity_col]], sep = '_')  
  cat(names(df_pt), '\n')
  rownames(df_pt) =  df_pt$Patient_ID
  save(df_pt, file = file.path(data_path, name_pt_file))
  return(df_pt)
}



select_points_by_name <- function(x, name_list = NULL){
  return(x[match(name_list, rownames(x)), ]) 
} 

open_data_fun <- function(str_desc = NULL, fn = NULL){
  df_data = read.csv(file = fn, 
                     header = TRUE, sep = '\t')
  if (!is.null(var_with_rownames[[str_desc]])){
    rownames(df_data) = df_data[[var_with_rownames[[str_desc]]]]
    df_data = df_data[, !(names(df_data) %in% c(var_with_rownames[[str_desc]],drop_columns[[str_desc]]))]
  }
  return(df_data)
}






embedded_data_types <<- list()
explain_data_types <<- list()


main <- function(data_types = NULL){  
  
  source('./filter_rows_by_p.R')
  source('./clean_data.R')
  source('./estimate_ID_twoNN.R')
  source('./remove_correlated_par.R')
  source('./cluster_features_par.R')
  source('./reduce_data_dimensionality.R')
  source('./clusterize_with_PINS.R')
  fn_pt = file.path(dati_rinominati, 'Pheno_Discovery_CoreCohort_FindingMS.txt')
  df_pt <<- create_pt_df(fn_pt)
  
  if (is.null(data_types)){ 
    data_types = c('miRna', 'RnaSeq', 'SNP') 
  }
  
  normalize_data_fun = list(miRna = min_max_norm, RnaSeq =  min_max_norm, SNP = NULL, methy =  min_max_norm)
  ID_col = list(pt = 'Patient_ID', miRna = 'ID_miRNASeq', 
                RnaSeq = 'ID_RNASeq', 
                SNP = 'ID_Genotyping', 
                methy = 'ID_Methylation')
  
  dim_split_corr = list(miRna = 505, RnaSeq = 1000, SNP = 1000, methy = 1000)
  cutoff = list(miRna = 0.8, RnaSeq = 0.8, SNP = 0.8, methy = 0.8)
  maxit_corr = list(miRna = 3, RnaSeq = 3, SNP = 3, methy = 3)
  method_corr = list(miRna = 'pearson', RnaSeq = 'pearson', SNP = 'pearson', methy = 'pearson')
  maxD = list(miRna = 300, RnaSeq = 10000, SNP = 100000, methy = NA)
  minD = list(miRna = NA, RnaSeq = NA, SNP = NA, methy = NA)
  dim_split_feat_clustering = list(miRna = 500, RnaSeq = 1000, SNP = 1000, methy = 1000)
  maxit_feature_clustering = list(miRna = 3, RnaSeq = 3, SNP = 1, methy = 3)

  
  fagg_twonn = 'median'
  dist_fun_twoNN = 'canberra'
  perc_points_ID_estimate = 0.95
  maxit_id_estimate = 11
  best_red = 'rpca'
  
  
  in_all_df = NULL
  
  embedded_data_types_fn = file.path(data_path, 'embedded_data_types.Rda')
  explain_data_types_fn = file.path(data_path, 'explain_data_types.Rda')
    
  
    for (str_desc in data_types){
      cat('\n\n*****************---------********************\n')
      
      ID_twonn = NULL
      cat('processing ', str_desc, ' data type\n')
      str_par = paste('_', fracVar, sep = '')
      task = 'clean_data'
      file_data_fn = file.path(data_path, paste(str_desc, str_par,  '_', task,'.Rda' , sep ='' ))
      ID_fn = file.path(data_path, paste(str_desc, str_par, '_', task,'_IDtwonn.Rda', sep =''))
      cat('-------------------------------------------------------\n')
      cat(task, '\n')
      if (file.exists(file_data_fn)){
        load(file = file_data_fn)
        load(file = ID_fn)
        cat('loading ', task,  ' data with dim', dim(mat_data), '\n')
        
      }else{
        df_data = open_data_fun(str_desc = str_desc, fn = fn_orig_data[[str_desc]])
        
        if (is.finite(fracVar)){
          fun_p = open_p_fun[[str_desc]]
          df_p_vals = fun_p(file_p[[str_desc]])
          df_data = filter_rows_by_p(df_mat = df_data, df_p_vals = df_p_vals)
        }
        
        mat_data = clean_data(df_data = df_data, str_desc = str_desc, 
                              normalize_f = normalize_data_fun[[str_desc]], 
                              ID_col = ID_col)
        cat('mat_data of ', str_desc, ' after ', task,
            ' = ', dim(mat_data), '\n')
        cat('saving in path ', file_data_fn, '\n')
        save(mat_data , file = file_data_fn)
      
        # ------------- id estimation
        cat('ID after ', task, '\n')
        ID_twonn[task] = estimate_ID_twoNN(mat_data, fagg_twonn = fagg_twonn,  dist_fun_twoNN, 
                                           maxit = maxit_id_estimate, perc_points = perc_points_ID_estimate)
        cat('estimated ids on ', str_desc, ' = ', ID_twonn, '\n')
        save(ID_twonn, file = ID_fn)
        #-----------------------------
        
      }
      
  ###########################
      str_par = paste(str_par, '_', cutoff[[str_desc]], sep = '')
      task = 'correlation_filter'
      file_data_fn = file.path(data_path, paste(str_desc, str_par,  '_', task,'.Rda' , sep ='' ))
      ID_fn = file.path(data_path, paste(str_desc, str_par, '_', task,'_IDtwonn.Rda', sep =''))
      cat('-------------------------------------------------------\n')
      cat(task, '\n')
      if (file.exists(file_data_fn)){
        load(file = file_data_fn)
        load(file = ID_fn)
        cat('loading ', task,  ' data with dim', dim(mat_data), '\n')
        
      }else{
        
        cat(task,  ' parameter ', cutoff[[str_desc]],  '\n')
        
        mat_data = remove_correlated_par(mat_data, dim_split = dim_split_corr[[str_desc]], 
                                         maxiter = maxit_corr[[str_desc]], 
                                         method = method_corr[[str_desc]], cutoff = cutoff[[str_desc]])
        cat('mat_data of ', str_desc, ' after ', task,
            ' = ', dim(mat_data), '\n')
        cat('saving in path ', file_data_fn, '\n')
        save(mat_data , file = file_data_fn)
  
        # ------------- id estimation
        cat('ID after ', task, '\n')
        ID_twonn[task] = estimate_ID_twoNN(mat_data, fagg_twonn = fagg_twonn,  dist_fun_twoNN, 
                                           maxit = maxit_id_estimate, perc_points = perc_points_ID_estimate)
        cat('estimated ids on ', str_desc, ' = ', ID_twonn, '\n')
        save(ID_twonn, file = ID_fn)
        #-----------------------------
        
        
      }
      
  ############################
      if (is.na(maxD[[str_desc]])){
        maxD[[str_desc]] = round(mean(ID_twonn))*10
        str_par = paste(str_par, '_IDx10', sep = '')
      }else{
        str_par = paste(str_par, '_', maxD[[str_desc]], sep = '')
      }
      task = 'feature_clustering'
      file_data_fn = file.path(data_path, paste(str_desc, str_par,  '_', task,'.Rda' , sep ='' ))
      ID_fn = file.path(data_path, paste(str_desc, str_par, '_', task,'_IDtwonn.Rda', sep =''))
      cat('-------------------------------------------------------\n')
      cat(task, '\n')
      if (file.exists(file_data_fn)){
        load(file = file_data_fn)
        load(file = ID_fn)
        cat('loading ', task,  ' data with dim', dim(mat_data), '\n')
        
      }else{
        
        cat(task,  'until dimension is greater than ', maxD[[str_desc]],  '\n')
        
        mat_data = cluster_features_par(mat_data, dim_split = dim_split_feat_clustering[[str_desc]], 
                                              maxD = maxD[[str_desc]], maxit = maxit_feature_clustering[[str_desc]],
                                              perturbation_based_clustering = FALSE, method = method_corr[[str_desc]])
        
        cat('mat_data of ', str_desc, ' after ', task,
            ' = ', dim(mat_data), '\n')
        cat('saving in path ', file_data_fn, '\n')
        save(mat_data , file = file_data_fn)
  
        # ------------- id estimation
        cat('ID after ', task, '\n')
        ID_twonn[task] = estimate_ID_twoNN(mat_data, fagg_twonn = fagg_twonn,  dist_fun_twoNN, 
                                           maxit = maxit_id_estimate, perc_points = perc_points_ID_estimate)
        cat('estimated ids on ', str_desc, ' = ', ID_twonn, '\n')
        save(ID_twonn, file = ID_fn)
        #-----------------------------
        
      }
      
      explain_data_types[[str_desc]] = mat_data
      
      
  #####################
      if (is.na(minD[[str_desc]])){
        minD[[str_desc]] = round(ID_twonn['correlation_filter']*1.5)
        cat('minD = ', minD[[str_desc]], '\n')
        str_par = paste(str_par, '_ID', '_', best_red, sep = '')
      }else{
        str_par = paste(str_par, '_', minD[[str_desc]], '_', best_red, sep = '')
      }
      best_red = NULL
      task = 'dimensionality_reduction'
      file_data_fn = file.path(data_path, paste(str_desc, str_par,  '_', task,'.Rda' , sep ='' ))
      ID_fn = file.path(data_path, paste(str_desc, str_par, '_', task,'_IDtwonn.Rda', sep =''))
      cat('-------------------------------------------------------\n')
      cat(task, '\n')
      if (file.exists(file_data_fn)){
        load(file = file_data_fn)
        load(file = ID_fn)
        cat('loading ', task,  ' data with dim', dim(mat_data), '\n')
      }else{
        cat(task,  'until dimension is greater than ', minD[[str_desc]],  '\n')
        
        mat_data = reduce_data_dimensionality(mat_data = mat_data, str_desc = str_desc, 
                                                  data_path = data_path,
                                                  dim_red = dim_red,
                                                  minD = minD[[str_desc]], 
                                                  fagg_twonn = fagg_twonn,
                                                  dist_fun_twoNN = dist_fun_twoNN,
                                                  best_red = best_red, 
                                                  ID_col = ID_col, str_par = str_par,
                                                  perc_points = perc_points_ID_estimate)
        cat('mat_data of ', str_desc, ' after ', task,
            ' = ', dim(mat_data), '\n')
        cat('saving in path ', file_data_fn, '\n')
        save(mat_data , file = file_data_fn)
      
        # ------------- id estimation
        cat('ID after ', task, '\n')
        ID_twonn[task] = estimate_ID_twoNN(mat_data, fagg_twonn = fagg_twonn,  dist_fun_twoNN, 
                                           maxit = maxit_id_estimate, perc_points = perc_points_ID_estimate)
        cat('estimated ids on ', str_desc, ' = ', ID_twonn, '\n')
        save(ID_twonn, file = ID_fn)
        #-----------------------------    
        
      }
      
      
      embedded_data_types[[str_desc]] = mat_data
      
  
      if (is.null(in_all_df)){ 
        in_all_df = row.names(mat_data)
      }else{
        in_all_df = intersect(in_all_df, row.names(mat_data))
      }
    }#end for (str_desc in data_types)
    embedded_data_types= lapply(embedded_data_types, select_points_by_name, name_list = in_all_df) 
    explain_data_types = lapply(explain_data_types, select_points_by_name, name_list = in_all_df) 
    
    save(embedded_data_types, file = embedded_data_types_fn)
    save(explain_data_types, file = explain_data_types_fn )
    
    cat('(embedded and clustered) data types saved \n')
    
  
  res_list = clusterize_with_PINS(embedded_data_types, explain_data_types, 
                       data_path = data_path, 
                       data_types = data_types, in_all_df = in_all_df)
  
  clusters = res_list[['clusters']]
  
  explain_data_list = res_list[['explain_data_list']]
  
  CLDs = NULL
  for (str_desc in names(explain_data_list)){ 
     CLDs = rbind(CLDs, explain_clusters(data = explain_data_list[[str_desc]], 
                                         clusters = as.factor(clusters$cluster2), 
                                         str_desc=str_desc))  
  }
   cat('computed CLDs:\n')
   print(CLDs)
   return(CLDs)
}      


