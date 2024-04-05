library(ape)
library(phangorn)
library(adephylo)
library(pbmcapply)

path_measure<-function(t){
  opt<-lapply(c("patristic","nNodes"), function(e) distRoot(t, method=e))
  names(opt)<- c("path_length","no_nodes")
  opt<-as.data.frame(opt)
  opt$cell_id<-row.names(opt)
  return(opt)
}
tumor_id_extract<-function(l){
  return(l[3])
}
###############################################################################################
myargs = commandArgs(trailingOnly=TRUE)

no_cores=50
tree_id<-myargs[1]

tree_file<- paste0("./no_root_l1000_lineagegrp_",tree_id,".nexus.trees")
no_root_tree<-  read.nexus(tree_file)
no_root_tree_summary_stats<-pbmclapply(no_root_tree,path_measure,mc.cores=no_cores,ignore.interactive = T)
no_root_tree_summary_stats<-data.table::rbindlist(no_root_tree_summary_stats, idcol = TRUE)
no_root_tree_summary_stats$tumor_id<-lapply(strsplit(no_root_tree_summary_stats$cell_id, "_"),tumor_id_extract)
no_root_tree_summary_stats$tumor_id<-as.character(no_root_tree_summary_stats$tumor_id)
node_bl_output_file <- paste0("path_length_vs_no_nodes_lineagegrp_",tree_id,".txt")
write.table(noquote(no_root_tree_summary_stats),node_bl_output_file,row.names = FALSE,sep="\t",quote = FALSE)

no_root_tree_summary_stats<-subset(no_root_tree_summary_stats, select = -c(.id))
no_root_tree_summary_stats_median<-aggregate(no_root_tree_summary_stats[, 1:2], list(no_root_tree_summary_stats$cell_id), median)
no_root_tree_summary_stats_median$tumor_id<-lapply(strsplit(no_root_tree_summary_stats_median[,1], "_"),tumor_id_extract)
no_root_tree_summary_stats_median$tumor_id<-as.character(no_root_tree_summary_stats_median$tumor_id)
no_root_tree_summary_stats_mean<-aggregate(no_root_tree_summary_stats[, 1:2], list(no_root_tree_summary_stats$cell_id), mean)
no_root_tree_summary_stats_mean$tumor_id<-lapply(strsplit(no_root_tree_summary_stats_mean[,1], "_"),tumor_id_extract)
no_root_tree_summary_stats_mean$tumor_id<-as.character(no_root_tree_summary_stats_mean$tumor_id)

node_bl_output_file <- paste0("median_path_length_vs_no_nodes_lineagegrp_",tree_id,".txt")
write.table(noquote(no_root_tree_summary_stats_median),node_bl_output_file,row.names = FALSE,sep="\t",quote = FALSE)
node_bl_output_file <- paste0("mean_path_length_vs_no_nodes_lineagegrp_",tree_id,".txt")
write.table(noquote(no_root_tree_summary_stats_mean),node_bl_output_file,row.names = FALSE,sep="\t",quote = FALSE)

