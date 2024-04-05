library(ape)
library(pbmcapply)

resolve_polytomies<-function(t){
  return(di2multi(t,0.000001))
}

rm(list = setdiff(ls(), lsf.str()))

myargs = commandArgs(trailingOnly=TRUE)
id<- myargs[1]

dir.create("./ml_trees", showWarnings = FALSE)
dir.create("./ml_trees_combined", showWarnings = FALSE)

node_pl_file<-read.delim(paste0("./path_length_vs_no_nodes_lineagegrp_",id,".txt"))
no_root_tree_file<-read.nexus(paste0("./no_root_l1000_lineagegrp_",id,".nexus.trees"))
resolved_tree<-lapply(no_root_tree_file,resolve_polytomies)

t_ids<-unique(node_pl_file$.id)
for ( i in 1:length(t_ids)){
  target_tree_id<-t_ids[i]
  target_tree<-resolved_tree[target_tree_id][[1]]
  target_data<-node_pl_file[node_pl_file$.id == target_tree_id,]
  target_data<-subset(target_data, select = -c(.id))
  dir.create(paste0("./ml_trees/tree_",i), showWarnings = FALSE)
  dummy_node_pl_file<-paste0("./ml_trees/tree_",i,"/path_length_vs_no_nodes_lineagegrp_",id,".txt")
  single_tree_file<-paste0("./ml_trees/tree_",i,"/no_root_tree",i,".nexus.trees")
  write.nexus(target_tree,file=single_tree_file)
  tmp_df<-target_data[c("cell_id","path_length","no_nodes","tumor_id")]
  dummy_pl_file<-tmp_df[,-c(ncol(tmp_df))]
  write.table(noquote(dummy_pl_file),dummy_node_pl_file,row.names = FALSE,sep="\t",quote = FALSE)
}

parallel_imp<-function(i){
  dummy_node_pl_file<-paste0("./ml_trees/tree_",i,"/path_length_vs_no_nodes_lineagegrp_",id,".txt")
  single_tree_file<-paste0("./ml_trees/tree_",i,"/no_root_tree",i,".nexus.trees")
  system(paste0("OMP_NUM_THREADS=1 ./BayesTraitsV4 ", single_tree_file," ",dummy_node_pl_file, " < ./regression_BT_run_file.txt > /dev/null"))
  data_start_line<-system(paste0("echo $(cat ./ml_trees/tree_",i,"/*Log.txt | grep -n Lh | cut -d : -f 1)"), intern = TRUE)
  system(paste0("tail -n +",data_start_line," ./ml_trees/tree_",i,"/*Log.txt > ./ml_trees/tree_",i,"/dataframe_tree_",i,"_path_length_vs_no_nodes_lineagegrp_",id,".txt"))
  system(paste0("mv ./ml_trees/tree_",i,"/dataframe_tree_",i,"_path_length_vs_no_nodes_lineagegrp_",id,".txt ./ml_trees_combined"))
  df<-read.delim(paste0("./ml_trees_combined/dataframe_tree_",i,"_path_length_vs_no_nodes_lineagegrp_",id,".txt"))
  df$tree_num<-i
  write.table(noquote(df),paste0("./ml_trees_combined/dataframe_tree_",i,"_path_length_vs_no_nodes_lineagegrp_",id,".txt"),row.names = FALSE,sep="\t",quote = FALSE)
  return(i)
}
no_cores<-40
opt<-pbmclapply(1:length(t_ids),parallel_imp,mc.cores=no_cores,ignore.interactive = T)

system(paste0("nawk 'FNR==1 && NR!=1{next;}{print}' ./ml_trees_combined/*.txt > ./ml_trees_combined/combined_path_length_vs_no_nodes_lineagegrp_",id,".txt"))
system("rm ./ml_trees_combined/dataframe*.txt")
