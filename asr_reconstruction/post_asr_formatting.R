library(ape)
library(pbmcapply)
library(tidytree)
options(warn=-1)

resolve_polytomies<-function(t){
  return(di2multi(t,0.000001))
}

label_nodes<-function(t){
  return(makeNodeLabel(t,"md5sum"))
}

most_likely_state_ML<-function(df,no_states){
  df[]<-pbmclapply(df,as.numeric,mc.cores=no_cores,ignore.interactive = T)
  no_trees<-ncol(df)
  no_nodes<-(nrow(df)/no_states)
  df$no_trees<-(ncol(df)-rowSums(is.na(df)))
  df$average<-matrixStats::rowMedians(as.matrix(df[,1:no_trees]), na.rm=T)
  mls<-c()
  for (i in seq(1,nrow(df),no_states)){
    print(i/nrow(df))
    comp_opt<-which.max(df[i:(i+(no_states-1)),ncol(df)])
    mls<-append(mls,(i+(comp_opt-1)))
  }
  return(df[mls,])
}

rec_node_format<-function(l){
  return(cbind(l[1],gsub("[\\(\\)]", "", regmatches(l, gregexpr("\\(.*?\\)", l)))[2]))
}

string_extraction<-function(df,no_states){
  tmp<-strsplit(row.names(df), " ",fixed=TRUE)
  opt<-as.data.frame(t(sapply(tmp,rec_node_format)))
  colnames(opt)<-c("rec_node_id","anc_state")
  df<-cbind(df,opt)
  if (any(df$average <= (1 - (1/no_states))) == TRUE){
    df[df$average <= (1 - (1/no_states)),]$anc_state<-NA
  }
  return(df)
}

asr_formatting<-function(df,no_states){
  tmp<-most_likely_state_ML(df,no_states)
  tmp<-string_extraction(tmp,no_states)
  tmp<-tmp[-c(1),]
  tmp$num_ord<-1:nrow(tmp)
  return(tmp)
}

list_formatting<-function(l){
  v1<-l[1]
  v2<-strsplit(v1[[1]],"_",fixed=TRUE)[[1]][3]
  return(v2)
}

update_internal_nodes<-function(tree,df){
  tree$node.label<-df$anc_state[match(interaction(tree$node.label),interaction(df$md5_node))]
  return(tree)
}

tips_per_node<-function(tree){
  no_tips<-length(tree$tip.label)
  node_list<-c((no_tips+1):(no_tips+tree$Nnode))
  holder=array(data=NA,c(0,no_tips+1))
  for ( i in 1:length(node_list)){
    desc_tips<-t(na.omit(tree$tip.label[phytools:::getDescendants(tree, node_list[i], curr=NULL)]))
    if (length(desc_tips) < no_tips){
      desc_tips<-cbind(desc_tips,t(rep("NA",(no_tips - length(desc_tips)))))
    }
    c<-cbind(tree$node.label[i],desc_tips)
    holder = rbind(holder,c)
  }
  return(as.data.frame(holder))
}

transition_branch_formating<-function(labeled_data,asr_tree,labeled_tree,tip){
  asr_data<-as_tibble(asr_tree)
  a<-rbind(ancestor(labeled_data,tip),labeled_data[labeled_data$label == tip,])
  a$order_column<-1:nrow(a)
  b<-rbind(ancestor(asr_data,tip),asr_data[asr_data$label == tip,])
  m<-merge(as.data.frame(a),as.data.frame(b),by.x=c("parent","node","branch.length"),by.y=c("parent","node","branch.length"))
  m<-m[order(m$order_column),]
  m<-m[ , !(names(m) %in% c("order_column"))]
  colnames(m)[c(2,4,6)]<-c("descendant","md5sum_descendant","location_descendant")
  m[grepl("_",m$location_descendant),]$location_descendant<-gsub('[[:digit:]]+', '',strsplit(m[grepl("_",m$location_descendant),]$location_descendant,"_",fixed=TRUE)[[1]][3])
  m[m=="NA"]<-NA  
  m$location_parent<-NA
  for (i in 1:nrow(m)){
    if (nrow(m[m$md5sum_parent == m$md5sum_descendant[i],]) != 0){
      m[m$md5sum_parent == m$md5sum_descendant[i],]$location_parent<-m$location_descendant[i]
    }
  }
  if (nrow(m[(m$location_descendant != "T") & (!is.na(m$location_descendant)),]) > 0){
    m[(m$location_descendant != "T") & (!is.na(m$location_descendant)),]$location_descendant<-"M"
  }
  m$transition_branch<-m$location_descendant != m$location_parent
  return(m)
}


tree_to_df<-function(idx){
  tree<-asr_tree[idx][[1]]
  labeled_tree<-label_nodes(tree)
  df<-as_tibble(labeled_tree)
  df$md5sum_parent<-NA
  for(i in 1:nrow(df)){
    df[df$parent == df$node[i],]$md5sum_parent<-df$label[i]
  }
  tips<-labeled_tree$tip.label
  for (j in 1:length(tips)){
    if (j == 1){
      transition_data<-transition_branch_formating(df,tree,labeled_tree,tips[j])
    }
    if (j != 1){
      transition_data<-rbind(transition_data,transition_branch_formating(df,tree,labeled_tree,tips[j]))
    }
  }
  pl_data<-partioned_path_length(transition_data)
  pl_data$tree_id<-names(asr_tree)[idx]
  return(pl_data)
}


partioned_path_length<-function(data){
  df<-as.data.frame(data)
  tips<-df$md5sum_descendant[grepl("_",df$md5sum_descendant)]
  tip_idx<-df[df$md5sum_descendant %in% tips,]
  tip_idx$idx<-row.names(tip_idx)
  tip_idx$idx<-as.numeric(tip_idx$idx)
  df$transition_path_length<-NA
  df$constant_path_length<-NA
  df$na_path_length<-NA
  df$no_nodes<-NA
  df$no_transitions<-NA
  df$tumor_id<-NA
  for (i in 1:length(tip_idx$idx)){
    if (i == 1){
      tmp<-df[c(1:tip_idx$idx[i]),]
    }
    if (i != 1){
      tmp<-df[c((tip_idx$idx[(i-1)]+1):tip_idx$idx[i]),]
    }
    df[tip_idx$idx[i],]$transition_path_length<-sum(tmp[!is.na(tmp$transition_branch) & tmp$transition_branch == "TRUE",]$branch.length,na.rm=TRUE)
    df[tip_idx$idx[i],]$constant_path_length<-sum(tmp[!is.na(tmp$transition_branch) & tmp$transition_branch == "FALSE",]$branch.length,na.rm=TRUE)
    df[tip_idx$idx[i],]$na_path_length<-sum(tmp[is.na(tmp$transition_branch),]$branch.length,na.rm=TRUE)
    df[tip_idx$idx[i],]$no_nodes<-(nrow(tmp)-1)
    df[tip_idx$idx[i],]$no_transitions<-nrow(tmp[!is.na(tmp$transition_branch) & tmp$transition_branch == "TRUE",])
    df[tip_idx$idx[i],]$tumor_id<-tail(tmp$location_descendant,1)
  }
  return(df)
}

matched_replacement_met_ids<-function(t){
  tip_ids<-t[[1]]$tip.label
  gt_met_locations<-unique(as.character(lapply(strsplit(tip_ids,".",fixed=TRUE),list_formatting)))
  gt_met_organ<-gsub('[[:digit:]]+', '',gt_met_locations)
  tmp<-as.data.frame(cbind(gt_met_locations,gt_met_organ))
  tmp$alpha_met_location<-"B"
  tmp[tmp$gt_met_organ == "T",]$alpha_met_location<-"A"
  return(tmp)
}


########################################################################################################
myargs = commandArgs(trailingOnly=TRUE)

no_cores=10
tree_id<-myargs[1]
no_root_tree<- read.nexus(paste0("./no_root_l1000_lineagegrp_",tree_id,".nexus.trees"))

matched_location_ids<-matched_replacement_met_ids(no_root_tree)
tip_ids<-no_root_tree[[1]]$tip.label
location_ids<-as.character(lapply(strsplit(tip_ids,".",fixed=TRUE),list_formatting))
matching_data_file<-as.data.frame(cbind(tip_ids,location_ids))
matching_data_file$location_ids<-gsub('[[:digit:]]+', '', matching_data_file$location_ids)
for ( i in 1:nrow(matching_data_file)){
  matching_data_file$location_ids[i]<-matched_location_ids[matched_location_ids$gt_met_organ == matching_data_file$location_ids[i],]$alpha_met_location
}


resolved_tree<-lapply(no_root_tree,resolve_polytomies)
labeled_tree<-pbmclapply(resolved_tree,label_nodes,mc.cores=no_cores,ignore.interactive = T)
##########################################################
df_tips_per_node<-data.table::fread(paste0("./df_tips_per_node_lineagegrp_",tree_id,".txt"))
btraits_formated<-df_tips_per_node[,-c(3)]
btraits_formated[is.na(btraits_formated)]<-""

tag_names<-cbind(rep("AddNode",nrow(btraits_formated)),paste0("RecNode",1:nrow(btraits_formated)),btraits_formated[,2])
colnames(tag_names)<-c("bayes_format","rec_node_id","tag_node_id")
tag_names<-as.data.frame(tag_names)
rm(btraits_formated)
##############################################################################
output<-data.table::fread(paste0("./asr_lineagegrp_",tree_id,"_dataframe.txt"),header=FALSE)
qmat<-t(output[3:4,])
colnames(qmat)<-qmat[1,]
write.table(qmat[-1,],paste0("./qmat_lineagegrp_",tree_id,".txt"),row.names = FALSE,sep="\t",quote = FALSE)

output<-as.data.frame(output[-c(1:5,nrow(output)),])
row.names(output)<-output[,1]
output<-output[,-c(1)]
print(dim(output))
asr_per_node<-asr_formatting(output,length(unique(matching_data_file$location_ids)))

asr_per_node<-merge(asr_per_node, tag_names[, c("rec_node_id", "tag_node_id")], by="rec_node_id")
asr_per_node<-merge(asr_per_node, df_tips_per_node[, c("tag_node_id","md5_node")], by="tag_node_id")
asr_per_node<-asr_per_node[order(asr_per_node$num_ord),]
matched_location_ids[matched_location_ids$gt_met_organ != "T",]$gt_met_organ<-"M"
asr_per_node$anc_state<-matched_location_ids$gt_met_organ[match(interaction(asr_per_node$anc_state), interaction(matched_location_ids$alpha_met_location))]

write.table(asr_per_node,paste0("./asr_lineagegrp_",tree_id,"_results_full.txt"),row.names = FALSE,sep="\t",quote = FALSE)
write.table(asr_per_node[,c(1,2,(ncol(asr_per_node)-4):ncol(asr_per_node))],paste0("./asr_lineagegrp_",tree_id,"_results_filtered.txt"),row.names = FALSE,sep="\t",quote = FALSE)

asr_tree<-pbmclapply(labeled_tree,update_internal_nodes,df=asr_per_node,mc.cores=no_cores,ignore.interactive = T)
write.nexus(asr_tree, file = paste0("asr_no_root_l1000_lineagegrp_",tree_id,".nexus.trees"),translate = TRUE)
##############################################################################
branch_transitions<-data.table::rbindlist(pbmclapply(1:length(asr_tree), tree_to_df,mc.cores=no_cores,ignore.interactive = T))
write.table(branch_transitions,paste0("./lineagegrp_",tree_id,"_branch_transitions.txt"),row.names = FALSE,sep="\t",quote = FALSE)

path_length_transitions<-branch_transitions[grepl("_",branch_transitions$md5sum_descendant),]
path_length_transitions<-path_length_transitions[,c(4,9:15)]
path_length_transitions$total_path_length<-path_length_transitions$transition_path_length + path_length_transitions$constant_path_length + path_length_transitions$na_path_length
path_length_transitions<-path_length_transitions[!((path_length_transitions$tumor_id != "T") & (path_length_transitions$no_transitions == 0)),]
write.table(path_length_transitions,paste0("./lineagegrp_",tree_id,"_path_length_transitions.txt"),row.names = FALSE,sep="\t",quote = FALSE)

median_path_length_transitions<-aggregate(path_length_transitions[, c(2:6,9)], list(path_length_transitions$md5sum_descendant), median)
median_path_length_transitions$tumor_id<-as.character(lapply(median_path_length_transitions$Group.1,list_formatting))
write.table(median_path_length_transitions,paste0("./lineagegrp_",tree_id,"_median_path_length_transitions.txt"),row.names = FALSE,sep="\t",quote = FALSE)




