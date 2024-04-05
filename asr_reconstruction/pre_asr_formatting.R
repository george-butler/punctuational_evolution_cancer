library(ape)
library(pbmcapply)
options(warn=-1)

resolve_polytomies<-function(t){
  return(di2multi(t,0.000001))
}

label_nodes<-function(t){
  return(makeNodeLabel(t,"md5sum"))
}

list_formatting<-function(l){
  v1<-l[1]
  v2<-strsplit(v1[[1]],"_",fixed=TRUE)[[1]][3]
  return(v2)
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
write.table(matching_data_file,paste0("./lineagegrp_",tree_id,"_data_file.txt"),row.names = FALSE,sep="\t",quote = FALSE,col.names=FALSE)


resolved_tree<-lapply(no_root_tree,resolve_polytomies)
write.nexus(resolved_tree, file = paste0("resolved_no_root_l1000_lineagegrp_",tree_id,".nexus.trees"),translate = TRUE)
labeled_tree<-pbmclapply(resolved_tree,label_nodes,mc.cores=no_cores,ignore.interactive = T)

df_tips_per_node<-data.table::rbindlist(pbmclapply(labeled_tree,tips_per_node,mc.cores=no_cores,ignore.interactive = T),fill=TRUE)
df_tips_per_node<-unique(df_tips_per_node,by = "V1")
df_tips_per_node[df_tips_per_node=="NA"]<-NA 
df_tips_per_node<-cbind(rep("AddTag",nrow(df_tips_per_node)),paste0("TagNode",1:nrow(df_tips_per_node)),df_tips_per_node)
colnames(df_tips_per_node)[2:3]<-c("tag_node_id","md5_node")
write.table(df_tips_per_node,paste0("./df_tips_per_node_lineagegrp_",tree_id,".txt"),row.names = FALSE,sep="\t",quote = FALSE)
##########################################################
btraits_formated<-df_tips_per_node[,-c(3)]
btraits_formated[is.na(btraits_formated)]<-""

system("echo 1 >> ./tmp1_asr_input_file.txt")
system("echo 1 >> ./tmp1_asr_input_file.txt")
system("echo ScaleTrees >> ./tmp1_asr_input_file.txt")
system("echo Kappa >> ./tmp1_asr_input_file.txt")

write.table(noquote(btraits_formated),"./tmp2_asr_input_file.txt",row.names = FALSE,sep="\t",quote = FALSE,col.names=FALSE)
fossil_val<-matched_location_ids[grepl("T",matched_location_ids$gt_met_organ),]$alpha_met_location
system(paste0("echo Fossil RecNode1 TagNode1 ", fossil_val, " >> ./tmp2_asr_input_file.txt"))
tag_names<-cbind(rep("AddNode",nrow(btraits_formated)),paste0("RecNode",1:nrow(btraits_formated)),btraits_formated[,2])
colnames(tag_names)<-c("bayes_format","rec_node_id","tag_node_id")
tag_names<-as.data.frame(tag_names)
write.table(tag_names[-c(1),],"./tmp3_asr_input_file.txt",row.names = FALSE,sep="\t",quote = FALSE,col.names=FALSE)
system(paste0("cat ./tmp1_asr_input_file.txt ./tmp2_asr_input_file.txt ./tmp3_asr_input_file.txt > asr_lineagegrp_",tree_id,"_input_file.txt"))
system(paste0("echo run >> asr_lineagegrp_",tree_id,"_input_file.txt"))
system("rm ./tmp*")
