clustering=function(unique.eu, diam){
 
  
  #tree.file <- read.table(paste(path.output,tree.name,"_states.tree",sep=""))
  tree.file <- read.table(paste0("Network_",diam,".net_states.tree"),sep="")
  colnames(tree.file) <- c("path", "flow", "name", "state_id", "node_id", "layer_id")
  par <-strsplit(as.character(tree.file[,1]),":")
  par <-lapply(par,function(x)x[-length(x)])
  for(i in 1:max(unlist(lapply(par,length)))){
    tree.file[,ncol(tree.file)+1]<-unlist(lapply(par,function(x)x[i]))
    colnames(tree.file)[ncol(tree.file)]<-paste("Level",i,sep="_")
  }
  
  loc.net=tree.file[which(grepl("_",tree.file[,"name"])),]
  levels=colnames(loc.net)[grepl( "Level_" , names(loc.net))]
  #unique.eu[, levels]=NA
  

  for( i in 1:nrow(unique.eu)){
    for(j in 1:nrow(loc.net)){
      if(unique.eu$tID[i]==loc.net$name[j]){
        #print(c(i,j))
        #unique.eu$level_4[i]=loc.net$Level_4[j]
        #unique.eu$level_3[i]=loc.net$Level_3[j]
        #unique.eu$level_2[i]=loc.net$Level_2[j]
        #unique.eu$level_1[i]=loc.net$Level_1[j]
        unique.eu[i, levels]= loc.net[j,  levels]
        #unique.eu[i, colnames(loc.net)[grepl( "Level_" , names(loc.net))]]= loc.net[j,  colnames(loc.net)[grepl( "Level_" , names(loc.net))]]
      }
    }
  }
  #colnames(unique.eu)
  
  eu.sub=unique.eu[,c("LIDNUM","NAME","LAT","LONG","TIMEUNIT","Genus_Species","level","cid","cID","tID",levels)]
  sub.sub.eu=unique(eu.sub)
  #colnames(sub.sub.eu)
  
  
  cluster_ID=unique.eu[,c("LIDNUM",levels)]
  #print(gg)
  #return(list(map=gg,cluster.ID=cluster_ID, loc.net=loc.net))
  return(list(sub.eu=sub.sub.eu,cluster.ID=cluster_ID))
}
