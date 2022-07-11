randomization<-function(europe){
  source("curveball abundance.R")
  library(dplyr)
  
  europe.reduced=europe[,c("LIDNUM","Genus_Species", "TIMEUNIT")]
  locations.reduced=unique(europe.reduced$LIDNUM)
  locations=unique(europe$LIDNUM)
  
  X<-split(europe.reduced,europe.reduced$TIMEUNIT)
  
  
  occurences=c()
  for(i in 1:length(X)){
    xx=X[[i]]
    #View(xx %>% group_by(Genus_Species))
    xx.g=xx %>% group_by(Genus_Species) %>% summarize(location = paste(sort(LIDNUM),collapse=", "))
    xx.g$location=strsplit(xx.g$location, split = ", ")
    
    incidence=c()
    for(j in 1:nrow(xx.g)){
      res <- table(factor(locations.reduced, levels=locations.reduced)[match(unlist(xx.g$location[j]),locations.reduced, nomatch=0)])
      #print(res[which(res>0)])
      incidence=rbind(incidence,cbind(xx.g$Genus_Species[j],t(res)))
    }
    species=incidence[,1]
    incidence=incidence[,-1]
    incidence=apply(incidence,2, as.numeric)
    rownames(incidence)=species
    col.margins=colSums(incidence)
    species.margins=rowSums(incidence)
    tot=sum(incidence)
    incidence=rbind(incidence,col.margins)
    species.margins=c(species.margins,tot)
    incidence=cbind(incidence,species.margins)
    occurences[[i]]=incidence
    names(occurences)[[i]]=names(X)[[i]]
  }
  
  rand.mat=c()
  
  for(i in 1:length(occurences)){
    cat("\n\n",i)
    m=occurences[[i]]
    
    l_hp=nrow(occurences[[i]])
    m0=curveball_abundance(m,5000)
    mm0=m0$rand.mat
    
    if(identical(m,mm0)==TRUE){
      cat("\n\nThe randomized and orginal matrices are identical \n")
    }
    
    col.margins=colSums(mm0)
    species.margins=rowSums(mm0)
    tot=sum(mm0)
    mm0=rbind(mm0,col.margins)
    species.margins=c(species.margins,tot)
    mm0=cbind(mm0,species.margins)
    rand.mat[[i]]=mm0
    if(any(is.na(mm0))==TRUE){
      break
    }
    column.margins=m["col.margins",]
    row.margins=m[,"species.margins"]
    col.margins=mm0["col.margins",]
    species.margins=mm0[,"species.margins"]
    cat("\n\nAre the randomized and orginal column margins identical?")
    print(table(col.margins==column.margins))
    cat("\n\nAre the randomized and orginal row margins identical?")
    print(table(row.margins==species.margins))
    
    
    if (length(which(row.margins!=species.margins))!=0 | length(which(col.margins!=column.margins))!=0) {
      print(which(row.margins!=species.margins))
      cat("\nOriginal row values",sum(m[which(row.margins!=species.margins)[1],]))
      cat("\nAfter randomization row values",sum(mm0[which(row.margins!=species.margins)[1],]),"\n\n")
      print(which(col.margins!=column.margins))
      cat("\nOriginal column values",sum(m[,which(col.margins!=column.margins)[1]]))
      cat("\nAfter randomization column values",sum(mm0[,which(col.margins!=column.margins)[1]]),"\n\n")
      
      cat("Matrix number ",i,"\n\n")
    }
    
    
  }
  names(rand.mat)=names(occurences)
  return(rand.mat)
}



findings<-function(rand.mat, candidate.loc500,first_occ_no_indet,lab){
  
  un.species=c()
  for(i in 1:length(rand.mat)){
    un.species=c(un.species,rownames(rand.mat[[i]]))
  }
  un.species=unique(un.species[-which(grepl("margins",un.species))])
  un.loc=c()
  for(i in 1:length(rand.mat)){
    un.loc=c(un.loc,colnames(rand.mat[[i]]))
  }
  un.loc=unique(un.loc[-which(grepl("margins",un.loc))])
  
  
  speciesXMN=c()
  
  for(i in 1:length(candidate.loc500)){ #For every location
    
    Mat <- matrix(0, nrow = length(un.species), ncol = length(rand.mat), dimnames = list(un.species, names(rand.mat)))
    for(j in 1:length(rand.mat)){ #For each time unit
      ds=rand.mat[[j]]
      ds=ds[-which(grepl("margins",rownames(ds))),-which(grepl("margins",colnames(ds)))]
      #subset randomized matrix keeping only the location in the radius and sum the occurrences per each species
      sub.rand=as.matrix(ds[,(colnames(ds) %in% candidate.loc500[[i]])])
      if(nrow(sub.rand)>0){ #If there is more than one location
        temp=rowSums(sub.rand)
        which(temp!=0,arr.ind = T)
      }
      #Match it by row names
      if(colnames(Mat)[j] == names(rand.mat)[j]){
        Mat[names(temp), j] <- temp
      }
      
    }
    dim(Mat)
    Mat=as.data.frame(Mat)
    Mat=Mat[rowSums(Mat)>0,]
    #Remove rows with only 0

    
    if(nrow(Mat)>=1){ #If the non-empty rows are more than 1
      #Create a capture history 
      ch=c()
      for(k in 1:nrow(Mat)){
        res=Mat[k,]
        y<-(res>0)*1
        ch[k]=paste(y,sep="",collapse = "") #Create a string of 0 for non appearances and 1 for appearances in a given time unit
      }
      Mat=data.frame(list(Mat,ch=ch),check.names = FALSE)
    }
    speciesXMN[[i]]=Mat
  }
  
  names(speciesXMN)=names(candidate.loc500)
  
  new.sp_rand500=c()
  for(j in 1:length(speciesXMN)){
    ds=as.data.frame(speciesXMN[[j]])
    dim(ds)
    
    nc=which(!grepl("MN",colnames(ds)))
    nc
    if(length(which(grepl("(indet)|(sp)",rownames(ds))))>=1){
      ds=ds[-which(grepl("(indet)|(sp)",rownames(ds))),]
    }
    dim(ds)
    if(nrow(ds)>=1){
      origin=c()
      for(i in 1:nrow(ds)){
        sp=first_occ_no_indet[which(first_occ_no_indet$Genus_Species==rownames(ds)[i]),]
        sp
        min.time=min(which(ds[i,-nc] != 0))#Find earlier time which has more than 0 occurrences
        origin[i]=(sp$Origination.time==lab[min.time])
        if(sp$Origination.time>lab[min.time]){ 
          stop("ERROR: First occurrence is later than the earliest occurrence in this dataset")
        }
      }
      if(length(which(origin==T))<=1){
        new.species=ds[which(origin==T),]
        orig.time=lab[min(which(new.species[,-nc] !=0))]
      }
      if(length(which(origin==T))>1){
        new.species=ds[which(origin==T),]
        orig.time=as.vector(lab[apply(new.species[,-nc] != 0, 1, function(x){min(which(x[-nc] !=0))})]) #Species originated in this location group
      }
      new.species=cbind(new.species,orig.time )
      new.sp_rand500[[j]]=new.species
    } 

    if(nrow(ds)<1){
      print(paste0("No species that is not indet or sp in location radius ",names(speciesXMN)[j]))
      new.sp_rand500[[j]]=paste0("No species that is not indet or sp in location radius ",names(speciesXMN)[j])
    }
    #Now we have to find out the origination time for each  
  }
  names(new.sp_rand500)=names(speciesXMN)
  
  
  #New species per time unit per location
  time.or_rand=c()
  for(i in 1:length(new.sp_rand500)){
    tv=table(new.sp_rand500[[i]]["orig.time"]) #Looks like the assignation of labels gives always Pre MN even when it is not the first occurrence
    tv
    time.or_rand=rbind(time.or_rand,tv[match(lab, names(tv))])
  }
  colnames(time.or_rand)=lab
  rownames(time.or_rand)=names(speciesXMN)
  time.or_rand[is.na(time.or_rand)] <- 0
  
  #Number of species per time unit per location
  time.sp_rand=c()
  for(i in 1:length(speciesXMN)){
    ds=speciesXMN[[i]]
    if(!is.null(nrow(ds))){ 
      tv=colSums(ds[,which(grepl("MN",colnames(ds)))])
      time.sp_rand=rbind(time.sp_rand,tv[match(lab, names(tv))])
    } else {
      
    }
  }
  colnames(time.sp_rand)=lab
  rownames(time.sp_rand)=names(speciesXMN)
  time.sp_rand[is.na(time.sp_rand)] <- 0
  
  #Origination fraction
  time.or.fraction_rand=time.or_rand/time.sp_rand
  time.or.fraction_rand[is.nan(time.or.fraction_rand)] <- 0
  
  return(list(New_species=time.or_rand,Total_species=time.sp_rand,Origination_fraction=time.or.fraction_rand,speciesXMN=speciesXMN))
}


probabilities<-function(speciesXMN,lab){
  #install.packages("RMark")
  library(RMark)

  p.time_rand500=c()
  Phi.dot=list(formula=~1)
  p.dot=list(formula=~1)
  Phi.time=list(formula=~time)
  p.time=list(formula=~time)
  
  for(j in 1:length(speciesXMN)){
    ds=as.data.frame(speciesXMN[[j]])
    ds$Species=rownames(ds)
    ds=ds[,c("Species","ch")]
    ds.process=process.data(ds,model="CJS")
    ds.ddl=make.design.data(ds.process)
    typeof(ds)
    dim(ds)
    if(nrow(ds)>0){ 
      mod=mod=mark(ds.process,ds.ddl,model.parameters=list(Phi=Phi.dot,p=p.time),output = FALSE)
      phi=as.numeric(mod$results$real$estimate[1])
      p=as.numeric(mod$results$real$estimate[-1])
      
      p.time_rand500=rbind(p.time_rand500,c(phi,p))
      cat("Estimate of phi =",phi, "\nEstimate of p = ",p)
      mod=c()
    } else {
      cat("Empty dataset for location ",names(speciesXMN)[j])
    }
    cleanup(lx = NULL,ask=FALSE, prefix = "mark")
  }
  rownames(p.time_rand500)=names(speciesXMN)
  colnames(p.time_rand500)=c("Phi",paste0("p_",lab)[-1])
  return(p.time_rand500)

}


characteristics<-function(find,prob,data){
  time.or.fraction=find$Origination_fraction
  time.prob=prob[,-1]
  time.or.frac=time.or.fraction[,-1]
  time.adj.frac=time.or.frac*(1/time.prob)
  time.adj.frac=as.data.frame(cbind(LIDNUM=rownames(time.adj.frac),time.adj.frac))
  names=rownames(time.adj.frac)
  time.adj.frac[,2:ncol(time.adj.frac)]=sapply(time.adj.frac[2:ncol(time.adj.frac)], as.numeric)
  rownames(time.adj.frac)=names
  bb=data[which(data$LIDNUM %in% time.adj.frac$LIDNUM),]
  bb=bb[!duplicated(bb$LIDNUM),]
  bb=bb[,c("LIDNUM","NAME","LAT","LONG")]
  adj.frac.coord=merge(time.adj.frac,bb,by= "LIDNUM")

  
#Total findings  
  lids=c()
  time.sp=find$Total_species[,-1]
  nfind.sp.fac=c()
  for(i in 2:(ncol(adj.frac.coord)-3)){
    adj.frac.coord[,i]=as.numeric(adj.frac.coord[,i])
    tr=as.numeric(quantile(adj.frac.coord[,i], 0.95))
    idx=which(adj.frac.coord[,i]>tr)
    ab.loc=adj.frac.coord[idx,c(1,i)]
    lids[[i-1]]=adj.frac.coord[idx,c("LIDNUM","NAME")]
    names(lids)[i-1]=colnames(adj.frac.coord)[i]
    
    loc=time.sp[rownames(time.sp) %in% ab.loc$LIDNUM,colnames(time.sp) %in% colnames(ab.loc)[2]]
    nfind.sp.fac[[i-1]]=loc
    names(nfind.sp.fac)[i-1]=colnames(adj.frac.coord)[i]
  }

  
  
#Different species
  tot.sp=c()
  no.indet500=find$speciesXMN
  for(i in 1:length(no.indet500)){
    ds=no.indet500[[i]]
    occ=colSums(ds[,2:(ncol(ds)-1)] != 0)
    tot.sp=rbind(tot.sp,occ)
  }
  rownames(tot.sp)=names(no.indet500)
  
  tot.sp.fac=c()
  for(i in 1:length(lids)){
    idx=unlist(lids[[i]]["LIDNUM"])
    adj.frac.coord[,names(lids)[i]]=as.numeric(adj.frac.coord[,names(lids)[i]])
    #idx.col=c("LIDNUM",names(lids)[i])
    #time=names(lids)[i]
    ab.loc=adj.frac.coord[adj.frac.coord$LIDNUM %in% idx,c("LIDNUM",names(lids)[i])]
    loc=tot.sp[rownames(tot.sp) %in% ab.loc$LIDNUM,colnames(tot.sp) %in% colnames(ab.loc)[2]]
    tot.sp.fac[[i]]=loc    
    names(tot.sp.fac)[i]=colnames(adj.frac.coord)[i]
  }
  
#Different species/Total findings  
  
  frac.dife=tot.sp/time.sp
  frac.dife[is.na(frac.dife)] <- 0
  
  frac.sp.fac=c()
  for(i in 1:length(lids)){
    idx=unlist(lids[[i]]["LIDNUM"])
    adj.frac.coord[,names(lids)[i]]=as.numeric(adj.frac.coord[,names(lids)[i]])
    #idx.col=c("LIDNUM",names(lids)[i])
    #time=names(lids)[i]
    ab.loc=adj.frac.coord[adj.frac.coord$LIDNUM %in% idx,c("LIDNUM",names(lids)[i])]
    loc=frac.dife[rownames(frac.dife) %in% ab.loc$LIDNUM,colnames(frac.dife) %in% colnames(ab.loc)[2]]
    frac.sp.fac[[i]]=loc
 }
  names(frac.sp.fac)=names(lids)
  
  
  
return(list(Sp.fac=lids,Total_findings_sp.fac=nfind.sp.fac, Total_findings_tot=time.sp, Different_species_sp.fac=tot.sp.fac,Different_species_tot=tot.sp, Fra_sp.fac=frac.sp.fac,Frac_tot=frac.dife))
}
