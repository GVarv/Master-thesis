curveball_abundance<-function(m,n){
  library(vecsets)
  #set.seed(23)
  #m=occurences[[4]]
  
  column.margins=m["col.margins",1:(ncol(m)-1)]
  row.margins=m[1:(nrow(m)-1),"species.margins"]
  m=m[-which(grepl("margins",rownames(m))),-which(grepl("margins",colnames(m)))]

  
  RC=dim(m)
  R=RC[1]  #Row number
  C=RC[2]  #Column number 
  hp=list() 
  na=c()
  k=0
  for (row in 1:R) {
    try=as.matrix(m[row,which(m[row,]!=0),drop=F]) #In which column the row is not empty?
    #try=m[row,which(m[row,]!=0)]
    hp[[row]]=unlist(lapply(colnames(try),function(x){rep(x,try[1,x])}))
    
  }
  names(hp)=rownames(m)
  hp2=hp
  hp3=hp
  #l_hp=length(hp)
  for (rep in 1:n){
    AB=sample(1:length(hp),2)
    #AB=c("Mesopithecus monspessulanus","Sus arvernensis")
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    
     #a=hp[["Sus arvernensis"]]
     #b=hp[["Dolomys nehringi"]]
     a
     b
    if(length(intersect(a,b))!=0){
        ab=vintersect(a,b) #This intersect count the location twice if they are in the intersection twice
    } else {
      ab=intersect(a,b)
    }
     #ab
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    
    #if(l_ab!=0){
     # print(AB)
       #print(hp[[AB[1]]])
      #cat(hp[[AB[2]]],"\n\n")
    #}
    if ((l_ab %in% c(l_a,l_b))==F){
      if(length(unlist (lapply (hp[AB], function (x) which (is.na (x)))))==0){
        
        c_ab=c(a,b)
        c_ab
        #sort(as.numeric(c_ab))
        #sort(as.numeric(vsetdiff(c_ab,ab)))
        tot=vsetdiff(vsetdiff(c_ab,ab,multiple = TRUE),ab) 
        tot 
        l_tot=length(tot)
  
        if(length(c_ab)==(l_tot+2*l_ab)){
          tot=sample(tot, l_tot, replace = FALSE, prob = NULL) #Shuffle the elements in c(a,b)
          tot
          L=l_a-l_ab
          n_a=c(ab,tot[1:L])
          #n_a
          if(length(n_a)==length(a)){
            hp[[AB[1]]] = n_a
          }
          n_b=c(ab,tot[(L+1):l_tot])
          #n_b
          if(length(n_b)==length(b)){
            hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])
          }
          
        
          if(length(unlist (lapply (hp[AB], function (x) which (is.na (x)))))>0){
          
            k=k+1
            cat("\n \nAfter analysis \n")
            print(hp[AB])
            cat("\nBefore analysis \n")
            print(hp2[AB])
          
            na[[k]]=names(hp2[AB])
            hp2=hp
        
          }
        }
      }
    }  
    #hp_l=lapply(hp,table)
    #hp_ll=lapply(hp_l,as.data.frame)
    #rm=matrix(0,R,C)
    #colnames(rm)=colnames(m)
    #rownames(rm)=rownames(m)
    #if(length(na)==0){
     # na=0
    #}
    #for (row in 1:R){
     # rm[row,as.character(unlist(hp_ll[[row]]["Var1"]))]=as.numeric(unlist(hp_ll[[row]]["Freq"]))
    #}
    #col.margins=colSums(rm)
    #species.margins=rowSums(rm)
    #if (length(which(row.margins!=species.margins))!=0 | length(which(col.margins!=column.margins))!=0) {
     # cat("\n \nAfter analysis \n")
      #print(hp[AB])
      #cat("\nBefore analysis \n")
      #print(hp2[AB])
      #break
    #}
  }
  hp_l=lapply(hp,table)
  hp_ll=lapply(hp_l,as.data.frame)
  rm=matrix(0,R,C)
  colnames(rm)=colnames(m)
  rownames(rm)=rownames(m)
  if(length(na)==0){
     na=0
  }
  for (row in 1:R){
    rm[row,as.character(unlist(hp_ll[[row]]["Var1"]))]=as.numeric(unlist(hp_ll[[row]]["Freq"]))
  }
  return(list(rand.mat=rm,na=na,hp2=hp2,orig.hp=hp3))
}
