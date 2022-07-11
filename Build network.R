build_network=function(n, filename){
  if (length(grep(pattern = "w32", x = version["os"]))) {
    eol <- "\n"
  }
  else {
    eol <- "\r\n"
  }

  
  nodes_1=rownames(n)
  nodes_2=colnames(n)
  
  cat(paste("*Vertices", sum(dim(n))), eol, 
      file = filename)
  cat(paste(1:dim(n)[1], " \"", nodes_1, "\"", eol, sep = ""), 
      file = filename, append = TRUE)
  cat(paste(seq(dim(n)[1] + 1, length = dim(n)[2]), " \"", 
            nodes_2, "\"", eol, sep = ""), file = filename, append = TRUE)
  cat("*Intra", eol, file = filename, append = TRUE)
  #Function melt to flatten matrix
  
  for (i in 1:dim(n)[1]) {
    for (j in 1:dim(n)[2]) {
      if (n[i, j] != 0) {
        cat(paste(n[i, j],i, j + dim(n)[1], 1,  eol), 
            file = filename, append = TRUE)
      }
    }
  }
}
