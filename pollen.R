
# row names are species/pollen taxa, colnames are sample ids

spectra_to_target_sum <- function(table, rand, summ) { 
  sample_ids <- colnames(table)
  stan=min(colSums(table)) # this line finds the lowest sum
  if(stan>summ) {stan <- summ}
  print(paste("all samples will be resampled to:", stan,"grains"))
  ssv <- data.frame(site=sample_ids, effort=NA,  q50=NA)
  retab <- rownames(table)
  i=7
  for(i in 1:NROW(sample_ids)){

    sit <- data.frame(rownames(table[table[,i]>0 ,]) , table[(table[,i]>0) ,i], stringsAsFactors = F)
    colnames(sit) <- c("taxa", "count") 
    rh <- character()
    e=1
    for(e in 1:nrow(sit)) { 
      rh <- c(rh, rep(as.character(sit[e,1]), sit[e,2]))
      #print(sit[e,2])
    }
    rhh <- data.frame(rh,1, stringsAsFactors = F)
    
    ## reshake
    poct <- integer()
    for (a in 1:rand){
      rhh <- rhh[sample(nrow(rhh)),]
      ddd <- aggregate(rhh[1:stan, 2], by = list(rhh[1:stan, 1]), "sum")
      sit <- merge(sit, ddd, by=1, all=T)
      poct <- c(poct,length(unique(rhh[1:stan, 1])))
    }
    print(paste(  "sample:" , sample_ids[i], "has:",median(poct), "pollen taxa"))  
    run <-  sit[,c(1,(which(poct==round(median(poct)))[1])+2)]
    retab <- merge(retab,run, by=1, all=T)
    
  }
  colnames(retab) <- c("taxa",sample_ids)
  return(retab)
}


age_to_tw <- function(df, agecolname, twdf, youngcol, oldcol)
{tw <- rep(NA, nrow(df))
for (i in seq(1,nrow(twdf),1))
{tw[which(df[,agecolname]>=twdf[i,youngcol] & df[,agecolname]<twdf[i,oldcol])] <- twdf[i,youngcol]+((twdf[i,oldcol]-twdf[i,youngcol])/2)}
return(cbind(df,tw))} 

cat("Here we go! We are useful functions to deal pollen assemblages and pollen-vegetation relationship\n\n")

