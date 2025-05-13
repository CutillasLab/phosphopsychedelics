
# Author: Pedro R. Cutillas
# Feb 2025


compare.by.ttest <- function(df.to.compare,
                             log.data=FALSE,
                             control.columns="",
                             test.columns="",
                             is.paired=FALSE)
{


  # Compare by ttest

  library(foreach)
  library(doParallel)

  nc <- ncol(df.to.compare)
  nr <- nrow(df.to.compare)

  if (log.data){
    df.s <- data.frame(prot=df.to.compare[,1], scale(log2(df.to.compare[,2:nc])))
  }else{
    df.s <- df.to.compare
  }

  if (length(control.columns)<2){
    print("Give me the names of control.columns and test.columns.")
    return()
  }
  if (length(test.columns)<2){
    print("Give me the names of control.columns and test.columns.")
    return()
  }

  registerDoParallel(cores = parallel::detectCores()-1)
  # print(paste("running on", parallel::detectCores()-1, "cores"))
  t1 <- Sys.time()

  modeled.values <- foreach(r = 1:nr, .combine = "rbind")%dopar%
    {

      cont <- unlist(df.s[r,control.columns])
      test <- unlist(df.s[r,test.columns])



      if (length(na.omit(cont))>2 & length(na.omit(test))>2){
        pval <- 1
        tryCatch(
          {pval <- t.test( cont,test, paired=is.paired )$p.value},error=function(e){}
        )
        fold <- sum(test,na.rm = T)-sum(cont,na.rm = T)

        result <- c(df.s[r,1],fold,pval)

        #ct <- cor.test(cont,test)

        #result <- c(df.s[r,1],fold,pval,ct$estimate,ct$p.value)
      }else{
        result <- c(df.s[r,1],0,1,0,1)
      }
      result
    }
  t2 <- Sys.time()
  print(t2-t1)
  stopImplicitCluster()

  df <- data.frame(modeled.values)
  colnames(df) <- c("protein","difference.test.vs.control","pvalues")#,"rvalue.pearson","pvalue.pearson")
  df$fold <- as.numeric(df$difference.test.vs.control)
  df$pvalues <- as.numeric(df$pvalues)
  df$FDR <- p.adjust(df$pvalues,method = "fdr")
  #df$abs_rvalue <- abs(as.numeric(df$rvalue.pearson))
  rownames(df) <- make.names(df$protein,unique = T)

  return(df)
}



enrichment.from.list.function <- function(
    list.of.peptides,
    background.list,
    prot_db, is.ksea=F){

  # Returns enrichment of phosphosites, proteins or transcripts in pathways ontologies or kinase-phosphosite relationsips
  #
  # list of peptides is the list of proteins, phosphosites or transcripts to be searched (transcripts needt o be converted to Uniprot names)
  # Accepts Uniprot names or phosphosite names in the form gene(xy), were gene = uniprot gene name; x=S, T or Y; y= amino acid position of phosphorylation
  # background.list is the background list of of proteins, phosphosites or transcripts
  # prot_db = any suitable ontology or kinase-substrate database
  # is.ksea = change to TRUE for analysing phosphoproteomics data against kinase-substrates
  # draw.tables.and.plots = if TRUE draws volcano plots of the enrichment data, requires ggplot2
  #
  # Enrichment = (q/k)/(j/m), where:
  #   q = peptides/proteins in list.of.peptides with a match in the ontology/K-S relationship dataset
  #   k = number of peptides/proteins in list.of.peptides
  #   j = peptides/proteins in background.list with a match in the ontology/K-S relationship dataset
  #   m = number of peptides/proteins in background.list
  #
  # p-values are calculated using the hypergeometric function and adjusted by FDR method
  #

  list.of.peptides <- gsub("..",");", list.of.peptides,fixed = T)
  background.list <- gsub("..",");", background.list,fixed = T)
  list.of.peptides <- gsub(".","(", list.of.peptides,fixed = T)
  background.list <- gsub(".","(", background.list,fixed = T)

  ll <- list.of.peptides[!grepl("_no_mod",list.of.peptides,fixed = T)]
  ll <- ll[!grepl("_Phospho (",ll,fixed=T)]
  list.of.peptides <- strsplit(ll,";")
  original.list.of.peptides <- list.of.peptides
  background.list <- strsplit(background.list,";")
  column.with.priors <- "genes"

  df.peptides <- data.frame(peptides=unlist(original.list.of.peptides),
                            proteins=unlist(list.of.peptides))

  if (is.ksea==TRUE){
    column.with.priors <- 3
  }else{

    #Convert phosphopeptides to genes for phosphoproteomics data when matching to proteomics sets
    x1 <-  grepl(")", as.character(list.of.peptides),fixed=T)
    if (TRUE %in% x1){ # True if phosphoproteomics data, false if proteomics data
      list.of.peptides <- unlist(lapply(df.peptides$proteins,
                                        function(x) strsplit(as.character(x),"(",fixed = TRUE)[[1]][1]))
      df.peptides$proteins <- list.of.peptides

      background.list <- unlist(lapply(background.list,
                                       function(x) strsplit(as.character(x),"(",fixed = TRUE)[[1]][1]))
    }else{
      column.with.priors <- "genes"
    }
  }

  list.of.peptides <- unique(trimws (unlist(list.of.peptides)))
  background.list <- unique(trimws(  unlist(background.list)))
  ####################



  df.ks <- protools2::protein_and_ks_sets[[prot_db]]

  #df.ks <- df.ks[order(-df.ks[,2]),]
  nr <- nrow(df.ks)




  pathway <- character(nr)
  pvalues <- numeric(nr)
  enrichment <- numeric(nr)
  counts <- numeric(nr)
  counts.bg <- numeric(nr)
  FDR <- numeric(nr)
  proteins <- character(nr)
  genes <- character(nr)

  m <- length(background.list)
  k <- length(list.of.peptides)

  bg.size <- numeric(nr)
  data.size <- numeric(nr)
  r=1

  for (r in 1:nr) {
    mym <- df.ks[r,2]
    kinase <- df.ks[r,1]
    if (is.na(mym) == F){
      if(mym>2){
        substrates <- as.character(df.ks[r,column.with.priors])
        ss <- c(unlist(strsplit(substrates,";")))
        #if (is.ksea==TRUE | is.ksea==T){
        #  ss <- paste(ss,";",sep = "")
        #}
        start.time <- Sys.time()

        q <- length(intersect(ss,list.of.peptides))
        j <- length(intersect(ss,background.list))
        n <- m-j

        if (j>0){
          prots <- intersect(ss,list.of.peptides)
          peptides <- df.peptides[df.peptides$proteins %in% prots,"peptides"]
          pvalue <- 1-phyper(q-1,j,m-j,k,lower.tail = TRUE, log.p = F)
          pathway[r] <- as.character(kinase)
          pvalues[r] <- pvalue
          enrichment[r] <- round((q/k)/(j/m), digits = 2)
          counts[r] <- q
          data.size[r] <- k
          counts.bg[r] <- j
          bg.size[r] <- m
          genes[r] <- paste(peptides,collapse = ";")
          proteins[r] <- paste(prots,collapse = ";")
        }
      }
      r=r+1
    }
  }
  results <- data.frame(pathway, pvalues,enrichment,counts, data.size,counts.bg,bg.size, proteins,genes)
  #results <- na.omit(results)
  if (nrow(results)>0){

    results <- subset(results, enrichment!=0 & counts>0)
    results <- results[order(-results$enrichment),]
    results$FDR <- p.adjust(results$pvalues, method = "BH")
  }

  return(results)
}

