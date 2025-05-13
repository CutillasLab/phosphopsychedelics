# Author: Pedro R. Cutillas
# Feb 2025


ks.enrichment <- function (dfx, ks_db, is.ksea = TRUE)
{
  dfx <- data.frame(dfx)
  rownames(dfx) <- make.names(dfx[, 1], unique = T)
  if (is.ksea == TRUE) {
    column.with.priors <- 3
  }else {
    x1 <- grepl("(", as.character(dfx[, 1]), fixed = T)
    if (TRUE %in% x1) {
      dfx[, 1] <- unlist(lapply(dfx[, 1], function(x) strsplit(as.character(x),
                                                               "(", fixed = TRUE)[[1]][1]))

      dfx[,1] <- paste0(dfx[,1],";")

      column.with.priors <- "genes"
    }
    x1 <- grepl(".", as.character(dfx[, 1]), fixed = T)
    if (TRUE %in% x1) {
      dfx[, 1] <- unlist(lapply(dfx[, 1], function(x) strsplit(as.character(x),
                                                               ".", fixed = TRUE)[[1]][1]))

      dfx[,1] <- paste0(dfx[,1],";")

      column.with.priors <- "genes"
    }



  }
  dfx[, 1] <- gsub("..", ");", dfx[, 1], fixed = T)
  dfx[, 1] <- gsub(".", "(", dfx[, 1], fixed = T)
  dfx[is.na(dfx)] <- 0
  nc <- ncol(dfx)
  df.ks <- protools2::protein_and_ks_sets[[ks_db]]
  nr <- nrow(df.ks)
  results.pvalues <- numeric(nr)
  results.zscores <- numeric(nr)
  results.q <- numeric(nr)
  results.m <- numeric(nr)
  results.distance <- numeric(nr)
  results.sites <- character(nr)
  kinases <- character(nr)
  r = 4
  for (r in 1:nr) {
    mym <- df.ks[r, 2]
    kinase <- df.ks[r, 1]
    if (is.na(mym) == F) {
      if (mym > 4) {
        substrates <- as.character(df.ks[r, column.with.priors])
        ss <- c(unlist(strsplit(substrates, ";")))
        start.time <- Sys.time()
        df.xx <- subset(dfx, dfx[, 1] %in% paste(ss,
                                                 ";", sep = ""))
        sites <- paste(unlist(rownames(df.xx)), collapse = ";")
        if (nrow(df.xx) > 2) {
          df.xx$prot.group <- unlist(kinase)
          sites.x <- data.frame(site = df.xx[, 1])
          sites.x$kinase <- unlist(kinase)
          c = 2
          ds <- numeric(nc - 1)
          pvals <- numeric(nc - 1)
          zscores <- numeric(nc - 1)
          ms <- numeric(nc - 1)
          qs <- numeric(nc - 1)
          values.all <- na.omit(as.numeric(subset(dfx[,
                                                      c], dfx[, c] != 0)))
          myvalues <- na.omit(as.numeric(subset(df.xx[,
                                                      c], df.xx[, c] != 0)))
          pval <- 1
          tryCatch({
            myks <- ks.test(values.all, myvalues)
            pval <- myks$p.value
          }, error = function(e) {
          })
          m <- 0
          q <- 0
          m <- nrow(df.xx)
          mysd <- sd(values.all, na.rm = T)
          mymedian <- median(myvalues, na.rm = T)
          mymedian.all <- median(values.all, na.rm = T)
          sd.all <- sd(values.all)
          results.distance[r] <- mymedian - mymedian.all
          results.pvalues[r] <- pval
          results.zscores[r] <- ((mymedian - mymedian.all) *
                                   ((sqrt(m)))/mysd)
          results.m[r] <- m
          results.sites[r] <- as.character(sites)
          kinases[r] <- unlist(kinase)
        }
      }
    }
  }
  xx <- data.frame(kinases, zscores = results.zscores, pvalues = results.pvalues,
                   m = results.m, distance = results.distance, sites = results.sites,
                   kinase_dataset = paste0(kinases, "_", ks_db))
  xx <- subset(xx, xx$m > 1)
  xx$qvalue <- p.adjust(xx$pvalues, method = "fdr")
  return(xx)
}
