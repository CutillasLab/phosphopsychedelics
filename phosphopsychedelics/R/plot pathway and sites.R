#Author: Pedro Rodriguez Cutillas
# Feb 2025

plot.selected.pathway <- function(pathway, x_factor="time_min", y_factor="site", value="mean",
                                  title="Abundance vs DMSO"){

  genes <- unique(unlist(strsplit(df.process[df.process$ontology==pathway,]$genes,";")))
  selected_sites <- foo[foo$gene %in% genes,]$peptide
  df <- data.frame()
  for (s in selected_sites){
    x <- plot_one_site(s)
    if (max(abs(x$data$mean))>0){
      df <- rbind.data.frame(x$data,df)
    }
  }
  df$drug <- df$treatment_1
  dfx <- reshape2::dcast(df,site~drug,value.var = "mean", fun.aggregate = mean)
  dfx$average <- apply(dfx[,2:ncol(dfx)],1,mean)
  dfx <- dfx[order(dfx$average),]



  kplot <- ggplot(df
                  ,aes_string(x=x_factor,y=y_factor,fill=value, size=abs(df$mean)))+
    geom_point(shape=21,color="grey")+
    scale_fill_gradient2(low="royalblue",high = "seagreen",mid = "white") +
    scale_y_discrete(limits=dfx$site)+
    scale_x_discrete(limits=c("0","5","15","30"))+
    theme_minimal()+
    # theme_pubr()+
    # theme_pubclean()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "right")+
    facet_wrap(drug_type~drug, nrow = 1)+
    labs(title = title,x="Treatment time (min)")


  return(kplot)
}

plot.selected.pathway_v2 <- function(pathway, x_factor="time_min", y_factor="site", value="mean",
                                     title="Abundance vs DMSO"){

  genes <- unique(unlist(strsplit(df.process[df.process$ontology==pathway,]$genes,";")))
  selected_sites <- foo[foo$gene %in% genes,]$peptide
  df <- data.frame()
  for (s in selected_sites){
    x <- plot_one_site(s)
    if (max(abs(x$data$mean))>1){
      df <- rbind.data.frame(x$data,df)
    }
  }
  df$drug <- df$treatment_1
  dfx <- reshape2::dcast(df,site~drug,value.var = "mean", fun.aggregate = mean)
  dfx$average <- apply(dfx[,2:ncol(dfx)],1,mean)
  dfx <- dfx[order(dfx$average),]



  kplot <- ggplot(df
                  ,aes_string(x=x_factor,y=y_factor,fill=value, size=abs(df$mean)))+
    geom_point(shape=21,color="grey")+
    scale_fill_gradient2(low="royalblue",high = "seagreen",mid = "white") +
    scale_y_discrete(limits=dfx$site)+
    scale_x_discrete(limits=c("0","5","15","30"))+
    theme_minimal()+
    # theme_pubr()+
    # theme_pubclean()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "right")+
    facet_wrap(drug_type~drug, nrow = 1)+
    labs(title = title,x="Treatment time (min)")


  return(kplot)
}




plot.selected.sites <- function(selected_sites, x_factor="time_min", y_factor="site", value="mean",
                                title="Abundance vs DMSO", max.fold=1){


  df <- data.frame()
  for (s in selected_sites){
    x <- plot_one_site(s)
    if (max(abs(x$data$mean))>max.fold){
      df <- rbind.data.frame(x$data,df)
    }
  }
  df$drug <- df$treatment_1
  dfx <- reshape2::dcast(df,site~drug,value.var = "mean", fun.aggregate = mean)
  dfx$average <- apply(dfx[,2:ncol(dfx)],1,mean)
  dfx <- dfx[order(dfx$average),]



  kplot <- ggplot(df
                  ,aes_string(x=x_factor,y=y_factor,fill=value, size=abs(df$mean)))+
    geom_point(shape=21,color="grey")+
    scale_fill_gradient2(low="royalblue",high = "seagreen",mid = "white") +
    scale_y_discrete(limits=dfx$site)+
    scale_x_discrete(limits=c("0","5","15","30"))+
    theme_minimal()+
    # theme_pubr()+
    # theme_pubclean()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "right")+
    facet_wrap(drug_type~drug, nrow = 1)+
    labs(title = title,x="Treatment time (min)")


  return(kplot)
}
