# Author: Pedro R. Cutillas
# Feb 2025


get_cluster_n_value <- function(df){

  library(factoextra)
  library(cluster)
  # compute gap statistic
  set.seed(123)
  gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
                      K.max = 10, B = 50)
  fviz_gap_stat(gap_stat)


}


plot_one_site <- function(site){

  c23 <- protools2::mycolors()$C23

  library(ggplot2)
  library(ggpubr )

  df.relative.to.dmso <- phosphopsychedelics::df.folds.relative.to.dmso

  treatments <-c("Ariadne" , "BNDF"   ,  "KETAMINE", "LSD"   ,
                 "PSILOCIN" ,"DOI"    ,  "DMT"   ,   "LISURIDE" ,"BrLSB"  )
  xx <- unique(df.design[,c("drug_type","treatment_1")])
  xx <- xx[xx$treatment_1 %in% treatments,]

  drug.types <- unique(df.design$drug_type)

  drug.types <- c("Hallucinogenic_serotoninergic_psychedelic" ,"non_hallucinogenic_serotoninergic_psychedelic","nonserotoninergic_psychedelic" ,
                  "Trkb_ligand")

  dfdmso <- data.frame(time_1="0MIN",
                       drug_type=xx$drug_type,
                       treatment_1=xx$treatment_1,
                       mean=0,
                       sd=0,
                       sem=0,
                       n=4)
  dfdmso <- dfdmso[dfdmso$drug_type %in% drug.types,]
  dfdmso$drug_class <- NA

  df.design$x <- paste(df.design$condition,df.design$biological_replicate,sep="_")

    if(site %in% rownames(df.relative.to.dmso)){

      name <- df.combi[grepl(site, df.combi$sites),"protein"][1]
      #$protein[1]


      df.pp.s <- data.frame(t(df.relative.to.dmso[site,2:ncol(df.relative.to.dmso)]))
      df.pp.s <- data.frame(x=rownames(df.pp.s),df.pp.s)
      df.pp.s<- merge.data.frame(df.pp.s,df.design,by="x")
      colnames(df.pp.s)[2] <- "signal"
      dfm <- miscscripts::data_summary(df.pp.s,varname = "signal",groupnames = c("time_1","drug_type", "drug_class", "treatment_1"))

      dfm <- dfm[dfm$treatment_1!="DMSO",]

      dfm <- rbind.data.frame(dfm,dfdmso)
      dfm$drug_type_f <- factor(dfm$drug_type,levels = drug.types)
      dfm$site <- site
      dfm$protein_name <- name

      if (nrow(na.omit(dfm))>0){


          dfm$time_min <- gsub("MIN","", dfm$time_1,fixed=T)
          # plot time vs mean phosphorylation
          pp <- ggplot(dfm,aes(x=time_min,y=mean,group=treatment_1,color=as.character(drug_type)))+
            #geom_boxplot()+
            geom_line()+
            geom_point()+
            geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                          position=position_dodge(.9))+
            scale_color_manual(values = c23)+
            theme_bw()+
            geom_hline(yintercept = 0)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "right")+
            facet_wrap(drug_type_f~treatment_1,  nrow = 1)+
            scale_x_discrete(limits=c("0", "5","15","30"))+
            labs(title = paste(site,"\n",name), color="Drug class",
                 x="treatment time (min)",y= "Mean Log2 fold over DMSO")+
            geom_hline(yintercept = 0.5, linetype=3)+
            geom_hline(yintercept = - 0.5, linetype=3)

           return(list(plot=pp,data=dfm))
        }
      }
    }







plot_each_site <- function(selected.sites, outputDir, outputFolder, tolerance = 0.5){

  c23 <- protools2::mycolors()$C23

  library(ggplot2)
  library(ggpubr )

  outputFolder <- paste0("./",outputFolder,"/")

  dir.create(file.path(outputDir, outputFolder), showWarnings = FALSE)

  treatments <-c("Ariadne" , "BNDF"   ,  "KETAMINE", "LSD"   ,
                 "PSILOCIN" ,"DOI"    ,  "DMT"   ,   "LISURIDE" ,"BrLSB"  )
  xx <- unique(df.design[,c("drug_type","treatment_1")])
  xx <- xx[xx$treatment_1 %in% treatments,]

  drug.types <- unique(df.design$drug_type)

  drug.types <- c("Hallucinogenic_serotoninergic_psychedelic" ,"nonserotoninergic_psychedelic" ,"non_hallucinogenic_serotoninergic_psychedelic",
                  "Trkb_ligand")

  dfdmso <- data.frame(time_1="0MIN",
                       drug_type=xx$drug_type,
                       treatment_1=xx$treatment_1,
                       mean=0,
                       sd=0,
                       sem=0,
                       n=4)
  dfdmso <- dfdmso[dfdmso$drug_type %in% drug.types,]
  dfdmso$drug_class <- NA

  for (site in selected.sites){

    if(site %in% rownames(df.relative.to.dmso)){

      name <- df.combi[grepl(site, df.combi$sites),"protein"]$protein[1]


      df.pp.s <- data.frame(t(df.relative.to.dmso[site,2:ncol(df.relative.to.dmso)]))
      df.pp.s <- data.frame(x=rownames(df.pp.s),df.pp.s)
      df.pp.s<- merge.data.frame(df.pp.s,df.design,by="x")
      colnames(df.pp.s)[2] <- "signal"
      dfm <- miscscripts::data_summary(df.pp.s,varname = "signal",groupnames = c("time_1","drug_type", "drug_class", "treatment_1"))

      dfm <- dfm[dfm$treatment_1!="DMSO",]

      dfm <- rbind.data.frame(dfm,dfdmso)
      dfm$drug_type_f <- factor(dfm$drug_type,levels = drug.types)

      if (nrow(na.omit(dfm))>0){
        if (max(abs(dfm$mean))>tolerance){

          dfm$time_min <- gsub("MIN","", dfm$time_1,fixed=T)
          # plot time vs mean phosphorylation
          pp <- ggplot(dfm,aes(x=time_min,y=mean,group=treatment_1,color=as.character(drug_type)))+
            #geom_boxplot()+
            geom_line()+
            geom_point()+
            geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                          position=position_dodge(.9))+
            scale_color_manual(values = c23)+
            theme_bw()+
            geom_hline(yintercept = 0)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "right")+
            facet_wrap(drug_type_f~treatment_1,  nrow = 2)+
            scale_x_discrete(limits=c("0", "5","15","30"))+
            labs(title = paste(site,"\n",name), color="Drug class",
                 x="treatment time (min)",y= "Mean Log2 fold over DMSO")+
            geom_hline(yintercept = 0.5, linetype=3)+
            geom_hline(yintercept = - 0.5, linetype=3)
          # plot(pp)
          # areas
          aas <- data.frame(treatment_1=treatments,t(df.areas[df.areas$site==site,treatments]))
          colnames(aas)[2] <- "area"
          aas <- unique(merge.data.frame(aas,df.design[,c("drug_type","drug_class", "treatment_1")],by="treatment_1", all.x = T,all.y = F))
          pp
          pp2 <- ggplot(aas,aes(x=drug_type,y=area,color=drug_type))+
            geom_boxplot()+
            geom_jitter(aes(shape=treatment_1))+
            scale_shape_manual(values=c(1:10))+
            scale_color_manual(values = c23)+
            stat_compare_means(comparisons = list(c("Hallucinogenic_serotoninergic_psychedelic" ,"non_hallucinogenic_serotoninergic_psychedelic")),
                               method = "t.test")+
            theme_bw()+
            theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position = "none",
              axis.text.x = element_blank(),
              #axis.title.x=element_blank()
            )+
            geom_hline(yintercept = 0, linetype=2)+
            ylim(c(NA,max(aas$area)+0.5))+
            labs(title = paste(site,"\n","area"), color="Drug class",y= "Area")+
            facet_wrap(~.)
          pp2


          pp3 <- ggplot(aas,aes(x=drug_class,y=area,color=drug_class))+
            geom_boxplot()+
            geom_jitter(aes(shape=treatment_1))+
            scale_shape_manual(values=c(1:10))+
            scale_color_manual(values = c23[10:25])+
            stat_compare_means(label.y = max(aas$area)+0.2)+
            theme_bw()+
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              legend.position = "none",
              #axis.text.x = element_blank(),
              #axis.title.x=element_blank()
            )+
            geom_hline(yintercept = 0, linetype=2)+
            ylim(c(NA,max(aas$area)+0.5))+
            labs(title = paste(site,"\n","area"), color="Drug class",y= "Area")+
            facet_wrap(~.)
          pp3





          #5 minutes

          xx5 <- dfm[dfm$time_1=="5MIN",]
          xx15<- dfm[dfm$time_1=="15MIN",]
          xx30<- dfm[dfm$time_1=="30MIN",]
          #xx5 <- data.frame(treatment_1=treatments,t(df.5min[df.5min$sites==site,treatments]))
          #colnames(xx5)[2] <- "Fold"
          #xx5 <- unique(merge.data.frame(xx5,df.design[,c("drug_type","treatment_1")],by="treatment_1", all.x = T,all.y = F))
          #xx5$treatment_time <- "5MIN"
          # 15 minutes
          #xx15 <- data.frame(treatment_1=treatments,t(df.15min[df.15min$sites==site,treatments]))
          #colnames(xx15)[2] <- "Fold"
          #xx15$treatment_time <- "15MIN"
          #xx15 <- unique(merge.data.frame(xx15,df.design[,c("drug_type","treatment_1")],by="treatment_1", all.x = T,all.y = F))
          #30 minutes
          #xx30 <- data.frame(treatment_1=treatments,t(df.30min[df.30min$sites==site,treatments]))
          #colnames(xx30)[2] <- "Fold"
          #xx30 <- unique(merge.data.frame(xx30,df.design[,c("drug_type","treatment_1")],by="treatment_1", all.x = T,all.y = F))
          #xx30$treatment_time <- "30MIN"
          # merge
          #aas$treatment_time <- "Area"
          xxx <- rbind.data.frame(xx5,xx15,xx30)
          xxx$time_1 <- factor(xxx$time_1,levels = c("5MIN","15MIN","30MIN"))

          if (max(abs(xxx$mean))>tolerance){

            # boxplots at differen times
            pp.times <- ggplot(xxx,aes(x=drug_type,y=mean,color=drug_type))+
              geom_boxplot()+
              geom_jitter(aes(shape=treatment_1))+
              scale_shape_manual(values=c(1:10))+
              scale_color_manual(values = c23)+
              stat_compare_means(comparisons = list(c("Hallucinogenic_serotoninergic_psychedelic" ,"non_hallucinogenic_serotoninergic_psychedelic")),
                                 method = "t.test")+
              theme_bw()+
              theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                legend.position = "none",
                legend.box = "vertical",
                axis.text.x = element_blank(),
                #axis.title.x=element_blank()
              )+
              geom_hline(yintercept = 0, linetype=2)+
              labs(title = paste(site,"\n","Fold over DMSO at diff times by drug type"), color="Drug class",
                   y= "Mean Log2 fold over DMSO")+
              ylim(c(NA,max(xxx$mean)+0.5))+
              facet_wrap(~time_1)

            pp.times
            pp.times2 <- ggplot(xxx,aes(x=drug_class,y=mean,color=drug_class))+
              geom_boxplot()+
              geom_jitter(aes(shape=treatment_1))+
              scale_shape_manual(values=c(1:10))+
              scale_color_manual(values = c23[10:25])+
              stat_compare_means(label.y = max(xxx$mean)+0.2)+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                legend.position = "none",
                legend.box = "vertical",
                #axis.text.x = element_blank(),
                #axis.title.x=element_blank()
              )+
              geom_hline(yintercept = 0, linetype=2)+
              labs(title = paste(site,"\n","Fold over DMSO at diff times by chemical class"), color="Drug class",
                   y= "Mean Log2 fold over DMSO")+
              ylim(c(NA,max(xxx$mean)+0.5))+
              facet_wrap(~time_1)

            pp.times2





            px1 <- cowplot::plot_grid(pp.times,pp2, rel_widths = c(1,0.4))
            px2 <- cowplot::plot_grid(pp.times2,pp3, rel_widths = c(1,0.4))

            pp.all <- cowplot::plot_grid(pp,px1,px2,nrow =3, rel_heights = c(1.2,.7,0.9) )
            pp.all
            #plot(pp.all)


            ggsave(filename = paste0(outputFolder,site,"plot.pdf"),plot = pp.all,width = 10, height = 10)
            ggsave(filename = paste0(outputFolder,site,"plot.png"),plot = pp.all,width = 10, height = 10)

          }
        }
      }
    }
  }
}



plot_medians_of_clusters <- function(df,kk,k){

  library(ggplot2)
  clusters <- c(1:k)

  cluster <- 1

  df.median.custers <- data.frame()

  for(c in clusters){
    print(c)
    s <- kk[kk$cluster==c,]$site
  #  print(s)

    x <- data.frame( cluster=c,median=apply(df[s,],2,mean,na.rm=T),
                     sd=apply(df[s,],2,sd,na.rm=T)
    )
    x$condition <- rownames(x)

    df.median.custers <-  rbind.data.frame(df.median.custers, x)

  }

  df.median.custers[is.na(df.median.custers)] <- 0

  df.median.custers <- merge.data.frame(df.median.custers,df.design,by="condition")


  x <- df.median.custers[df.median.custers$drug_type=="Hallucinogenic_serotoninergic_psychedelic",
                         c("cluster","median","sd","treatment_1","time_1")]

  treatments <- unique(x$treatment_1)

  xx <- data.frame(cluster= rep( c(1:k), times=length(treatments)),
                   median=0,
                   sd=0,
                   treatment_1=rep(treatments,each=k),
                   time_1="0MIN"
  )

  xx[order(xx$cluster),]


  x1 <- rbind.data.frame(x,xx)
  colnames(x)
  #f.median.custers$c <- as.character(df.median.custers$c)
  pp1 <- ggplot(x1,
                aes(x=time_1,y=median, group=treatment_1))+
    geom_point(color="red")+
    geom_line(color="red")+
    #scale_color_gradient2(low="blue",mid = "white",high = "red",midpoint = 0)+
    geom_hline(yintercept = 0,linetype=2)+
    scale_x_discrete(limits= c("0MIN", "5MIN","15MIN","30MIN"))+
    theme_bw()+
    theme(legend.position = "right")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)#,
          #axis.text.y = element_blank(),
          #axis.title.y=element_blank()
    )+
    geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.2,
                  position=position_dodge(.9),
                  color="red")+
    #ylim(c(-1,1))+
    facet_grid(cluster~treatment_1)
  pp1


  x <- df.median.custers[df.median.custers$drug_type=="non_hallucinogenic_serotoninergic_psychedelic",
                         c("cluster","median","sd","treatment_1","time_1")]
  treatments <- unique(x$treatment_1)

  xx <- data.frame(cluster= rep( c(1:k), times=length(treatments)),
                   median=0,
                   sd=0,
                   treatment_1=rep(treatments,each=k),
                   time_1="0MIN"
  )
  x2 <- rbind.data.frame(x,xx)

  pp2 <- ggplot(x2,
                aes(x=time_1,y=median, group=treatment_1))+
    geom_point(color="blue")+
    geom_line(color="blue")+
    #scale_color_gradient2(low="blue",mid = "white",high = "red",midpoint = 0)+
    geom_hline(yintercept = 0,linetype=2)+
    scale_x_discrete(limits= c("0MIN", "5MIN","15MIN","30MIN"))+
    theme_bw()+
    theme(legend.position = "right")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)#,
          #axis.text.y = element_blank(),
          #axis.title.y=element_blank()
    )+
    geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.2,
                  position=position_dodge(.9),
                  color="blue")+
    #ylim(c(-1,1))+
    facet_grid(cluster~treatment_1)

  x <- df.median.custers[df.median.custers$drug_type %in% c("Trkb_ligand","nonserotoninergic_psychedelic"),
                         c("cluster","median","sd","treatment_1","time_1")]
  treatments <- unique(x$treatment_1)

  xx <- data.frame(cluster= rep( c(1:k), times=length(treatments)),
                   median=0,
                   sd=0,
                   treatment_1=rep(treatments,each=k),
                   time_1="0MIN"
  )
  x3 <- rbind.data.frame(x,xx)

  pp3 <- ggplot(x3,
                aes(x=time_1,y=median, group=treatment_1))+
    geom_point(color="seagreen")+
    geom_line(color="seagreen")+
    #scale_color_gradient2(low="blue",mid = "white",high = "red",midpoint = 0)+
    geom_hline(yintercept = 0,linetype=2)+
    scale_x_discrete(limits= c("0MIN", "5MIN","15MIN","30MIN"))+
    theme_bw()+
    theme(legend.position = "right")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)#,
          #axis.text.y = element_blank(),
          #axis.title.y=element_blank()
    )+
   # ylim(c(-1,1))+
    geom_errorbar(aes(ymin=median-sd, ymax=median+sd), width=.2,
                  position=position_dodge(.9),
                  color="seagreen")+

    facet_grid(cluster~treatment_1)

  pp4 <- cowplot::plot_grid(pp1,pp2,pp3,nrow=1, rel_widths = c(1,0.9,0.8))



  return(list(plot=pp4,data=rbind(x1,x2,x3)))

}


plot_kmeans <- function(df,k){

  library(ggplot2)

  df.areas.selected.m <- reshape2::melt(as.matrix(  df))
  colnames(df.areas.selected.m) <- c("sites","condition","value")


  kk <- kmeans(df,centers=k)
  kk <- data.frame(kk$cluster)

  kk <- data.frame(sites=rownames(kk),kk)
  kk <- kk[order(kk$kk.cluster),]

  #diff.sig.sites <- diff.sig.sites[order(diff.sig.sites$diff),]
  #diff.sig.sites$x <- "F_Ha.v.NH"

  df.areas.selected.m <- merge.data.frame(df.areas.selected.m,kk,by="sites")

  foo <- data.frame(condition=df.areas.selected.m$condition,do.call('rbind', strsplit(as.character(df.areas.selected.m$condition),'_',fixed=TRUE)))
  colnames(foo)[2:3] <- c("drug","time")

  df.areas.selected.m <- merge.data.frame(foo,df.areas.selected.m,by="condition")
  df.areas.selected.m <- merge.data.frame(df.areas.selected.m,df.design,by="condition")


  pplot1 <- ggplot(df.areas.selected.m[df.areas.selected.m$drug_type=="Hallucinogenic_serotoninergic_psychedelic",],
                   aes(x=time,y=sites,color=value))+
    geom_point(size=3,shape=15)+
    scale_color_gradient2(low="blue",mid = "white",high = "red",midpoint = 0)+
    #scale_y_discrete(limits=kk$sites)+
    scale_x_discrete(limits= c("5MIN","15MIN","30MIN"))+
    theme_bw()+
    theme(legend.position = "right")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)#,
          #axis.text.y = element_blank(),
          #axis.title.y=element_blank()
    )+
    facet_grid(kk.cluster~drug, scales = "free_y")+
    labs(title = "Hallucinogenic_serotoninergic_psychedelic" )
  pplot1

  x <- df.areas.selected.m[df.areas.selected.m$drug_type=="non_hallucinogenic_serotoninergic_psychedelic",]
  #x <- x[x$drug_type!="Vehicle",]

  pplot2 <- ggplot(x,
                   aes(x=time,y=sites,color=value))+
    geom_point(size=3,shape=15)+
    scale_color_gradient2(low="blue",mid = "white",high = "red",midpoint = 0)+
    #scale_y_discrete(limits=kk$sites)+
    scale_x_discrete(limits= c("5MIN","15MIN","30MIN"))+
    theme_bw()+
    theme(legend.position = "right")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)#,
          #axis.text.y = element_blank(),
          #axis.title.y=element_blank()
    )+
    facet_grid(kk.cluster~drug, scales = "free_y")+
    labs(title = "Non_hallucinogenic_serotoninergic_psychedelic" )
  pplot2

  x <- df.areas.selected.m[df.areas.selected.m$drug_type %in% c("Trkb_ligand","nonserotoninergic_psychedelic"),]

  pplot3 <- ggplot(x,
                   aes(x=time,y=sites,color=value))+
    geom_point(size=3,shape=15)+
    scale_color_gradient2(low="blue",mid = "white",high = "red",midpoint = 0)+
    #scale_y_discrete(limits=kk$sites)+
    scale_x_discrete(limits= c("5MIN","15MIN","30MIN"))+
    theme_bw()+
    theme(legend.position = "right")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)#,
          #axis.text.y = element_blank(),
          #axis.title.y=element_blank()
    )+
    facet_grid(kk.cluster~drug, scales = "free_y")+
    labs(title = "Other" )
  pplot3



  pplot4 <- cowplot::plot_grid(pplot1,pplot2,pplot3,nrow=1, rel_widths = c(1,0.9,0.8))

  pplot4

  return(list(plot=pplot4, kk=kk))
  #ggsave(filename = "heatmapSigSites_KmeansSigSitesOnly.pdf",plot = pplot4, width = 17, height = 10 )
  ###################
}
