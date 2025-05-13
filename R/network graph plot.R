#Author: Pedro Rodriguez Cutillas
# May 2025


make_network_graph_bipartite <- function(df.network,graph.title=""){



  library(igraph)
  library(ggplot2)
  library(ggrepel)
#df.network <- x
#E(g)

  #a <- as.character(df.network$node1)
 # a <-ifelse(nchar(a) > 25, paste0(strtrim(a, 25), '.x'), a)
 # df.network$node1 <- a

 # a <- as.character(df.network$node1)
 # a <-ifelse(nchar(a) > 25, paste0(strtrim(a, 25), '.x'), a)
 # df.network$node2 <- a


  df.network$node1 <- as.factor(df.network$node1)
  df.network$node2 <- as.factor(df.network$node2)

g <- graph.data.frame(df.network, directed = F)
#if (is_bipartite(g)==F){
  g <- g %>%
    set_vertex_attr(name = "type", value = names(V(g)) %in% df.network$node1)
  g


  #}
#g <- make_bipartite_graph(edges=E(g),types = )
is_bipartite(g)
  #V(g)$name <- V(g)$name
 # V(g)$type <- bipartite_mapping(g)$type
  #V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
  #V(g)$shape <- ifelse(V(g)$type, "circle", "square")
  #E(g)$color <- "lightgray"
  V(g)$label.color <- ifelse(V(g)$type, "black", "purple4")
  ## V(g)$label.font <-  2
  #V(g)$label.cex <- ifelse(V(g)$type, 0.7, 1.2)
  ## V(g)$label.dist <-0
  #V(g)$frame.color <-  "gray"
  #V(g)$size <- ifelse(V(g)$type,3, 8)
  #V(g)$size <- degree(g)

  #fr.all <- layout_with_graphopt(g)

  #fr.all <- layout.circle(g)
  fr.all <- layout_as_bipartite(g)

  df <- data.frame(name= V(g)$name,
                   color=V(g)$label.color,
                   type=ifelse(V(g)$type,1,2)
                   #,
                   #font.size=V(g)$label.cex,
                   #degree=degree(g, mode="all")
  )

  fr.all.df <- as.data.frame(fr.all)  ## convert the layout to a data.frame
  fr.all.df$name <- V(g)$name# c(df.network$pathway,df.network$protein)# colnames(xx[,2:ncol(xx)])  ## add in the species codes

  gx <- get.data.frame(g)
  gx$from.x <- fr.all.df$V1[match(gx$from, fr.all.df$name)]  #  match the from locations from the node data.frame we previously connected
  gx$from.y <- fr.all.df$V2[match(gx$from, fr.all.df$name)]
  gx$to.x <- fr.all.df$V1[match(gx$to, fr.all.df$name)]  #  match the to locations from the node data.frame we previously connected
  gx$to.y <- fr.all.df$V2[match(gx$to, fr.all.df$name)]
  gx$name <- fr.all.df$name[match(gx$from, fr.all.df$name)]
  xxx <- merge.data.frame(fr.all.df ,df, by="name")
  xxx <- merge.data.frame(xxx, gx,by.x = "name",by.y = "from")

  pplot <- ggplot() +
    geom_segment(data=gx,
                 aes(x=from.x,xend = to.x, y=from.y,yend = to.y,
                     color=as.numeric(rvalue),
                     #size=-log10(as.numeric(qvalue)),
                     #alpha=as.numeric(rvalue),
                     alpha=-log10(as.numeric(qvalue))
                     ),size=1.5
                 )+
    scale_color_gradient2(high = "black",low="white", mid ="grey",
                          midpoint = median(as.numeric(gx$rvalue)))+ # add line type
    #scale_color_grey()+
    #geom_point(data=fr.all.df,aes(x=V1,y=V2),size=21,colour="black") +  # adds a black border around the nodes

    # geom_point(data=xxx,aes(x=V1,y=V2,shape=as.character(type),
    #color=zscore, size=zscore)) +
    #scale_linewidth(range = c(0,0.1))+
    geom_point(data=xxx,aes(x=V1,y=V2), size=4)+#,shape=as.character(type))) +

    #geom_label_repel(data=xxx,aes(x=V1,y=V2,label=ifelse(type==2, name.y,"")),
    #                        colour="purple4",size=3) + # add the node labels


    geom_text(data=xxx,aes(x=V1,y=V2,label=name),
               colour="purple4",size=3) + # add the node labels

    geom_text(data=xxx,aes(x=to.x,y=to.y,label=to),
              colour="seagreen",size=3)+

    #unique(xxx$type)


    #scale_colour_manual(values=c("2"="purple4","1"="black"))+  # add colour scaling for group membership
    #scale_size_manual(values=c(0.5,5))+
    #scale_color_gradient2(low="blue",high = "purple",mid = "white")+
    #scale_linetype_manual(values=c("0"="dashed","1"="solid"))+ # add linteyp scaling for within and between groups
    scale_x_continuous(expand=c(0,1))+  # expand the x limits
    scale_y_continuous(expand=c(0,1))+ # expand the y limits
    theme_bw()+  # use the ggplot black and white theme
    theme(
      axis.text.x = element_blank(),  # remove x-axis text
      axis.text.y = element_blank(), # remove y-axis text
      axis.ticks = element_blank(),  # remove axis ticks
      axis.title.x = element_blank(), # remove x-axis labels
      axis.title.y = element_blank(), # remove y-axis labels
      panel.background = element_blank(),
      panel.border =element_blank(),
      panel.grid.major = element_blank(),  #remove major-grid labels
      panel.grid.minor = element_blank(),  #remove minor-grid labels
      plot.background = element_blank()#,
      #legend.position = "none"
    )+
    labs(title=graph.title, color="R-value",alpha="-log10(qvalue)")+
    coord_flip()

  plot(pplot)

  return(pplot)


}
