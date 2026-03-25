
# load useful packages
pacman::p_load(here, tidyverse, sf, terra, patchwork)

# set CRS parameters
CRS_utm <- '+init=EPSG:20935'

# FUNCTIONS

# pretty Tukey plot
GGTukey<-function(Tukey){
  A<-require("tidyverse")
  if(A==TRUE){
    library(tidyverse)
  } else {
    install.packages("tidyverse")
    library(tidyverse)
  }
  B<-as.data.frame(Tukey[1])
  colnames(B)[2:3]<-c("min",
                      "max")
  C<-data.frame(id=row.names(B),
                min=B$min,
                max=B$max)
  D<-C%>%
    ggplot(aes(id))+
    geom_errorbar(aes(ymin=min,
                      ymax=max),
                  width = 0.2)+
    geom_hline(yintercept=0,
               color="red")+
    labs(x=NULL)+
    coord_flip()+
    theme(title=element_text(color="black",size=15),
          axis.text = element_text(color="black",size=10),
          axis.title = element_text(color="black",size=10),
          panel.grid=element_line(color="grey75"),
          axis.line=element_blank(),
          plot.background=element_rect(fill="white",color="white"),
          panel.background=element_rect(fill="white"),
          panel.border = element_rect(colour = "black", fill = NA,size=0.59),
          legend.key= element_rect(color="white",fill="white")
    )
  return(D)
}