#########################################
#### Working with STRUCTURE outuput SEE: https://www.royfrancis.com/pophelper/articles/index.html#input-files-1
#########################################
setwd("//wsl.localhost/Ubuntu/home/sksjds/structure_analysis/results/")

setwd("C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/structure_output_100000")
# install dependencies
#install.packages(c("ggplot2","gridExtra","label.switching","tidyr","remotes"),repos="https://cloud.r-project.org")
#remotes::install_github('royfrancis/pophelper')

library(pophelper)
library(tcltk)
library(gridExtra)


# read some structure files with confidence intervals
directory_path <- "C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/structure_output_100000"
sfiles1 <- list.files(path = directory_path, 
                      pattern = "output", #the information on the name of the files
                      full.names=TRUE)
slist1 <- readQ(files=sfiles1, readci=F) #there no interval confidence (ic)

# just checking the information content in the list
names(attributes(slist1[[1]])) 

#We can grab the individual labels 
slist <- readQ(files = sfiles1, indlabfromfile = T)
head(slist[[1]])

#A tabulated table of STRUCTURE runs
head(tabulateQ(slist))
# basic usage
tr1 <- tabulateQ(qlist=slist)
sr1 <- summariseQ(tr1)

#Evanno Method
#basic usage
evannoMethodStructure(data=sr1)
#another usage
em <- evannoMethodStructure(summariseQ(tabulateQ(slist)))

# to export a plot
evannoMethodStructure(data=sr1,exportplot=T,exportpath=getwd())
evannoMethodStructure(
  data = sr1,
  exportplot = T,
  exportpath = "C:/Users/givan/OneDrive/Área de Trabalho/Repositorios/diversity_LMH/"
)

# do not compute plot, only return results as table
em <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=F)

# to export plot and table
evannoMethodStructure(data=sr1,exportplot=T,writetable=T,na.rm=T,exportpath=getwd())

# returns both data and plot
p <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=T)

# to return only plot and save it
p <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=T,returndata=F)

r1 <- summariseQ(tabulateQ(slist))
p <- evannoMethodStructure(data=sr1,exportplot=F,returnplot=T,returndata=F,basesize=12,linesize=0.7)
grid.arrange(p)

#plotting the graph
slist <- readQ(sfiles1, indlabfromfile=T)

slist1 <- alignK(slist[c(4, #K=4
                         #5, #K=3
                         7, #K=2
                        # 9, #K=5
                         #10,
                         11)]) #K=10
p1 <- plotQ(slist1,imgoutput="join", # "sep" if only one K plot / "join" if two or more
            returnplot=T,
            exportplot=F,
            basesize=11, 
            ordergrp = T, 
            sortind = "all", 
            sharedindlab = F, 
            showindlab = F, 
            showlegend = F)
grid.arrange(p1$plot[[1]])
