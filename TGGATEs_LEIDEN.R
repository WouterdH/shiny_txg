## TO DO
  # general
  	# 1) 	change DOSE_LEVEL & TIME radioButtons tot multiple select
  	# 2)  Histograms 

  # Module review tab
  	# 1) module enrichment
  	# 2) search by gene name

  # Gene tab
  	# 1) DOSE_LEVEL * TIME plots for all compounds for selected gene

  # Module correlation tab
    # 1) Add pearson R/R2

  # Experiment correlation tab
    # 1) Add pearson R/R2




## Settings
options(stringsAsFactors = FALSE)

library(shiny)
library(DT)

library(ape)    
library(igraph)
library(pheatmap)

library(ggplot2)
library(gplots)

library(grid)
library(RColorBrewer)

library(scales)


## Data
load('/Users/Wouter/STACK/ownCloud/TOX/gitProjects/shiny/TGGATEs_LEIDEN/modulesLeidenParsed_with_mergedModules.RData')

compounds <<- unique(datLeiden$human_hep_experiments_information$COMPOUND_NAME)
modules   <<- unique(datLeiden$module_definitions$module)
modules   <<- modules[order(as.numeric(do.call(rbind, strsplit(modules, ':'))[, 2]))]
genes     <<- sort(datLeiden$module_definitions$human.gene.symbol)


## launch webUI
  setwd('/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/')

  library(shiny)

  runApp('TGGATEs_LEIDEN', port = 8080)
