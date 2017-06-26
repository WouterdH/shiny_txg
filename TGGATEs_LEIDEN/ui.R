library(shiny)

## Source seperate tabs
source('tabTXG.R')
source('tabEXPERIMENTREVIEW.R')
source('tabMODULES.R')
source('tabGENES.R')
source('tabCOMP_COR.R')
source('tabMODULE_COR.R')


## Data
#load('/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs/modulesLeidenParsed_with_mergedModules.RData')
#compounds <<- unique(datLeiden$human_hep_experiments_information$COMPOUND_NAME)
#  names(compounds) <<- compounds

#modules <<- unique(datLeiden$module_definitions$module)
#modules <<- modules[order(as.numeric(do.call(rbind, strsplit(modules, ':'))[, 2]))]
#  names(modules) <<- modules


## UI
shinyUI(fluidPage( 
  # Title
  titlePanel("TG-GATEs WGCNA Modules - Leiden - 0.1"),

 # 


  # Main 
  column(12,
    tabsetPanel( 
      tabTXG,
      tabEXPERIMENTREVIEW,
      tabMODULES,
      tabGENES,
      tabCOMP_COR,
      tabMODULE_COR
    )
  )#,
  
  #imageOutput('lacdr_image')

  # ,

  #column(2,
  #  imageOutput('lacdr_image', height = 128)
  #)

))
