## Perform GO term enrichment on parsed shiny data
## modulesLeidenParsed_with_mergedModules.RData

## Settings
  options(stringsAsFactors = FALSE)

  library(topGO)
  library(hgu133plus2hsentrezg.db)
  # library(parallel)
  
  
## Data
  load('/data/wouter/WGCNA/modulesLeidenParsed_with_mergedModules.RData')

  geneDF <- data.frame(module = datLeiden$module_definitions$module,
                       probe  = datLeiden$module_definitions$human.entrez.gene.id,
                       entrez = sapply(strsplit(datLeiden$module_definitions$human.entrez.gene.id, '_'), function(x) x[1]))


## Functions
  runGO <- function(dat, type = c('BP', 'MF', 'CC'), runFisherTest = TRUE, fixedModule = NULL, .debug = FALSE) {
    if(.debug) {
      dat    <- geneDF
      type   <- 'BP'
    }

    if(is.null(fixedModule)) {
      modules <- sort(unique(dat$module))

      tmp <- lapply(modules, function(.module) {
        cat(.module, '\n')

        entrezNames <- dat$probe
        myInterestingEntrez <- dat$probe[which(dat$module == .module)]
        
        entrezList <- factor(as.integer(entrezNames %in% myInterestingEntrez))
          names(entrezList) <- entrezNames

        GOdata <- new('topGOdata', 
                  ontology = type, 
                  allGenes = entrezList, 
                  annot = annFUN.db, affyLib = 'hgu133plus2hsentrezg.db',
                  nodeSize = 5)

        if(runFisherTest) resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") else resultFisher <- NA

        rslt <- list(GOdata = GOdata, resultFisher = resultFisher, entrezList = entrezList)

        return(rslt)
      })

      names(tmp) <- modules

    } else {
      .module <- fixedModule
      cat(.module, '\n')
      entrezNames <- dat$probe
      myInterestingEntrez <- dat$probe[which(dat$module == .module)]
        
      entrezList <- factor(as.integer(entrezNames %in% myInterestingEntrez))
        names(entrezList) <- entrezNames

      GOdata <- new('topGOdata', 
                ontology = type, 
                allGenes = entrezList, 
                annot = annFUN.db, affyLib = 'hgu133plus2hsentrezg.db',
                nodeSize = 5)

      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")  

      tmp <- list(GOdata = GOdata, resultFisher = resultFisher, entrezList = entrezList)
    }

    return(tmp)
  }


## Enrichizzle
  # BP
    library(GO.db)

    GO <- as.list(GOTERM)

    BP <- runGO(dat = geneDF, type = 'BP', runFisherTest = TRUE)
    GOannotation <- data.frame(t(sapply(GO[names(score(BP[[1]]$resultFisher))], function(x) c(GOID(x), Term(x), Ontology(x), Definition(x)) )))
      colnames(GOannotation) <- c('goID', 'goTERM', 'goONTOLOGY', 'goDEFINITION')
      for(i in 1:length(BP)) BP[[i]]$GOannotation <- GOannotation
     

    save(BP, file = '/data/wouter/WGCNA/BP_enrichment_20Jun2017.RData')

    BP_results <- lapply(BP, function(x) {
      tmp <- x[c('entrezList', 'resultFisher', 'GOannotation')]
      tmp$resultFisher <- score(tmp$resultFisher)

      tmp
    } )

    save(BP_results, file = '/data/wouter/WGCNA/BP_enrichment_20Jun2017_results.RData')

  # MF
    library(GO.db)

    GO <- as.list(GOTERM)

    MF <- runGO(dat = geneDF, type = 'MF', runFisherTest = TRUE)
    GOannotation <- data.frame(t(sapply(GO[names(score(MF[[1]]$resultFisher))], function(x) c(GOID(x), Term(x), Ontology(x), Definition(x)) )))
      colnames(GOannotation) <- c('goID', 'goTERM', 'goONTOLOGY', 'goDEFINITION')
      for(i in 1:length(MF)) MF[[i]]$GOannotation <- GOannotation
     

    save(MF, file = '/data/wouter/WGCNA/MF_enrichment_20Jun2017.RData')

    MF_results <- lapply(MF, function(x) {
      tmp <- x[c('entrezList', 'resultFisher', 'GOannotation')]
      tmp$resultFisher <- score(tmp$resultFisher)

      tmp
    } )

    save(MF_results, file = '/data/wouter/WGCNA/MF_enrichment_20Jun2017_results.RData')

  # CC
    library(GO.db)

    GO <- as.list(GOTERM)

    CC <- runGO(dat = geneDF, type = 'CC', runFisherTest = TRUE)
    GOannotation <- data.frame(t(sapply(GO[names(score(CC[[1]]$resultFisher))], function(x) c(GOID(x), Term(x), Ontology(x), Definition(x)) )))
      colnames(GOannotation) <- c('goID', 'goTERM', 'goONTOLOGY', 'goDEFINITION')
      for(i in 1:length(CC)) CC[[i]]$GOannotation <- GOannotation
     

    save(CC, file = '/data/wouter/WGCNA/CC_enrichment_20Jun2017.RData')

    CC_results <- lapply(CC, function(x) {
      tmp <- x[c('entrezList', 'resultFisher', 'GOannotation')]
      tmp$resultFisher <- score(tmp$resultFisher)

      tmp
    } )

    save(CC_results, file = '/data/wouter/WGCNA/CC_enrichment_20Jun2017_results.RData')


  
  



## Save data
