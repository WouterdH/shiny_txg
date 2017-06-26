## Parse Leiden modules into shiny format

## Settings
  options(stringsAsFactors = FALSE)

  library(affy)
  library(reshape2)


## WGCNA modules Leiden
  load('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/getTG-GATES/moduleBuilds/MERGED_FC_ThresholdFDR_1e-03_customCDF_onlyEntrezProbes_2017-05-24.RData') # merged modules
  load('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/getTG-GATES/human/raw_rmaData_customCDF.RData')

## WGCNA parsed modules for shiny - Lilly template
  load('/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs/modulesLillyParsed.RData')

## Parse
  # human_hep_experiments_information
  human_hep_experiments_information <- pData(rmaData)[which(!duplicated(pData(rmaData)$COMPOUND_TIME_DOSE) & !grepl('Control', pData(rmaData)$COMPOUND_TIME_DOSE)), ]
    human_hep_experiments_information$COMPOUND_NAME <- toupper(human_hep_experiments_information$COMPOUND_NAME)
      human_hep_experiments_information$COMPOUND_NAME[grep('TNF', human_hep_experiments_information$COMPOUND_NAME)] <- 'TNF-ALPHA, RAT'
    human_hep_experiments_information$DOSE <- as.numeric(human_hep_experiments_information$DOSE) 
    human_hep_experiments_information$DOSE_LEVEL <- gsub('High', 'HI', gsub('Middle', 'MED', gsub('Low', 'LO', human_hep_experiments_information$DOSE_LEVEL)))
    human_hep_experiments_information$TIME <- as.numeric(unlist(strsplit(human_hep_experiments_information$SACRI_PERIOD, ' hr')))/24

  human_hep_experiments_information <- human_hep_experiments_information[, colnames(human_hep_experiments_information)[c(1:21, 67, 68, 69)]]
    rownames(human_hep_experiments_information) <- human_hep_experiments_information$COMPOUND_TIME_DOSE



  # module_definitions
  module_definitions <- data.frame(derivation.system     = 'TG-GATEs HHCYT',
                     horvath_module_merged = mergedModules$mergedModules$colors, 
                     human.entrez.gene.id  = colnames(mergedModules$datExpr),
                     human.gene.symbol     = sapply(colnames(mergedModules$datExpr), function(.prb) symbol$external_gene_name[which(symbol$entrezGene_probes == .prb)] ))

  moduleSizes <- as.data.frame(sort(table(module_definitions$horvath_module_merged), decreasing = TRUE))
    colnames(moduleSizes) <- c('horvath_module_merged', 'size')
  moduleSizes$module <- paste0('WGCNA|HHCYT|TGGATEs|LEIDEN:', rownames(moduleSizes))

  module_definitions$module <- sapply(module_definitions$horvath_module_merged, function(.module) moduleSizes$module[which(moduleSizes$horvath_module_merged == .module)] )
  
  tmp_expr <- mergedModules$datExpr_zScored
  tmp_egs  <- mergedModules$mergedModules$newMEs_scaled
    rownames(tmp_egs) <- rownames(tmp_expr)

  ProbeModuleCorMatrix <- cor(data.frame(tmp_expr, tmp_egs))

  module_definitions$correlationEG <- sapply(1:nrow(module_definitions), function(i) {
    .horvath_module_merged <- paste0('ME', module_definitions$horvath_module_merged[i])
    .prb <- paste0('X', module_definitions$human.entrez.gene.id[i])

    ProbeModuleCorMatrix[.horvath_module_merged, .prb]
  })


  # module_hubgenes
  module_hubgenes <- data.frame(module = sort(unique(module_definitions$module)))
  module_hubgene_cor <- do.call(rbind, lapply(module_hubgenes$module, function(.mod) {
        .horvath_module_merged <- paste0('ME', as.character(moduleSizes$horvath_module_merged[which(moduleSizes$module == .mod)]))
        
        .egs  <- mergedModules$mergedModules$newMEs_scaled[, .horvath_module_merged]
        .expr <- mergedModules$datExpr_zScored[, which(paste0('ME', mergedModules$mergedModules$colors) == .horvath_module_merged)]

        .hubgene <- sort(apply(.expr, 2, function(.gene_expr) cor(.gene_expr, .egs )), decreasing = TRUE)[1]
        
        return(data.frame(probe = names(.hubgene), eigengene_correl = .hubgene))
      }))

  module_hubgenes <- data.frame(module_hubgenes, module_hubgene_cor)
    module_hubgenes$symbol <- sapply(module_hubgenes$probe, function(.prb) symbol$external_gene_name[which(symbol$entrezGene_probes == .prb)] )

    rownames(module_hubgenes) <- 1:nrow(module_hubgenes)


  # TODO: module_preservation


  # scored_experiments_absEG  
  egs <- data.frame(mergedModules$mergedModules$newMEs_scaled)
    colnames(egs) <- sapply(gsub('ME', '', colnames(egs)), function(.horvath_module_merged) moduleSizes$module[which(moduleSizes$horvath_module_merged == .horvath_module_merged)] )
    egs$experiment <- rownames(mergedModules$datExpr_zScored)

  egs <- melt(egs, id.vars = 'experiment')
    colnames(egs) <- c('COMPOUND_TIME_DOSE', 'module', 'pc1_score')

  egs_meta <- data.frame(do.call(rbind, strsplit(egs$COMPOUND_TIME_DOSE, '_')))
    colnames(egs_meta) <- c('COMPOUND', 'TIME', 'DOSE_LEVEL')
    egs_meta$DOSE_LEVEL <- gsub('High', 'HI', gsub('Middle', 'MED', gsub('Low', 'LO', egs_meta$DOSE_LEVEL)))
    egs_meta$TIME <- as.numeric(gsub('hr', '', egs_meta$TIME))
    egs_meta$COMPOUND_NAME <- sapply(egs_meta$COMPOUND, function(.compound) unique(human_hep_experiments_information$COMPOUND_NAME[which(human_hep_experiments_information$COMPOUND.Abbr. == .compound)]) )

  scored_experiments_absEG <- data.frame(egs_meta, egs)


  # Module enrichment
  load('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/getTG-GATES/moduleBuilds/BP_enrichment_20Jun2017_results.RData')
  load('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/getTG-GATES/moduleBuilds/MF_enrichment_20Jun2017_results.RData')
  load('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/getTG-GATES/moduleBuilds/CC_enrichment_20Jun2017_results.RData')

  getModuleEnrichment <- function(module) {
    bp <- BP_results[[module]]
    mf <- MF_results[[module]]
    cc <- CC_results[[module]]

    bp$GOannotation <- bp$GOannotation[names(which(bp$resultFisher < 0.05)), ]
    bp$GOannotation$pval   <- bp$resultFisher[rownames(bp$GOannotation)]
    
    mf$GOannotation <- mf$GOannotation[names(which(mf$resultFisher < 0.05)), ]
    mf$GOannotation$pval   <- mf$resultFisher[rownames(mf$GOannotation)]

    cc$GOannotation <- cc$GOannotation[names(which(cc$resultFisher < 0.05)), ]
    cc$GOannotation$pval   <- cc$resultFisher[rownames(cc$GOannotation)]

    go <- rbind(bp$GOannotation, mf$GOannotation, cc$GOannotation)

    if(nrow(go) == 0) {
      go[1, ]   <- NA
      go$module <- NA
    } else {
      go$module <- module
    }
   
    return(go)
  }

  if(FALSE) load('/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs_LEIDEN/modulesLeidenParsed_with_mergedModules.RData')

  module_enrichment <- do.call(rbind, lapply(sort(unique(module_definitions$module)), getModuleEnrichment))
  module_enrichment <- module_enrichment[order(module_enrichment$pval), ]


## Save data
  datLeiden <- list(human_hep_experiments_information = human_hep_experiments_information,
            module_definitions = module_definitions,
            module_hubgenes = module_hubgenes,
            module_enrichment = module_enrichment,
            scored_experiments_absEG = scored_experiments_absEG)

  save(datLeiden, 
     file = '/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs_LEIDEN/modulesLeidenParsed.RData')

## Save data with all requirements for downstream goals.
  mergedModules <- mergedModules[c('mergedModules', 'datExpr', 'datExpr_zScored')]
  
  load('/Users/Wouter/stack/ownCloud/TOX/ExperimentData/getTG-GATES/human/limmaFits_customCDF_onlyEntrezProbes.RData')
  limmaFits_list <- lapply(DEGlist_allEntrez, function(x) x[sort(rownames(DEGlist_allEntrez[[1]])), ])

  limmaFits <- data.frame(limmaFits_list)
    colnames(limmaFits) <- gsub('w\\.', 'w__', gsub('e\\.', 'e__', gsub('h\\.', 'h__', colnames(limmaFits))))
  
  # TXG map phylogram
    library(igraph)
    library(ape)

    newME <- mergedModules$mergedModules$newMEs_scaled
      
    newMEd <- as.dist(1 - abs(cor(newME)))
    newMEcl <- hclust(newMEd, method = 'ward')
    newMEph <- as.phylo(newMEcl)

    if(FALSE) .newMEph <- newMEph 

    newMEel <- graph.edgelist(newMEph$edge)

      sampledEdgesForLengthening <- data.frame(start = c( 405, 442, 439, 421, 421, 458, 406, 438, 422, 414, 429, 430, 526, 417, 418, 447, 474, 448, 404, 409, 425, 426, 446, 408, 427, 428, 449, 450, 486, 415, 433, 462, 412, 419, 478, 423, 484, 432, 436, 456, 424, 444),
                                               stop  = c( 441, 510, 464, 440, 439, 521, 422, 482, 437, 429, 472, 453, 565, 459, 447, 473, 494, 452, 408, 425, 488, 446, 491, 410, 470, 449, 495, 485, 552, 433, 499, 475, 416, 467, 513, 483, 501, 436, 443, 489, 431, 456))
      sampledEdgesForLengthening <- apply(sampledEdgesForLengthening, 1, function(x) which(newMEph$edge[,1] == x[1] & newMEph$edge[,2] == x[2] ) )
   
    newMEph$edge.length[sampledEdgesForLengthening] <- newMEph$edge.length[sampledEdgesForLengthening] * 3
    
    newMEph$tip.label <- sapply(gsub('ME', '', newMEph$tip.label), function(.horvath_module_merged) datLeiden$module_definitions$module[which(datLeiden$module_definitions$horvath_module_merged == .horvath_module_merged)[1]])

    txgMap <- list(newMEd = newMEd, newMEcl = newMEcl, newMEph = newMEph, newMEel = newMEel)
   
  save(datLeiden, mergedModules, limmaFits, symbol, txgMap,
     file = '/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs_LEIDEN/modulesLeidenParsed_with_mergedModules.RData')

  









