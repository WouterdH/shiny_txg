

## Data
#load('/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs_LEIDEN/modulesLeidenParsed_with_mergedModules.RData')
#compounds <<- unique(datLeiden$human_hep_experiments_information$COMPOUND_NAME)
#  names(compounds) <<- compounds

#modules <<- unique(datLeiden$module_definitions$module)
#modules <<- modules[order(as.numeric(do.call(rbind, strsplit(modules, ':'))[, 2]))]
#  names(modules) <<- modules

## Functions
  # tabEXPERIMENTREVIEW.R
    plotDoseResponse <- function(dat) {
      axis_limits <- max(abs(dat$pc1_score))
      axis_limits <- c(-1 * axis_limits, axis_limits)

      P <- ggplot(dat, aes(x = DOSE_LEVEL, y = pc1_score)) +
                 theme(legend.position = 'none', 
                       axis.title = element_text(size = 14),
                       axis.text = element_text(size = 12),
                       strip.text.x = element_text(as.character(unique(dat$TIME)))) +
                 geom_hline(yintercept = 0, colour = 'grey70') +
                 xlab('') + ylab('Eigengene score') +
                 ylim(axis_limits) +
                 geom_line(aes(group = module, colour = module), size = 1.25, alpha = 0.4) +
                 geom_point(aes(fill = module), colour = 'black', pch = 21, size = 3)

      if(nrow(dat) > 0) P <- P + scale_x_discrete(limits = c('LO', 'MED', 'HI'))

      print(P)
    }  

    plotTimeResponse <- function(dat) {
      axis_limits <- max(abs(dat$pc1_score))
      axis_limits <- c(-1 * axis_limits, axis_limits)

      P <- ggplot(dat, aes(x = TIME, y = pc1_score)) +
                 theme(legend.position = 'none', 
                       axis.title = element_text(size = 14),
                       axis.text = element_text(size = 12),
                       strip.text.x = element_text()) +
                 geom_hline(yintercept = 0, colour = 'grey70') +
                 xlab('') + ylab('Eigengene score') +
                 ylim(axis_limits) +
                 geom_line(aes(group = module, colour = module), size = 1.25, alpha = 0.4) +
                 geom_point(aes(fill = module), colour = 'black', pch = 21, size = 3)

      if(nrow(dat) > 0) P <- P + scale_x_discrete(limits = c('2', '8', '24'))

      print(P)
    }  

    getEGS <- function(dat, mod, time, dose) {
      if(FALSE) {
        dat <- EXPERIMENTREVIEW_doseResponse_dat
        mod <- 'WGCNA|HHCYT|TGGates:301'
        time <- '2'
        dose <- 'HI'
      }

      ind <- which(dat$module == mod & dat$TIME == time & dat$DOSE_LEVEL == dose)

      if(length(ind) == 0) {
        return(NA)
      } else {
        return(dat$pc1_score[ind])
      }
    }

  # tabMODULES.R
    plotModuleCircleGraph <- function(fc_dat) {
      if(FALSE) {
        fc_dat <- MODULES_fc_dat

        input <- list(MODULES_module = 'WGCNA|HHCYT|TGGATEs|LEIDEN:139',
                              MODULES_time = 24,
                              MODULES_dose = 'HI',
                              MODULES_compound = 'DIETHYL MALEATE')
      }

      fc_color_range <- c(-1 * max(abs(as.numeric(fc_dat$logFC))), max(abs(as.numeric(fc_dat$logFC))))
      .graphVertexColorsPalette <- data.frame(logFC = rev(seq(fc_color_range[1], fc_color_range[2], length.out = 10)),
                                       color = colorRampPalette(brewer.pal(n = 11, name = "RdYlGn"))(10))

      fc_dat$color <- .graphVertexColorsPalette[sapply(as.numeric(fc_dat$logFC), function(.fc) which(.graphVertexColorsPalette$logFC < .fc)[1] ), 'color']
      .vertexColors <- fc_dat$color
        names(.vertexColors) <- fc_dat$probe

      .prbs <- fc_dat$probe
      .prbCorMatrix <- cor(mergedModules$datExpr_zScored[, .prbs])
        diag(.prbCorMatrix) <- 0

      .graph <- graph.adjacency(.prbCorMatrix, weighted = TRUE, mode = "lower")
        V(.graph)$color        <- fc_dat$color
        V(.graph)$size         <- 50 * abs(as.numeric(fc_dat$corEG))
        V(.graph)$frame.color  <- 'white'
        V(.graph)$label        <- sapply(names(V(.graph)), function(.prb) fc_dat$gene[which(fc_dat$probe == .prb)])
        V(.graph)$label.family <- 'mono'
        V(.graph)$label.font   <- 2
        V(.graph)$label.color  <- 'black'
        V(.graph)$label.dist   <- 2.25
        V(.graph)$label.degree <- 0


        E(.graph)$width <- 10 * abs(E(.graph)$weight)
        E(.graph)$color <- 'grey80'


        plot(.graph, 
             layout = layout_in_circle(.graph))
    }

## Define server logic required to plot data
shinyServer(function(input, output) {
  ## tabSetPanel - ui.R
    output$lacdr_image <- renderImage({
      list(src = '/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs_LEIDEN/images/logoLACDR.JPG',
         #contentType = 'image/png',
         width = 180,
         height = 55,
         alt = "LACDR")
    }, deleteFile = FALSE)

  ## UI - tabCOMP_COR.R
    output$selectedCompund_x <- renderText(input$compcor_compound_x)
    output$selectedCompund_y <- renderText(input$compcor_compound_y)

    output$compound_cor_plot <- renderPlot({
        if(FALSE) input <- list(compound_x = 'ACETAMINOPHEN', compound_y = 'ACARBOSE', 
                                x_time = 24, y_time = 24,
                                x_dose = 'HI', y_dose = 'HI')

        # Filter based on radioSelection
          x_EGS <- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$COMPOUND_NAME == input$compcor_compound_x & 
                                                           datLeiden$scored_experiments_absEG$DOSE_LEVEL == input$compcor_x_dose &
                                                           datLeiden$scored_experiments_absEG$TIME == input$compcor_x_time), ]
            rownames(x_EGS) <- x_EGS$module
          y_EGS <- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$COMPOUND_NAME == input$compcor_compound_y & 
                                                           datLeiden$scored_experiments_absEG$DOSE_LEVEL == input$compcor_y_dose &
                                                           datLeiden$scored_experiments_absEG$TIME == input$compcor_y_time), ]
            rownames(y_EGS) <- y_EGS$module

          overlappingModules <- intersect(rownames(x_EGS), rownames(y_EGS))

          compound_cor_plot_dat <- data.frame(module  = x_EGS[overlappingModules, 'module'],
                                                 TIME  = x_EGS[overlappingModules, 'TIME'],
                                                 DOSE  = x_EGS[overlappingModules, 'DOSE_LEVEL'],
                                                 x_EGS = x_EGS[overlappingModules, 'pc1_score'],
                                                 y_EGS = y_EGS[overlappingModules, 'pc1_score'])

          compound_cor_plot_dat$mean_EGS <- rowMeans(compound_cor_plot_dat[, c('x_EGS', 'y_EGS')])
          
          # egs_color_range <- c(-1.1 * max(abs(compound_cor_plot_dat$mean_EGS)), 1.1 * max(abs(compound_cor_plot_dat$mean_EGS)))
          # egs_color_palette <- data.frame(EGS = rev(seq(egs_color_range[1], egs_color_range[2], length.out = 100)),
          #                                 color = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(100))

          # compound_cor_plot_dat$color <- egs_color_palette[sapply(compound_cor_plot_dat$mean_EGS, function(.egs) which(egs_color_palette$EGS < .egs)[1] ), 'color']

        # Aesthetics
          axis_limits <- max(abs(c(compound_cor_plot_dat$x_EGS, compound_cor_plot_dat$y_EGS)))
          axis_limits <- c(-1 * axis_limits, axis_limits)

          myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
          sc <- scale_colour_gradientn(colours = myPalette(399), limits = axis_limits)
          sf <- scale_fill_gradientn(colours = myPalette(399), limits = axis_limits)

          compound_cor_plot_dat <<- compound_cor_plot_dat

    	    ggplot(compound_cor_plot_dat, aes(x = x_EGS, y = y_EGS)) + 
    	    	xlab(input$compcor_compound_x) + ylab(input$compcor_compound_y) +
    	    	xlim(axis_limits) + ylim(axis_limits) +
            stat_smooth(method = 'lm', fullrange = TRUE) + 
            geom_point(aes(fill = mean_EGS), colour = 'black', pch = 21, size = 5, alpha = 0.8) +
    	    	theme(legend.position = 'none', 
                  axis.title = element_text(size = 14),
                  axis.text = element_text(size = 12)) +
           sc + sf

    })

    output$brush_info_compcor <- renderDataTable({
      COMP_COR_brushedPoints <- brushedPoints(compound_cor_plot_dat[, c('module', 'x_EGS', 'y_EGS')], 
                                              input$compound_cor_plot_brush)

      COMP_COR_brushedPoints$x_EGS <- round(COMP_COR_brushedPoints$x_EGS, digits = 2)
      COMP_COR_brushedPoints$y_EGS <- round(COMP_COR_brushedPoints$y_EGS, digits = 2)

      x_EGS_colname <- datLeiden$human_hep_experiments_information$COMPOUND.Abbr.[which(datLeiden$human_hep_experiments_information$COMPOUND_NAME == input$compcor_compound_x)[1]]
      y_EGS_colname <- datLeiden$human_hep_experiments_information$COMPOUND.Abbr.[which(datLeiden$human_hep_experiments_information$COMPOUND_NAME == input$compcor_compound_y)[1]]


      colnames(COMP_COR_brushedPoints) <- c('module', x_EGS_colname, y_EGS_colname)

      COMP_COR_brushedPoints

    }, rownames = FALSE)

  ## UI - tabEXPERIMENTREVIEW.R
    output$EXPERIMENTREVIEW_doseResponse <- renderPlot({
      EXPERIMENTREVIEW_doseResponse_dat <<- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$COMPOUND_NAME == input$EXPERIMENTREVIEW_compound), ]
        EXPERIMENTREVIEW_doseResponse_dat$DOSE_LEVEL <<- factor(EXPERIMENTREVIEW_doseResponse_dat$DOSE_LEVEL, ordered = TRUE, levels = c('LO', 'MED', 'HI'))
        EXPERIMENTREVIEW_doseResponse_dat$TIME <<- factor(EXPERIMENTREVIEW_doseResponse_dat$TIME, ordered = TRUE, levels = c('2', '8', '24'))

      axis_limits <- max(abs(EXPERIMENTREVIEW_doseResponse_dat$pc1_score))
      axis_limits <- c(-1 * axis_limits, axis_limits)

      EXPERIMENTREVIEW_doseResponse_dat$selected <<- EXPERIMENTREVIEW_doseResponse_dat$module %in% unique(EXPERIMENTREVIEW_doseResponse_dat$module[which(EXPERIMENTREVIEW_doseResponse_dat$pc1_score <= input$EXPERIMENTREVIEW_egs_thresholds[2] & EXPERIMENTREVIEW_doseResponse_dat$pc1_score >= input$EXPERIMENTREVIEW_egs_thresholds[1])])
      EXPERIMENTREVIEW_doseResponse_dat <<- EXPERIMENTREVIEW_doseResponse_dat[EXPERIMENTREVIEW_doseResponse_dat$selected, ]

      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
      sc <- scale_colour_gradientn(colours = myPalette(399), limits = axis_limits)
      sf <- scale_fill_gradientn(colours = myPalette(399), limits = axis_limits)
      
      P <- ggplot(EXPERIMENTREVIEW_doseResponse_dat, aes(x = DOSE_LEVEL, y = pc1_score)) +
                 theme(legend.position = 'none', 
                       axis.title = element_text(size = 14),
                       axis.text = element_text(size = 12)) +
                 geom_hline(yintercept = 0, colour = 'grey70') +
                 xlab('') + ylab('Eigengene score') +
                 ylim(axis_limits) +
                 geom_line(aes(group = module, colour = pc1_score), size = 2, alpha = 0.4) +
                 geom_point(aes(fill = pc1_score), colour = 'black', pch = 21, size = 5, alpha = 0.8) +
                 facet_grid(. ~ TIME) +
                 sc + sf

      if(nrow(EXPERIMENTREVIEW_doseResponse_dat) > 0) P <- P + scale_x_discrete(limits = c('LO', 'MED', 'HI'))

      P

    })

    output$EXPERIMENTREVIEW_timeResponse <- renderPlot({
      EXPERIMENTREVIEW_doseResponse_dat <<- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$COMPOUND_NAME == input$EXPERIMENTREVIEW_compound), ]
        EXPERIMENTREVIEW_doseResponse_dat$DOSE_LEVEL <<- factor(EXPERIMENTREVIEW_doseResponse_dat$DOSE_LEVEL, ordered = TRUE, levels = c('LO', 'MED', 'HI'))
        EXPERIMENTREVIEW_doseResponse_dat$TIME <<- factor(EXPERIMENTREVIEW_doseResponse_dat$TIME, ordered = TRUE, levels = c('2', '8', '24'))

      axis_limits <- max(abs(EXPERIMENTREVIEW_doseResponse_dat$pc1_score))
      axis_limits <- c(-1 * axis_limits, axis_limits)

      EXPERIMENTREVIEW_doseResponse_dat$selected <<- EXPERIMENTREVIEW_doseResponse_dat$module %in% unique(EXPERIMENTREVIEW_doseResponse_dat$module[which(EXPERIMENTREVIEW_doseResponse_dat$pc1_score <= input$EXPERIMENTREVIEW_egs_thresholds[2] & EXPERIMENTREVIEW_doseResponse_dat$pc1_score >= input$EXPERIMENTREVIEW_egs_thresholds[1])])
      EXPERIMENTREVIEW_doseResponse_dat <<- EXPERIMENTREVIEW_doseResponse_dat[EXPERIMENTREVIEW_doseResponse_dat$selected, ]

      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
      sc <- scale_colour_gradientn(colours = myPalette(399), limits = axis_limits)
      sf <- scale_fill_gradientn(colours = myPalette(399), limits = axis_limits)

      P <- ggplot(EXPERIMENTREVIEW_doseResponse_dat, aes(x = TIME, y = pc1_score)) +
                 theme(legend.position = 'none', 
                       axis.title = element_text(size = 14),
                       axis.text = element_text(size = 12)) +
                 geom_hline(yintercept = 0, colour = 'grey70') +
                 xlab('') + ylab('Eigengene score') +
                 ylim(axis_limits) +
                 geom_line(aes(group = module, colour = pc1_score), size = 2, alpha = 0.4) +
                 geom_point(aes(fill = pc1_score), colour = 'black', pch = 21, size = 5) +
                 facet_grid(. ~ DOSE_LEVEL) +
                 sc + sf

      if(nrow(EXPERIMENTREVIEW_doseResponse_dat) > 0) P <- P + scale_x_discrete(limits = c('2', '8', '24'))

      P

    })

    output$EXPERIMENTREVIEW_brush_info <- renderDataTable({
      EXPERIMENTREVIEW_brushedPoints <- brushedPoints(EXPERIMENTREVIEW_doseResponse_dat[, c('module', 'pc1_score', 'DOSE_LEVEL', 'TIME')], 
                                              input$EXPERIMENTREVIEW_doseResponse_brush)


      EXPERIMENTREVIEW_brushedPoints$pc1_score <- round(EXPERIMENTREVIEW_brushedPoints$pc1_score, digits = 2)
     
      EXPERIMENTREVIEW_brushedPoints <- EXPERIMENTREVIEW_brushedPoints[, c('module', 'DOSE_LEVEL', 'TIME', 'pc1_score')]
        colnames(EXPERIMENTREVIEW_brushedPoints) <- c('module', 'dose', 'time', 'EGS')

      EXPERIMENTREVIEW_brushedPoints
      
    }, rownames = FALSE)

    output$EXPERIMENTREVIEW_heatmap <- renderPlot({
      if(length(input$EXPERIMENTREVIEW_doseResponse_brush) > 0) {
        EXPERIMENTREVIEW_brushedPoints <- brushedPoints(EXPERIMENTREVIEW_doseResponse_dat[, c('module', 'pc1_score', 'DOSE_LEVEL', 'TIME')], 
                                              input$EXPERIMENTREVIEW_doseResponse_brush)


        EXPERIMENTREVIEW_brushedPoints$pc1_score <- round(EXPERIMENTREVIEW_brushedPoints$pc1_score, digits = 2)
       
        EXPERIMENTREVIEW_brushedPoints <- EXPERIMENTREVIEW_brushedPoints[, c('module', 'DOSE_LEVEL', 'TIME', 'pc1_score')]
          colnames(EXPERIMENTREVIEW_brushedPoints) <- c('module', 'dose', 'time', 'EGS')

        .selectedExperiment_time <- as.character(EXPERIMENTREVIEW_brushedPoints$time[1])
        .selectedExperiment_dose <- gsub('LO', 'Low', gsub('MED', 'Middle', gsub('HI', 'High', EXPERIMENTREVIEW_brushedPoints$dose[1])))
        .selectedModule_definitions <- datLeiden$module_definitions[unlist(sapply(as.character(EXPERIMENTREVIEW_brushedPoints$module), function(.mod) which(datLeiden$module_definitions$module == .mod) )), c('human.entrez.gene.id', 'human.gene.symbol', 'module')]
        
        .selected_probes      <- unique(.selectedModule_definitions$human.entrez.gene.id)
        .selected_fcCols      <- grep('logFC', grep(unique(datLeiden$human_hep_experiments_information$COMPOUND.Abbr.[which(datLeiden$human_hep_experiments_information$COMPOUND_NAME == input$EXPERIMENTREVIEW_compound)]), colnames(limmaFits) , value = TRUE), value = TRUE)
         
        .selectedExperiment_heatmap_dat <<- limmaFits[.selected_probes, .selected_fcCols]                                      

        .heatmap_color_breaks <- seq(-1 * ceiling(max(abs(range(.selectedExperiment_heatmap_dat)))),
                                      1 * ceiling(max(abs(range(.selectedExperiment_heatmap_dat)))),
                                      length.out = 101)

        .heatmap_row_annotation <- .selectedModule_definitions[, 'module', drop = FALSE]
        .heatmap_row_annotation$module <- do.call(rbind, strsplit(.heatmap_row_annotation$module, '\\|'))[, 4]


        .heatmap_col_labels <- toupper(apply(do.call(rbind, strsplit(colnames(.selectedExperiment_heatmap_dat), '_'))[, 1:3], 1, function(x) paste(x, collapse = ' - ')))

        ## Pheatmap draw colnames adjustment
        draw_colnames_45 <<- function (coln, gaps, ...) {
          coord = pheatmap:::find_coordinates(length(coln), gaps)
          x = coord$coord - 0.5 * coord$size
          res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 2, hjust = 1, rot = 45, gp = gpar(...))
          return(res)
        }

        assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))

        if(nrow(.selectedExperiment_heatmap_dat) > 30) .heatmap_show_rownames <- FALSE else .heatmap_show_rownames <- TRUE
        if(.heatmap_show_rownames) {
          .heatmap_row_labels <- datLeiden$module_definitions[.selected_probes, 'human.gene.symbol']
        } else {
          .heatmap_row_labels <- rep('', nrow(.selectedExperiment_heatmap_dat))
        }

        pheatmap(.selectedExperiment_heatmap_dat, 
                 clustering_distance_rows = 'euclidean',
                 clustering_distance_cols = 'correlation',
                 breaks = .heatmap_color_breaks,
                 scale = 'none',
                 annotation_row = .heatmap_row_annotation, 
                 color = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(101),
                 show_rownames = .heatmap_show_rownames,
                 labels_col = .heatmap_col_labels,
                 labels_row = .heatmap_row_labels) 
      } else {
        return()
      }
    })

  ## UI - tabMODULE_COR.R
    output$selectedModule_x <- renderText(input$modcor_module_x)
    output$selectedModule_y <- renderText(input$modcor_module_y)

    output$module_cor_plot <- renderPlot({
        if(FALSE) input <- list(module_x = 'WGCNA|HHCYT|TGGATEs|LEIDEN:375', module_y = 'WGCNA|HHCYT|TGGATEs|LEIDEN:275', 
                                x_time = 24, y_time = 24,
                                x_dose = 'HI', y_dose = 'MED')

        # Filter based on radioSelection
          x_EGS <- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$module == input$modcor_module_x & 
                                                           datLeiden$scored_experiments_absEG$DOSE_LEVEL == input$modcor_x_dose &
                                                           datLeiden$scored_experiments_absEG$TIME == input$modcor_x_time), ]
            rownames(x_EGS) <- x_EGS$COMPOUND
          y_EGS <- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$module == input$modcor_module_y & 
                                                           datLeiden$scored_experiments_absEG$DOSE_LEVEL == input$modcor_y_dose &
                                                           datLeiden$scored_experiments_absEG$TIME == input$modcor_y_time), ]
            rownames(y_EGS) <- y_EGS$COMPOUND

          overlappingCompounds <- intersect(rownames(x_EGS), rownames(y_EGS))

          module_cor_plot_dat <- data.frame(compound  = x_EGS[overlappingCompounds, 'COMPOUND_NAME'],
                                             TIME  = x_EGS[overlappingCompounds, 'TIME'],
                                             DOSE  = x_EGS[overlappingCompounds, 'DOSE_LEVEL'],
                                             x_EGS = x_EGS[overlappingCompounds, 'pc1_score'],
                                             y_EGS = y_EGS[overlappingCompounds, 'pc1_score'])

          module_cor_plot_dat$mean_EGS <- rowMeans(module_cor_plot_dat[, c('x_EGS', 'y_EGS')])

        # Aesthetics
          axis_limits <- max(abs(c(module_cor_plot_dat$x_EGS, module_cor_plot_dat$y_EGS)))
          axis_limits <- c(-1 * axis_limits, axis_limits)

          myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
          sc <- scale_colour_gradientn(colours = myPalette(160), limits = axis_limits)
          sf <- scale_fill_gradientn(colours = myPalette(160), limits = axis_limits)

          module_cor_plot_dat <<- module_cor_plot_dat

          ggplot(module_cor_plot_dat, aes(x = x_EGS, y = y_EGS)) + 
            xlab(input$modcor_module_x) + ylab(input$modcor_module_y) +
            xlim(axis_limits) + ylim(axis_limits) +
            stat_smooth(method = 'lm', fullrange = TRUE) + 
            geom_point(aes(fill = mean_EGS), colour = 'black', pch = 21, size = 5, alpha = 0.8) +
            theme(legend.position = 'none', 
                  axis.title = element_text(size = 14),
                  axis.text = element_text(size = 12))  +
            sc + sf
    })

    output$brush_info_modcor <- renderDataTable({
      MODCOR_brushedPoints <- brushedPoints(module_cor_plot_dat[, c('compound', 'x_EGS', 'y_EGS')], 
                                  input$module_cor_plot_brush)

      MODCOR_brushedPoints$x_EGS <- round(MODCOR_brushedPoints$x_EGS, digits = 2)
      MODCOR_brushedPoints$y_EGS <- round(MODCOR_brushedPoints$y_EGS, digits = 2)

      x_EGS_colname <- gsub('WGCNA\\|', '', input$modcor_module_x)
      y_EGS_colname <- gsub('WGCNA\\|', '', input$modcor_module_y)


      colnames(MODCOR_brushedPoints) <- c('compound', x_EGS_colname, y_EGS_colname)

      MODCOR_brushedPoints

      
    }, rownames = FALSE)

  ## UI - tabMODULES.R
    output$MODULES_doseResponse_plot <- renderPlot({
      .modules_module   <- input$MODULES_module

      .modules_dat <- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$module == .modules_module), ]
        .modules_dat$DOSE_LEVEL <- factor(.modules_dat$DOSE_LEVEL, ordered = TRUE, levels = c('LO', 'MED', 'HI'))
        .modules_dat$TIME <- factor(.modules_dat$TIME, ordered = TRUE, levels = c('2', '8', '24'))

      axis_limits <- max(abs(.modules_dat$pc1_score))
      axis_limits <- c(-1 * axis_limits, axis_limits)

      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
      sc <- scale_colour_gradientn(colours = myPalette(158), limits = axis_limits)
      sf <- scale_fill_gradientn(colours = myPalette(158), limits = axis_limits)

      .modules_dat$selected <- .modules_dat$COMPOUND_NAME %in% unique(.modules_dat$COMPOUND_NAME[which(.modules_dat$pc1_score >= input$MODULES_egs_thresholds[1] & .modules_dat$pc1_score <= input$MODULES_egs_thresholds[2])])
      MODULES_doseResponse_dat <<- .modules_dat[.modules_dat$selected, ]

      P <- ggplot(MODULES_doseResponse_dat, aes(x = DOSE_LEVEL, y = pc1_score)) +
                 theme(legend.position = 'none', 
                       axis.title = element_text(size = 14),
                       axis.text = element_text(size = 12)) +
                 geom_hline(yintercept = 0, colour = 'grey70') +
                 xlab('') + ylab('Eigengene score') +
                 ylim(axis_limits) +
                 geom_line(aes(group = COMPOUND_NAME, colour = pc1_score), size = 2, alpha = 0.4) +
                 geom_point(aes(fill = pc1_score), colour = 'black', pch = 21, size = 5, alpha = 0.8) +
                 facet_grid(. ~ TIME) +
                 sc + sf

      P           
    })

    output$MODULES_timeResponse_plot <- renderPlot({
      .modules_module   <- input$MODULES_module

      .modules_dat <- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$module == .modules_module), ]
        .modules_dat$DOSE_LEVEL <- factor(.modules_dat$DOSE_LEVEL, ordered = TRUE, levels = c('LO', 'MED', 'HI'))
        .modules_dat$TIME <- factor(.modules_dat$TIME, ordered = TRUE, levels = c('2', '8', '24'))

      axis_limits <- max(abs(.modules_dat$pc1_score))
      axis_limits <- c(-1 * axis_limits, axis_limits)

      myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
      sc <- scale_colour_gradientn(colours = myPalette(158), limits = axis_limits)
      sf <- scale_fill_gradientn(colours = myPalette(158), limits = axis_limits)

      .modules_dat$selected <- .modules_dat$COMPOUND_NAME %in% unique(.modules_dat$COMPOUND_NAME[which(.modules_dat$pc1_score >= input$MODULES_egs_thresholds[1] & .modules_dat$pc1_score <= input$MODULES_egs_thresholds[2])])
      MODULES_doseResponse_dat <<- .modules_dat[.modules_dat$selected, ]

      P <- ggplot(MODULES_doseResponse_dat, aes(x = TIME, y = pc1_score)) +
                 theme(legend.position = 'none', 
                       axis.title = element_text(size = 14),
                       axis.text = element_text(size = 12)) +
                 geom_hline(yintercept = 0, colour = 'grey70') +
                 xlab('') + ylab('Eigengene score') +
                 ylim(axis_limits) +
                 geom_line(aes(group = COMPOUND_NAME, colour = pc1_score), size = 2, alpha = 0.4) +
                 geom_point(aes(fill = pc1_score), colour = 'black', pch = 21, size = 5, alpha = 0.8) +
                 facet_grid(. ~ DOSE_LEVEL) +
                 sc + sf

      P           
    })

    output$MODULES_brush_info <- renderDataTable({
      MODULES_brushedPoints <- brushedPoints(MODULES_doseResponse_dat[, c('COMPOUND_NAME', 'pc1_score', 'DOSE_LEVEL', 'TIME')], 
                                              input$MODULES_doseResponse_plot_brush)


      MODULES_brushedPoints$pc1_score <- round(MODULES_brushedPoints$pc1_score, digits = 2)
     
      MODULES_brushedPoints <- MODULES_brushedPoints[, c('COMPOUND_NAME', 'DOSE_LEVEL', 'TIME', 'pc1_score')]
        colnames(MODULES_brushedPoints) <- c('compound', 'dose', 'time', 'EGS')

      MODULES_brushedPoints <<- MODULES_brushedPoints
      
    }, rownames = FALSE, selection = 'single')

    output$MODULES_module_dat <- renderDataTable({
      cat('Parsing MODULES_module_dat ...\n')

      if(length(input$MODULES_brush_info_rows_selected) > 0) {
          .modules_module   <- input$MODULES_module

          .modules_selectedRows <- input$MODULES_brush_info_rows_selected

          .modules_compound <- MODULES_brushedPoints[.modules_selectedRows, 'compound']
          .modules_time     <- MODULES_brushedPoints[.modules_selectedRows, 'time'] 
          .modules_dose     <- MODULES_brushedPoints[.modules_selectedRows, 'dose']

          .tmp_human_hep_experiments_information <- datLeiden$human_hep_experiments_information[, c('COMPOUND.Abbr.', 'COMPOUND_NAME', 'SACRI_PERIOD', 'DOSE_LEVEL', 'COMPOUND_TIME_DOSE')]
          .tmp_module_definitions <- datLeiden$module_definitions[, c('module', 'human.entrez.gene.id', 'human.gene.symbol', 'correlationEG')]

          .modules_module_genes <-  .tmp_module_definitions[which(.tmp_module_definitions$module == .modules_module), c('human.entrez.gene.id', 'human.gene.symbol')]
            colnames(.modules_module_genes) <- c('probe', 'gene')
          .modules_compound_abbr <- unique(.tmp_human_hep_experiments_information$COMPOUND.Abbr.[which(.tmp_human_hep_experiments_information$COMPOUND_NAME == .modules_compound)])
          .modules_time_char <- paste0(.modules_time, ' hr')

          expInd <- which(.tmp_human_hep_experiments_information$SACRI_PERIOD == .modules_time_char &
                           .tmp_human_hep_experiments_information$COMPOUND.Abbr == .modules_compound_abbr &
                           .tmp_human_hep_experiments_information$DOSE_LEVEL == .modules_dose)

          expID <- .tmp_human_hep_experiments_information$COMPOUND_TIME_DOSE[expInd]

          .modules_fc_dat <- data.frame(.modules_module_genes, limmaFits[.modules_module_genes$probe, grep(expID, colnames(limmaFits))])
          .modules_fc_dat <- .modules_fc_dat[, c('probe', 'gene',
                                                  grep('logFC', colnames(.modules_fc_dat), value = TRUE), 
                                                  grep('P.Value', colnames(.modules_fc_dat), value = TRUE),
                                                  grep('adj', colnames(.modules_fc_dat), value = TRUE))]
            colnames(.modules_fc_dat) <- c('probe', 'gene', 'logFC', 'P', 'FDR')

          .modules_fc_dat <- .modules_fc_dat[order(.modules_fc_dat$P), ]
            rownames(.modules_fc_dat) <- 1:nrow(.modules_fc_dat)

          .modules_fc_dat$corEGS <- round(.tmp_module_definitions[.modules_fc_dat$probe, 'correlationEG'], digits = 2)
          .modules_fc_dat$logFC <- round(.modules_fc_dat$logFC, digits = 2)

          #.modules_fc_dat$P <- format(.modules_fc_dat$P, digits = 3, scientific = TRUE)
          #.modules_fc_dat$FDR <- format(.modules_fc_dat$FDR, digits = 3, scientific = TRUE)

          .modules_fc_dat <- .modules_fc_dat[, -which(colnames(.modules_fc_dat) == 'probe')]

          datatable(.modules_fc_dat, rownames = FALSE) %>% formatCurrency(columns = c('logFC'), currency = '') %>% formatSignif(columns = c('P', 'FDR'), digits = 2)
      } else {
        return()
      }
     
    })

    output$MODULES_gene_plot <- renderPlot({
      cat('Parsing MODULES_gene_plot ...\n')

      if(length(input$MODULES_brush_info_rows_selected) > 0) {
        .modules_module   <- input$MODULES_module

        .modules_selectedRows <- input$MODULES_brush_info_rows_selected

        .modules_compound <- MODULES_brushedPoints[.modules_selectedRows, 'compound']
        .modules_time     <- MODULES_brushedPoints[.modules_selectedRows, 'time'] 
        .modules_dose     <- MODULES_brushedPoints[.modules_selectedRows, 'dose']

        .tmp_human_hep_experiments_information <- datLeiden$human_hep_experiments_information[, c('COMPOUND.Abbr.', 'COMPOUND_NAME', 'SACRI_PERIOD', 'DOSE_LEVEL', 'COMPOUND_TIME_DOSE')]
        .tmp_module_definitions <- datLeiden$module_definitions[, c('module', 'human.entrez.gene.id', 'human.gene.symbol', 'correlationEG')]

        .modules_module_genes <-  .tmp_module_definitions[which(.tmp_module_definitions$module == .modules_module), c('human.entrez.gene.id', 'human.gene.symbol')]
          colnames(.modules_module_genes) <- c('probe', 'gene')
        .modules_compound_abbr <- unique(.tmp_human_hep_experiments_information$COMPOUND.Abbr.[which(.tmp_human_hep_experiments_information$COMPOUND_NAME == .modules_compound)])
        .modules_time_char <- paste0(.modules_time, ' hr')

        expInd <- which(.tmp_human_hep_experiments_information$SACRI_PERIOD == .modules_time_char &
                         .tmp_human_hep_experiments_information$COMPOUND.Abbr == .modules_compound_abbr &
                         .tmp_human_hep_experiments_information$DOSE_LEVEL == .modules_dose)

        expID <- .tmp_human_hep_experiments_information$COMPOUND_TIME_DOSE[expInd]

        .modules_fc_dat <- data.frame(.modules_module_genes, limmaFits[.modules_module_genes$probe, grep(expID, colnames(limmaFits))])
        .modules_fc_dat <- .modules_fc_dat[, c('probe', 'gene',
                                                grep('logFC', colnames(.modules_fc_dat), value = TRUE), 
                                                grep('P.Value', colnames(.modules_fc_dat), value = TRUE),
                                                grep('adj', colnames(.modules_fc_dat), value = TRUE))]
          colnames(.modules_fc_dat) <- c('probe', 'gene', 'logFC', 'P', 'FDR')

        .modules_fc_dat <- .modules_fc_dat[order(.modules_fc_dat$P), ]
          rownames(.modules_fc_dat) <- 1:nrow(.modules_fc_dat)

        .modules_fc_dat$corEG <- as.character(round(.tmp_module_definitions[.modules_fc_dat$probe, 'correlationEG'], digits = 2))
        .modules_fc_dat$logFC <- as.character(round(.modules_fc_dat$logFC, digits = 2))

        MODULES_fc_dat <<- .modules_fc_dat

        fc_color_range <- c(-1.1 * max(abs(as.numeric(MODULES_fc_dat$logFC))), 1.1 * max(abs(as.numeric(MODULES_fc_dat$logFC))))
        .graphVertexColorsPalette <- data.frame(logFC = rev(seq(fc_color_range[1], fc_color_range[2], length.out = 100)),
                                         color = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(100))

        MODULES_fc_dat$color <- .graphVertexColorsPalette[sapply(as.numeric(MODULES_fc_dat$logFC), function(.fc) which(.graphVertexColorsPalette$logFC < .fc)[1] ), 'color']
        .vertexColors <- MODULES_fc_dat$color
          names(.vertexColors) <- MODULES_fc_dat$probe

        MODULES_fc_dat$shape <- c('circle', 'square')[(1 + MODULES_fc_dat$probe %in% datLeiden$module_hubgenes$probe)]


        MODULES_selectedExperiment_modEGS <- datLeiden$scored_experiments_absEG[which(datLeiden$scored_experiments_absEG$COMPOUND_TIME_DOSE == expID & 
                                                                                      datLeiden$scored_experiments_absEG$module == input$MODULES_module), 'pc1_score']

        plotTitle <- paste0('module EGS: ', round(MODULES_selectedExperiment_modEGS, digits = 3))

        .prbs <- MODULES_fc_dat$probe
        .prbCorMatrix <- cor(mergedModules$datExpr_zScored[, .prbs])
          diag(.prbCorMatrix) <- 0

        .graph <- graph.adjacency(.prbCorMatrix, weighted = TRUE, mode = "lower")
          V(.graph)$color        <- MODULES_fc_dat$color
          V(.graph)$shape        <- MODULES_fc_dat$shape
          V(.graph)$size         <- 50 * abs(as.numeric(MODULES_fc_dat$corEG))
          V(.graph)$frame.color  <- MODULES_fc_dat$color
          V(.graph)$label        <- sapply(names(V(.graph)), function(.prb) MODULES_fc_dat$gene[which(MODULES_fc_dat$probe == .prb)])
          V(.graph)$label.family <- 'mono'
          V(.graph)$label.font   <- 2
          V(.graph)$label.color  <- 'black'
          V(.graph)$label.dist   <- 0
          V(.graph)$label.degree <- 0


          E(.graph)$width <- 10 * abs(E(.graph)$weight)
          E(.graph)$color <- 'grey80'


          plot(.graph, main = plotTitle,
               layout = layout_in_circle(.graph))
        } else {
          .modules_module <- input$MODULES_module
          .modules_probes <- rownames(datLeiden$module_definitions)[which(datLeiden$module_definitions$module == .modules_module)]

          .prbCorMatrix <- cor(mergedModules$datExpr_zScored[, .modules_probes])
            diag(.prbCorMatrix) <- 0

          .graph <- graph.adjacency(.prbCorMatrix, weighted = TRUE, mode = "lower")
            V(.graph)$color        <- 'grey80'
            V(.graph)$shape        <- c('circle', 'square')[1 + as.numeric(.modules_probes %in% datLeiden$module_hubgenes$probe)]
            V(.graph)$size         <- 30
            V(.graph)$frame.color  <- 'grey80'
            V(.graph)$label        <-  sapply(names(V(.graph)), function(.prb) datLeiden$module_definitions$human.gene.symbol[which(datLeiden$module_definitions$human.entrez.gene.id == .prb)])
            V(.graph)$label.family <- 'mono'
            V(.graph)$label.font   <- 2
            V(.graph)$label.color  <- 'black'
            V(.graph)$label.dist   <- 0
            V(.graph)$label.degree <- 0



            E(.graph)$width <- 10 * abs(E(.graph)$weight)
            E(.graph)$color <- 'grey80'

          plot(.graph, main = '',
               layout = layout_in_circle(.graph))
        }
     
    })

    output$MODULES_module_enrichment <- renderDataTable({
      .enrichment <- datLeiden$module_enrichment[which(datLeiden$module_enrichment$module == input$MODULES_module), c('goID', 'goTERM', 'goONTOLOGY', 'pval')]
        colnames(.enrichment) <- c('GO ID', 'GO Term', 'GO Ontology', 'Enrichment Pval')

      datatable(.enrichment, rownames = FALSE) %>% formatSignif(columns = c('Enrichment Pval'), digits = 2)
     
    })

  ## UI - tabTXG.R
    output$TXG_map <- renderPlot({
      if(FALSE) {
        input <- list(TXG_compound = 'DOXORUBICIN',
                      TXG_time = 24,
                      TXG_dose = 'HI',
                      TXG_egs_thresholds = c(-25, 2))
      }

      plot(txgMap[[3]], type = 'fan',  # phylogram cladogram fan unrooted radial
           cex = 0.25, 
           label.offset = 0.25, 
           edge.width = 3, 
           edge.color = 'gray60',
           tip.color = 'gray80',
           no.margin = TRUE,
           use.edge.length = TRUE,
           show.tip.label = FALSE,
           rotate.tree = 40)

      .TXG_map_plot <- get('last_plot.phylo', envir = .PlotPhyloEnv)

      .txg_nodeX <- .TXG_map_plot$xx[1:.TXG_map_plot$Ntip]
      .txg_nodeY <- .TXG_map_plot$yy[1:.TXG_map_plot$Ntip]

      .txg_compound <- input$TXG_compound
      .txg_time     <- input$TXG_time 
      .txg_dose     <- input$TXG_dose

      .txg_EGS_ind <- which(datLeiden$scored_experiments_absEG$COMPOUND_NAME == .txg_compound &
                            datLeiden$scored_experiments_absEG$TIME == .txg_time &
                            datLeiden$scored_experiments_absEG$DOSE_LEVEL == .txg_dose)

      .txg_EGS_dat <- datLeiden$scored_experiments_absEG[.txg_EGS_ind, c('module', 'pc1_score')]
        rownames(.txg_EGS_dat) <- .txg_EGS_dat$module
      .txg_EGS_dat <- .txg_EGS_dat[txgMap[[3]]$tip.label, ]

      .txg_EGS_modulePlotData <- data.frame(x = .txg_nodeX, 
                                            y = .txg_nodeY, 
                                            EGS = .txg_EGS_dat$pc1_score)

      fc_color_range <- c(-1.1 * max(abs(.txg_EGS_modulePlotData$EGS)), 1.1 * max(abs(.txg_EGS_modulePlotData$EGS)))
      .graphVertexColorsPalette <- data.frame(EGS = rev(seq(fc_color_range[1], fc_color_range[2], length.out = 100)),
                                              color = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(100))

      .txg_EGS_modulePlotData$color <- .graphVertexColorsPalette[sapply(.txg_EGS_modulePlotData$EGS, function(.egs) which(.graphVertexColorsPalette$EGS < .egs)[1] ), 'color']
      .txg_EGS_modulePlotData <- .txg_EGS_modulePlotData[order(abs(.txg_EGS_modulePlotData$EGS), decreasing = TRUE), ]

      .txg_EGS_modulePlotData <-  .txg_EGS_modulePlotData[which( .txg_EGS_modulePlotData$EGS > input$TXG_egs_thresholds[1] & .txg_EGS_modulePlotData$EGS < input$TXG_egs_thresholds[2]), ]


      points(x = 1.025 * .txg_EGS_modulePlotData$x, 
             y = 1.025 * .txg_EGS_modulePlotData$y, 
            cex = 1 + abs(.txg_EGS_modulePlotData$EGS),
            pch = 21,
            # col = .txg_EGS_modulePlotData$color,
            col = 'gray60',
            bg = .txg_EGS_modulePlotData$color)

      # Highlight selected modules from TXG_map_data
      selectedRows <- input$TXG_map_data_rows_selected

      if(length(selectedRows) > 0) {
        for(i in seq(24, 0, -2)) {
          points(x   = 1.025 * .txg_EGS_modulePlotData$x[selectedRows], 
               y   = 1.025 * .txg_EGS_modulePlotData$y[selectedRows], 
               cex = 1 + abs(.txg_EGS_modulePlotData$EGS[selectedRows]),
               pch = 21,
               col = alpha('black', 0.05),
               bg  = alpha('black', 0.05),
               lwd = i)
        }

        for(i in seq(12, 0, -2)) {
          points(x   = 1.025 * .txg_EGS_modulePlotData$x[selectedRows], 
               y   = 1.025 * .txg_EGS_modulePlotData$y[selectedRows], 
               cex = 1 + abs(.txg_EGS_modulePlotData$EGS[selectedRows]),
               pch = 21,
               col = alpha('black', 0.1),
               bg  = alpha('black', 0.1),
               lwd = i)
        }

        for(i in seq(8, 0, -1)) {
          points(x   = 1.025 * .txg_EGS_modulePlotData$x[selectedRows], 
               y   = 1.025 * .txg_EGS_modulePlotData$y[selectedRows], 
               cex = 1 + abs(.txg_EGS_modulePlotData$EGS[selectedRows]),
               pch = 21,
               col = alpha('black', 0.1),
               bg  = alpha('black', 0.1),
               lwd = i)
        }

        points(x   = 1.025 * .txg_EGS_modulePlotData$x[selectedRows], 
               y   = 1.025 * .txg_EGS_modulePlotData$y[selectedRows], 
               cex = 1 + abs(.txg_EGS_modulePlotData$EGS[selectedRows]),
               pch = 21,
               col = 'black',
               bg  = 'black',
               lwd = 1)

        points(x   = 1.025 * .txg_EGS_modulePlotData$x[selectedRows], 
               y   = 1.025 * .txg_EGS_modulePlotData$y[selectedRows], 
               cex = 1 + abs(.txg_EGS_modulePlotData$EGS[selectedRows]),
               pch = 21,
               col = .txg_EGS_modulePlotData$color[selectedRows],
               bg  = .txg_EGS_modulePlotData$color[selectedRows],
               lwd = 0)




      }
    })

    output$TXG_map_data <- renderDataTable({
      .TXG_map_plot <- get('last_plot.phylo', envir = .PlotPhyloEnv)

      .txg_nodeX <- .TXG_map_plot$xx[1:.TXG_map_plot$Ntip]
      .txg_nodeY <- .TXG_map_plot$yy[1:.TXG_map_plot$Ntip]

      .txg_compound <- input$TXG_compound
      .txg_time     <- input$TXG_time 
      .txg_dose     <- input$TXG_dose

      .txg_EGS_ind <- which(datLeiden$scored_experiments_absEG$COMPOUND_NAME == .txg_compound &
                            datLeiden$scored_experiments_absEG$TIME == .txg_time &
                            datLeiden$scored_experiments_absEG$DOSE_LEVEL == .txg_dose)

      .txg_EGS_dat <- datLeiden$scored_experiments_absEG[.txg_EGS_ind, c('module', 'pc1_score')]
        rownames(.txg_EGS_dat) <- .txg_EGS_dat$module
      .txg_EGS_dat <- .txg_EGS_dat[txgMap[[3]]$tip.label, ]

      .txg_EGS_modulePlotData <- data.frame(x = .txg_nodeX, 
                                            y =.txg_nodeY, 
                                            EGS = .txg_EGS_dat$pc1_score,
                                            module = .txg_EGS_dat$module)

      .txg_EGS_modulePlotData <- .txg_EGS_modulePlotData[order(abs(.txg_EGS_modulePlotData$EGS), decreasing = TRUE), ]
      .txg_EGS_modulePlotData <- .txg_EGS_modulePlotData[which( .txg_EGS_modulePlotData$EGS > input$TXG_egs_thresholds[1] & .txg_EGS_modulePlotData$EGS < input$TXG_egs_thresholds[2]), ]

      .txg_EGS_modulePlotData$EGS <- round(.txg_EGS_modulePlotData$EGS, digits = 2)

      TXG_EGS_all <<- .txg_EGS_modulePlotData  
      
      .txg_EGS_modulePlotData[, c('module', 'EGS')]  

    }, options = list(lengthMenu = c(5, 10, 25, 50), pageLength = 25,
                      paging = TRUE, searching = TRUE),
       rownames = FALSE)


    output$TXG_map_brush_data <- renderDataTable({
      TXG_map_brush_input <<- input$TXG_map_brush[c('xmin', 'xmax', 'ymin', 'ymax')]

      selectedInd <- which(TXG_EGS_all$x > TXG_map_brush_input$xmin & TXG_EGS_all$x < TXG_map_brush_input$xmax &
                           TXG_EGS_all$y > TXG_map_brush_input$ymin & TXG_EGS_all$y < TXG_map_brush_input$ymax)

      TXG_EGS_all[selectedInd, c('module', 'EGS')]

     # TXG_map_brushedPoints <- brushedPoints(TXG_EGS_all, 
     #                                         input$TXG_map_brush)

     # TXG_map_brushedPoints
    },  options = list(lengthMenu = c(5, 10, 25), pageLength = 10,
                      paging = TRUE, searching = TRUE),
        rownames = FALSE)

  ## UI - tabGENES.R
    output$GENES_gene <- renderText(input$GENES_gene)

    output$GENES_filteredExperiments <- renderDataTable({
      if(FALSE) {
        input <- list(GENES_gene = 'SRXN1',
                      GENES_fcThreshold = c(-3, 1),
                      GENES_dose = c('LO', 'HI'),
                      GENES_time = c('2 hr', '8 hr', '24 hr'))
      }

      .genes_probe <- datLeiden$module_definitions$human.entrez.gene.id[which(datLeiden$module_definitions$human.gene.symbol == input$GENES_gene)]
      .genes_expID <- datLeiden$human_hep_experiments_information$COMPOUND_TIME_DOSE[which(datLeiden$human_hep_experiments_information$DOSE_LEVEL %in% input$GENES_dose & 
                                                                                           datLeiden$human_hep_experiments_information$SACRI_PERIOD %in% input$GENES_time)]

      .genes_geneExpData <- datLeiden$human_hep_experiments_information[.genes_expID, c('COMPOUND_NAME', 'DOSE_LEVEL', 'SACRI_PERIOD')]
        colnames(.genes_geneExpData) <- c('compound', 'dose', 'time')
        .genes_geneExpData$time <- as.numeric(gsub(' hr', '', .genes_geneExpData$time))

      # Annoying parsing issue; colnames(limmaFits) doesn't allow columns starting with numbers
        rownames(.genes_geneExpData) <- gsub('2NF', 'X2NF', rownames(.genes_geneExpData))

      .genes_geneExpData$logFC <- unlist(limmaFits[.genes_probe, paste0(rownames(.genes_geneExpData), '__logFC')])
      .genes_geneExpData$FC <- 2^.genes_geneExpData$logFC
      .genes_geneExpData$Pval     <- unlist(limmaFits[.genes_probe, paste0(rownames(.genes_geneExpData), '__P.Value')])
      .genes_geneExpData$FDR   <- unlist(limmaFits[.genes_probe, paste0(rownames(.genes_geneExpData), '__adj.P.Val')])

      .genes_geneExpData <- .genes_geneExpData[which(.genes_geneExpData$logFC > input$GENES_fcThreshold[1] & .genes_geneExpData$logFC < input$GENES_fcThreshold[2]), ]
      .genes_geneExpData <- .genes_geneExpData[order(abs(.genes_geneExpData$logFC), decreasing = TRUE), ]

      rownames(.genes_geneExpData) <- 1:nrow(.genes_geneExpData)

      #.genes_geneExpData$logFC <- round(.genes_geneExpData$logFC, digits = 2)
      #.genes_geneExpData$absFC <- round(.genes_geneExpData$absFC, digits = 2)
      #.genes_geneExpData$P   <- format(.genes_geneExpData$P, digits = 3, scientific = TRUE)
      #.genes_geneExpData$FDR <- format(.genes_geneExpData$FDR, digits = 3, scientific = TRUE) 

      datatable(.genes_geneExpData, rownames = FALSE) %>% formatCurrency(columns = c('logFC', 'FC'), currency = '') %>% formatSignif(columns = c('Pval', 'FDR'), digits = 2)

    })

    


})






















