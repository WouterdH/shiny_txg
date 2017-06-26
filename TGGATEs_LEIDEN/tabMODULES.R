 
library(shiny)

#load('/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs/modulesLeidenParsed_with_mergedModules.RData')
#compounds <- unique(datLeiden$human_hep_experiments_information$COMPOUND_NAME)
#  names(compounds) <- compounds

tabContent <- fluidPage(
				fluidRow(
				  column(width = 5, offset = 0, br(), 
				  	fluidRow(
				  	 	column(12,
				  	 		selectInput('MODULES_module', 'Module:', selected = 'WGCNA|HHCYT|TGGATEs|LEIDEN:144', 
			      		        		modules, width = '100%')
				  	 	)
				  	),

				  	#fluidRow(
				  	# 	column(12, 
				  	# 		selectInput('MODULES_compound', 'Compound:', selected = 'DIETHYL MALEATE', 
			      	#	        		compounds, width = '100%')
				  	# 	)
				  	#),

				  	#fluidRow(
			      	#    column(5, offset = 2, 
			      	#      radioButtons('MODULES_time', label = '',
			      	#      choices = list('2 hour' = 2, '8 hour' = 8, '24 hour' = 24),
			      	#      selected = 24)
			      	#    ),
  
			      	#    column(5, offset = 0, 
			      	#      radioButtons('MODULES_dose', label = '',
			      	#      choices = list('Low' = 'LO', 'Medium' = 'MED', 'High' = 'HI'),
			      	#      selected = 'HI')  
			      	#    )
			      	#),

				  	fluidRow(	
				  	  column(width = 12, offset = 0,
						
						plotOutput('MODULES_gene_plot' , height = 512, width = 512)
				  	  )
				  	),

				  	fluidRow(
				  	  column(12, offset = 0, 
				  	  	dataTableOutput('MODULES_module_enrichment')
				  	  )
				  	)
				  ),

				  column(width = 7, offset = 0, br(),
				  	fluidRow(column(width = 10, offset = 1, 
				  	  sliderInput("MODULES_egs_thresholds", label = '',
        						min = -17, max = 25, value = c(-17, 25), width = '100%')
				  	)),

				  	fluidRow(	
				  	  column(width = 12, offset = 0,
						
						plotOutput('MODULES_doseResponse_plot', height = 256, #width = 768,
							click = 'MODULES_doseResponse_plot_click',
				  			brush = brushOpts(id = 'MODULES_doseResponse_plot_brush'))
				  	  )
				  	),

				  	fluidRow(	
				  	  column(width = 12, offset = 0,
						
						plotOutput('MODULES_timeResponse_plot', height = 256, #width = 768,
							click = 'MODULES_doseResponse_plot_click',
				  			brush = brushOpts(id = 'MODULES_doseResponse_plot_brush'))
				  	  )
				  	),

				  	fluidRow(
				  	  column(width = 12, offset = 0, 
				  	  	dataTableOutput('MODULES_brush_info')
				  	  )
				  	),

			      	fluidRow(
			      	 column(12, offset = 0, br(),
			      	  	# h4('Foldchange data:'),
			      	  	# div(tableOutput('MODULES_module_dat_bak'), style = 'font-size:85%')
			      	  	dataTableOutput('MODULES_module_dat')
			      	  )
			      	)

				  )
				)
			  )

tabMODULES <- tabPanel('Modules', tabContent, id = 'MODULES')









