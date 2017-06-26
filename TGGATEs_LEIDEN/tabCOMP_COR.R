 
library(shiny)

#load('/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs/modulesLeidenParsed_with_mergedModules.RData')
#compounds <- unique(datLeiden$human_hep_experiments_information$COMPOUND_NAME)
#  names(compounds) <- compounds


tabContent <- fluidPage(
				fluidRow(
				  column(width = 6, offset = 0,
				    column(6, offset = 0, br(),
			      	  selectInput('compcor_compound_x', 'Compound 1:', 
			      		        compounds),

			      	  fluidRow(
			      	    column(5, offset = 1, 
			      	      radioButtons('compcor_x_time', label = '',
			      	      choices = list('2 hour' = 2, '8 hour' = 8, '24 hour' = 24),
			      	      selected = 24)
			      	    ),
  
			      	    column(5, offset = 1, 
			      	      radioButtons('compcor_x_dose', label = '',
			      	      choices = list('Low' = 'LO', 'Medium' = 'MED', 'High' = 'HI'),
			      	      selected = 'HI')  
			      	    )
			      	  )

			        ),

			        column(6, offset = 0, br(),
			      	  selectInput('compcor_compound_y', 'Compound 2:', 
			      		        compounds),

			      	  fluidRow(
			      	    column(5, offset = 1, 
			      	      radioButtons('compcor_y_time', label = '',
			      	      choices = list('2 hour' = 2, '8 hour' = 8, '24 hour' = 24),
			      	      selected = 8)
			      	    ),
  
			      	    column(5, offset = 1, 
			      	      radioButtons('compcor_y_dose', label = '',
			      	      choices = list('Low' = 'LO', 'Medium' = 'MED', 'High' = 'HI'),
			      	      selected = 'HI')  
			      	    )
			      	  )

			        ),

			        fluidRow(
					  column(offset = 0, width = 12, br(),
					  	# h4('Selected modules:'),
					  	dataTableOutput('brush_info_compcor')
					  )
					)

				  ),

				  column(width = 6, offset = 0,
				  	plotOutput('compound_cor_plot', height = 580, width = 540,
				  		       click = 'compound_cor_plot_click',
				  			   brush = brushOpts(id = 'compound_cor_plot_brush')
				  	)
				  )
				)
			  )

tabCOMP_COR <- tabPanel('Compound correlation', tabContent, id = 'COMPS')









