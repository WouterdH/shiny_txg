 
library(shiny)

#load('/Users/Wouter/stack/ownCloud/TOX/gitProjects/shiny/TGGATEs/modulesLeidenParsed_with_mergedModules.RData')
#modules <- unique(datLeiden$module_definitions$module)
#modules <- modules[order(as.numeric(do.call(rbind, strsplit(modules, ':'))[, 2]))]
#  names(modules) <- modules


tabContent <- fluidPage(
				fluidRow(
				  column(width = 6, offset = 0,
				    column(6, offset = 0, br(),
			      	  selectInput('modcor_module_x', 'Module 1:', 
			      		        modules, width = '100%'),

			      	  fluidRow(
			      	    column(5, offset = 1, 
			      	      radioButtons('modcor_x_time', label = '',
			      	      choices = list('2 hour' = 2, '8 hour' = 8, '24 hour' = 24),
			      	      selected = 24)
			      	    ),
  
			      	    column(5, offset = 1, 
			      	      radioButtons('modcor_x_dose', label = '',
			      	      choices = list('Low' = 'LO', 'Medium' = 'MED', 'High' = 'HI'),
			      	      selected = 'HI')  
			      	    )
			      	  )

			        ),

			        column(6, offset = 0, br(),
			      	  selectInput('modcor_module_y', 'Module 2:', 
			      		        modules, width = '100%'),

			      	  fluidRow(
			      	    column(5, offset = 1, 
			      	      radioButtons('modcor_y_time', label = '',
			      	      choices = list('2 hour' = 2, '8 hour' = 8, '24 hour' = 24),
			      	      selected = 8)
			      	    ),
  
			      	    column(5, offset = 1, 
			      	      radioButtons('modcor_y_dose', label = '',
			      	      choices = list('Low' = 'LO', 'Medium' = 'MED', 'High' = 'HI'),
			      	      selected = 'HI')  
			      	    )
			      	  )

			        ),

			        fluidRow(
					  column(offset = 0, width = 12, br(),
					  	#h4('Selected compounds:'),
					  	dataTableOutput('brush_info_modcor')
					  )
					)

				  ),

				  column(width = 6, offset = 0,
				  	plotOutput('module_cor_plot', height = 580, width = 540,
				  		       click = 'module_cor_plot_click',
				  			   brush = brushOpts(id = 'module_cor_plot_brush')
				  	)
				  )
				)
			  )

tabMODULE_COR <- tabPanel('Module correlation', tabContent, id = 'MODULES')









