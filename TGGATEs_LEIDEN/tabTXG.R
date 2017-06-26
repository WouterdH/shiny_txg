 
library(shiny)

# tabTXG   <- tabPanel('TXG map', 'You found the TXG tab!\n', id = 'TXG')
tabContent <- fluidPage(
				fluidRow(
				  column(width = 4, offset = 0, br(), 
				  	selectInput('TXG_compound', 'Compound:', selected = 'DIETHYL MALEATE', 
			      		        compounds, width = '100%'),

				  	fluidRow(
			      	    column(5, offset = 2, 
			      	      radioButtons('TXG_time', label = '',
			      	      choices = list('2 hour' = 2, '8 hour' = 8, '24 hour' = 24),
			      	      selected = 24)
			      	    ),
  
			      	    column(5, offset = 0, 
			      	      radioButtons('TXG_dose', label = '',
			      	      choices = list('Low' = 'LO', 'Medium' = 'MED', 'High' = 'HI'),
			      	      selected = 'HI')  
			      	    )
			      	),

				  	sliderInput("TXG_egs_thresholds", label = h3("Eigengene score:"),
        						min = -25, max = 25, value = c(-25, 25), width = '100%'),
    				
			      	fluidRow(
			      	 column(12, offset = 0, br(),
			      	  	dataTableOutput('TXG_map_data')
			      	  )
			      	)
				  ),

				  column(width = 8, offset = 0, br(),
				  	fluidRow(	
				  	  column(width = 12, offset = 1,
						
						plotOutput('TXG_map', height = 720, width = 720,
								   click = 'TXG_map_click',
								   brush = brushOpts(id = 'TXG_map_brush')
							
				  		),

				  		dataTableOutput('TXG_map_brush_data', width = 720)
				  	  )
				  	)
				  )
				)
			  )

tabTXG <- tabPanel('TXG Map', tabContent, id = 'TXG')









