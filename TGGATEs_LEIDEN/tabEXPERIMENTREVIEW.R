 
library(shiny)

tabContent <- fluidPage(
				fluidRow(
				  column(width = 5, offset = 0, br(), 
				  	selectInput('EXPERIMENTREVIEW_compound', 'Compound:', selected = 'DOXORUBICIN', 
			      		        compounds, width = '100%'),

				  	sliderInput("EXPERIMENTREVIEW_egs_thresholds", label = h3("Eigengene score:"),
        						min = -17, max = 25, value = c(-17, 25),
        						width = '100%'),
    				
			      	fluidRow(
			      	 column(12, offset = 0, br(),
			      	  	dataTableOutput('EXPERIMENTREVIEW_brush_info')
			      	  )
			      	)
				  ),

				  column(width = 7, offset = 0, br(),
				  	fluidRow(	
				  	  column(width = 12, offset = 0,
						
						plotOutput('EXPERIMENTREVIEW_doseResponse', height = 256, #width = 768,
							click = 'EXPERIMENTREVIEW_doseResponse_click',
				  			brush = brushOpts(id = 'EXPERIMENTREVIEW_doseResponse_brush'))
				  	  )
				  	),

				  	fluidRow(	
				  	  column(width = 12, offset = 0,
						
						plotOutput('EXPERIMENTREVIEW_timeResponse', height = 256, #width = 768,
							click = 'EXPERIMENTREVIEW_doseResponse_click',
				  			brush = brushOpts(id = 'EXPERIMENTREVIEW_doseResponse_brush'))
				  	  )
				  	),

				  	fluidRow(
				  	  column(width = 12, offset = 0,
				  	  	#h5('Tentative heatmap of selected modules & compound.')
				  	  	plotOutput('EXPERIMENTREVIEW_heatmap', height = 512)
				  	  	
				  	  )
				  	)

				  )
				)
			  )

tabEXPERIMENTREVIEW <- tabPanel('Compounds', tabContent, id = 'EXPERIMENTREVIEW')









