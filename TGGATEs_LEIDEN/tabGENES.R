 
library(shiny)


tabContent <- fluidPage(
				fluidRow(
				  column(width = 5, offset = 0, br(), 
				  	fluidRow(
				  	 	column(12,
				  	 		textInput('GENES_gene', '', value = 'SRXN1', 
			      		        	  width = '100%')
				  	 	)
				  	),

				  	br(),

				  	fluidRow(
				  		sliderInput("GENES_fcThreshold", label = h5("Gene logFC cutoff"),
        						min = -7, max = 9, value = c(-7, 9),
        						width = '100%')
			      	),

			      	fluidRow(
				  	 	column(5, offset = 1, 
				  	 		checkboxGroupInput('GENES_dose', '', 
				  	 							choices  = c('Low' = 'LO', 'Medium' = 'MED', 'High' = 'HI'),
				  	 							selected = c('LO', 'MED', 'HI'), 
			      		        				width = '100%')
				  	 	),

				  	 	column(5, offset = 1, 
				  	 		checkboxGroupInput('GENES_time', '', 
				  	 							choices = c('2 hour' = '2 hr', '8 hour' = '8 hr', '24 hour' = '24 hr'),
				  	 							selected = c('2 hr', '8 hr', '24 hr'), 
			      		        				width = '100%')
				  	 	)


				  	),


			      	fluidRow(
			      	 column(12, offset = 0, br(),
			      	  	div(dataTableOutput('GENES_filteredExperiments', width = '100%'), style = 'font-size:80%')
			      	  )
			      	)
				  ),

				  column(width = 7, offset = 0, br(),
				  	fluidRow(	
				  	  column(width = 12, offset = 0, br()
						
						
				  	  )
				  	)

				  )
				),

				br(), br(), br(), br(),

				verbatimTextOutput('GENES_gene')
			  )

tabGENES <- tabPanel('Genes', tabContent, id = 'GENES')









