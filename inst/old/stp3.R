

stateProf = function() {
message("NOTE: takes some time to initialize for 40000 promoter regions")
ver = R.version
 if (ver$minor <= "1.3") {
 rowRanges = rowData
 assign("rowRanges<-", get("rowData<-"), .GlobalEnv)
}
#library(shiny)
#library(ggplot2)
#library(GenomicRanges)
#library(Homo.sapiens)
#library(erma)
ermaset = makeErmaSet()[,4:17]
availSyms = keys(Homo.sapiens, keytype="SYMBOL")
availSyms = sort(availSyms[grep("^[A-Z]", availSyms)])

ui = fluidPage(titlePanel("ChromImpute state profiles"),
#
# should be selectize
#
   fluidRow(column(3,selectizeInput('sym', 'Gene symbol', choices=NULL,
           selected="CD28", options=list(placeholder="[symbol]"),
           multiple=FALSE)),
   column(3,numericInput('scope', 'Scope', value=50000))),
   fluidRow( plotOutput('plot') )
   )

server = function(input, output, session) {
   INI = reactive({
      if (is.null(input$sym) | nchar(input$sym)==0) return("CD28")
      else input$sym
      })
   output$plot = renderPlot(
      stateProfile(ermaset, INI(), ctsize=12, width=input$scope)
      )
   updateSelectizeInput(session, 'sym', choices = availSyms, server = TRUE)
}

shinyApp(ui=ui, server=server)
}
