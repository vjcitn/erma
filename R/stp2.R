
stateProfOLD = function() {
message("NOTE: takes some time to initialize for 40000 promoter regions")
ver = R.version
 if (ver$minor <= "1.3") {
 rowRanges = rowData
 assign("rowRanges<-", get("rowData<-"), .GlobalEnv)
}
library(shiny)
library(ggplot2)
library(GenomicRanges)
library(Homo.sapiens)
library(erma)
ermaset = makeErmaSet()[,4:17]
availSyms = keys(Homo.sapiens, keytype="SYMBOL")
availSyms = availSyms[grep("^[A-Z]", availSyms)]

ui = fluidPage(
#
# should be selectize
#
   fluidRow(column(3,selectInput('sym', 'Gene symbol', choices=sort(availSyms),
          multiple=FALSE, selected="CD28")),
   column(3,numericInput('scope', 'Scope', value=50000))),
   fluidRow( plotOutput('plot') )
   )

server = function(input, output) {
   output$plot = renderPlot(
      stateProfile(ermaset, input$sym, ctsize=12, width=input$scope)
   )
}
shinyApp(ui=ui, server=server)
}
