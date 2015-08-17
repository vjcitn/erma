stdvname = function() "Standardized.Epigenome.name"

getShort = function()
 get(load(dir(
    system.file("data",package="erma"),full.names=TRUE, pattern="short")))

selectByShort = function(gf, shortnames) {
 short = getShort()
 stdn = colData(gf)[,"Standardized.Epigenome.name"]
 inds = match(shortnames, short)
 gf[, which(stdn %in% names(short)[inds])]
}

setGeneric("shortNames", function(x) standardGeneric("shortNames"))
setMethod("shortNames", "ErmaSet", function(x) {
  getShort()[colData(x)[,stdvname()]]
})
 
 
stateProf = function(ermaset, shortCellType=TRUE, ctsize=10) {
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
#ermaset = makeErmaSet()
short_celltype = get(load(dir(
    system.file("data",package="erma"),full.names=TRUE, pattern="short")))
availSyms = keys(Homo.sapiens, keytype="SYMBOL")
availSyms = sort(availSyms[grep("^[A-Z]", availSyms)])

ui = fluidPage(titlePanel("ChromImpute state profiles"),
#
# should be selectize
#
   fluidRow(column(2,selectizeInput('sym', 'Gene symbol', choices=NULL,
           selected="CD28", options=list(placeholder="CD28"),
           multiple=FALSE)),
   column(2,numericInput('upstr', 'Upstream', value=2000)),
   column(2,numericInput('downstr', 'Downstream', value=200)),
   column(2,selectInput('short', 'Celltype format?', selected="short",
              choices=c("short", "long")))),
   fluidRow( plotOutput('plot') )
   )

server = function(input, output, session) {
   INI = reactive({
      if (is.null(input$sym) | nchar(input$sym)==0) return("CD28")
      else input$sym
      })
   output$plot = renderPlot({
      shortCellType = (input$short == "short")
      stateProfile(ermaset, INI(), ctsize=ctsize, upstream=input$upstr,
              downstream=input$downstr, shortCellType=shortCellType)
      })
   updateSelectizeInput(session, 'sym', choices = availSyms, server = TRUE)
}

shinyApp(ui=ui, server=server)
}
