#here is the code for the app

#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.15/bioc"))
library(shiny)
library(dplyr)
library(ape)
library(BiocManager)
library(ggtree)

#read in dataframes
D04_LSN <- read.csv(here::here("D04_LSN.csv"))
D04_LSR <- read.csv(here::here("D04_LSR.csv"))
D04_RBT <- read.csv(here::here("D04_RBT.csv"))
D06_BCC <- read.csv(here::here("D06_BCC.csv"))
D06_MIS <- read.csv(here::here("D06_MIS.csv"))
D07_LAC <- read.csv(here::here("D07_LAC.csv"))
D12_RBK <- read.csv(here::here("D12_RBK.csv"))
D13_MIS <- read.csv(here::here("D13_MIS.csv"))

possibledat<- c("D04_LSN","D04_LSR", "D04_RBT",
               "D06_BCC","D06_MIS","D07_LAC",
               "D12_RBK","D13_MIS")
playouts <- c("rectangular", "circular", "unrooted")

#### Writing the Helper Function ####
startplot <- function(dat, hmap = FALSE, phylolayout = 'rectangular'){
  #use PC
  rownames(dat) <- dat$ID
  distmatrix<-dist(dat %>% 
                     select(-c(ID,tot)))
  
  dathclust<-hclust(distmatrix)
  
  if (hmap) {
    heatmap(as.matrix(distmatrix))
    
  } else {
    #plot(dathclust)
    phylodat<-as.phylo(dathclust)
    ggtree(phylodat, layout = phylolayout, branch.length = 'none') +
      geom_tiplab()
      
  }
  #return(p)
}

#### Building the App page ####

ui <- fluidPage(
  headerPanel("Heat Map and Dendrogram of Mutations Among Similar Regions"),
  selectInput('site', 'Sample', possibledat),
  selectInput('ptype', 'Plot Type',
              c("Phylogeny","Heatmap"),
              selected = "Phylogeny"),
  uiOutput('phylayout'),
  plotOutput('plot1')
  
)

server <- function(input, output) {
  dat <- reactive({
    get(input$site)})
  
  output$phylayout <- renderUI({
    if(input$ptype != 'Heatmap'){
      selectInput('phylayout', "Phylogeny Layout",
                  playouts,
                  selected = "rectangular")
    }
  })
  
  output$plot1 <- renderPlot({
    if (input$phylayout == "unrooted"){
      phylotype <- 'daylight'
    }else{
      phylotype <- input$phylayout
    }
    if(input$ptype != "Heatmap"){
      startplot(dat(), hmap = FALSE, phylolayout = phylotype)
    }else{
      startplot(dat(), hmap = TRUE)
    }
  })
}

shinyApp(ui = ui, server = server)