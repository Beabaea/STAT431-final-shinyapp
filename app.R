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



#will move into sever, lets have it be reactive to 
Mutations <- NULL
for (i in 1:nrow(dat)){ #can't figure this out without a loop but luckily df are small :,)
  subdat<-dat[i,]%>%select(-c(ID,tot))
  muts<-vecmutations(subdat)
  Mutations <- c(Mutations, muts)
}
table <- data.frame('ID' = dat$ID, 'Total number of mutations' = dat$tot, 'Mutations' = Mutations)


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
      geom_tiplab() #+ 
      #ggplot2::xlim(0, 100)
      
      
  }
}

#very small helper function to create vector of mutations in one sample
vecmutations <- function(dat){
  toString(colnames((dat[,colSums(dat) > 0])))
}

#### Building the App page ####

ui <- fluidPage(
  headerPanel("Heat Map and Dendrogram of Mutations Among Similar Regions"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput('site', 'Sample', possibledat),
      radioButtons('ptype', 'Plot Type',
                   c("Phylogeny","Heatmap"),
                   selected = "Phylogeny"),
      uiOutput('phylayout')
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput('plot1')),
        tabPanel("Data Table",tableOutput('infotab')),
        tabPanel("General Information",
                 HTML("This project is by Beatrice Weier. <br> 
                 These plots take data from 
                 <a 'https://www.nature.com/articles/s41586-020-2785-8'> this paper </a> ;
                 an article with data regarding mutations in melanocytes from 
                 different subject. The R code takes this data and creates a 
                 purely objective scaling of mutations with equal distances for 
                 all mutations. A big disclaimer is that the distance are 
                 calculated from the gene name, not the specific mutation within 
                 the gene. I hope to eventually better this app and make 
                 phylogenies accountign for biologicalvariance and probabilites 
                of certain mutations"))
      )
    )
  )
  
)

server <- function(input, output) {
  dat <- reactive({
    get(input$site)})
  
  output$phylayout <- renderUI({
    if(input$ptype != 'Heatmap'){
      radioButtons('phylayout', "Phylogeny Layout",
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
  
  output$infotab <- renderTable({
    Mutations <- NULL
    
    for (i in 1:nrow(dat())){ #can't figure this out without a loop but luckily df are small :,)
      subdat<-dat()[i,]%>%select(-c(ID,tot))
      muts<-vecmutations(subdat)
      Mutations <- c(Mutations, muts)
    }
    
    data.frame('ID' = dat()$ID, 'Number_of_Mutations' = dat()$tot, 'Genes' = Mutations)
    
  })
  
  #output$txt <- renderText(
   # "This project is by Beatrice Weier. These plots take data from https://www.nature.com/articles/s41586-020-2785-8 ;
    #an article with data regarding mutations in melanocytes from different subject. 
    #The R code takes this data and creates a purely objective scaling of mutations 
    #with equal distances for all mutations. A big disclaimer is that the distance 
    #are calculated from the gene name, not the specific mutation within the gene.
    #I hope to eventually better this app and make phylogenies accountign for biological
    #variance and probabilites of certain mutations"
  #)
}

shinyApp(ui = ui, server = server)