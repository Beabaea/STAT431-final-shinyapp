#here is the code for the app
#install_version("nloptr", version = "2.0.1", repos = "http://cran.us.r-project.org")
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.14/bioc"))
library(shiny)
library(dplyr)
library(ape)
library(BiocManager)
library(ggtree)
library(ggplot2)

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
      geom_tiplab() +
      scale_x_continuous(expand = c(.25, .25))+
      scale_y_continuous(expand = c(.25, .25))
  }
}

#very small helper function to create vector of mutations in one sample
vecmutations <- function(dat){
  toString(colnames((dat[,colSums(dat) > 0])))
}

#function for bootstrapping
bootci <- function(data,uno,dos){
  subdat<-data %>% 
    filter(ID %in% c(uno,dos)) %>%
    select(-c(ID,tot))
  
  num<-100 #number of bootstraps
  boot<-sjstats::bootstrap(subdat,n=num)
  tibboots<-lapply(boot$strap,FUN=as_tibble)
  preddist<-sapply(tibboots,FUN=dist)
  obs <- dist(subdat)
  
  #CI code found from https://www.programmingr.com/statistics/confidence-interval-in-r/
  center <- mean(preddist, na.rm = TRUE)
  stddev <- sd(preddist, na.rm = TRUE)
  error <- qnorm(0.95)*stddev/sqrt(num)
  
  paste("The 95% confidence interval places the distance between these points as",
        round(center-error,3), "to",
        round(center+error,3), sep=" ")
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
        tabPanel("Plot", 
                 plotOutput('plot1', width = "90%"),
                 uiOutput('CI'),
                 textOutput('confint')),
        
        tabPanel("Data Table",tableOutput('infotab')),
        tabPanel("General Information",
                 HTML("This project is by Beatrice Weier. <br> 
                 These plots take data from the Nature paper
                 <a 'https://www.nature.com/articles/s41586-020-2785-8'> 
                 The genomic landscapes of individual melanocytes from human skin </a>;
                 an article with data regarding mutations in melanocytes from 
                 different subject. The data is depicted by the patient (P#), the
                 site of cell (e.g. RBK is right Back), and a number to indicate
                 individual cell.The data used included specific regarding the 
                 mutation. <br>
                 The R code takes this data and creates a 
                 purely objective scaling of mutations with equal distances for 
                 all mutations/genes. A big disclaimer is that the distance are 
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
  
  #render new buttons
  output$phylayout <- renderUI({
    if(input$ptype != 'Heatmap'){
      radioButtons('phylayout', "Phylogeny Layout",
                  playouts,
                  selected = "rectangular")
    }
  })
  #Render selector for specific data CI
  output$CI <- renderUI({
    selectInput('sub',"Select Two For Bootstrap analysis",
                dat()$ID, multiple=TRUE)
  })
  
  output$confint <- renderText(
    if (length(input$sub)==2){
      bootci(dat(),input$sub[1],input$sub[2])
    }else{
      warning("Select Only Two")
    }

    )
  #for plot
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
  
  #this is for data table tab
  output$infotab <- renderTable({
    Mutations <- NULL
    
    for (i in 1:nrow(dat())){ #can't figure this out without a loop but luckily df are small :,)
      subdat<-dat()[i,]%>%select(-c(ID,tot))
      muts<-vecmutations(subdat)
      Mutations <- c(Mutations, muts)
    }
    
    data.frame('ID' = dat()$ID, 'Number_of_Mutations' = dat()$tot, 'Genes' = Mutations)
  })
}

shinyApp(ui = ui, server = server)