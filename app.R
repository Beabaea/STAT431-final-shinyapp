#here is the code for the app

library(shiny)

#read in dataframes
D04_LSN <- read.csv("C:/Users/Beatrice Weier/Documents/D04_LSN.csv")
D04_LSR <- read.csv("C:/Users/Beatrice Weier/Documents/D04_LSR.csv")
D04_RBT <- read.csv("C:/Users/Beatrice Weier/Documents/D04_RBT.csv")
D06_BCC <- read.csv("C:/Users/Beatrice Weier/Documents/D06_BCC.csv")
D06_MIS <- read.csv("C:/Users/Beatrice Weier/Documents/D06_MIS.csv")
D07_LAC <- read.csv("C:/Users/Beatrice Weier/Documents/D07_LAC.csv")
D12_RBK <- read.csv("C:/Users/Beatrice Weier/Documents/D12_RBK.csv")
D13_MIS <- read.csv("C:/Users/Beatrice Weier/Documents/D13_MIS.csv")

possibledat<-c("D04_LSN","D04_LSR", "D04_RBT",
               "D06_BCC","D06_MIS","D07_LAC",
               "D12_RBK","D13_MIS")

#### Writing the Helper Function ####
startplot <- function(dat, hmap = FALSE){
  #use PC
  rownames(dat) <- dat$ID
  distmatrix<-dist(dat[,-1])
  
  dathclust<-hclust(distmatrix)
  
  #plots typical dendrogram or heatmap
  # ifelse(hmap == FALSE, plot(dathclust), 
  #              heatmap(as.matrix(distmatrix)))
  
  if (hmap) {
    heatmap(as.matrix(distmatrix))
    
  } else {
    plot(dathclust)
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
  plotOutput('plot1')
  
)

server <- function(input, output) {
  dat <- reactive({
    get(input$site)})
  
  output$plot1 <- renderPlot({
    if(input$ptype != "Heatmap"){
      startplot(dat(), hmap = FALSE)
    }else{
      startplot(dat(), hmap = TRUE)
    }
  })
}

shinyApp(ui = ui, server = server)