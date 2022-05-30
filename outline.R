#this page is the outline of what the code should follow

library(dplyr)
library(tidyr)
library(tidyverse)

alldat <- readxl::read_xlsx(here::here("HumanMelDat.xlsx")) %>% 
  select(GENE,SAMPLE) %>%
  filter(SAMPLE != 'D04_LAC05')

#create new data sets to be used for shiny
#makes a big list with each obs a sample
listWgenes<- split(alldat, alldat$SAMPLE)

for (i in 1:length(listWgenes)){
  name <- names(listWgenes)[i]
  x <- as.data.frame(table(listWgenes[[i]]$GENE)) %>%
    pivot_wider(names_from = Var1, 
                values_from = Freq) %>%
    mutate(ID = name) %>%
    select(ID, everything())
  assign(name, x)
}

#making new lists with df to then merge
#note only those with 3+ df are used
D04_LSN <- list(D04_LSN01, D04_LSN02, D04_LSN03,
              D04_LSN04, D04_LSN05) %>% 
            reduce(full_join, by='ID') %>%
            replace(is.na(.), 0) %>%
            as.data.frame()
  

D04_LSR <- list(D04_LSR01, D04_LSR02, D04_LSR03,
                D04_LSR04, D04_LSR05, D04_LSR06,
                D04_LSR09)%>% 
            reduce(full_join, by='ID') %>%
            replace(is.na(.), 0) %>%
            as.data.frame()

D04_RBT <- list(D04_RBT02, D04_RBT06, D04_RBT12,
                D04_RBT20)%>% 
          reduce(full_join, by='ID') %>%
          replace(is.na(.), 0) %>%
          as.data.frame()

D06_BCC <- list(D06_BCC01, D06_BCC03, D06_BCC04,
                D06_BCC05, D06_BCC09)%>% 
          reduce(full_join, by='ID') %>%
          replace(is.na(.), 0) %>%
          as.data.frame()

D06_MIS <- list(D06_MIS03, D06_MIS04, D06_MIS06,
                D06_MIS07, D06_MIS08, D06_MIS10,
                D06_MIS11, D06_MIS12)%>% 
          reduce(full_join, by='ID') %>%
          replace(is.na(.), 0) %>% 
          as.data.frame()

D07_LAC <- list(D07_LAC03, D07_LAC04, D07_LAC05)%>% 
          reduce(full_join, by='ID') %>%
          replace(is.na(.), 0) %>%
          as.data.frame()

D12_RBK <- list(D12_RBK01, D12_RBK02, D12_RBK03)%>% 
          reduce(full_join, by='ID') %>%
          replace(is.na(.), 0) %>%
          as.data.frame()
 
D13_MIS <- list(D13_MIS07, D13_MIS08, D13_MIS17,
                D13_MIS19, D13_MIS20, D13_MIS21)%>% 
          reduce(full_join, by='ID') %>%
          replace(is.na(.), 0) %>%
          as.data.frame()

#write a helper function that will be primary use for app
startplot <- function(dat,heatmap=F){
  #use PC
  rownames(dat) <- dat$ID
  distmatrix<-dist(dat[,-1])
  
  dathclust<-hclust(distmatrix)
  dathclust$labels <- dat$ID
  
  #plots typical dendrogram or heatmap
  p <- ifelse(heatmap==F, yes = plot(dathclust), 
              no = heatmap(as.matrix(distmatrix)))
  return(p)
}

startplot(D04_LSR, heatmap = TRUE)
