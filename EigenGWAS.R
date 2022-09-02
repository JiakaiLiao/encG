library(shiny)
library(ggplot2)
library(ggrepel) # geom_text_repel
library(reshape2)
library(jcolors)
library(grid)
library(data.table)
library(dplyr)
library(DT)
library(tidyr)
library(shinycssloaders)
library(shinythemes)

#Import data
file_name = list() 
files <- dir(path = "data/test/",
             full.names = F,
             pattern = ".frq")
for (i in seq_along(files)){
  X = strsplit(files[i], "\\.")[[1]][1]
  file_name[[i]] = X
}

f1 <-function(x, y){
  cohortName <- c(x, y)
  wcName <- c(x)
  ctlName <- c(y)
  pos_list = list()
  for (i in seq_along(wcName)){
    X = read.table(paste0("data/test/", wcName[i], ".frq"), header = TRUE)
    pos_list[[i]] = X$SNP
  }
  
  for (j in seq_along(ctlName)){
    X = read.table(paste0("data/reference/", ctlName[j], ".frq"), header = TRUE)
    pos_list[[length(wcName) + j]] = X$SNP
  }
  position <- Reduce(intersect, pos_list)
  
  ctlX = read.table("data/reference/CHB.frq", header = TRUE)
  ctlX = subset(ctlX, ctlX$SNP %in% position)
  
  M = nrow(ctlX)
  Frq = matrix(NA, M, length(cohortName))
  for(i in seq_along(wcName)){
    X = read.table(paste0("data/test/", wcName[i],".frq"),header = TRUE)
    A = subset(X, X$SNP %in% position)
    dt = cbind(ctlX, A)
    # Change major and minor
    changeA1 = which(dt[,3]!=dt[,9])
    newMaf = dt[,11]
    newMaf[changeA1]=1-dt[changeA1,11]
    Frq[,i] = newMaf
  }
  
  for(i in seq_along(ctlName)){
    X = read.table(paste0("data/reference/", ctlName[i],".frq"),header = TRUE)
    A = subset(X, X$SNP %in% position)
    dt = cbind(ctlX, A)
    # Change major and minor
    changeA1 = which(dt[,3]!=dt[,9])
    newMaf = dt[,11]
    newMaf[changeA1]=1-dt[changeA1,11]
    Frq[,length(wcName) + i] = newMaf
  }
  #return(Frq)
  ## Scale and PCA
  Gs = scale(t(Frq))#将Frq矩阵转置，再做标准化
  Sig = tcrossprod(Gs)#计算Y乘Y的转置
  eigL = eigen(Sig)
  eigVec = eigL$vectors

  df <- data.frame(PC1 = eigVec[,1], PC2 = eigVec[,2], flag = cohortName)
  df$flag = factor(df$flag, levels = cohortName)
  return(df)
  }
f2 <-function(x, y){
  pos_list = list()
  cohortName <- x
  cohortN <- vector("integer",length(cohortName))
  ctlName <- y
  ctlN <- vector("integer",length(ctlName))
  for (i in seq_along(cohortName)){
    X = read.table(paste0("data/test/", cohortName[i], ".frq"), header = TRUE)
    pos_list[[i]] = X$SNP
    cohortN[i] = ceiling(mean(X$NCHROBS))
  }
  
  for (j in seq_along(ctlName)){
    Y = read.table(paste0("data/reference/", ctlName[j], ".frq"), header = TRUE)
    pos_list[[length(cohortName)+ j]] = Y$SNP
    ctlN[j] = ceiling(mean(Y$NCHROBS))
  }
  position <- Reduce(intersect, pos_list)
  ctlX = list()
  
  for(j in seq_along(ctlName)){
    Y = read.table(paste0("data/reference/", ctlName[j], ".frq"), header = TRUE)
    A = subset(Y, Y$SNP %in% position)
    ctlX[[j]] = A
  }
  
  M = nrow(ctlX[[j]])
  Fst = matrix(NA, length(cohortName), length(ctlName))
  for (i in seq_along(cohortName)){
    X = read.table(paste0("data/test/", cohortName[i], ".frq"), header = TRUE)
    A = subset(X, X$SNP %in% position)
    for(j in seq_along(ctlName)){
      dt = cbind(ctlX[[j]], A)
      n1 = ctlN[j]
      n2 = cohortN[i]
      n = n1+n2
      p1 = dt[,5]
      p2 = dt[,11]
      # Change major and minor in p2
      changeA1 = which(dt[,3]!=dt[,9])
      p2[changeA1]=1-dt[changeA1,11]
      # average p
      pMean = n1*p1/n+n2*p2/n
      Fst[i,j] = mean( 2*( n1*(p1-pMean)^2 + n2*(p2-pMean)^2 ) /n/pMean/(1-pMean) )
    }
  }
  Fst_sc = t(apply(1/Fst, 1, function(x) {x/sum(x)}))
  rownames(Fst_sc) = cohortName
  colnames(Fst_sc) = ctlName
  df = data.frame(cohortName, Fst_sc)
  colnames(df)[-c(1)] = ctlName
  df = reshape2::melt(df, id.vars = "cohortName")
  return(df)
}

button_color_css <- "
#DivCompClear, #FinderClear, #EnterTimes{
/* Change the background color of the update button
to blue. */
background: DodgerBlue;
/* Change the text size to 15 pixels. */
font-size: 15px;
}"

#Define UI
ui <- fluidPage(
  
#Navbar structure for UI  
  navbarPage("EigenGWAS", theme = shinytheme("united"),
             tabPanel("Genetic map", icon = icon("chart-bar"),
                      tags$style(button_color_css),
                      sidebarLayout(
                        sidebarPanel(
                          
                          titlePanel("Genetic map"),
                          #shinythemes::themeSelector(),
                          checkboxGroupInput("area1", h3("Within-Chinese"),
                                             choices = file_name,
                                             inline = TRUE,
                                             width = '200px',
                                             selected = c("CAS1","CAS2","Fudan","SBWCH","WBBC",
                                                          "YiKon1","YiKon2")),
                          
                          checkboxGroupInput("area2", h3("Reference area"),
                                              choices = list("ACB", "ASW", "BEB", "CEU","CHB",
                                                             "CHS", "CDX", "CLM", "ESN", "FIN", 
                                                             "GBR", "GIH", "GWD", "IBS", "ITU", 
                                                             "JPT", "KHV", "LWK", "MSL", "MXL", 
                                                             "PEL", "PJL", "PUR", "STU", "TSI", 
                                                             "YRI"),
                                              inline = TRUE ,
                                              width = '200px' ,
                                              selected = c("CHB", "CHS", "CDX")),
      
                          actionButton("Btn",label = "Run")),

                       mainPanel(
                         fluidRow(
                            column(12,
                               withSpinner(plotOutput(outputId = "scatterplotFinder"))),
                         
                           hr(),
                           br(),
                         
                           column(12,
                                withSpinner(plotOutput(outputId = "structure")))))
      )
    ),
   
    tabPanel("Data explorer",icon = icon("table"),
             fileInput("file", 
                       h4("Choose file from directory"),
                       multiple = TRUE,
                       accept =  c("text/csv",
                                   "text/comma-separated-values",
                                   "text/plain",
                                   ".csv")),
             
             checkboxInput("header", "Header", TRUE),
             
             tags$hr(),
             
             DTOutput("dataTable", width = "100%")
    ), 
  )
)
server <- function(input, output){
  
  #size = x*1024*2, x代表输入文件最大为 x MB.
  options(shiny.maxRequestSize=100*1024^2)
  
  data1  <- reactive({
    f1(input$area1, input$area2)})

  output$scatterplotFinder <- renderPlot({
    if(input$Btn == 0)
      return()
    
    df_map <- data1()
    ggplot(df_map, aes(x = PC1, y = PC2, color = flag, label = flag))+
      geom_point(cex = 2, shape = 17)+
      geom_text_repel(max.overlaps = 20)+
      scale_color_manual(values = as.vector(jcolors("pal8")[c(1 : length(df_map$flag))]))+
      # scale_x_reverse()+
      xlab(label = "PC 1")+
      ylab(label = "PC 2")+
      theme_bw()+
      theme(legend.position = "none",
            axis.line = element_line(colour = "black"),
            panel.grid = element_blank(),
            panel.border = element_blank())
  })

  data2 <- reactive({
    f2(input$area1, input$area2)})
  output$structure <- renderPlot({

    if(input$Btn == 0)
      return()

    df <- data2()
    ggplot(df, aes(x=cohortName, y=value, color=variable, fill = variable))+
      geom_bar(stat = 'identity',position = "stack")+
      scale_color_jcolors(palette = "pal3")+
      scale_fill_jcolors(palette = "pal3")+
      labs(x = "Cohort", y = "Genetic composition",color = "", fill = "")+
      theme_bw()+
      theme(axis.line = element_line(colour = "black"),
            panel.border = element_blank(),
            panel.grid = element_blank())
  })
##Data Explorer##########
  output$dataTable <- renderDT({
    req(input$file)
    
    tryCatch({
        inFile <- input$file
        if (is.null(inFile)) return(NULL)
        df1 <- fread(inFile$datapath, header = input$header)},
        
        error = function(e){
          stop(safeError(e))
        }
      )
    datatable(df1)  
    })
}

shinyApp(ui = ui, server = server)