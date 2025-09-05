############## Authors:
############## Charuka D. Wickramasinghe, Karmanos Cancer Institute, School of Medicine, Wayne State University, Detroit, MI, USA,
############## NSSM. Hapuhinna,  Department of Mathematics and Statistics, Northern Kentucky University, KY, USA
############## Date: 09/03/2025 

################################################################################
############################# Call Libraries ###################################
library(shiny)
library(shinythemes)
library(ggplot2)
library(qqplotr)
library(ggpubr)
library(cowplot)
library(deSolve)
library(ggrepel)
library(shinyalert)
library(tidyverse)
library(tidyr)
library(gridExtra)
library(knitr)
library(plotly)
library(shinydashboard)
library(shinydashboardPlus)
library(tableHTML)
library(DT)
library(formattable)
library(dplyr)
library(ggpubr)
library(pROC)
library(Metrics)
library(DescTools)
library(data.table)
library(reshape2)
library(fresh)
library(shinyWidgets)
library(htmltools)
library(bslib)
library(statmod)
library(pracma) 
library("readxl")
library(ggformula)
library(ggtext)
library(DEoptim)
library(numDeriv)
library(MASS)
library(optimx)
library(qpNCA)
library(DescTools)
library(PKNCA)

################################################################################
################################################################################
################################################################################
################################################################################
############# Define system of ODE  ############################################

mb4.model <- function (t, x, params,artconc) {
    Cbb = x[1]
    Cbm = x[2]
    Cccsf = x[3]
    Cscsf = x[4]
    
    Vbb = params[1]
    Vbm = params[2]
    Vccsf = params[3]
    Vscsf = params[4]
    
    Qbrain = params[5]
    Qbulkbc = params[6]
    Qbulkcb = params[7]
    QCsink = params[8]
    QSsink = params[9]
    QSin = params[10]
    QSout = params[11]
    Qglybm = params[12]
    Qglyccsf = params[13]
    Qglyscsf = params[14]
    
    PSB = params[15]
    PSC = params[16]
    PSE = params[17]
    CLeffbm = params[18]
    CLupbm = params[19]
    CLeffccsf = params[20]
    CLupccsf = params[21]
    
    lambb = params[22]
    lambm = params[23]
    lamccsf = params[24]
    fubb = params[25]
    fubm = params[26]
    fuccsf = params[27]
    
    Cart = artconc(t)
    
    dCbbdt <- (1/Vbb)*(
        Qbrain*(Cart-Cbb)
        + PSB*(lambm*fubm*Cbm - lambb*fubb*Cbb)
        + CLeffbm*fubm*Cbm
        - CLupbm*fubb*Cbb
        + PSC*(lamccsf*fuccsf*Cccsf - lambb*fubb*Cbb)
        + CLeffccsf*fuccsf*Cccsf # modified on 09142022
        - CLupccsf*fubb*Cbb # modified on 09142022
        + QCsink*Cccsf
        + QSsink*Cscsf
    )
    dCbmdt <- (1/Vbm)*(
        PSB*(lambb*fubb*Cbb - lambm*fubm*Cbm)
        - CLeffbm*fubm*Cbm
        + CLupbm*fubb*Cbb
        + PSE*(lamccsf*fuccsf*Cccsf - lambm*fubm*Cbm)
        + Qbulkcb*Cccsf # modified on 09142022
        - Qbulkbc*fubm*Cbm # modified on 09142022
        - Qglybm*fubm*Cbm # modified on 09142022
    )
    dCccsfdt <- (1/Vccsf)*(
        PSC*(lambb*fubb*Cbb - lamccsf*fuccsf*Cccsf)
        + CLupccsf*fubb*Cbb # modified on 09142022
        - CLeffccsf*fuccsf*Cccsf # modified on 09142022
        + PSE*(lambm*fubm*Cbm - lamccsf*fuccsf*Cccsf)
        + Qbulkbc*fubm*Cbm # modified on 09142022
        - Qbulkcb*Cccsf # modified on 09142022
        + QSout*Cscsf - QSin*Cccsf - QCsink*Cccsf
        - Qglyccsf*Cccsf # modified on 09142022
    )
    dCscsfdt <- (1/Vscsf)*(
        QSin*Cccsf - QSout*Cscsf - QSsink*Cscsf
        - Qglyscsf*Cscsf # modified on 09142022
    )
    dxdt <- c(dCbbdt,dCbmdt,dCccsfdt,dCscsfdt)
    list(dxdt)
}

run.sim4 <- function(time,conc,param){
    ttimes0 = time
    tconcs0 = conc
    bparms = as.numeric(param)
    
    artconc = approxfun(x=ttimes0,y=tconcs0,method="linear",rule=2)
    
    btimes = ttimes0

    bstart <- c(Cbb=0,Cbm=0,Cccsf=0,Cscsf=0)
    
    taconc = ode(
        func=mb4.model
        ,y=bstart
        ,times=btimes
        ,parms=bparms
        ,artconc=artconc
    )
    
    tsimd = as.data.frame(cbind(taconc,tconcs0))
    dimnames(tsimd)[[2]] = c("Time","Cbb","Cbm","Cccsf","Cscsf","Cart")
    tsimd
}   


################################################################################
################################################################################
# This function is for the sensitivity analysis

sen4 <- function(time,conc,param){
  ttimes0 = time
  tconcs0 = conc
  bparms = as.numeric(param)
  
  artconc = approxfun(x=ttimes0,y=tconcs0,method="linear",rule=2)
  
  btimes = ttimes0
  
  bstart <- c(Cbb=0,Cbm=0,Cccsf=0,Cscsf=0 )
  
  taconc = ode(
    func=mb4.model
    ,y=bstart
    ,times=btimes
    ,parms=bparms
    ,artconc=artconc
  )

  tsimd = as.data.frame(cbind(taconc,tconcs0))
  dimnames(tsimd)[[2]] = c("Time","Cbb","Cbm","Cccsf","Cscsf", "Cart")
  

  tsimd
  
}   




################################################################################
############# Adding Losss function to paramater estimation ####################

Loss.4comp <- function(param,conc,time,Obs_Cbb,Obs_Cbm,Obs_Cccsf,Obs_Cscsf){
  ttimes0 = time
  tconcs0 = conc
  bparms = as.numeric(param)
  artconc = approxfun(x=ttimes0,y=tconcs0,method="linear",rule=2)
  btimes = ttimes0
  bstart <- c(Cbb=0,Cbm=0,Cccsf=0,Cscsf=0)
  solve4comp = ode(
    func=mb4.model
    ,y=bstart
    ,times= btimes
    ,parms=bparms
    ,artconc=artconc
    ,method = lsoda
  )
  
  Simdata = as.data.frame(cbind(solve4comp,tconcs0))
  dimnames(Simdata)[[2]] = c("Time","Cbb","Cbm","Cccsf","Cscsf","Cart")
  
  
  # loss <- 0
  loss <- c()
  for (i in btimes) {
    Cbb_Obs <- Obs_Cbb[i]
    Cbb_Sim <- Simdata[i,2]
    
    Cbm_Obs <- Obs_Cbm[i]
    Cbm_Sim <- Simdata[i,3]
    
    Cccsf_Obs <- Obs_Cccsf[i]
    Cccsf_Sim <- Simdata[i,4]
    
    Cscsf_Obs <- Obs_Cscsf[i]
    Cscsf_Sim <- Simdata[i,5]
    
    residuals_brain4 <- ((Cbb_Obs - Cbb_Sim)^2 + 
                           (Cbm_Obs - Cbm_Sim)^2 +
                           (Cccsf_Obs -  Cccsf_Sim)^2 + 
                           (Cscsf_Obs - Cscsf_Sim)^2)
    
    
    loss[i] <-  residuals_brain4
    
  }
  
  Total_Loss <- sum(loss)
  Total_Loss

}  



################################################################################
################################################################################
################################################################################
################################################################################


############################ Define user interface #############################

ui <- fluidPage(
  tags$head(tags$style(HTML('.progress-bar {background-color: blue;}'))),
  tags$head(tags$style(type="text/css",".shiny-output-error { visibility: hidden;}", 
                       ".shiny-output-error:before { visibility: hidden; }")),
  theme = shinytheme("cerulean"),br(),
  titlePanel(div("TumorPred: Version 1.0 (2025)", style = "color: green")),br(),
  
########################### Define sidebar panel.###############################

sidebarPanel(width = 4,tags$figure(class = "top",tags$img(src = "lg.png",
                                  width = 400, alt = "CBM"),tags$figcaption("")),
               
conditionalPanel('input.method === "Upload Data"',br(),br(),fileInput('file1', 
                                      'Choose csv File',accept=c('text/csv',
                  'text/comma-separated-values,text/plain', '.csv')), checkboxInput('header', 'Header', TRUE), radioButtons('sep', 'Separator',
                  c(Comma=',',Semicolon=';',Tab='\t'),','), radioButtons('quote', 'Quote',c(None='','Double Quote'='"','Single Quote'="'"),'"')),
    
  
  conditionalPanel('input.method === "Module 1: Model Simulations"', 
                   
                   br(),
                   br(),
                   br(),  
                   tags$hr(style="border-top: 5px solid black;"),
                   br(),
                   
                   h4( tags$b(" Please select inputs for model simulation:",style = "color: maroon")),
                   br(), 
                   
                   
                   selectInput("user_time", 
                               label = "Time",
                               ""
                   ),
                   selectInput("user_cart", 
                               label = "Plasma concentration",
                               ""
                   ),
                   
                   
                   selectInput("user_param", 
                               label = "Model parameters",
                               ""
                   ),
                   
                   
                   numericInput("user_PKmintime", 
                                label = "Minimum time point for Cmax, Tmax and AUC ",value = 0,min = 0,max = 100000
                                
                   ), 
                   
                   numericInput("user_PKmaxtime", 
                                label = "Maximum time point for Cmax, Tmax and AUC",value = 48,min = 0,max = 100000
                                
                   ),
                   

                   
                   actionButton("Resinput", 
                                
                                "Submit for Results", 
                                class = "btn-success"), 
                   tags$hr(style="border-top: 5px solid black;"),
                   br()),
  
      
  conditionalPanel('input.method === "Module 3: Parameter Estimation"',br(),  tags$hr(style="border-top: 5px solid black;"),
                  h4( tags$b("Please select inputs for parameter estimation:",style = "color: maroon" )),br(), 
                  selectInput("isim4_timep", label = "Time", ""),
                  selectInput("isim4_cartp", 
                      label = "Plasma concentration",
                      ""
          ),
          
          numericInput("maxiter", 
                       label = "Number of iterations",value = 1,min = 1,max = 1e10
                       
          ),         
          
          numericInput("rel_tol", 
                       label = "Relative tolarance",value = 1e-2,min = 0,max = 1e10
                       
          ), 
          
          numericInput("VTR_Value", 
                       label = "VTR",value = 0.1,min = 0,max = 1e10
                       
          ),

          fluidRow(
            column(6,
                   numericInput("Vbbmin", 
                                label = "Vbb min",value = 0.06,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Vbmmin", 
                                label = "Vbm min",value = 1.1,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Vccsfmin", 
                                label = "Vccsf min",value = 0.1,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Vscsfmin", 
                                label = "Vscsf min",value = 0.02,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Qbrainmin", 
                                label = "Qbrain min",value = 37.5,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QbulkBCmin", 
                                label = "QbulkBC min",value = 0.005,min = 0,max =1e10
                                
                   ),
                   
                   numericInput("QbulkCBmin", 
                                label = "QbulkCB min",value = 0.005,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QCsinkmin", 
                                label = "QCsink min",value = 0.01,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QSsinkmin", 
                                label = "QSsink min",value = 0.007,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QSinmin", 
                                label = "QSin min",value = 0.01,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QSoutmin", 
                                label = "QSout min",value = 0.007,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Qglybmmin", 
                                label = "Qglybm min",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QglyCCSFmin", 
                                label = "QglyCCSF min",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QglySCSFmin", 
                                label = "QglySCSF min",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("PSBmin", 
                                label = "PSB min",value = 134.5,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("PSCmin", 
                                label = "PSC min",value = 67,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("PSEmin", 
                                label = "PSE min",value = 299,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("CLeffbmmin", 
                                label = "CLeffbm min",value = 109,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("CLupbmmin", 
                                label = "CLupbm min",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("CLeffccsfmin", 
                                label = "CLeffccsf min",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("CLupccsfmin", 
                                label = "CLupccsf min",value = 11.5,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("lamdabbmin", 
                                label = "lamdabb min",value = 0.03,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("lamdabmmin", 
                                label = "lamdabm min",value = 0.01,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("lamdaccsfmin", 
                                label = "lamdaccsf min",value = 0.02,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("fubbmin", 
                                label = "fubb min",value = 0.1,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("fubmmin", 
                                label = "fubm min",value = 0.04,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("fuccsfmin", 
                                label = "fuccsf min",value = 0.9,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("fuscsfmin", 
                                label = "fuscsf min",value = 0.9,min = 0,max = 1e10
                                
                   ),
                   
            ),
            column(6,
                   numericInput("Vbbmax", 
                                label = "Vbb max",value = 0.064952435,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Vbmmax", 
                                label = "Vbm max",value = 1.104115461,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Vccsfmax", 
                                label = "Vccsf max",value = 0.103984624,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Vscsfmax", 
                                label = "Vscsf max",value = 0.025996156,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Qbrainmax", 
                                label = "Qbrain max",value = 38.5,min = 0,max =1e10
                                
                   ),
                   
                   numericInput("QbulkBCmax", 
                                label = "QbulkBC max",value = 0.005164106,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QbulkCBmax", 
                                label = "QbulkCB max",value = 0.005164106,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QCsinkmax", 
                                label = "QCsink max",value = 0.01277633,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QSsinkmax", 
                                label = "QSsink max",value = 0.007761342,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QSinmax", 
                                label = "QSin max",value = 0.015251337,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QSoutmax", 
                                label = "QSout max",value = 0.007489995,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("Qglybmmax", 
                                label = "Qglybm max",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QglyCCSFmax", 
                                label = "QglyCCSF max",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("QglySCSFmax", 
                                label = "QglySCSF max",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("PSBmax", 
                                label = "PSB max",value = 135.5,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("PSCmax", 
                                label = "PSC max",value = 68,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("PSEmax", 
                                label = "PSE max",value = 301,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("CLeffbmmax", 
                                label = "CLeffbm max",value = 111,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("CLupbmmax", 
                                label = "CLupbm max",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("CLeffccsfmax", 
                                label = "CLeffccsf max",value = 0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("CLupccsfmax", 
                                label = "CLupccsf max",value = 12.5,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("lamdabbmax", 
                                label = "lamdabb max",value = 0.033,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("lamdabmmax", 
                                label = "lamdabm max",value = 0.017,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("lamdaccsfmax", 
                                label = "lamdaccsf max",value = 0.03,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("fubbmax", 
                                label = "fubb max",value = 0.13,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("fubmmax", 
                                label = "fubm max",value = 0.05,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("fuccsfmax", 
                                label = "fuccsf max",value = 1.0,min = 0,max = 1e10
                                
                   ),
                   
                   numericInput("fuscsfmax", 
                                label = "fuscsf max",value = 1.0,min = 0,max = 1e10
                                
                   ),
            )
          ),
          
          
          numericInput("PKmintime", 
                       label = "Minimum time point for Cmax, Tmax and AUC ",value = 0,min = 0,max = 1e10
                       
          ), 
          
          numericInput("PKmaxtime", 
                       label = "Maximum time point for Cmax, Tmax and AUC",value = 48,min = 0,max = 1e10
                       
          ),
          
          br(),
          actionButton("subPtable", 
                       "Submit for parameter estimation", 
                       class = "btn-success"), 
          br(),
          br(),
          tags$hr(style="border-top: 5px solid black;"),
          br(),
          
        ),     
        
        
  
        
#################################################################################################  
# Begin Adding B4_IIV
conditionalPanel(
  'input.method === "Module 2: Sensitivity Analysis"',

br(),  
tags$hr(style="border-top: 5px solid black;"),
br(),

h4( tags$b(" Please select inputs for parameter sensitivity analysis:",style = "color: maroon")),
br(), 
 

selectInput("Sen_time", 
            label = "Time",
            ""
),
selectInput("Sen_cart", 
            label = "Plasma concentration",
            ""
),

 
selectInput("isim4_param", 
            label = "Model parameters",
            ""
),


  selectInput("Pickpara4", 
              label = "Select a parameter: ",
              list("Vbb","Vbm","Vccsf","Vscsf"
                   ,"Qbrain","Qbulk,BC","Qbulk,CB"
                   ,"QCsink","QSsink","QSin","QSout"
                   ,"Qgly,bm","Qgly,CCSF","Qgly,SCSF","PSB","PSC","PSE"
                   ,"CLeff,bm","CLup,bm","CLeff,ccsf","CLup,ccsf"
                   ,"lambdabb","lambdabm","lambdaccsf"
                   ,"fubb","fubm","fuccsf")
  ),
  
  
  numericInput("sim4.Pickminpara", 
               label = "The minimum value of the selected parameter",value = 0.01,min = 0,max = 1e10
               
  ), 
  
  numericInput("sim4.Pickmaxpara", 
               label = "The maximum value of the selected parameter", value = 10,min = 0,max = 1e10
               
  ),
  

  numericInput("sim4.Pickparanum", 
               label = "The number of parameter points of the selected parameter ", value = 3,min = -1,max = 1e10
               
  ), 
  
  

 numericInput("PKtimemin", 
              label = "Minimum time point for Cmax, Tmax and AUC ",value = 0,min = 0,max = 1e10
              
 ), 
 
 numericInput("PKtimemax", 
              label = "Maximum time point for Cmax, Tmax and AUC",value = 48,min = 0,max = 1e10
              
 ), 
 

br(),
  actionButton("sim4.Sensub", 
               "Submit for parameter sensitivity analysis", 
               class = "btn-success"),


br(),
# tags$hr(style="border-color: black;"),




tags$hr(style="border-top: 5px solid black;"), 
 
  
),
        


),  
########################### end of side bar panel. #############################
################################################################################
################################################################################
################################################################################


tags$hr(style="border-top: 5px solid red;"),
##### Begin main panel
  



########################### Define main panel.##################################

mainPanel(width = 8,
tabsetPanel(type =  "pills",id = 'method',
            
            tabPanel('About TumorPred ',br(), 
            tags$hr(style="border-top: 5px solid red;"),
            h2("Welcome to TumorPred 1.0 !"),
            
            
            
            
            tags$h3("TumorPred is a web-based platform (built with R/Shiny) 
            for running model simulations, 
                    conducting sensitivity analysis, and performing parameter 
                    estimation in compartmental modeling, 
                    specifically for the human central nervous system and brain 
                    tumors.", style = "color: maroon;"),
            tags$hr(style="border-top: 5px solid red;"),    
            h3("Schematic Illustration of the brain model structure"),  
                 fluidRow(tags$figure(class = "top",tags$img(src = "B4M.png",width = 600,alt = "CBM"),tags$figcaption("")),),
                 tags$hr(style="border-top: 5px solid red;"),),
            
            tabPanel("Upload Data", tags$hr(style="border-top: 5px solid red;"),
                    navbarPage(title = 'DataTable Options',
                    tabPanel('Data',     DT::dataTableOutput('tab.data'))
                   # tabPanel('Summary',verbatimTextOutput('summ.data'))
                    )
                    ),
            

            tabPanel('Module 1: Model Simulations',     
                     br(),
                     tags$hr(style="border-top: 5px solid red;"),
                     
                     tabsetPanel( type =  "pills",       #tabset panel "A" 
                                  
                                  
                                  
                                  tabPanel("View :: Parameters : Concentration profiles : Pharmacokinetic Parameters ",  
                                           tags$hr(style="border-top: 5px solid red;"),
                                           conditionalPanel('input.method === "Module 1: Model Simulations"',
                                                            #condition = "(input.isim4_timep !== 'undefined' && input.isim4_timep.length > 0 &&
                                                            # input.isim4_cartp !== 'undefined' && input.isim4_cartp.length > 0 
                                                            #)",
                                                            
                                                            #  tags$hr(style="border-top: 5px solid red;"),
                                                            navbarPage(
                                                              title = '',
                                                              tabPanel('User input parameters', tableOutput("userparameters")),
                                                              tabPanel('Concentrations table', tableOutput("userConctab")),
                                                              tabPanel('Concentration plots',plotOutput('userConplot',width = 900, height = 600 )),
                                                              tabPanel('Concentration log plots',plotOutput('userconlogplot',width = 900, height = 600 )),
                                                              tabPanel('Pharmacokinetic parameters (Cmax, Tmax, AUC)', tableOutput("userPK")),
                                                              
                                                            ), # end of navbarPage
                                                            
                                                            
                                           ) # end of conditional panel
                                           
                                  ), # end of tab panel
                                  
                                  
                                  tabPanel('Download :: Results',     
                                           tags$hr(style="border-top: 5px solid red;"),
                                           conditionalPanel( 'input.method === "Module 1: Model Simulations"',
                                                             #   tags$hr(),
                                                             
                                                             # navbarPage(
                                                             #   title = '',
                                                             
                                                             br(),
                                                             
                                                             downloadButton("userparm", "User parameters",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             # tags$hr(),
                                                             br(),
                                                             
                                                             downloadButton("userconct", "Concentration table",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             # tags$hr(),
                                                             br(),
                                                             downloadButton("usercontpltdwn", "Concentration plots",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             # tags$hr(),
                                                             br(),
                                                             downloadButton("userconclotplt", "Concentration log plots",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             #  tags$hr(),
                                                             
                                                             br(),
                                                             downloadButton("userpkparm", "Pharmacokinetic parameters (Cmax, Tmax, AUC)",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             
                                                             #  ), # end of navbarPage
                                                             
                                                             
                                           ) # end of conditional panel
                                  ), # end of tab panel
                                  
                                  

                                  
                                  
                     ), #end of tabset panel "A"
                     
                     
                     
                     
            ), 
            
          
            tabPanel('Module 2: Sensitivity Analysis',
                     
                     br(),
                     tags$hr(style="border-top: 5px solid red;"),
                     
                     #tags$h3("Result 3: Sensitivity analysis", style = "color: maroon;"),   
                     tabsetPanel( type =  "pills",  
                                  
                                  
                                  tabPanel('View :: Parameter Sensitivity Simulation Results',     
                                           tags$hr(style="border-top: 5px solid red;"),
                                           conditionalPanel( 'input.method === "Module 2: Sensitivity Analysis"',
                                                             # condition = "(input.isim4_timep !== 'undefined' && input.isim4_timep.length > 0 &&
                                                             #          input.isim4_cartp !== 'undefined' && input.isim4_cartp.length > 0 
                                                             #         )",
                                                             
                                                             #  tags$hr(style="border-top: 5px solid red;"),
                                                             navbarPage(
                                                               title = '',
                                                               tabPanel('Concentration table', tableOutput("sim4.sendata4a")),
                                                               tabPanel('Concentration plots',plotOutput('sim4.senplot4b',width = 900, height = 1000)),
                                                               tabPanel('Concentration log plots',plotOutput('sim4.senplot4blog',width = 900, height = 1000)),
                                                               tabPanel('Pharmacokinetic parameters (Cmax, Tmax, AUC)', tableOutput("sim4.tmax4")),
                                                               tabPanel('AUC scatter plots',plotOutput('sim4.aucplot',width = 900, height = 800)),
                                                               # tabPanel('AUC scatter table', tableOutput("sim4.scatble")),
                                                             ), # end of navbarPage
                                                             
                                                             
                                           ) # end of conditional panel
                                  ), # end of tab panel
                                  
                                  
                                  
                                  tabPanel('Download :: Parameter Sensitivity Simulation Results',     
                                           tags$hr(style="border-top: 5px solid red;"),
                                           conditionalPanel( 'input.method === "Module 2: Sensitivity Analysis"',
                                                             #   tags$hr(),
                                                             
                                                             # navbarPage(
                                                             #   title = '',
                                                             
                                                             br(),
                                                             
                                                             downloadButton("sim4.sencondata", "Concentration table",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             # tags$hr(),
                                                             br(),
                                                             downloadButton("sim4.senconplot", "Concentration plots",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             # tags$hr(),
                                                             br(),
                                                             downloadButton("sim4.senconlogplot", "Concentration log plots",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             #  tags$hr(),
                                                             
                                                             br(),
                                                             downloadButton("sim4.parmako4", "Pharmacokinetic parameters (Cmax, Tmax, AUC)",class = "btn-block"),
                                                             tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                                             
                                                             #  ), # end of navbarPage
                                                             
                                                             
                                           ) # end of conditional panel
                                  ), # end of tab panel
                                  
                                  
                     ), # end of third tabset panel
                     
                     #), # end of tabset panel
                     
                     
                     
                     
            ), # end of tab panel
            
            
            
            

            tabPanel('Module 3: Parameter Estimation',     
                     br(),
                     tags$hr(style="border-top: 5px solid red;"),
            
                     tabsetPanel( type =  "pills",       #tabset panel "A" 
                     
                     
                     
            tabPanel("View :: Estimated Parameters : Concentration profiles : Pharmacokinetic Parameters ",  
                     tags$hr(style="border-top: 5px solid red;"),
                     conditionalPanel('input.method === "Module 3: Parameter Estimation"',
                                      #condition = "(input.isim4_timep !== 'undefined' && input.isim4_timep.length > 0 &&
                                     # input.isim4_cartp !== 'undefined' && input.isim4_cartp.length > 0 
                                      #)",
                                      
                                      #  tags$hr(style="border-top: 5px solid red;"),
                                      navbarPage(
                                        title = '',
                                        tabPanel('Estimated parameters', tableOutput("para_table")),
                                        tabPanel('Concentrations table', tableOutput("singlerun.sim4pesti")),
                                        tabPanel('Concentration plots',plotOutput('sinplot.sim4pesti',width = 900, height = 600 )),
                                        tabPanel('Concentration log plots',plotOutput('sinlogplot.sim4pesti',width = 900, height = 600 )),
                                        tabPanel('Pharmacokinetic parameters (Cmax, Tmax, AUC)', tableOutput("pkparameters.sim4pesti")),
                                        
                                      ), # end of navbarPage
                                      
                                      
                     ) # end of conditional panel
                     
            ), # end of tab panel

            
            tabPanel('Download :: Results',     
                     tags$hr(style="border-top: 5px solid red;"),
                     conditionalPanel( 'input.method === "Module 3: Parameter Estimation"',
                                       #   tags$hr(),
                                       
                                       # navbarPage(
                                       #   title = '',
                                       
                                       br(),
                                       
                                       downloadButton("paraestidata", "Estimated parameters",class = "btn-block"),
                                       tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                       # tags$hr(),
                                       br(),
                                       
                                       downloadButton("paracondata", "Concentration table",class = "btn-block"),
                                       tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                       # tags$hr(),
                                       br(),
                                       downloadButton("paraconplot", "Concentration plots",class = "btn-block"),
                                       tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                       # tags$hr(),
                                       br(),
                                       downloadButton("paraconlogplot", "Concentration log plots",class = "btn-block"),
                                       tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                       #  tags$hr(),
                                       
                                       br(),
                                       downloadButton("parmako4", "Pharmacokinetic parameters (Cmax, Tmax, AUC)",class = "btn-block"),
                                       tags$head(tags$style(".btn-block{background:white;} .btn-block{color: blue;}")),
                                       
                                       #  ), # end of navbarPage
                                       
                                       
                     ) # end of conditional panel
            ), # end of tab panel
            
                     ), #end of tabset panel "A"
            ), 
            
#####################################################################################################################            
# Begin Adding  B4_IIV






#### Download sample input files
tabPanel("A Sample Input File", 
         
         tags$hr(style="border-top: 5px solid maroon;"),
         br(), 
         h2("A sample input file for 4 compartment brain model"),
         br(),  
         downloadButton("B4InputData", "Download Input File for Brain 4 Model"),
         br(), 
         
         br(), 
         h2(""),
         
         br(),
         br(), 
         h4("The input file contains a number of columns, named as Time (sampling times),
         Plasma (drug plasma concentrations), Parameters (system- and drug-specific parameter names), 
         Cbb (Observed data for brain blood concentration), Cbm (Observed data for brain mass concentration), 
         Cccsf (Observed data for cranial CSF concentration),
         and Cscsf (Observed data for spinal CSF concentration."),
            
         br(), 
         br(),
         
         h4("Important Note 1: User must keep the same colum names as shown in the 
         sample input file for different simulations for different drugs."),
         
         br(),
         
         h4("Important Note 2: Estimating 27 parameters for a give data set may take a while. 
             Thus it is reccomended always to run for one iteration first to see the compelxity of the analysis.
             Then user can perform analysis for higher numebr of iterations based on the accuracy needed.")
      
),


tabPanel("Help", 
         
         tags$hr(style="border-top: 5px solid red;"),
         
         h2("Contact"),
         br(),   
         br(), 
         tags$h4("Please direct your questions related to TumorPred version 1.0 to the following email address, queries 
            concerning pharmacology, statistics, mathematics, and technical matters can also be directed there.",style = "color: green;"),
         br(), 
         h4("gi6036@wayne.edu"),
         
         h4("hapuhinnan1@nku.edu"),
         br(),  
         
         h2("Share your feedback"),
         br(),
         tags$h4("We value your opinion and strive to align our vision and concepts with your needs and expectations. 
         Your feedback is crucial in ensuring that you get the most out of our software. Please share your feedback and opnion to:
            ",style = "color: green;"),
         br(),
         h4("wickramasinghec@karmanos.org"),
         
),



       ) # end of tabset panel of main panel 

    )   # end of main panel



########################### end main panel.#####################################
################################################################################
################################################################################
################################################################################
    
) # end of ui fluid page


########################### end ui fluid page. #################################
################################################################################
################################################################################
################################################################################





########################### Define server  #####################################

server <- function(input, output, session) {
    data <- reactive({
        study.data <- NULL
        inFile <- input$file1 
        if (is.null(inFile)){
            return(NULL)
        } 
        study.data <- read.csv(inFile$datapath, header=input$header, sep=input$sep, 
                               quote=input$quote)
        study.data
    })

  
    observe({
      updateSelectInput(
        session,
        "isim4_param",
        choices=c("",names(data())[grep("sim",names(data()),ignore.case=T)]))
    })
    
    observe({
      updateSelectInput(
        session,
        "isim4_legpos",
        choices=c("",names(data())))
    })
    

 
 
  ##### update parameter input codes
    
    observe({
      updateSelectInput(
        session,
        "isim4_timep",
        choices=c("",names(data())[grep("time",names(data()),ignore.case=T)]))
    })
    
    
    observe({
      updateSelectInput(
        session,
        "isim4_cartp",
        choices=c("",names(data())[grep("plasma",names(data()),ignore.case=T)]))
    })
    
    
    observe({
      updateSelectInput(
        session,
        "Sen_time",
        choices=c("",names(data())[grep("time",names(data()),ignore.case=T)]))
    })
    
    
    observe({
      updateSelectInput(
        session,
        "Sen_cart",
        choices=c("",names(data())[grep("plasma",names(data()),ignore.case=T)]))
    })
    
    
    
    
    observe({
      updateSelectInput(
        session,
        "user_time",
        choices=c("",names(data())[grep("time",names(data()),ignore.case=T)]))
    })
    
    
    observe({
      updateSelectInput(
        session,
        "user_cart",
        choices=c("",names(data())[grep("plasma",names(data()),ignore.case=T)]))
    })
    
    
    
    observe({
      updateSelectInput(
        session,
        "user_param",
        choices=c("",names(data())[grep("sim",names(data()),ignore.case=T)]))
    })
    
    
  
    
# End Adding B4_IIV
####################################################################################################
    
 
  
 #### Input the user data file   
    output$tab.data <- DT::renderDataTable(
        DT::datatable(data(), options = list(pageLength = 10),editable=FALSE)
    )
    output$summ.data <- renderPrint({
        tdata = data()
        totalobs.pos = grep("total",names(tdata),ignore.case = T)
        unboundobs.pos = grep("unbound",names(tdata),ignore.case = T)
        
        if(length(totalobs.pos)>0){
            tobb = tdata[,c(totalobs.pos:(totalobs.pos+3))]
            thead = paste0(as.character(tobb[1,]),"_total")
            tvbb = tobb[-1,]

            tnbb = length(na.omit(as.numeric(as.character(tvbb[,1]))))
            if(tnbb==0){
                tdata[,c(totalobs.pos:(totalobs.pos+3))] = NA
            }else{
                tvbb2 = tvbb[1:tnbb,]
                totald = matrix(as.numeric(as.character(unlist(tvbb2))),dim(tvbb2)[1],dim(tvbb2)[2])
                tdata[,c(totalobs.pos:(totalobs.pos+3))] = NA
                tdata[1:dim(totald)[1],c(totalobs.pos:(totalobs.pos+3))] = totald
            }
            
            names(tdata)[c(totalobs.pos:(totalobs.pos+3))] = thead
        }
        
        if(length(unboundobs.pos)>0){
            uobb = tdata[,c(unboundobs.pos:(unboundobs.pos+3))]
            uhead = paste0(as.character(uobb[1,]),"_unbound")
            uvbb = uobb[-1,]
            
            unbb = length(na.omit(as.numeric(as.character(uvbb[,1]))))
            if(unbb==0){
                tdata[,c(unboundobs.pos:(unboundobs.pos+3))] = NA
            }else{
                uvbb2 = uvbb[1:unbb,]
                unboundd = matrix(as.numeric(as.character(unlist(uvbb2))),dim(uvbb2)[1],dim(uvbb2)[2])
                tdata[,c(unboundobs.pos:(unboundobs.pos+3))] = NA
                tdata[1:dim(unboundd)[1],c(unboundobs.pos:(unboundobs.pos+3))] = unboundd
            }
            
            names(tdata)[c(unboundobs.pos:(unboundobs.pos+3))] = uhead
        }

        Hmisc::describe(tdata)
    })
    
    output$desc.data <- renderText(
        extract_help(input$rname.data, input$file1, to="html")
    )

    
    
    
################################## Output results of user input file##################
#####################################################################################
    
    
    #To get the parameter list from data file
    user_paramtable4 <- reactive({
      if(input$user_param==""){
        shinyalert("Oops!", "Please select a list of parameters", type = "error")
        return(NULL)
      }else{
        # Get the parameter column and remove empty/NA values
        iparam <- na.omit(as.character(data()[,input$user_param]))
        
        # Filter out any empty strings and keep only the first 27 values
        iparam <- iparam[iparam != ""]
        
        if(length(iparam) > 27) {
          iparam <- iparam[1:27]  # Take only the first 27 parameters
        }
      }
      
      if(length(iparam) != 27){
        shinyalert("Oops!", "Some parameters are missing", type = "error")
        return(NULL)
      }
      
      # Convert to numeric
      as.numeric(iparam)
      
      
      
      
      
    })
    
    #To get parameter list with the parameter names 
    User_iparamout4 <- reactive({
      userpara <- user_paramtable4()
      
      if(is.null(userpara)) {
        return(NULL)
      }
      
      panames = c(
        "Vbb","Vbm","Vccsf","Vscsf"
        ,"Qbrain","Qbulk,BC","Qbulk,CB"
        ,"QCsink","QSsink","QSin","QSout"
        ,"Qgly,bm","Qgly,CCSF","Qgly,SCSF","PSB","PSC","PSE"
        ,"CLeff,bm","CLup,bm","CLeff,ccsf","CLup,ccsf"
        ,"lambdabb","lambdabm","lambdaccsf"
        ,"fubb","fubm","fuccsf"
      )
      
      # Get the original parameter values as strings from the data
      original_params <- na.omit(as.character(data()[,input$user_param]))
      original_params <- original_params[original_params != ""]
      
      if(length(original_params) > 27) {
        original_params <- original_params[1:27]
      }
      
      tpatab = data.frame(Parameter = panames, Value = original_params)
      tpatab
      
      
    })
    
    
    
    
    
    #To output a list of user selected  parameter 
    output$userparameters <- renderTable({
      if (input$Resinput>0) { 
        withProgress(User_iparamout4() ,message = 'Calculation in progress...', 
                      value = 0.3 )
      }
    })  
    
    
    
    #To downolad user selected  parameter  
    output$userparm <- downloadHandler(
      filename = "user_selected_parameters.csv",
      content = function(file) {
        if (input$Resinput>0) { 
          
          withProgress(message = 'Downloading in progress',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          write.csv(User_iparamout4(), file, row.names = FALSE)
        }
      }
    )   
    
    
    
    #To output concentration for  of user selected  parameter    
    Conc_for_user_param <- reactive({
      itime <- na.omit(data()[,input$user_time])
      iconc <- na.omit(data()[,input$user_cart])
      iparam <- na.omit(data()[,input$user_param])
      sinrun <- run.sim4(time=itime,conc=iconc,param=iparam)
      sinrun
    }) 
    
    #To output the concentrations table to main panel
    output$userConctab <- renderTable({
      if (input$Resinput>0) { 
        withProgress(Conc_for_user_param()  ,message = 'Calculation in progress...', 
                     value = 0.3 )
      }
    })
    
    
    #To download the concentrations table 
    output$userconct <- downloadHandler(
      filename = "ConcentationTable.csv",
      content = function(file) {
        if (input$Resinput>0) { 
          
          withProgress(message = 'Downloading in progress',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          write.csv(Conc_for_user_param(), file, row.names = FALSE)
        }
      }
    )   
    
    
   
################### Concentration plots#########################################
    #Reactive function to get concentrations profiles to main panel
    
    User_Concentrations <- reactive({
      singleConc <- Conc_for_user_param()
      ObsData <- data()
      colors <- c("Predicted" = "blue", "Observed" = "red")
      
      plot_Cbbsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,2]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cbb, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbb") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cbmsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,3]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cbm, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbm") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      
      plot_Cccsfsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,4]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cccsf, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cccsf") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cscsfsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,5]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cscsf, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cscsf") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16))      
      
      
      singlist <- list(plot_Cbbsin, plot_Cbmsin,plot_Cccsfsin,plot_Cscsfsin)
      
      grid.arrange(grobs = singlist, ncol = 2, nrow =2)
      #AllSingle <- grid.arrange(grobs = singlist)
      # print(AllSingle)
    })
    
    
    #To get concentration profiles to main panel
    output$userConplot <- renderPlot({ 
      if (input$Resinput>0) { 
        withProgress(User_Concentrations()  ,message = 'Calculation in progress...', 
                     value = 0.3 )
      }
    }) 
    
    
    #Reactive function to download concentration profiles
    
    user_Conc_down <- reactive({
      singleConc <- Conc_for_user_param()
      
      ObsData <- data()
      colors <- c("Predicted" = "blue", "Observed" = "red")
      #singleConc[,2]
      plot_Cbbsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,2]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cbb, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbb") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "bottom",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cbmsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,3]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cbm, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbm") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "bottom",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      
      plot_Cccsfsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,4]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cccsf, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cccsf") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "bottom",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cscsfsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,5]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cscsf, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cscsf") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "bottom",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16))      
      
      
      singlist <- list(plot_Cbbsin, plot_Cbmsin,plot_Cccsfsin,plot_Cscsfsin)
      
      #AllSingle <- grid.arrange(grobs = singlist)
      
      #dev.off()
      #grid.arrange(grobs = singlist)

      singlist <- list(print(plot_Cbbsin),print(plot_Cbmsin),print(plot_Cccsfsin),print(plot_Cscsfsin))
      dev.off()
      grid.arrange(grobs = singlist)
    })
    
    
    #To  download the concentration profiles
    output$usercontpltdwn <- downloadHandler(
      filename <- function(){ "Concentration_plots.pdf"} ,
      content = function(file) {
        if (input$Resinput>0) { 
          withProgress(message = 'Downloading in progress...',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          pdf(file,paper="default")
          plot(user_Conc_down())
          dev.off()
        }
      }
    )     
    
    
 
    #Reactive function to get log concentration plots  to main panel for user parameters    
    user_singlelogrunplot_para <- reactive({
      singleConc <- Conc_for_user_param()
      
      ObsData <- data()
      colors <- c("Predicted" = "blue", "Observed" = "red")
      
      plot_Cbbsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,2])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cbb), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbb") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cbmsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,3])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cbm), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbm") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      
      plot_Cccsfsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,4])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cccsf), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cccsf") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cscsfsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,5])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cscsf), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cscsf") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16))      
      
      singloglist <- list(plot_Cbbsin, plot_Cbmsin,plot_Cccsfsin,plot_Cscsfsin)
      
      #   AlllogSingle <- grid.arrange(grobs = singloglist ,ncol=2)
      #   print(AlllogSingle)
      grid.arrange(grobs = singloglist,ncol=2, nrow =2)
      
      
      
      
      
    })  
    
    
    
    
    #To get concentration log profiles to main panel for user parameters
    output$userconlogplot <- renderPlot({ 
      if (input$Resinput>0) { 
        withProgress(user_singlelogrunplot_para()  ,message = 'Calculation in progress...', 
                     detail= 'This may take a while based on the number of iterations provided...',
                     value = 0.3 )
      }
    })  
    
    
    
    #Reactive function to download log concentration plots       
    user_singlelogrunplot_paralog <- reactive({
      singleConc <- Conc_for_user_param()
      
      ObsData <- data()
      colors <- c("Predicted" = "blue", "Observed" = "red")
      
      plot_Cbbsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,2])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cbb), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbb") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cbmsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,3])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cbm), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbm") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      
      plot_Cccsfsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,4])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cccsf), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cccsf") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cscsfsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,5])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cscsf), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cscsf") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16))      
      
      
      singloglist <- list(print(plot_Cbbsin), print(plot_Cbmsin), print(plot_Cccsfsin), print(plot_Cscsfsin))
      
      # AlllogSingle <- grid.arrange(grobs = singloglist ,ncol=2)
      # print(AlllogSingle)
      
      dev.off()
      grid.arrange(grobs = singloglist)
      
    })
    
    
    #To download log concentration plots 
    output$userconclotplt <- downloadHandler(
      filename <- function(){ "Concentration_plots.pdf"} ,
      content = function(file) {
        if (input$Resinput>0) { 
          withProgress(message = 'Downloading in progress...',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          pdf(file,paper="default")
          plot(user_singlelogrunplot_paralog())
          dev.off()
        }
      }
    )       
    
    
    
    
    #Reactive function to generate PK parameters for user given parameters
    user_PKdata1_para <- reactive({
      PkData_2nd <- Conc_for_user_param()
      PkData <- subset(PkData_2nd,  PkData_2nd$Time >= input$user_PKmintime & PkData_2nd$Time <=input$user_PKmaxtime)
      
      Cmax <- PkData %>% summarise_all(max)
      Cmax<- Cmax[,-1]
      
      t_max<- apply(PkData[,-1], 2, function(x) 
        PkData$Time[x %in% max(x)])
      
      
      maxt <- as.data.frame(t_max)
      Cbb<- maxt$t_max[1]
      Cbm<- maxt$t_max[2]
      Cccsf<- maxt$t_max[3]
      Cscsf<- maxt$t_max[4]
      Cart<- maxt$t_max[5]
      
      Tmax <- data.frame(Cbb,Cbm,Cccsf,Cscsf,Cart)
      
      Cbb_AUC<- AUC(x=PkData$Time, y=PkData$Cbb,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      Cbm_AUC<- AUC(x=PkData$Time, y=PkData$Cbm,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      Cccsf_AUC<- AUC(x=PkData$Time, y=PkData$Cccsf,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      Cscsf_AUC<- AUC(x=PkData$Time, y=PkData$Cscsf,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      Cart_AUC<- AUC(x=PkData$Time, y=PkData$Cart,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      AUC1 <- data.frame(Cbb_AUC,Cbm_AUC,Cccsf_AUC,Cscsf_AUC,Cart_AUC)
      
      names(AUC1)[1] <- "Cbb"
      names(AUC1)[2] <- "Cbm"
      names(AUC1)[3] <- "Cccsf"
      names(AUC1)[4] <- "Cscsf"
      names(AUC1)[5] <- "Cart"
      
      PK_parameter <- c("C_max","T_max","AUC")
      ParaData <- rbind(Cmax,Tmax,AUC1)
      
      PK_Data <- cbind(PK_parameter,ParaData)
      PK_Data
      
    })
    
    
    #To output the PK parameter table to main panel
    output$userPK <- renderTable({
      if (input$Resinput>0) { 
        withProgress( user_PKdata1_para() ,message = 'Calculation in progress...', 
                      detail= 'This may take a while based on the number of iterations provided...',
                      value = 0.3 )
      }
    })  
    
    
    
    #To downolad the PK parameter table 
    output$userpkparm <- downloadHandler(
      filename = "pkparameters.csv",
      content = function(file) {
        if (input$Resinput>0) { 
          
          withProgress(message = 'Downloading in progress',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          write.csv(user_PKdata1_para(), file, row.names = FALSE)
        }
      }
    )   
    
################################################################################        
################################################################################    
################################################################################   
################################################################################
##### Parameter estimation    
##### This reactive function is to define the loss function to estimate parameters
    
    Parameter_estimationB4 <- reactive({
      
      Loss.4comp <- function(param,conc,time,Obs_Cbb,Obs_Cbm,Obs_Cccsf,Obs_Cscsf){
        ttimes0 = time
        tconcs0 = conc
        bparms = as.numeric(param)
        artconc = approxfun(x=ttimes0,y=tconcs0,method="linear",rule=2)
        btimes = ttimes0
        bstart <- c(Cbb=0,Cbm=0,Cccsf=0,Cscsf=0)
        solve4comp = ode(
          func=mb4.model
          ,y=bstart
          ,times= btimes
          ,parms=bparms
          ,artconc=artconc
          ,method = lsoda
        )
        
        Simdata = as.data.frame(cbind(solve4comp,tconcs0))
        dimnames(Simdata)[[2]] = c("Time","Cbb","Cbm","Cccsf","Cscsf","Cart")
        
        
        # loss <- 0
        loss <- c()
        for (i in btimes) {
          Cbb_Obs <- Obs_Cbb[i]
          Cbb_Sim <- Simdata[i,2]
          
          Cbm_Obs <- Obs_Cbm[i]
          Cbm_Sim <- Simdata[i,3]
          
          Cccsf_Obs <- Obs_Cccsf[i]
          Cccsf_Sim <- Simdata[i,4]
          
          Cscsf_Obs <- Obs_Cscsf[i]
          Cscsf_Sim <- Simdata[i,5]
          
          residuals_brain4 <- (  (Cbb_Obs - Cbb_Sim)^2 + 
                                   (Cbm_Obs - Cbm_Sim)^2 +
                                   (Cccsf_Obs -  Cccsf_Sim)^2 + 
                                   (Cscsf_Obs - Cscsf_Sim)^2)
          
          
          loss[i] <-  residuals_brain4
          
        }
        
        Total_Loss <- sum(loss)
        Total_Loss
        
      }  
      
      
      
      ptime <- na.omit(data()[,input$isim4_timep])
      pconc <- na.omit(data()[,input$isim4_cartp])  
      
      Pdata <- data()
      
      lower1 <- c(input$Vbbmin,
                  input$Vbmmin,
                  input$Vccsfmin,
                  input$Vscsfmin,
                  input$Qbrainmin,
                  input$QbulkBCmin,
                  input$QbulkCBmin,
                  input$QCsinkmin,
                  input$QSsinkmin,
                  input$QSinmin,
                  
                  input$QSoutmin,
                  input$Qglybmmin,
                  input$QglyCCSFmin,
                  input$QglySCSFmin,
                  input$PSBmin,
                  input$PSCmin,
                  input$PSEmin,
                  input$CLeffbmmin,
                  input$CLupbmmin,
                  input$CLeffccsfmin,
                  
                  input$CLupccsfmin,
                  input$lamdabbmin,
                  input$lamdabmmin,
                  input$lamdaccsfmin,
                  input$fubbmin,
                  input$fubmmin,
                  input$fuccsfmin,
                  input$fuscsfmin
      )
      
      
      
      upper1 <- c(input$Vbbmax,
                  input$Vbmmax,
                  input$Vccsfmax,
                  input$Vscsfmax,
                  input$Qbrainmax,
                  input$QbulkBCmax,
                  input$QbulkCBmax,
                  input$QCsinkmax,
                  input$QSsinkmax,
                  input$QSinmax,
                  
                  input$QSoutmax,
                  input$Qglybmmax,
                  input$QglyCCSFmax,
                  input$QglySCSFmax,
                  input$PSBmax,
                  input$PSCmax,
                  input$PSEmax,
                  input$CLeffbmmax,
                  input$CLupbmmax,
                  input$CLeffccsfmax,
                  
                  input$CLupccsfmax,
                  input$lamdabbmax,
                  input$lamdabmmax,
                  input$lamdaccsfmax,
                  input$fubbmax,
                  input$fubmmax,
                  input$fuccsfmax,
                  input$fuscsfmax
      )
      
      
      TotalnoiceCbb <- Pdata$CbbS
      TotalnoiceCbm <- Pdata$CbmS
      TotalnoiceCccsf <- Pdata$CccsfS
      TotalnoiceCscsf <- Pdata$CscsfS
      
      Max_Iteration <- input$maxiter
      RelTol <- input$rel_tol
      VTRVal <- input$VTR_Value
      
      set.seed(1234)
      myestimations <- DEoptim(fn=Loss.4comp,
                               lower = lower1,
                               upper = upper1,
                               control = DEoptim.control(VTR = VTRVal, itermax = Max_Iteration,
                                                         reltol=RelTol),
                               conc=pconc,
                               time=ptime,
                               Obs_Cbb=TotalnoiceCbb,
                               Obs_Cbm=TotalnoiceCbm,
                               Obs_Cccsf=TotalnoiceCccsf,
                               Obs_Cscsf=TotalnoiceCscsf,
                               fnMap=NULL)
       myestimations$optim$bestmem
    }) 
    
    
# Generate parameter table
    Parameter_table <- reactive({
      panames = c(
        "Vbb","Vbm","Vccsf","Vscsf"
        ,"Qbrain","Qbulk,BC","Qbulk,CB"
        ,"QCsink","QSsink","QSin","QSout"
        ,"Qgly,bm","Qgly,CCSF","Qgly,SCSF","PSB","PSC","PSE"
        ,"CLeff,bm","CLup,bm","CLeff,ccsf","CLup,ccsf"
        ,"lambdabb","lambdabm","lambdaccsf"
        ,"fubb","fubm","fuccsf","fuscsf"
      )
      plist <- Parameter_estimationB4()
      Ptablefinal = data.frame(Parameter = panames,Estimated_value = plist)
      Ptablefinal 
    }) 
    
  
#To output the estimated parameter table to main panel
    output$para_table <- renderTable({
      if (input$subPtable > 0) { 
        
        withProgress(Parameter_table()  ,message = 'Calculation in progress...', 
                      detail= 'This may take a while based on the number of iterations provided...',
                      value = 0.3 )
      }
    })  
    
#To download the estimated parameter table 
    output$paraestidata <- downloadHandler(
      filename = "Estimatedparameters.csv",
      content = function(file) {
        if (input$subPtable>0) { 
          
          withProgress(message = 'Downloading in progress',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          write.csv(Parameter_table(), file, row.names = FALSE)
        }
      }
    )
    

    
#To get concentration table based on the estimated parameters
    sim4singlerun_para <- reactive({
      itime <- na.omit(data()[,input$isim4_timep])
      iconc <- na.omit(data()[,input$isim4_cartp])
      iparam <- Parameter_estimationB4()
      sinrun <- run.sim4(time=itime,conc=iconc,param=iparam)
      sinrun 
    }) 
    
#To output the concentrations table to main panel
    output$singlerun.sim4pesti <- renderTable({
      if (input$subPtable>0) { 
        withProgress(sim4singlerun_para()  ,message = 'Calculation in progress...', 
                     detail= 'This may take a while based on the number of iterations provided...',
                     value = 0.3 )
      }
    })

 
#To download the concentrations table 
    output$paracondata <- downloadHandler(
      filename = "ConcentationTable.csv",
      content = function(file) {
        if (input$subPtable>0) { 
          
          withProgress(message = 'Downloading in progress',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          write.csv(sim4singlerun_para(), file, row.names = FALSE)
        }
      }
    )   
    
    
#Reactive function to get concentrations profiles to main panel
    
    singlerunplot_para <- reactive({
      singleConc <- sim4singlerun_para()
      ObsData <- data()
      colors <- c("Predicted" = "blue", "Observed" = "red")
      
      plot_Cbbsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,2]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cbb, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbb") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cbmsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,3]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cbm, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbm") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      
      plot_Cccsfsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,4]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cccsf, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cccsf") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cscsfsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,5]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cscsf, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cscsf") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16))      
      
      
      singlist <- list(plot_Cbbsin, plot_Cbmsin,plot_Cccsfsin,plot_Cscsfsin)
      
      grid.arrange(grobs = singlist, ncol = 2, nrow=2)
      #AllSingle <- grid.arrange(grobs = singlist)
      # print(AllSingle)
    })
    
    
#To get concentration profiles to main panel
    output$sinplot.sim4pesti <- renderPlot({ 
      if (input$subPtable>0) { 
        withProgress(singlerunplot_para()  ,message = 'Calculation in progress...', 
                     detail= 'This may take a while based on the number of iterations provided...',
                     value = 0.3 )
      }
    }) 
    

#Reactive function to download concentration profiles
    
    singlerunplotdown <- reactive({
      singleConc <- sim4singlerun_para()
      
      ObsData <- data()
      colors <- c("Predicted" = "blue", "Observed" = "red")
      
      plot_Cbbsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,2]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cbb, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbb") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "bottom",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cbmsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,3]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cbm, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbm") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "bottom",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      
      plot_Cccsfsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,4]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cccsf, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cccsf") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "bottom",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cscsfsin <- ggplot(singleConc, aes(x = Time,y = singleConc[,5]) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = Cscsf, color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cscsf") + 
        theme_light() + labs(x = "Time",y = "Concentration",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "bottom",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16))      
      
      
      singlist <- list(plot_Cbbsin, plot_Cbmsin,plot_Cccsfsin,plot_Cscsfsin)
      
      #AllSingle <- grid.arrange(grobs = singlist)
      
      #dev.off()
      #grid.arrange(grobs = singlist)
      
      
      
      
      singlist <- list(print(plot_Cbbsin),print(plot_Cbmsin),print(plot_Cccsfsin),print(plot_Cscsfsin))
      
      dev.off()
      
      grid.arrange(grobs = singlist)
      
    })
    
    
#To  download the concentration profiles
    output$paraconplot <- downloadHandler(
      filename <- function(){ "Concentration_plots.pdf"} ,
      content = function(file) {
        if (input$subPtable>0) { 
          withProgress(message = 'Downloading in progress...',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          pdf(file,paper="default")
          plot(singlerunplotdown())
          dev.off()
        }
      }
    )     
    
 
    
    
#Reactive function to get log concentration plots  to main panel     
    singlelogrunplot_para <- reactive({
      singleConc <- sim4singlerun_para()
      
      ObsData <- data()
      colors <- c("Predicted" = "blue", "Observed" = "red")
      
      plot_Cbbsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,2])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cbb), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbb") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cbmsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,3])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cbm), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbm") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      
      plot_Cccsfsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,4])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cccsf), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cccsf") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cscsfsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,5])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cscsf), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cscsf") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16))      
      
      singloglist <- list(plot_Cbbsin, plot_Cbmsin,plot_Cccsfsin,plot_Cscsfsin)
      
   #   AlllogSingle <- grid.arrange(grobs = singloglist ,ncol=2)
   #   print(AlllogSingle)
      grid.arrange(grobs = singloglist,ncol=2, nrow =2)
      
      

      

    })  
    
    
    
    
    #To get concentration log profiles to main panel
    output$sinlogplot.sim4pesti <- renderPlot({ 
      if (input$subPtable>0) { 
        withProgress(singlelogrunplot_para()  ,message = 'Calculation in progress...', 
                     detail= 'This may take a while based on the number of iterations provided...',
                     value = 0.3 )
      }
    })  
    
    
    
#Reactive function to download log concentration plots       
    singlelogrunplot_paralog <- reactive({
      singleConc <- sim4singlerun_para()
      
      ObsData <- data()
      colors <- c("Predicted" = "blue", "Observed" = "red")
      
      plot_Cbbsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,2])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cbb), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbb") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cbmsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,3])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cbm), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cbm") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      
      plot_Cccsfsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,4])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cccsf), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cccsf") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16)) 
      
      plot_Cscsfsin <- ggplot(singleConc, aes(x = Time,y = log10(singleConc[,5])) ) + geom_line(size=0.5,aes(color = "Predicted"))+
        geom_point(ObsData,mapping=aes(x = Time ,y = log10(Cscsf), color = "Observed") ,shape =1, fill = "white", size = 2, stroke = 0.5)+
        ggtitle("Cscsf") + 
        theme_light() + labs(x = "Time",y = "Log10(Concentration)",color = "Concentrations") + scale_color_manual(values = colors) +
        theme(legend.position = "right",plot.title = element_text(size = 16, face = "bold"),
              legend.title=element_text(size=16), text = element_text(size=16),
              legend.text=element_text(size=16))      
      
      
      singloglist <- list(print(plot_Cbbsin), print(plot_Cbmsin), print(plot_Cccsfsin), print(plot_Cscsfsin))
      
     # AlllogSingle <- grid.arrange(grobs = singloglist ,ncol=2)
     # print(AlllogSingle)
      
      dev.off()
      grid.arrange(grobs = singloglist)
      
    })
    
    
#To download log concentration plots 
    output$paraconlogplot <- downloadHandler(
      filename <- function(){ "Concentration_plots.pdf"} ,
      content = function(file) {
        if (input$subPtable>0) { 
          withProgress(message = 'Downloading in progress...',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          pdf(file,paper="default")
          plot(singlelogrunplot_paralog())
          dev.off()
        }
      }
    )     
    
    
    
#Reactive function to generate PK parameters
    PKdata1_para <- reactive({
      PkData_2nd <- sim4singlerun_para()
      PkData <- subset(PkData_2nd,  PkData_2nd$Time >= input$PKmintime & PkData_2nd$Time <=input$PKmaxtime)
      
      Cmax <- PkData %>% summarise_all(max)
      Cmax<- Cmax[,-1]
      
      t_max<- apply(PkData[,-1], 2, function(x) 
        PkData$Time[x %in% max(x)])
      
      
      maxt <- as.data.frame(t_max)
      Cbb<- maxt$t_max[1]
      Cbm<- maxt$t_max[2]
      Cccsf<- maxt$t_max[3]
      Cscsf<- maxt$t_max[4]
      Cart<- maxt$t_max[5]
      
      Tmax <- data.frame(Cbb,Cbm,Cccsf,Cscsf,Cart)
      
      Cbb_AUC<- AUC(x=PkData$Time, y=PkData$Cbb,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      Cbm_AUC<- AUC(x=PkData$Time, y=PkData$Cbm,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      Cccsf_AUC<- AUC(x=PkData$Time, y=PkData$Cccsf,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      Cscsf_AUC<- AUC(x=PkData$Time, y=PkData$Cscsf,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      Cart_AUC<- AUC(x=PkData$Time, y=PkData$Cart,from = min(PkData$Time, na.rm = TRUE), to = max(x=PkData$Time, na.rm = TRUE),method = "trapezoid")
      
      AUC1 <- data.frame(Cbb_AUC,Cbm_AUC,Cccsf_AUC,Cscsf_AUC,Cart_AUC)
      
      names(AUC1)[1] <- "Cbb"
      names(AUC1)[2] <- "Cbm"
      names(AUC1)[3] <- "Cccsf"
      names(AUC1)[4] <- "Cscsf"
      names(AUC1)[5] <- "Cart"
      
      PK_parameter <- c("C_max","T_max","AUC")
      ParaData <- rbind(Cmax,Tmax,AUC1)
      
      PK_Data <- cbind(PK_parameter,ParaData)
      PK_Data
      
    })
    
    
#To output the PK parameter table to main panel
    output$pkparameters.sim4pesti <- renderTable({
      if (input$subPtable>0) { 
        withProgress( PKdata1_para() ,message = 'Calculation in progress...', 
                      detail= 'This may take a while based on the number of iterations provided...',
                      value = 0.3 )
      }
    })  
    
    
    
#To downolad the PK parameter table 
    output$parmako4 <- downloadHandler(
      filename = "pkparameters.csv",
      content = function(file) {
        if (input$subPtable>0) { 
          
          withProgress(message = 'Downloading in progress',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          write.csv(PKdata1_para(), file, row.names = FALSE)
        }
      }
    )   
    
    

################################################################################
################################################################################
################################################################################
#######  sensitivity analysis
    
#To get the parameter list from data file
iparamtable4 <- reactive({
      if(input$isim4_param==""){

        shinyalert("Oops!", "Please select a list of parameters", type = "error")
        forceInput <- NULL
        
      }else{
        iparam <- na.omit(data()[,input$isim4_param])
        
      }
      

      if(length(iparam)!=27){
        shinyalert("Oops!", "Some parameters are missing", type = "error")
        forceInput <- NULL
      }
      
      iparam
    })
    
    
#To get parameter list with the parameter names 
    iparamout4 <- reactive({
      panames = c(
        "Vbb","Vbm","Vccsf","Vscsf"
        ,"Qbrain","Qbulk,BC","Qbulk,CB"
        ,"QCsink","QSsink","QSin","QSout"
        ,"Qgly,bm","Qgly,CCSF","Qgly,SCSF","PSB","PSC","PSE"
        ,"CLeff,bm","CLup,bm","CLeff,ccsf","CLup,ccsf"
        ,"lambdabb","lambdabm","lambdaccsf"
        ,"fubb","fubm","fuccsf"
      )
      tpatab = data.frame(Parameter=panames,Value = iparamtable4())
      tpatab
    })
    
    
    
    
    

    ## Following reactive function is to output a data table for concentrations for 
    ## individually selected parameters for sensitivity analysis for brain4  
    
    
    
    Indiparadata4 <- reactive({
      itime <- na.omit(data()[,input$Sen_time])
      iconc <- na.omit(data()[,input$Sen_cart])
      
      pa4list <- iparamout4()
      para4list<- as.data.frame(pa4list)
      
      mylist4 <- seq(input$sim4.Pickminpara, input$sim4.Pickmaxpara, length.out= input$sim4.Pickparanum)
      mylist4<- as.data.frame(mylist4)
      
      
      
      Cbb_con <- matrix(0,nrow = length(itime),ncol = input$sim4.Pickparanum+1)
      Cbm_con <- matrix(0,nrow = length(itime),ncol = input$sim4.Pickparanum+1)
      Cccsf_con <- matrix(0,nrow = length(itime),ncol = input$sim4.Pickparanum+1)
      Cscsf_con <- matrix(0,nrow = length(itime),ncol = input$sim4.Pickparanum+1)
      Cart_con <- matrix(0,nrow = length(itime),ncol = input$sim4.Pickparanum+1)
      
      
      cc<- seq(input$sim4.Pickminpara, input$sim4.Pickmaxpara, length.out= input$sim4.Pickparanum)
      colnames(Cbb_con) <- c("Time",cc )
      colnames(Cbm_con) <- c("Time",cc )
      colnames(Cccsf_con) <- c("Time",cc )
      colnames(Cscsf_con) <- c("Time",cc )
      colnames(Cart_con) <- c("Time",cc )
      
      
      
      ulist <- c(
        "Vbb","Vbm","Vccsf","Vscsf"
        ,"Qbrain","Qbulk,BC","Qbulk,CB"
        ,"QCsink","QSsink","QSin","QSout"
        ,"Qgly,bm","Qgly,CCSF","Qgly,SCSF","PSB","PSC","PSE"
        ,"CLeff,bm","CLup,bm","CLeff,ccsf","CLup,ccsf"
        ,"lambdabb","lambdabm","lambdaccsf"
        ,"fubb","fubm","fuccsf"
      )
      
      
      
      for (i in 1:input$sim4.Pickparanum) { 
        
        
        #  This if else is to assign a value to user input paremeter
        if ( input$Pickpara4 == ulist[1] ) {j=1} 
        else if (input$Pickpara4 == ulist[2]) {j=2} 
        else if ( input$Pickpara4 == ulist[3]) {j=3} 
        else if ( input$Pickpara4 == ulist[4]) {j=4} 
        else if ( input$Pickpara4 == ulist[5]) {j=5} 
        else if ( input$Pickpara4 == ulist[6]) {j=6} 
        else if ( input$Pickpara4 == ulist[7]) {j=7} 
        else if ( input$Pickpara4 == ulist[8]) {j=8} 
        else if ( input$Pickpara4 == ulist[9]) {j=9} 
        else if ( input$Pickpara4 == ulist[10]) {j=10} 
        else if ( input$Pickpara4 == ulist[11]) {j=11} 
        else if ( input$Pickpara4 == ulist[12]) {j=12} 
        else if ( input$Pickpara4 == ulist[13]) {j=13} 
        else if ( input$Pickpara4 == ulist[14]) {j=14} 
        else if ( input$Pickpara4 == ulist[15]) {j=15} 
        else if ( input$Pickpara4 == ulist[16]) {j=16} 
        else if ( input$Pickpara4 == ulist[17]) {j=17} 
        else if ( input$Pickpara4 == ulist[18]) {j=18} 
        else if ( input$Pickpara4 == ulist[19]) {j=19} 
        else if ( input$Pickpara4 == ulist[20]) {j=20} 
        else if ( input$Pickpara4 == ulist[21]) {j=21} 
        else if ( input$Pickpara4 == ulist[22]) {j=22} 
        else if ( input$Pickpara4 == ulist[23]) {j=23} 
        else if ( input$Pickpara4 == ulist[24]) {j=24} 
        else if ( input$Pickpara4 == ulist[25]) {j=25} 
        else if ( input$Pickpara4 == ulist[26]) {j=26} 
        else{j=27}
        
        
        
        
        newpa4 <- replace(para4list$Value, j ,mylist4[i,1])
       
        
        
        ConSen4i <- sen4(time = itime,conc= iconc,param=newpa4)
        
        #ConSen4i <- format(round(ConSen4i, 4), nsmall = 4)
        
        
        
        Cbbs <- ConSen4i[,2]
        Cbms <- ConSen4i[,3]
        Cccsfs  <- ConSen4i[,4]
        Cscsfs  <- ConSen4i[,5]
        Carts  <- ConSen4i[,6]
        
        
        
        
        Cbb_con[,1] <- ConSen4i[,1]   # This is to save a time column for cbb concentrations
        Cbb_con[,i+1] <- Cbbs         # This will save the cbb concentrations starting from the second column
        
        Cbm_con[,1] <- ConSen4i[,1]
        Cbm_con[,i+1] <- Cbms 

        Cccsf_con[,1] <- ConSen4i[,1]
        Cccsf_con[,i+1] <- Cccsfs   
        
        
        Cscsf_con[,1] <- ConSen4i[,1]
        Cscsf_con[,i+1] <- Cscsfs     
        
        Cart_con[,1] <- ConSen4i[,1]
        Cart_con[,i+1] <- Carts
        
        
      }
      
      
      
      Cbb_con<- as.data.frame(Cbb_con)
      Cbm_con<- as.data.frame(Cbm_con)
      Cccsf_con<- as.data.frame(Cccsf_con)
      Cscsf_con<- as.data.frame(Cscsf_con)
      Cart_con<- as.data.frame(Cart_con)
      
      
      
      
      Allcons <- list(Cbb_con,Cbm_con,Cccsf_con,Cscsf_con,Cart_con)
      Allcons
    })
    
    
    
    output$sim4.sendata4a <- renderTable({
      if (input$sim4.Sensub>0) { 
        
        withProgress(message = 'Calculation in progress...',
                     value = 0, {
                       for (i in 1:2) {
                         incProgress(1/2)
                         Sys.sleep(0.25)
                       }
                     })
        
        
        Indiparadata4() 
      }
    })
    

    output$sim4.sencondata <- downloadHandler(
      filename = "Sen_Conc_Data.csv",
      content = function(file) {
        
        if (input$sim4.Sensub>0) {
          
          
          withProgress(message = 'Download in progress...',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          
          
        write.csv(Indiparadata4(), file, row.names = FALSE)
        }
        
      }
    )
    

    
# Following reactive function is to generate T max and C max for user selected time range under sensitivity for barin 4
# CTAinfo means Cmax_Tmax_AUC  
  
  CTAinfo <- reactive({
    CTAtable4 <- Indiparadata4()
    
    
    #Cbb cmax tmax auc table
    
    Cbbinfo4 <- CTAtable4[[1]]
    
    Cbb4_user <- subset(Cbbinfo4,  Cbbinfo4$Time >= input$PKtimemin & Cbbinfo4$Time <=input$PKtimemax)

    Cbb_CMax4 <- Cbb4_user %>% summarise_all(max)
    Cbb_CMax4 <- Cbb_CMax4[,-1]
    Cbb_CMax4 <- as.data.frame(Cbb_CMax4)
    
    Tmax_Cbb <- apply(Cbb4_user[,-1], 2, function(x) Cbb4_user$Time[x %in% max(x)])
    Cbb_Tmax4 <- as.data.frame(t(Tmax_Cbb))
    
    AUC_CBB4 <- c()
    for (m in 1:input$sim4.Pickparanum) {
      AUC_CBB_Data4 <- AUC(x=Cbb4_user$Time, y=Cbb4_user[,m+1],from = min(Cbb4_user$Time, na.rm = TRUE), to = max(x=Cbb4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_CBB4[m] <- AUC_CBB_Data4
      
    }
    
    Cbb_AUC4 <- as.data.frame(t(AUC_CBB4))
    colnames(Cbb_AUC4) <- seq(input$sim4.Pickminpara, input$sim4.Pickmaxpara, length.out= input$sim4.Pickparanum)
    
    Cbb_PK_All <- rbind(Cbb_CMax4,Cbb_Tmax4,Cbb_AUC4)
    Cbb_PK_All %>% mutate_if(is.numeric, round,digits=3)
    Parameter_values <- c("Cmax of Cbb", "Tmax of Cbb", "AUC of Cbb")
    
    CBB_PK4 <- cbind(Parameter_values,Cbb_PK_All)
    rows4 <- nrow(CBB_PK4)
    CBB_PK4[rows4+1,1] <-NA

    
    #Cbm cmax tmax auc table

    Cbminfo4 <- CTAtable4[[2]]
    
    Cbm4_user <- subset(Cbminfo4,  Cbminfo4$Time >= input$PKtimemin & Cbminfo4$Time <=input$PKtimemax)
    
    Cbm_CMax4 <- Cbm4_user %>% summarise_all(max)
    Cbm_CMax4 <- Cbm_CMax4[,-1]
    Cbm_CMax4 <- as.data.frame(Cbm_CMax4)
    
    Tmax_Cbm <- apply(Cbm4_user[,-1], 2, function(x) Cbm4_user$Time[x %in% max(x)])
    Cbm_Tmax4 <- as.data.frame(t(Tmax_Cbm))
    
    AUC_CBM4 <- c()
    for (m in 1:input$sim4.Pickparanum) {
      AUC_CBM_Data4 <- AUC(x=Cbm4_user$Time, y=Cbm4_user[,m+1],from = min(Cbm4_user$Time, na.rm = TRUE), to = max(x=Cbm4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_CBM4[m] <- AUC_CBM_Data4
      
    }
    
    Cbm_AUC4 <- as.data.frame(t(AUC_CBM4))
    colnames(Cbm_AUC4) <- seq(input$sim4.Pickminpara, input$sim4.Pickmaxpara, length.out= input$sim4.Pickparanum)
    
    Cbm_PK_All <- rbind(Cbm_CMax4,Cbm_Tmax4,Cbm_AUC4)
    Cbm_PK_All %>% mutate_if(is.numeric, round,digits=3)
    Parameter_values <- c("Cmax of Cbm", "Tmax of Cbm", "AUC of Cbm")
    
    CBM_PK4 <- cbind(Parameter_values,Cbm_PK_All)
    rows4 <- nrow(CBM_PK4)
    CBM_PK4[rows4+1,1] <-NA
    
    
    #Cccsf cmax tmax auc table
    
    
    Cccsfinfo4 <- CTAtable4[[3]]
    
    Cccsf4_user <- subset(Cccsfinfo4,  Cccsfinfo4$Time >= input$PKtimemin & Cccsfinfo4$Time <=input$PKtimemax)
    
    Cccsf_CMax4 <- Cccsf4_user %>% summarise_all(max)
    Cccsf_CMax4 <- Cccsf_CMax4[,-1]
    Cccsf_CMax4 <- as.data.frame(Cccsf_CMax4)
    
    Tmax_Cccsf <- apply(Cccsf4_user[,-1], 2, function(x) Cccsf4_user$Time[x %in% max(x)])
    Cccsf_Tmax4 <- as.data.frame(t(Tmax_Cccsf))
    
    AUC_Cccsf4 <- c()
    for (m in 1:input$sim4.Pickparanum) {
      AUC_Cccsf_Data4 <- AUC(x=Cccsf4_user$Time, y=Cccsf4_user[,m+1],from = min(Cccsf4_user$Time, na.rm = TRUE), to = max(x=Cccsf4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_Cccsf4[m] <- AUC_Cccsf_Data4
      
    }
    
    Cccsf_AUC4 <- as.data.frame(t(AUC_Cccsf4))
    colnames(Cccsf_AUC4) <- seq(input$sim4.Pickminpara, input$sim4.Pickmaxpara, length.out= input$sim4.Pickparanum)
    
    Cccsf_PK_All <- rbind(Cccsf_CMax4,Cccsf_Tmax4,Cccsf_AUC4)
    Cccsf_PK_All %>% mutate_if(is.numeric, round,digits=3)
    Parameter_values <- c("Cmax of Cccsf", "Tmax of Cccsf", "AUC of Cccsf")
    
    Cccsf_PK4 <- cbind(Parameter_values,Cccsf_PK_All)
    rows4 <- nrow(Cccsf_PK4)
    Cccsf_PK4[rows4+1,1] <-NA
    
    
    
    
    
    #Cscsf cmax tmax auc table
    
    
    Cscsfinfo4 <- CTAtable4[[4]]
    
    Cscsf4_user <- subset(Cscsfinfo4,  Cscsfinfo4$Time >= input$PKtimemin & Cscsfinfo4$Time <=input$PKtimemax)
    
    Cscsf_CMax4 <- Cscsf4_user %>% summarise_all(max)
    Cscsf_CMax4 <- Cscsf_CMax4[,-1]
    Cscsf_CMax4 <- as.data.frame(Cscsf_CMax4)
    
    Tmax_Cscsf <- apply(Cscsf4_user[,-1], 2, function(x) Cscsf4_user$Time[x %in% max(x)])
    Cscsf_Tmax4 <- as.data.frame(t(Tmax_Cscsf))
    
    AUC_Cscsf4 <- c()
    for (m in 1:input$sim4.Pickparanum) {
      AUC_Cscsf_Data4 <- AUC(x=Cscsf4_user$Time, y=Cscsf4_user[,m+1],from = min(Cscsf4_user$Time, na.rm = TRUE), to = max(x=Cscsf4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_Cscsf4[m] <- AUC_Cscsf_Data4
      
    }
    
    Cscsf_AUC4 <- as.data.frame(t(AUC_Cscsf4))
    colnames(Cscsf_AUC4) <- seq(input$sim4.Pickminpara, input$sim4.Pickmaxpara, length.out= input$sim4.Pickparanum)
    
    Cscsf_PK_All <- rbind(Cscsf_CMax4,Cscsf_Tmax4,Cscsf_AUC4)
    Cscsf_PK_All %>% mutate_if(is.numeric, round,digits=3)
    Parameter_values <- c("Cmax of Cscsf", "Tmax of Cscsf", "AUC of Cscsf")
    
    Cscsf_PK4 <- cbind(Parameter_values,Cscsf_PK_All)
    rows4 <- nrow(Cscsf_PK4)
    Cscsf_PK4[rows4+1,1] <-NA
    
 
    
    
    
    #Cart cmax tmax auc table
    
    Cartinfo4 <- CTAtable4[[5]]
    
    Cart4_user <- subset(Cartinfo4,  Cartinfo4$Time >= input$PKtimemin & Cartinfo4$Time <=input$PKtimemax)
    
    Cart_CMax4 <- Cart4_user %>% summarise_all(max)
    Cart_CMax4 <- Cart_CMax4[,-1]
    Cart_CMax4 <- as.data.frame(Cart_CMax4)
    
    Tmax_Cart <- apply(Cart4_user[,-1], 2, function(x) Cart4_user$Time[x %in% max(x)])
    Cart_Tmax4 <- as.data.frame(t(Tmax_Cart))
    
    AUC_Cart4 <- c()
    for (m in 1:input$sim4.Pickparanum) {
      AUC_Cart_Data4 <- AUC(x=Cart4_user$Time, y=Cart4_user[,m+1],from = min(Cart4_user$Time, na.rm = TRUE), to = max(x=Cart4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_Cart4[m] <- AUC_Cart_Data4
      
    }
    
    Cart_AUC4 <- as.data.frame(t(AUC_Cart4))
    colnames(Cart_AUC4) <- seq(input$sim4.Pickminpara, input$sim4.Pickmaxpara, length.out= input$sim4.Pickparanum)
    
    Cart_PK_All <- rbind(Cart_CMax4,Cart_Tmax4,Cart_AUC4)
    Cart_PK_All %>% mutate_if(is.numeric, round,digits=3)
    Parameter_values <- c("Cmax of Cart", "Tmax of Cart", "AUC of Cart")
    
    Cart_PK4 <- cbind(Parameter_values,Cart_PK_All)
    rows4 <- nrow(Cart_PK4)
    Cart_PK4[rows4+1,1] <-NA 
    
    
   
    AllPKinfo4 <- rbind(CBB_PK4,CBM_PK4,Cccsf_PK4,Cscsf_PK4,Cart_PK4)
    AllPKinfo4
    
    
    

  })
    
    
    
  output$sim4.tmax4 <- renderTable({
    if (input$sim4.Sensub>0){
      
      withProgress(message = 'Calculation in progress...',
                   value = 0, {
                     for (i in 1:2) {
                       incProgress(1/2)
                       Sys.sleep(0.25)
                     }
                   })
      
      
      CTAinfo()
    }
    
    
    
  })  
    
    
  output$sim4.parmako4 <- downloadHandler(
    filename = "CTA_PK.csv",
    content = function(file) {
      
      if (input$sim4.Sensub>0) {
        
        withProgress(message = 'Download in progress...',
                     value = 0, {
                       for (i in 1:2) {
                         incProgress(1/2)
                         Sys.sleep(0.25)
                       }
                     })
        
        
        write.csv(CTAinfo(), file, row.names = FALSE)
      }
      
    }
  )
  
  
  
  
  # Following reactive function is to show the AUC scatter plots for brain4 model:
  
  
  ScatterAUC4 <- reactive({
    
    AUCtable4a <- Indiparadata4()
    
    Cbbinfo4 <- AUCtable4a[[1]]
    Cbb4_user <- subset(Cbbinfo4,  Cbbinfo4$Time >= input$PKtimemin & Cbbinfo4$Time <=input$PKtimemax)
    Cbb4_user %>% mutate_if(is.numeric, round, digits=2)
    
    Cbminfo4 <- AUCtable4a[[2]]
    Cbm4_user <- subset(Cbminfo4,  Cbminfo4$Time >= input$PKtimemin & Cbminfo4$Time <=input$PKtimemax)
    Cbm4_user %>% mutate_if(is.numeric, round, digits=2)
    
    Cccsfinfo4 <- AUCtable4a[[3]]
    Cccsf4_user <- subset(Cccsfinfo4,  Cccsfinfo4$Time >= input$PKtimemin & Cccsfinfo4$Time <=input$PKtimemax)
    Cccsf4_user %>% mutate_if(is.numeric, round, digits=2)
    
    
    Cscsfinfo4 <- AUCtable4a[[4]]
    Cscsf4_user <- subset(Cscsfinfo4,  Cscsfinfo4$Time >= input$PKtimemin & Cscsfinfo4$Time <=input$PKtimemax)
    Cscsf4_user %>% mutate_if(is.numeric, round, digits=2)
    
    
    Cartinfo4 <- AUCtable4a[[5]]
    Cart4_user <- subset(Cartinfo4,  Cartinfo4$Time >= input$PKtimemin & Cartinfo4$Time <=input$PKtimemax)
    Cart4_user %>% mutate_if(is.numeric, round, digits=2)
    
    
    AUC_CBB4 <- c()
    AUC_CBM4 <- c()
    AUC_Cccsf4 <- c()
    AUC_Cscsf4 <- c()
    AUC_Cart4 <- c()
    
    
    
    
    
    for (m in 1:input$sim4.Pickparanum) {
      
     # round(Cbb4_user,digits = 2)
      AUC_CBB_Data4 <- AUC(x=Cbb4_user$Time, y=Cbb4_user[,m+1],from = min(Cbb4_user$Time, na.rm = TRUE), to = max(x=Cbb4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_CBB4[m] <- AUC_CBB_Data4
     # AUC_CBB4[m] %>% mutate_if(is.numeric, round, digits=2)
      #round(AUC_CBB4,digits = 2)
      
      
      
      AUC_CBM_Data4 <- AUC(x=Cbm4_user$Time, y=Cbm4_user[,m+1],from = min(Cbm4_user$Time, na.rm = TRUE), to = max(x=Cbm4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_CBM4[m] <- AUC_CBM_Data4
      #AUC_CBM4[m] %>% mutate_if(is.numeric, round, digits=2)
      
      AUC_Cccsf_Data4 <- AUC(x=Cccsf4_user$Time, y=Cccsf4_user[,m+1],from = min(Cccsf4_user$Time, na.rm = TRUE), to = max(x=Cccsf4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_Cccsf4[m] <- AUC_Cccsf_Data4
     # AUC_Cccsf4[m] %>% mutate_if(is.numeric, round, digits=2)
      
      
      AUC_Cscsf_Data4 <- AUC(x=Cscsf4_user$Time, y=Cscsf4_user[,m+1],from = min(Cscsf4_user$Time, na.rm = TRUE), to = max(x=Cscsf4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_Cscsf4[m] <- AUC_Cscsf_Data4
    #  AUC_Cscsf4[m] %>% mutate_if(is.numeric, round, digits=2)
      
      AUC_Cart_Data4 <- AUC(x=Cart4_user$Time, y=Cart4_user[,m+1],from = min(Cart4_user$Time, na.rm = TRUE), to = max(x=Cart4_user$Time, na.rm = TRUE),method = "trapezoid")
      AUC_Cart4[m] <- AUC_Cart_Data4
    #  AUC_Cart4[m] %>% mutate_if(is.numeric, round, digits=2)
      
      
      
    }
    
    Cbb_AUC4 <- as.data.frame(AUC_CBB4)
    Cbb_AUC4 %>% mutate_if(is.numeric, round, digits=2)
    
    
    Cbm_AUC4 <- as.data.frame(AUC_CBM4)
    Cbm_AUC4 %>% mutate_if(is.numeric, round, digits=2)
    
    Cccsf_AUC4 <- as.data.frame(AUC_Cccsf4)
    Cccsf_AUC4 %>% mutate_if(is.numeric, round, digits=2)
    
    Cscsf_AUC4 <- as.data.frame(AUC_Cscsf4)
    Cscsf_AUC4 %>% mutate_if(is.numeric, round, digits=2)
    
    Cart_AUC4 <- as.data.frame(AUC_Cart4)
    Cart_AUC4 %>% mutate_if(is.numeric, round, digits=2)
    
    
    Parameter_Values <- seq(input$sim4.Pickminpara, input$sim4.Pickmaxpara, length.out= input$sim4.Pickparanum)
    Parameter_Values<- as.data.frame(Parameter_Values)
    
    AUC_Table <- data.frame( Parameter_Values, Cbb_AUC4, Cbm_AUC4,Cccsf_AUC4,Cscsf_AUC4,Cart_AUC4)
    #AUC_Table %>% mutate_if(is.numeric, round,digits=3)
   # AUC_Table %>% mutate(across(where(is.numeric), round, digits=3))
   # round_df(AUC_Table, digits=1, rf = "round")
   # round(AUC_Table, digits = 0)
    
    
    scaleFUN <- function(x) sprintf("%.2f", x) # This function is to round the y axis to two decimal places
    
    
    
    Cbb_Scatter <- ggplot(AUC_Table, aes(x = Parameter_Values, y = AUC_CBB4)) +
      geom_point(color='red',size = 3)+  theme_light() 
    
    Cbm_Scatter <- ggplot(AUC_Table, aes(x = Parameter_Values, y = AUC_CBM4)) + 
      geom_point(color='red',size = 3)+theme_light() 

    
    
    
    
    Cccsf_Scatter <- ggplot(AUC_Table, aes(x = Parameter_Values, y = AUC_Cccsf4)) +
      geom_point(color='red',size = 3)+theme_light() 
    
    Cscsf_Scatter <- ggplot(AUC_Table, aes(x = Parameter_Values, y = AUC_Cscsf4)) +
      geom_point(color='red',size = 3)+theme_light()     
    
    Cart_Scatter <- ggplot(AUC_Table, aes(x = Parameter_Values, y = AUC_Cart4)) +
      geom_point(color='red',size = 3)+theme_light()   
    
    
    senplotauc <- list(Cbb_Scatter,Cbm_Scatter,Cccsf_Scatter,Cscsf_Scatter,Cart_Scatter)
    #nelplot <- list(print(Cbb_Scatter),print(Cbm_Scatter),print(Cccsf_Scatter),print(Cscsf_Scatter),print(Cart_Scatter))
    
    
    All_Scatter<- grid.arrange(grobs = senplotauc,ncol=2)
 
    
  })
  

  
  output$sim4.aucplot <- renderPlot({
    if (input$sim4.Sensub>0) { 
      
      withProgress(message = 'Creating plots...',
                   value = 0, {
                     for (i in 1:2) {
                       incProgress(1/2)
                       Sys.sleep(0.25)
                     }
                   })
      
      
      ScatterAUC4() 
    }
  })
  
    
    # This reactive function is to "output plots" for concentration plots for individually 
    # selected parameters for sensitivity analysis into main panel
    
    
    
    Indiconcenplot4 <- reactive({
      
      ##### new idea to make the code efficient begin
      call_list4 <- Indiparadata4()
      #### new idea to make the code efficient end
      
    
      Cbbplot <- call_list4[[1]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cbbplota <- ggplot(Cbbplot, aes(x = Time, y = concentration, group = group)) +
        geom_line(aes(color =group),linewidth=0.8) + theme_light() + ggtitle("Cbb")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cbmplot <- call_list4[[2]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cbmplota <- ggplot(Cbmplot, aes(x = Time, y = concentration, group = group)) +
        geom_line(aes(color =group),linewidth=0.8)+ theme_light()  +ggtitle("Cbm")+ theme(legend.position = "bottom")+ labs(x = "Time")
      

      Cccsfplot <- call_list4[[3]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cccsfplota <- ggplot(Cccsfplot, aes(x = Time, y = concentration, group = group)) + 
        geom_line(aes(color =group),linewidth=0.8)  + theme_light() +ggtitle("Cccsf")+ theme(legend.position = "bottom")+ labs(x = "Time")
      

      Cscsfplot <- call_list4[[4]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cscsfplota <- ggplot(Cscsfplot, aes(x = Time, y = concentration, group = group)) + 
        geom_line(aes(color =group),linewidth=0.8)+ theme_light()  +ggtitle("Cscsf")  +  theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cartplot <- call_list4[[5]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cartplota <- ggplot(Cartplot, aes(x =Time, y = concentration, group = group)) + 
        geom_line(aes(color =group),linewidth=0.8) + theme_light() + ggtitle("Cart")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      grid.arrange(Cbbplota, Cbmplota,Cccsfplota,Cscsfplota, Cartplota,ncol=2)
      
    })
    
    
    
    output$sim4.senplot4b <- renderPlot({
      if (input$sim4.Sensub>0) { 
        
        withProgress(message = 'Creating plots...',
                     value = 0, {
                       for (i in 1:2) {
                         incProgress(1/2)
                         Sys.sleep(0.25)
                       }
                     })
        
        
        
        Indiconcenplot4() 
      }
    })  
    
   
   
    
    
    # This reactive function is to "download plots" for concentration plots for individually 
    # selected parameters for sensitivity analysis into main panel
    
    
    Indiconcenplot4down <- reactive({
      
      ##### new idea to make the code efficient begin
      call_list4 <- Indiparadata4()
      #### new idea to make the code efficient end
      
      
      Cbbplot <- call_list4[[1]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cbbplota <- ggplot(Cbbplot, aes(x = Time, y = concentration, group = group)) +
        geom_line(aes(color =group),linewidth=0.8) + theme_light() + ggtitle("Cbb")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cbmplot <- call_list4[[2]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cbmplota <- ggplot(Cbmplot, aes(x = Time, y = concentration, group = group)) +
        geom_line(aes(color =group),linewidth=0.8)+ theme_light()  +ggtitle("Cbm")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cccsfplot <- call_list4[[3]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cccsfplota <- ggplot(Cccsfplot, aes(x = Time, y = concentration, group = group)) + 
        geom_line(aes(color =group),linewidth=0.8)  + theme_light() +ggtitle("Cccsf")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cscsfplot <- call_list4[[4]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cscsfplota <- ggplot(Cscsfplot, aes(x = Time, y = concentration, group = group)) + 
        geom_line(aes(color =group),linewidth=0.8)+ theme_light()  +ggtitle("Cscsf")  +  theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cartplot <- call_list4[[5]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cartplota <- ggplot(Cartplot, aes(x =Time, y = concentration, group = group)) + 
        geom_line(aes(color =group),linewidth=0.8) + theme_light() + ggtitle("Cart")+ theme(legend.position = "bottom")+ labs(x = "Time")
    
      
      
      nelplot <- list(print(Cbbplota),print(Cbmplota),print(Cccsfplota),print(Cscsfplota),print(Cartplota))
      
      dev.off()
      
      grid.arrange(grobs = nelplot)
      
    })
    
    
    
    output$sim4.senconplot <- downloadHandler(
      filename <- function(){ "SenPlts4down.pdf"} ,
      content = function(file) {
        
        if (input$sim4.Sensub>0) {
          
          
          withProgress(message = 'Download in progress...',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          
          
        pdf(file,paper="default")
        plot(Indiconcenplot4down() )
        dev.off()
        }
        
      }
    )
    

    # This reactive function is to "output log plots" for concentration log plots for individually 
    # selected parameters for sensitivity analysis into main panel
    
    
    
    logconcenplot4 <- reactive({
      
      ##### new idea to make the code efficient begin
      call_list4 <- Indiparadata4()
      #### new idea to make the code efficient end
      
      
      Cbbplot <- call_list4[[1]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cbbplota <- ggplot(Cbbplot, aes(x = Time, y = log10(concentration), group = group)) +
        geom_line(aes(color =group),linewidth=0.8) + theme_light() + ggtitle("Cbb")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cbmplot <- call_list4[[2]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cbmplota <- ggplot(Cbmplot, aes(x = Time, y = log10(concentration), group = group)) +
        geom_line(aes(color =group),linewidth=0.8)+ theme_light()  +ggtitle("Cbm")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cccsfplot <- call_list4[[3]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cccsfplota <- ggplot(Cccsfplot, aes(x = Time, y = log10(concentration), group = group)) + 
        geom_line(aes(color =group),linewidth=0.8)  + theme_light() +ggtitle("Cccsf")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cscsfplot <- call_list4[[4]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cscsfplota <- ggplot(Cscsfplot, aes(x = Time, y = log10(concentration), group = group)) + 
        geom_line(aes(color =group),linewidth=0.8)+ theme_light()  +ggtitle("Cscsf")  +  theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cartplot <- call_list4[[5]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cartplota <- ggplot(Cartplot, aes(x =Time, y = log10(concentration), group = group)) + 
        geom_line(aes(color =group),linewidth=0.8) + theme_light() + ggtitle("Cart")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      grid.arrange(Cbbplota, Cbmplota,Cccsfplota,Cscsfplota, Cartplota,ncol=2)
      
    })
    
    
    
    output$sim4.senplot4blog <- renderPlot({
      if (input$sim4.Sensub>0) { 
        
        withProgress(message = 'Creating plots...',
                     value = 0, {
                       for (i in 1:2) {
                         incProgress(1/2)
                         Sys.sleep(0.25)
                       }
                     })
        
        logconcenplot4() 
      }
    })  
    
    
    
    
    # This reactive function is to "download log plots" for concentration log plots for individually 
    # selected parameters for sensitivity analysis into main panel
    
    
    logconcenplot4down <- reactive({
      
      ##### new idea to make the code efficient begin
      call_list4 <- Indiparadata4()
      #### new idea to make the code efficient end
      
      
      Cbbplot <- call_list4[[1]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cbbplota <- ggplot(Cbbplot, aes(x = Time, y = log10(concentration), group = group)) +
        geom_line(aes(color =group),linewidth=0.8) + theme_light() + ggtitle("Cbb")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cbmplot <- call_list4[[2]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cbmplota <- ggplot(Cbmplot, aes(x = Time, y = log10(concentration), group = group)) +
        geom_line(aes(color =group),linewidth=0.8)+ theme_light()  +ggtitle("Cbm")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cccsfplot <- call_list4[[3]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cccsfplota <- ggplot(Cccsfplot, aes(x = Time, y = log10(concentration), group = group)) + 
        geom_line(aes(color =group),linewidth=0.8)  + theme_light() +ggtitle("Cccsf")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cscsfplot <- call_list4[[4]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cscsfplota <- ggplot(Cscsfplot, aes(x = Time, y = log10(concentration), group = group)) + 
        geom_line(aes(color =group),linewidth=0.8)+ theme_light()  +ggtitle("Cscsf")  +  theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      Cartplot <- call_list4[[5]] %>% gather(key ="group", value = "concentration",  1:input$sim4.Pickparanum+1)
      Cartplota <- ggplot(Cartplot, aes(x =Time, y = log10(concentration), group = group)) + 
        geom_line(aes(color =group),linewidth=0.8) + theme_light() + ggtitle("Cart")+ theme(legend.position = "bottom")+ labs(x = "Time")
      
      
      #grid.arrange(Cbbplota, Cbmplota,Cccsfplota,Cscsfplota, Cartplota,ncol=2)
      
      
      nelplot <- list(print(Cbbplota),print(Cbmplota),print(Cccsfplota),print(Cscsfplota),print(Cartplota))
      
      dev.off()
      
      grid.arrange(grobs = nelplot)
      
    })
    
    
    
    output$sim4.senconlogplot <- downloadHandler(
      filename <- function(){ "logPlts4down.pdf"} ,
      content = function(file) {
        if (input$sim4.Sensub>0) {
          
          withProgress(message = 'Download in progress...',
                       value = 0, {
                         for (i in 1:2) {
                           incProgress(1/2)
                           Sys.sleep(0.25)
                         }
                       })
          
          
          
        pdf(file,paper="default")
        plot(logconcenplot4down() )
        dev.off()
        }
      }
    )
    

    
    
######## End  sensitivity analysis #############################################
################################################################################
################################################################################

########## Sample Input  files download
    
    # Define the download handler
    output$B4InputData <- downloadHandler(
      filename = function() {
        "Brain4_input.csv"  # Name of the file when it's downloaded
      },
      content = function(file) {
        # Copy the file from the www directory to the user's download location
        file.copy("www/SimcypData.csv", file)
      }
    )
    
    
    
    
}

########################### end server  ########################################
################################################################################
################################################################################
################################################################################



######## Combine ui and server to create a functioning web application #########

shinyApp(ui = ui, server = server)





################
# end of codes #
################