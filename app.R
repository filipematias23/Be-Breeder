######################
### Be-Breeder 2.0 ###
######################

######################################################################################
######################################################################################
###                                                                                ###
### Be-Breeder 2.0 is a R/Shiny application for statistical analysis and plant     ###
### breeding developed by Prof. Dr. Roberto Fritsche Neto, Dr. Filipe Inacio       ###
### Matias and PhD candidates Julia Silva Morosini and Fernando Garcia Espolador,  ###
### who integrate the Laboratory of Allogamous Plant Breeding group - Department   ###
### of Genetics, "Luiz de Queiroz"College of Agriculture, University of Sao Paulo, ###
### Brazil.                                                                        ###  
###                                                                                ###
### This app is free for online use, just the intellectual property is reserved in ###
### periodicals by scientific citations of Be-Breeder papers.                      ###
###                                                                                ###
###                                                                                ### 
### Contact:                                                                       ###
###                                                                                ###
### Avenida Padua Dias, 11 - Caixa Postal 83                                       ###
### CEP: 13418-900 - Piracicaba - Sao Paulo - Brasil                               ###
### Phone: 55 (19) 3429-4125 Line 45                                               ###
###                                                                                ###
### E-mail: roberto.neto@usp.br                                                    ###
###         filipematias23@usp.br                                                  ###
###         julia.morosini@usp.br                                                  ###
###                                                                                ###
### https://twitter.com/Alogamas_ESALQ                                             ###
### http://www.genetica.esalq.usp.br/alogamas/index2.html                          ###
###                                                                                ###
######################################################################################
######################################################################################

### Libraries ###

  library(shiny)
  library(shinydashboard)
  library(Matrix)
  library(lme4)
  library(lmerTest)
  library(agricolae)
  library(rrBLUP)
  library(adegenet)
  library(poppr)
  library(ape)
  library(matrixcalc)
  library(data.table)
  library(cluster)
  library(fpc)
  library(stringr)
  library(snpReady)
  library(corrgram)
  library(corrr)
  library(tableHTML)
  source("NCII.R")
  source("diallels.R")
  source("gardner.R")

##################
### Shiny - Ui ###
##################

# #####################################################################################################################################################################################
  
  header <- dashboardHeader(title = "Be-Breeder 2.0",
                            titleWidth = 250)
  
# #####################################################################################################################################################################################
  
  sidebar <- dashboardSidebar(
    tags$head(tags$style(HTML('.main-header .logo {
                              font-family: "Georgia", Times, "Times New Roman", serif;
                              font-weight: bold;
                              font-size: 25px;
                              }'))),
    width = 250,
    sidebarUserPanel("User Menu:"),
    sidebarMenu(
      id = "tabs",
      menuItem("Start", tabName = "start",icon = icon("cogs"), badgeLabel = "new", badgeColor = "green"),
      menuItem("Learning",icon = icon("graduation-cap"),
               menuSubItem("Inbreeding Effect", tabName = "Lear1", icon = icon("angle-double-right")),
               menuSubItem("Qualitative x Quantitative", tabName = "Lear3", icon = icon("angle-double-right")),
               menuSubItem("Progeny Size", tabName = "Lear4", icon = icon("angle-double-right")),
               menuSubItem("Selection Effect (HWE)", tabName = "Lear5", icon = icon("angle-double-right")),
               menuSubItem("Genetic Variance Components", tabName = "Lear6", icon = icon("angle-double-right")),
               menuSubItem("Synthetic Pop. Construction", tabName = "Lear7", icon = icon("angle-double-right")),
               menuItem("Recurrent Selection", icon = icon("angle-double-right"),
                        menuSubItem("Intrapopulation",tabName = "Lear8", icon = icon("angle-double-right")),
                        menuSubItem("Reciprocal",tabName = "Lear81", icon = icon("angle-double-right"))),
               menuItem("Hybrids (Jenkins)",icon = icon("angle-double-right"),
                        menuSubItem("Hybrid Prediction", tabName = "Lear91", icon = icon("angle-double-right")),
                        menuSubItem("Number of Hybrids", tabName = "Lear92", icon = icon("angle-double-right"))),
               menuSubItem("Genotype x Environment", tabName = "Lear10", icon = icon("angle-double-right")),
               menuSubItem("Heterosis", tabName = "Lear11", icon = icon("angle-double-right")),
               menuSubItem("Tester Effect", tabName = "Lear12", icon = icon("angle-double-right")),
               menuSubItem("Genetic Drift", tabName = "Lear13", icon = icon("angle-double-right")),
               menuSubItem("Indirect Selection", tabName = "Lear14", icon = icon("angle-double-right")),
               menuSubItem("Residual Effect on Selection", tabName = "Lear15", icon = icon("angle-double-right")),
               menuSubItem("Replicates vs Population size", tabName = "Lear16", icon = icon("angle-double-right")),
               menuSubItem("Economics in Plant Breeding", tabName = "Lear17", icon = icon("angle-double-right")),
               menuSubItem("Linkage Disequilibrium (LD)", tabName = "Lear18", icon = icon("angle-double-right")),
               menuSubItem("Cycles to Reduce LD", tabName = "Lear19", icon = icon("angle-double-right"))),
      menuItem("Phenotypic Breeding", icon = icon("leaf"), 
               menuItem("Experimental Analysis", icon = icon("angle-double-right"),
                        menuSubItem("", tabName = "subitem11", icon = icon("th")),
                        menuSubItem("Statistical Model", tabName = "subitem12")),
               menuItem("Diallel Analysis", icon = icon("venus-mars"),
                        menuItem("Griffing Design",icon = icon("angle-double-right"),
                                 menuSubItem("Dataset",tabName ="subitem1311"),
                                 menuSubItem("Analysis",tabName ="subitem1312")),
                        menuItem("Gardner and Eberhart Design",icon = icon("angle-double-right"),
                                 menuSubItem("Dataset",tabName ="subitem1321"),
                                 menuSubItem("Analysis",tabName ="subitem1322")),
                        menuItem("Factorial Design",icon = icon("angle-double-right"),
                                 menuSubItem("Dataset",tabName ="subitem1331"),
                                 menuSubItem("Analysis",tabName ="subitem1332"))),
               menuItem("Index Selection", icon = icon("angle-double-right"),
                        menuSubItem("Index File",tabName ="subitem141"),
                        menuSubItem("Index Analysis",tabName ="subitem142")),
               menuItem("Correlation Analysis",icon = icon("angle-double-right"),
                        menuSubItem("File Input",tabName ="subitem151"),
                        menuSubItem("Coefficients and graphs",tabName ="subitem152")),
               menuItem("Path Analysis",icon = icon("angle-double-right"),
                        menuSubItem("Trait File",tabName ="subitem161"),
                        menuSubItem("Path Analysis",tabName ="subitem162")),
               menuItem("Biplot Analysis",icon = icon("angle-double-right"),
                        menuSubItem("Dataset",tabName ="subitem171"),
                        menuSubItem("GE Biplot",tabName ="subitem172"),
                        menuSubItem("GE Cluster",tabName ="subitem174")),
               menuItem("Experiment Designs",icon = icon("angle-double-right"),
                        menuSubItem("Designs",tabName ="subitem18"))),
      menuItem("Molecular Breeding", icon = icon("sitemap"),
               menuItem("Genotyping Data",icon = icon("wrench"),
                        menuItem("Quality Control", icon = icon("angle-double-right"),
                        menuSubItem("Sample Set",tabName = "mol121", icon = icon("angle-double-right")),
                        menuSubItem("HapMap",tabName = "mol1212", icon = icon("angle-double-right")),
                        menuSubItem("Raw Data",tabName = "mol122", icon = icon("angle-double-right"))),
                        menuItem("Kinship Matrix", icon = icon("angle-double-right"),
                        menuSubItem("Z Matrix",tabName = "mol11", icon = icon("angle-double-right")),
                        menuSubItem("Kinship Matrix", tabName = "mol13", icon = icon("angle-double-right")))),
               menuItem("Genomic Selection (GS)",icon = icon("angle-double-right"),
                        menuItem("GS analysis",icon = icon("angle-double-right"),
                        menuSubItem("Phenotypic file", tabName = "mol21", icon = icon("angle-double-right")),
                        menuSubItem("Z Matrix", tabName = "mol22", icon = icon("angle-double-right")),
                        menuSubItem("rrBLUP", tabName = "mol23", icon = icon("angle-double-right"))),
                        menuItem("Prediction and Selection",icon = icon("angle-double-right"),
                                 menuSubItem("Marker effect", tabName = "mol241", icon = icon("angle-double-right")),
                                 menuSubItem("Z Matrix", tabName = "mol242", icon = icon("angle-double-right")),
                                 menuSubItem("Selection", tabName = "mol243", icon = icon("angle-double-right")))),
               menuItem("Genomic Association (GWAS)",icon = icon("angle-double-right"),
                                 menuSubItem("Phenotypic file", tabName = "mol31", icon = icon("angle-double-right")),
                                 menuSubItem("HapMap|t(Z)", tabName = "mol32", icon = icon("angle-double-right")),
                                 menuSubItem("Scores", tabName = "mol33", icon = icon("angle-double-right")),
                                 menuSubItem("Manhattan plot", tabName = "mol34", icon = icon("angle-double-right"))),
               menuItem("Diversity Analysis",icon = icon("angle-double-right"),
                        menuItem("Genetic Diversity",icon = icon("angle-double-right"),
                        menuSubItem("Dataset",tabName = "mol41",icon = icon("angle-double-right")),
                        menuSubItem("Diversity Summary",tabName = "mol42",icon = icon("angle-double-right")),
                        menuSubItem("Graphs",tabName = "mol43",icon = icon("angle-double-right"))),
                        menuItem("Discriminant Analysis",icon = icon("angle-double-right"),
                                 menuSubItem("Dataset",tabName = "mol61",icon = icon("angle-double-right")),
                                 menuSubItem("Graphs",tabName = "mol62",icon = icon("angle-double-right"))),
                        menuItem("Population Genetics",icon = icon("angle-double-right"),
                                 menuSubItem("Subpopulation Groups",tabName = "mol51",icon = icon("angle-double-right")),
                                 menuSubItem("Z Matrix",tabName = "mol52",icon = icon("angle-double-right")),
                                 menuSubItem("PopGen Analysis",tabName = "mol53",icon = icon("angle-double-right"))))),
      menuItem("Information", icon = icon("info-circle"),
               menuSubItem("About Be-Breeder 2.0", tabName = "info1", icon = icon("angle-double-right")),
               menuSubItem("Please cite us", tabName = "info2", icon = icon("angle-double-right")),
               menuSubItem("Team Contact", tabName = "info3", icon = icon("angle-double-right"))),
      menuItem(" ",icon = icon("angle-right"), badgeLabel = "Download Manual",badgeColor = "red",href = "http://vencovsky.esalq.usp.br:3838/shiny/manual/Manual_BeBreeder.pdf"),
      htmlTemplate("flagC.html",align = "center")
      ))
  
# #####################################################################################################################################################################################
  
  body <- dashboardBody(
    tabItems(
      
# #####################################################################################################################################################################################
      
      tabItem("start", box(status = "primary", height = 700, width = 600, align = "center",

                           fluidRow(
                             infoBox(hr(), "Home Page",icon=icon("home"), fill = T, color = "aqua",
                                     href = "http://www.genetica.esalq.usp.br/alogamas/index2.html"),
                             infoBox(hr(), "USP Page",icon=icon("users"), fill = T, color = "aqua",
                                     href = "http://www.en.esalq.usp.br"),
                             infoBox(hr(), "Twitter Page",icon=icon("twitter"), fill = T, color = "aqua",
                                     href = "https://twitter.com/Alogamas_ESALQ?ref_src=twsrc%5Etfw")),
                           br(),
                           tags$style(make_css(list('.box', 
                                                    c('text-align','font-size', 'font-family', 'color'), 
                                                    c('top','25px', 'arial', 'red')))),
                           strong("version 2.0: NEW FEATURES AVAILABLE!"),
                           br(),
                            br(),
                           a(id = "web_button", class = "btn action_button",
                             href = "http://vencovsky.esalq.usp.br:3838/shiny/manual/Manual_BeBreeder.pdf",
                             img(src="logobb19.jpg", height = 320, width = 950, align = "center"))
                           )),

# #####################################################################################################################################################################################

# Learning

# #####################################################################################################################################################################################

      tabItem("Lear1",
              fluidPage(
                titlePanel("Inbreeding Effect"),
                h4(em("In here, one can observe the fluctuation in additive and dominance genetic variances between and within populations as a function of the inbreeding coefficient 
                      (Wright's F statistics), and visualize how the genotype frequencies change according to the level of self-generation", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    selectInput("Lear1.1", "Self-Generation:",
                                choices = c("F1","F2","F3","F4","F5","F6","F7","F8","F9","F.inf"),width = 150),
                    sliderInput("Lear1.2","Wright's inbreeding coefficient (F)",0,1,0),
                    tableOutput("cont.Lear.1.0"),
                    helpText("ABOUT:"),
                    helpText('The results are given using the following equations:'),
                    img(src='Lear1.jpg',height = 200, width = 170)),
                  mainPanel(plotOutput("cont.Lear.1"))))),
      
# ######################################################################################################################################################################################
      
      tabItem("Lear3",
              fluidPage(
                titlePanel("Qualitative x Quantitative"),
                h4(em("In here, one can simulate the genetic structure of a trait, choosing the number of genes with dominance, partial dominance, and additive effects. The aim is to 
                      observe how the intra-allelic interactions influence the number of phenotypic classes and the frequency of each class in the population", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear3.1", "Number of dominant genes",0,30,15,step =1),
                    sliderInput("Lear3.2", "Number of partial dominant genes",0,30,10,step =1),
                    sliderInput("Lear3.3", "Number of additive genes",0,30,15,step =1)),
                  mainPanel(plotOutput("cont.Lear.3"))))),
      
# ######################################################################################################################################################################################
      
      tabItem("Lear4",
              fluidPage(
                titlePanel("Progeny Size"),
                h4(em("In here, one can estimate the number of individuals to be evaluated in the progeny in order to obtain a genotype carrying a trait of interest controled by n genes in a 
                      population that has undergone m inbreeding generations, with a certain probability (P)", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    numericInput("Lear4.1","Number of Genes (n)",min = 0,value = 0),
                    numericInput("Lear4.2","Self-Generation (m)",min = 0,value = 0),
                    sliderInput("Lear4.3", "Precision (P)",0,1,0.95),
                    helpText("ABOUT:"),
                    helpText("Progeny size is given by 'nº', calculated as:"),
                    img(src='Lear4.jpg',height = 150, width = 160)),
                  tableOutput("cont.Lear.4")))),
                  #verbatimTextOutput("cont.Lear.4")))),
      
# ######################################################################################################################################################################################
      
      tabItem("Lear5",
              fluidPage(
                titlePanel("Selection Effect (HWE)"),
                h4(em("In here, one can simulate different scenarios to understand how a locus with 2 alleles is affected by selection regarding the Hardy-Weinberg equilibrium (HWE) according to 
                      the number of genotypes (original and selected), their correspondent phenotypic value, the inbreeding rate, and the broad-sense heritability", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    numericInput("Lear5.2","Number of AA genotype",min = 0,value = 1),
                    numericInput("Lear5.3","Number of Aa genotype",min = 0,value = 0),
                    numericInput("Lear5.4","Number of aa genotype",min = 0,value = 0),
                    sliderInput("Lear5.41","Rate of Inbreeding",0,1,0),
                    numericInput("Lear5.5","Yield of AA",min = 0,value = 0),
                    numericInput("Lear5.6","Yield of Aa",min = 0,value = 0),
                    numericInput("Lear5.7","Yield of aa",min = 0,value = 0),
                    sliderInput("Lear5.11", "Heritability",0,1,1),
                    numericInput("Lear5.8","Number of AA selected",min = 0,value = 0),
                    numericInput("Lear5.9","Number of Aa selected",min = 0,value = 0),
                    numericInput("Lear5.10","Number of aa selected",min = 0,value = 0),
                    helpText("ABOUT:"),
                    helpText("Results are calculated using the following equations:"),
                    img(src='Lear5.jpg',height = 160, width = 230),
                    helpText("DS = selection differential"),
                    helpText("SG = selection gain (response to selection)")),
                  mainPanel(
                  #textOutput("con.Lear.5.01"),
                  htmlOutput("phr1"),
                  tableOutput("cont.Lear.5.0"),
                  htmlOutput("phr2"),
                  tableOutput("cont.Lear.5.1"),
                  htmlOutput("phr3"),
                  tableOutput("cont.Lear.5.2"),
                  htmlOutput("phr4"),
                  tableOutput("cont.Lear.5.3"),
                  htmlOutput("phr5"),
                  tableOutput("cont.Lear.5.4"),
                  htmlOutput("phr6"),
                  tableOutput("cont.Lear.5.5"))))),

                  #verbatimTextOutput("cont.Lear.5")))),
      
# ######################################################################################################################################################################################
      
      tabItem("Lear6",
              fluidPage(
                titlePanel("Genetic Variance Components"),
                h4(em("In here, one can set up the additive (a) and dominance (d) effects of a locus with two alleles (p and q) in order to observe how they impact on 
                      the magnitude and relationship among the genetic variance components", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear6.2","Frequency of allele (p)",0,1,0.6),
                    sliderInput("Lear6.1","Dominance Effect (d)",0,2.5,0.5),
                    sliderInput("Lear6.3","Additive Effect (a)",0,1,0.3),
                    #verbatimTextOutput("cont.Lear.6.1")),
                    helpText("ABOUT:"),
                    helpText("Results are calculated using the following equations:"),
                    img(src='Lear6.1.jpg',height = 200, width = 200),
                    helpText("add = average degree of dominance")),
                  mainPanel(plotOutput("cont.Lear.6"),
                            tableOutput("cont.Lear.6.1"))))),
      
# ######################################################################################################################################################################################
      
      tabItem("Lear7",
              fluidPage(
                titlePanel("Constructing Synthetic Populations"),
                h4(em("In here, the aim is to visualize how the number of alleles of a gene relates to genetic variability due to heterozigosity. The user provides the number of alleles
                      and their respective initial frequencies in the composition of the population.", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear7.1","Number of alleles",0,10,10,step = 1),
                    textAreaInput("Lear7.2","Frequency of each allele (Sum = 1.0)",
                                  value = "0.1,\n0.1,\n0.1,\n0.1,\n0.1,\n0.1,\n0.1,\n0.1,\n0.1,\n0.1",
                                  rows = 10, width = 80)),
                  mainPanel(
                  tableOutput("cont.Lear.7"),
                  helpText("ABOUT:"),
                  helpText("Results are calculated using the following equations:"),
                  img(src='Lear7.jpg',height = 150, width = 250))
                  ))),
                  #verbatimTextOutput("cont.Lear.7")))),
      
# ######################################################################################################################################################################################
      
      tabItem("Lear8",
              fluidPage(
                titlePanel("Intrapopulation Recurrent Selection"),
                h4(em("Considering the IRS breeding scheme, the user can simulate different selection scenarios defining the evaluation and recombination units and derived parameters to obtain 
                      the estimates of response to selection, effective size of the population (Ne), and inbreeding coefficient (Wright'F)", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    radioButtons('Lear8.1', '| Evaluation/Recombination | c | D1 | Ne |',
                                 c("| HS/HS | 1/4 | 0 |  4  |"='MI/MI',
                                   "| HS/S1 | 1/2 | 0 |  1  |"='MI/S1',
                                   "| FS/FS | 1/2 | 0 |  2  |"='IC/IC',
                                   "| FS/S1 | 1/2 | 0 |  1  |"='IC/S1',
                                   "| S1/S1 |  1  |1/2|  1  |"='S1/S1',
                                   "| S2/S2 | 3/2 |5/4| 2/3 |"='S2/S2'),
                                 'MI/MI'),
                    numericInput("Lear8.2","Number of Progenies Evaluated",min = 0,value = 0),
                    sliderInput("Lear8.3","Inbreeding Depression",0,1,0),
                    sliderInput("Lear8.4","Heritability",0,1,0),
                    sliderInput("Lear8.5","Selection Intensity",0,50,0,step = 1),
                    sliderInput("Lear8.6","Number of cycles",1,20,1,step = 1),
                    sliderInput("Lear8.7","Cycle Interval (year)",0.2,7,1,step = 0.1),
                    helpText("ABOUT:"),
                    helpText("Results are calculated using the following equation:"),
                    img(src='Lear8.jpg',width = 250,height = 110)),
                  mainPanel(tableOutput("cont.Lear.8a1"),
                            plotOutput("cont.Lear.8g"))))),

# ######################################################################################################################################################################################
      
      tabItem("Lear81",
              fluidPage(
                titlePanel("Reciprocal Recurrent Selection"),
                h4(em("Considering the RRS breeding scheme, the user can simulate different selection scenarios defining the evaluation and recombination units and derived parameters to obtain 
                      the estimates of response to selection, effective size of the population (Ne) for each group, and inbreeding coefficient (Wright'F) for each group", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    radioButtons('Lear81.1', '| Evaluation/Recombination | c | Ne |',
                                 c("| HS/HS | 1/8 | 4 |"='MI/MI',
                                   "| HS/S1 | 1/4 | 1 |"='MI/S1',
                                   "| TC/HS | 1/16| 4 |"='TC/MI',
                                   "| FS/S1 | 1/4 | 1 |"='IC/S1'),
                                 'MI/MI'),
                    numericInput("Lear81.21","Number of Progenies - Group 1",min = 0,value = 0),
                    numericInput("Lear81.22","Number of Progenies - Group 2",min = 0,value = 0),
                    sliderInput("Lear81.41","Heritability - Group 1",0,1,0),
                    sliderInput("Lear81.42","Heritability - Group 2",0,1,0),
                    sliderInput("Lear81.51","Selection Intensity (i1)",1,50,1,step = 1),
                    sliderInput("Lear81.52","Selection Intensity (i2)",1,50,1,step = 1),
                    sliderInput("Lear81.6","Number of cycles",1,20,1,step = 1),
                    sliderInput("Lear81.7","Cycle Interval (year)",0.2,7,1,step = 0.1),
                    helpText("ABOUT:"),
                    helpText("Results are calculated using the following equation:"),
                    img(src='Lear8.1.jpg',width = 250,height = 110)),
                  mainPanel(tableOutput("cont.Lear.8.1"),
                            plotOutput("cont.Lear.81g"))))),

# ######################################################################################################################################################################################
      
      tabItem("Lear91",
              fluidPage(
                titlePanel("Hybrid Prediction"),
                h4(em("The user can obtain the predicted value for three-way hybrids (TH) and double-cross hybrids (DH) through the input of a .txt file
                      containing the mean phenotypic values for the single-cross hybrids (SH)", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    fileInput('file.L9', 'Choose File:',
                              accept=c('text/csv', 
                                       'text/comma-separated-values,text/plain', 
                                       '.csv',".txt")),
                    tags$hr(),
                    checkboxInput('header.L9', 'Header', FALSE),
                    radioButtons('sep.L9', 'Separator',
                                 c(Comma=',',
                                   Semicolon=';',
                                   Tab='\t'),
                                 ','),
                    radioButtons('quote.L9', 'Quote',
                                 c(None='',
                                   'Double Quote'='"',
                                   'Single Quote'="'"),
                                 ''),
                    selectInput("Lear91.out", "Look Options:",
                                choices = c("Table",
                                            "Three-way Cross Hybrids Prediction",
                                            "Double Cross Hybrids Prediction")),
                    checkboxInput("ex.L.9","Example", FALSE),
                    helpText("ABOUT:"),
                    helpText("Results are calculated using the following equation:"),
                    img(src='Lear9.1.jpg',width = 300,height = 150)),
                  tableOutput("cont.Lear.9.1")))),

# ######################################################################################################################################################################################

      tabItem("Lear92",
              fluidPage(
                titlePanel("Number of Hybrids"),
                h4(em("The user can obtain the potential number of single-cross hybrids (SH), three-way hybrids (TH), and double-cross hybrids (DH), given the number of lines in a 
                      breeding population (either for one or two heterotic groups)", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    selectInput("Lear9.2.1", "Population Structure",
                                choices = c("One Heterotic Group",
                                            "Two Heterotic Groups")),
                    # uiOutput("res3.lear92"),
                    numericInput("Lear9.2.2","Number of Lines", min = 0,value = 0),
                    numericInput("Lear9.2.3","Number of Lines in Group 1", min = 0,value = 0),
                    numericInput("Lear9.2.4","Number of Lines in Group 2", min = 0,value = 0),
                    helpText("ABOUT:"),
                    helpText("Results are calculated using the following equations:"),
                    img(src='Lear9.2.jpg',height = 300, width = 220)),
                  tableOutput("cont.Lear.9.2")))),
      
# ######################################################################################################################################################################################
      
      tabItem("Lear10",
              fluidPage(
                titlePanel("Genotype x Environment"),
                h4(em("In here, one can simulate different degrees of interaction considering three genotypes and two environments in order to verify how the GxE interaction influences the variances
                      and the broad-sense heritability", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear10.7","Coefficient of Variation (CV)",5,20,5),
                    sliderInput("Lear10.8","Number of Repetitions (r)",1,4,1,step = 1),
                    sliderInput("Lear10.1","Yield of Genotype 1 in Environment 1",0,5,0),
                    sliderInput("Lear10.2","Yield of Genotype 1 in Environment 2",0,5,0),
                    sliderInput("Lear10.3","Yield of Genotype 2 in Environment 1",0,5,0),
                    sliderInput("Lear10.4","Yield of Genotype 2 in Environment 2",0,5,0),
                    sliderInput("Lear10.5","Yield of Genotype 3 in Environment 1",0,5,0),
                    sliderInput("Lear10.6","Yield of Genotype 3 in Environment 2",0,5,0),
                    helpText("ABOUT:"),
                    helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"))),
                  mainPanel(htmlOutput("phr10.1"),
                            tableOutput("cont.Lear.10a.1"),
                            htmlOutput("phr10.2"),
                            tableOutput("cont.Lear.10a.2"),
                            htmlOutput("phr10.3"),
                            tableOutput("cont.Lear.10a.3"),
                            plotOutput("cont.Lear.10.2"))))),
      
# ######################################################################################################################################################################################
      
      tabItem("Lear11",
              fluidPage(
                titlePanel("Heterosis"),
                h4(em("The aim in here is to understand how the heterosis fluctuates according to the number of genes, the dominance deviation, and the genetic divergence between
                      the parents (lines) for the simulated hybrids", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear11.1","Number of Genes",0,100,0, step = 10),
                    sliderInput("Lear11.2","Dominance Deviation (d)",0,1,0, step = 0.1),
                    sliderInput("Lear11.3","Genetic Divergence",0,1,0, step = 0.1),
                    helpText("ABOUT:"),
                    helpText("Results are calculated using the following equations:"),
                    img(src='Lear11.jpg',height = 200, width = 150)),
                  tableOutput("cont.Lear.11")))),

# #####################################################################################################################################################################################      
      
      tabItem("Lear12",
              fluidPage(
                titlePanel("Tester Effect"),
                h4(em("In here, the user can determine the effect of testers in a experimental condition in order to verify their impacts on the genetic variability", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear12.1","Frequency of allele (p)",0,1,0),
                    sliderInput("Lear12.2","r",0,1,0),
                    sliderInput("Lear12.3","Additive Effect (a)",0,1,0),
                    sliderInput("Lear12.4","Dominance Effect (d)",0,2,0, step = 0.1),
                    sliderInput("Lear12.5","Inbreeding coefficient (F)",0,1,0),
                    verbatimTextOutput("cont.Lear.12.1"),
                    helpText("ABOUT:"),
                    helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"))),
                  mainPanel(plotOutput("cont.Lear.12"))))),
                  # mainPanel(tableOutput("cont.Lear.12.1"),
                  #           plotOutput("cont.Lear.12"))))),

# #####################################################################################################################################################################################      

      tabItem("Lear13",
              fluidPage(
                titlePanel("Genetic Drift"),
                h4(em("In here, one can understand how the occurence of genetic drift is influenced by the frequency of the allele 'a' (q) and the selection intensity adopted. From that, 
                      is is possible to estimate how many generations are needed to eliminate this allele", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear13.1","Frequency of allele 'a' (q)",0,1,0),
                    sliderInput("Lear13.2","Selection Intensity",0,1,0.1),
                    verbatimTextOutput("cont.Lear.13.1"),
                    helpText("ABOUT:"),
                    helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"))),
                  mainPanel(plotOutput("cont.Lear.13"))))),

# #####################################################################################################################################################################################      

      tabItem("Lear14",
              fluidPage(
                titlePanel("Indirect Selection"),
                h4(em("In here, the user can verify the possibility of selecting a trait based on another, which is called Indirect selection. Providing the values for the parameters
                      requested, it is possible to calculate the response to the indirect selection", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear14.1","Heritability of trait 'x'",0,1,0),
                    sliderInput("Lear14.2","Correlation between 'x' and 'y'",0,1,0),
                    numericInput("Lear14.3","Number of individuals evaluated",min = 0,value = 10),
                    numericInput("Lear14.4","Number of individuals selected",min = 0,value = 0),
                    helpText("ABOUT:"),
                    helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"))),
                  mainPanel(tableOutput("cont.Lear.14.1"),
                            plotOutput("cont.Lear.14"))))),

# #####################################################################################################################################################################################

      tabItem("Lear15",
              fluidPage(
                titlePanel("Residual Effect on Selection"),
                h4(em("The aim in this tab is to analyze if the selection process in a breeding population is being reliable and capturing the best breeding values. Among other variables,
                      it is displayed the genetic variability (GV) for each scenario", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    numericInput("Lear15.1","Number of individuals evaluated",min = 0,value = 10),
                    sliderInput("Lear15.2","Selection Intensity (%)",0,20,0),
                    sliderInput("Lear15.3","Heritability",0,1,0),
                    helpText("ABOUT:"),
                    helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"))),
                    tableOutput("cont.Lear.15")))),

# #####################################################################################################################################################################################

      tabItem("Lear16",
              fluidPage(
                titlePanel("Replicates vs Population size"),
                h4(em("In here, the user here may infer what the impact of the number of replicates on the population size is, and especially, what is the increment of that 
                      on the response to selection (RS), simulating different scenarios in function of the parameters above", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear16.1","Number of replicates",1,5,1),
                    numericInput("Lear16.2","Number of plots available",min = 0,value = 0),
                    numericInput("Lear16.3","Number of genotypes selected",min = 0,value = 0),
                    sliderInput("Lear16.4","Trait Heritability",0,1,0),
                    verbatimTextOutput("cont.Lear.16.1"),
                    helpText("ABOUT:"),
                    helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"))),
                  mainPanel(plotOutput("cont.Lear.16"))))),

# #####################################################################################################################################################################################

      tabItem("Lear17",
              fluidPage(
                titlePanel("Economics in Plant Breeding"),
                h4(em("In here, the aim is to provide a broader framework of a breeding process, considering the cost spent per plot, the viability, and the expected genetic gain. By that,
                      the user can observe the influence of the nature of the trait and the experimental size on the total costs", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    numericInput("Lear17.1","Number of genotypes evaluated",min = 0,value = 10),
                    numericInput("Lear17.2","Number of genotypes selected",min = 0,value = 0),
                    sliderInput("Lear17.3","Trait Heritability",0,1,0),
                    sliderInput("Lear17.4","Number of replicates",1,5,1),
                    sliderInput("Lear17.5","Time (year/cycle)",0.2,7,0.2),
                    numericInput("Lear17.6","Cost to evaluate one plot (U$)",min = 0,value = 10),
                    verbatimTextOutput("cont.Lear.17.1"),
                    helpText("ABOUT:"),
                    helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"))),
                  mainPanel(plotOutput("cont.Lear.17"))))),

# ######################################################################################################################################################################################

      tabItem("Lear18",
              fluidPage(
                titlePanel("Linkage Disequilibrium (LD)"),
                h4(em("In here, one can simulate different conditions informing the frequencies of the haplotypes AB, Ab, aB, and ab, and the rates of the alleles A, a, B,
                      and b. The application will then estimate whether the loci are independent", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear18.1","Frequency of haplotype AB",0,1,0),
                    sliderInput("Lear18.2","Frequency of haplotype Ab",0,1,0),
                    sliderInput("Lear18.3","Frequency of haplotype aB",0,1,0),
                uiOutput('ui.lear18.1'),
                    sliderInput("Lear18.5","Frequency of A",0,1,0),
                uiOutput('ui.lear18.2'),
                    sliderInput("Lear18.7","Frequency of B",0,1,0),
                uiOutput('ui.lear18.3')),
                  mainPanel(tableOutput("cont.Lear.18"),
                            helpText("ABOUT:"),
                            helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank")))))),

# ######################################################################################################################################################################################

      tabItem("Lear19",
              fluidPage(
                titlePanel("Cycles to Reduce LD"),
                h4(em("In here, combining the loss of linkage disequilibrium (LD) aimed to the recombination ratio between loci in a population, it is possible to compute how many cycles are
                      needed to dissipate the LD", style = "color:gray")),
                br(),
                sidebarLayout(
                  sidebarPanel(
                    sliderInput("Lear19.1","Loss of LD (%)",0,100,5),
                    sliderInput("Lear19.2","Recombination fraction between loci",0,0.5,0.05),
                    verbatimTextOutput("cont.Lear.19.1"),
                    helpText("ABOUT:"),
                    helpText("Equations for this tab are detailed on GitHub:", a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"))),
                  mainPanel(plotOutput("cont.Lear.19"))))),

# ######################################################################################################################################################################################

# Phenotypic Breeding
      
# ######################################################################################################################################################################################
      
      tabItem("subitem11",
              fluidPage(
                titlePanel("Uploading Files"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput('file1', 'Choose File:',
                              accept=c('text/csv', 
                                       'text/comma-separated-values,text/plain', 
                                       '.csv',".txt")),
                    tags$hr(),
                    checkboxInput('header', 'Header', FALSE),
                    radioButtons('sep', 'Separator',
                                 c(Comma=',',
                                   Semicolon=';',
                                   Tab='\t'),
                                 ','),
                    radioButtons('quote', 'Quote',
                                 c(None='',
                                   'Double Quote'='"',
                                   'Single Quote'="'"),
                                 ''),
                    selectInput("table.1", "Look Options:",
                                choices = c("Table","Data Structure","Names"),width = 150),
                    checkboxInput("info.11",strong(code("HELP")), FALSE),
                    checkboxInput("ex.11","Example", FALSE)),
                  verbatimTextOutput("contents11")))),

# #####################################################################################################################################################################################

      tabItem("subitem12", 
              fluidPage(
                titlePanel("Statistical Model"),
                sidebarLayout(
                  sidebarPanel( 
                    checkboxInput('res.12.name', 'Sources of Variation', FALSE),
                    textInput("model","Type the Statistical Model", width = 700,
                              value = "y~Block+Site+(1|Genotype)+(1|Genotype:Site)"),
                    checkboxInput('ml', 'ML(Default = REML)', FALSE),
                    checkboxInput('res.12.fix', 'Genotype as Fixed (Default = Random)', FALSE),
                    actionButton("res.12.run","run",icon = icon("hand-o-up")),
                    selectInput("res.12.1", "Choose Results:",
                                choices = c("Mean",
                                            "Summary",
                                            "ANOVA",
                                            "Analysis of Deviance",
                                            "BLUE","Adjusted Means",
                                            "BLUP","Predicted Means"),width = 200),
                    sliderInput("res.12.2", "Selection Intensity",0,100,100),
                    textInput("name.txt.12","Type the File Name"),
                    downloadButton('downloadData12', 'Download'),
                    checkboxInput("info.12",strong(code("HELP")), FALSE)),
                  verbatimTextOutput("contents12")))),

# #####################################################################################################################################################################################

tabItem("subitem1311",
        fluidPage(
          titlePanel("Uploading Files"),
          sidebarLayout(
            sidebarPanel(
              fileInput('file1311', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('header1311', 'Header', FALSE),
              radioButtons('sep1311', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quote1311', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.1311",strong(code("HELP")), FALSE),
              checkboxInput("ex.1311","Example", FALSE)),
            verbatimTextOutput("contents1311")))),

# #####################################################################################################################################################################################

tabItem("subitem1321",
        fluidPage(
          titlePanel("Uploading Files"),
          sidebarLayout(
            sidebarPanel(
              fileInput('file1321', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('header1321', 'Header', FALSE),
              radioButtons('sep1321', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quote1321', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.1321",strong(code("HELP")), FALSE),
              checkboxInput("ex.1321","Example", FALSE)),
            verbatimTextOutput("contents1321")))),

# #####################################################################################################################################################################################

tabItem("subitem1331",
        fluidPage(
          titlePanel("Uploading Files"),
          sidebarLayout(
            sidebarPanel(
              fileInput('file1331', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('header1331', 'Header', FALSE),
              radioButtons('sep1331', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quote1331', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.1331",strong(code("HELP")), FALSE),
              checkboxInput("ex.1331","Example", FALSE)),
            verbatimTextOutput("contents1331")))),

# #####################################################################################################################################################################################

tabItem("subitem1312",               
        fluidPage(
          titlePanel("Griffing Design"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.1312", "Choose Results:",
                          choices = c("Output")),
              textInput("name.txt.1312","Type the Name File"),
              downloadButton('downloadData1312', 'Download'),
              checkboxInput("info.1312",strong(code("HELP")), FALSE)),
            verbatimTextOutput("contents1312")))),

# #####################################################################################################################################################################################

tabItem("subitem1322",               
        fluidPage(
          titlePanel("Gardner and Eberhart Design"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.1322", "Choose Results:",
                          choices = c("Output")),
              textInput("name.txt.1322","Type the Name File"),
              downloadButton('downloadData1322', 'Download'),
              checkboxInput("info.1322",strong(code("HELP")), FALSE)),
            verbatimTextOutput("contents1322")))),

# #####################################################################################################################################################################################

tabItem("subitem1332",               
        fluidPage(
          titlePanel("Factorial Design"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.1332", "Choose Results:",
                          choices = c("Output")),
              textInput("name.txt.1332","Type the Name File"),
              downloadButton('downloadData1332', 'Download'),
              checkboxInput("info.1332",strong(code("HELP")), FALSE)),
            verbatimTextOutput("contents1332")))),

# #####################################################################################################################################################################################

      tabItem("subitem141", 
              fluidPage(
                titlePanel("Index Selection"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput('file14', 'Choose File',
                              accept=c('text/csv', 
                                       'text/comma-separated-values,text/plain', 
                                       '.csv',".txt")),
                    tags$hr(),
                    checkboxInput('header.14', 'Header', FALSE),
                    radioButtons('sep.14', 'Separator',
                                 c(Comma=',',
                                   Semicolon=';',
                                   Tab='\t'),
                                 ','),
                    radioButtons('quote.14', 'Quote',
                                 c(None='',
                                   'Double Quote'='"',
                                   'Single Quote'="'"),
                                 ''),
                    checkboxInput("info.141",strong(code("HELP")), FALSE),
                    checkboxInput("ex.14","Example", FALSE)),
                  verbatimTextOutput("contents14.1")))),
      
# #####################################################################################################################################################################################

      tabItem("subitem142", 
              fluidPage(
                titlePanel("Index Analysis"),
                sidebarLayout(
                  sidebarPanel(
                    checkboxInput('res.14.name', 'Trait Names', FALSE),
                    textInput("vec.14", "Type a vector with weights (comma delimited)", "0.3,-0.5,0.1,-0.1"),
                    sliderInput("res.14", "Selection Intensity",0,100,100),
                    textInput("name.txt.14","Type the File Name"),
                    downloadButton('downloadData14', 'Download'),
                    checkboxInput("info.142",strong(code("HELP")), FALSE)),
                  verbatimTextOutput("contents14.2")))),

# #####################################################################################################################################################################################

tabItem("subitem151",
        fluidPage(
          titlePanel("Uploading Files"),
          sidebarLayout(
            sidebarPanel(
              fileInput('file15', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('header15', 'Header', FALSE),
              radioButtons('sep15', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quote15', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.151",strong(code("HELP")), FALSE),
              checkboxInput("ex.15.0","Example", FALSE)),
            verbatimTextOutput("contents151")))),

# #####################################################################################################################################################################################      

tabItem("subitem152", 
        fluidPage(
          titlePanel("Correlation Analysis"),
          sidebarLayout(
            sidebarPanel(
             selectInput("res.15.2", "Choose Results:",
                          choices = c("Genotype Correlation",
                                      "pValue")),
             numericInput("res.15.3","Minimum Value of Correlation",min = c(-1),max = 1,value = 0.1,step = 0.01),
              textInput("name.txt.15","Type the File Name"),
              downloadButton('downloadData15', 'Download'),
             checkboxInput("info.152",strong(code("HELP")), FALSE)),
            mainPanel(verbatimTextOutput("contents15.2"),
                      plotOutput("contents15.1"),
                      plotOutput("contents15.3"))))),

# #####################################################################################################################################################################################
      
      tabItem("subitem161",
              fluidPage(
                titlePanel("Uploading Files"),
                sidebarLayout(
                  sidebarPanel(
                    fileInput('file16', 'Choose File.txt',
                              accept=c('text/csv', 
                                       'text/comma-separated-values,text/plain', 
                                       '.csv',".txt")),
                    tags$hr(),
                    checkboxInput('header16', 'Header', FALSE),
                    radioButtons('sep16', 'Separator',
                                 c(Comma=',',
                                   Semicolon=';',
                                   Tab='\t'),
                                 ','),
                    radioButtons('quote16', 'Quote',
                                 c(None='',
                                   'Double Quote'='"',
                                   'Single Quote'="'"),
                                 ''),
                    checkboxInput("info.161",strong(code("HELP")), FALSE),
                    checkboxInput("ex.16","Example", FALSE)),
                  verbatimTextOutput("contents161")))),
      
# #####################################################################################################################################################################################      
      
      tabItem("subitem162", 
              fluidPage(
                titlePanel("Path Analysis"),
                sidebarLayout(
                  sidebarPanel(
                    checkboxInput('res.16.name', 'Trait Names', FALSE),
                    textInput("name.16", "Main Trait", "trait1"),
                    selectInput("res.16.2", "Choose Results:",
                                choices = c("Path Analysis",
                                            "Traits Correlation")),
                    textInput("name.txt.16","Type the File Name"),
                    downloadButton('downloadData16', 'Download'),
                    checkboxInput("info.162",strong(code("HELP")), FALSE)),
                  verbatimTextOutput("contents16.2")))),

# #####################################################################################################################################################################################      

tabItem("subitem171",
        fluidPage(
          titlePanel("Uploading Files"),
          sidebarLayout(
            sidebarPanel(
              fileInput('file17', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('header17', 'Header', FALSE),
              radioButtons('sep17', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quote17', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.171",strong(code("HELP")), FALSE),
              checkboxInput("ex.17","Example", FALSE)),
            verbatimTextOutput("contents171")))),

# #####################################################################################################################################################################################      

tabItem("subitem172", 
        fluidPage(
          titlePanel("GE Biplot Analysis"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.172", "Choose Results:",
                          choices = c("Summary","Scores of G","Scores of E"),width = 200),
              textInput("name.txt.172","Type the File Name"),
              downloadButton('downloadData172', 'Download'),
              checkboxInput("info.172",strong(code("HELP")), FALSE)),
            mainPanel(verbatimTextOutput("contents172"),
                      plotOutput("contents173"))))),


# #####################################################################################################################################################################################      

tabItem("subitem174", 
        fluidPage(
          titlePanel("GE Cluster Analysis"),
          sidebarLayout(
            sidebarPanel(
              numericInput("num.174","Determine number of clusters",min = 1,value = 3),
              selectInput("summary.174", "Choose Results:",
                          choices = c("Cluster information","Genotype groups","Cluster mean"),width = 200),
              checkboxInput("sel.175","Cluster Selection", FALSE),
              numericInput("num.175","Determine number of clusters",min = 1,value = 3),
              textInput("name.txt.174","Type the File Name"),
              downloadButton('downloadData174', 'Download')),
            mainPanel(verbatimTextOutput("contents174"),
                      plotOutput("contents175"))
          ))),

# #####################################################################################################################################################################################      

tabItem("subitem18", 
        fluidPage(
          titlePanel("Experiment Designs"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.18", "Choose Designs:",
                          choices = c("Completely randomized design",
                                      "Randomized complete block design",
                                      "Latin square design",
                                      "Balanced Incomplete Block Designs",
                                      "Lattice designs",
                                      "Alpha designs",
                                      "Augmented block designs",
                                      "Split plot designs",
                                      "Factorial"),width = 400),
              uiOutput('ui.18.1'),
              uiOutput('ui.18.3'),
              uiOutput('ui.18.2'),
              numericInput("num.18.4","Seed",value = 100),
              textInput("name.txt.18","Type the File Name"),
              downloadButton('downloadData18', 'Download'),
              checkboxInput("info.18",strong(code("HELP")), FALSE)),
            verbatimTextOutput("contents18")))),

# #####################################################################################################################################################################################      
      
# Genotypic Breeding

# #####################################################################################################################################################################################      

tabItem("mol121",
        fluidPage(
          titlePanel("Uploading File"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol121', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol121', 'Header', FALSE),
              radioButtons('sepmol121', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol121', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              radioButtons('na.stringsmol121', 'na.strings',
                           c(None='',
                             'NA'='NA',
                             '-'="-"),
                           ''),
              checkboxInput("info.mol121",strong(code("HELP")), FALSE),
              checkboxInput("exmol121","Example", FALSE)),
            verbatimTextOutput("contentsmol121")))),

# #####################################################################################################################################################################################      

tabItem("mol122", 
        fluidPage(
          titlePanel("Raw Data"),
          sidebarLayout(
            sidebarPanel(
              numericInput("res1.mol122","MAF:",value = 0.05,min = 0,max = 1),
              numericInput("res2.mol122","Call Rate:",value = 0.95,min = 0,max = 1),
              numericInput("res3.mol122","Sweep Sample:",value = 0.5,min = 0,max = 1),
              checkboxInput("res4.mol122","Input Data",TRUE),
              selectInput("res4.5.mol122", "Choose Input Method:",
                          choices = c( "wright","mean")),
              selectInput("res5.mol122", "Choose Data Frame:",
                          choices = c("long","wide")),
              selectInput("res6.mol122", "Choose Outfile:",
                          choices = c("012","-101","structure")),
              actionButton("res7.mol122","run",icon = icon("hand-o-up")),
              selectInput("res8.mol122", "Choose Output:",
                          choices = c("M.clean", "Hapmap", "Report")),
              downloadButton('downloadDatamol122', 'Download'),
              checkboxInput("info.mol122",strong(code("HELP")), FALSE)),
            verbatimTextOutput("contentsmol122")))),

# #####################################################################################################################################################################################      

tabItem("mol11",
        fluidPage(
          titlePanel("Uploading File"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol11', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol11', 'Header', FALSE),
              radioButtons('sepmol11', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol11', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol11",strong(code("HELP")), FALSE),
              checkboxInput("exmol11","Example", FALSE)),
            verbatimTextOutput("contentsmol11")))),

# #####################################################################################################################################################################################      

tabItem("mol1212",
        fluidPage(
          titlePanel("HapMap"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol1212', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol1212', 'Header', FALSE),
              radioButtons('sepmol1212', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol1212', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol1212",strong(code("HELP")), FALSE),
              checkboxInput("exmol1212","Example", FALSE)),
            verbatimTextOutput("contentsmol1212")))),

# #####################################################################################################################################################################################      

tabItem("mol13", 
        fluidPage(
          titlePanel("Kinship Matrix"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res1.mol13", "Choose method:",
                          choices = c("VanRaden", "UAR", "UARadj","GK")),
              selectInput("res2.mol13", "Choose frame:",
                          choices = c("long","wide")),
              uiOutput("res3.mol13"),
              actionButton("res5.mol13","run",icon = icon("hand-o-up")),
              textInput("name.txt.mol13","Type the File Name"),
              downloadButton('downloadDatamol13', 'Download'),
              checkboxInput("info.mol13",strong(code("HELP")), FALSE)),
            verbatimTextOutput("contentsmol13")))),

# #####################################################################################################################################################################################      

tabItem("mol21",
        fluidPage(
          titlePanel("Phenotypic file"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol21', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol21', 'Header', FALSE),
              radioButtons('sepmol21', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol21', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol21",strong(code("HELP")), FALSE),
              checkboxInput("exmol21","Example", FALSE)),
            verbatimTextOutput("contentsmol21")))),

# #####################################################################################################################################################################################      

tabItem("mol22",
        fluidPage(
          titlePanel("Z Matrix"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol22', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol22', 'Header', FALSE),
              radioButtons('sepmol22', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol22', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol22",strong(code("HELP")), FALSE),
              checkboxInput("exmol22","Example", FALSE)),
            verbatimTextOutput("contentsmol22")))),

# #####################################################################################################################################################################################      

tabItem("mol23", 
        fluidPage(
          titlePanel("Genomic selection (GS)"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.mol23", "Choose Results:",
                          choices = c("Marker effects",
                                      "Breeding values",
                                      "Accuracy")),
              textInput("name.txt.mol23","Type the File Name"),
              downloadButton('downloadDatamol23', 'Download')),
            verbatimTextOutput("contentsmol23")))),

# #####################################################################################################################################################################################      

tabItem("mol241",
        fluidPage(
          titlePanel("Marker effect"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol241', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol241', 'Header', FALSE),
              radioButtons('sepmol241', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol241', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol241",strong(code("HELP")), FALSE),
              checkboxInput("exmol241","Example", FALSE)),
            verbatimTextOutput("contentsmol241")))),

# #####################################################################################################################################################################################      

tabItem("mol242",
        fluidPage(
          titlePanel("Z Matrix"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol242', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol242', 'Header', FALSE),
              radioButtons('sepmol242', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol242', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol242",strong(code("HELP")), FALSE),
              checkboxInput("exmol242","Example", FALSE)),
            verbatimTextOutput("contentsmol242")))),

# #####################################################################################################################################################################################      

tabItem("mol243", 
        fluidPage(
          titlePanel("Genomic selection (GS)"),
          sidebarLayout(
            sidebarPanel(
              sliderInput("res.mol243", "Selection Intensity",0,100,100),
              textInput("name.txt.mol243","Type the File Name"),
              downloadButton('downloadDatamol243', 'Download')),
            verbatimTextOutput("contentsmol243")))),

# #####################################################################################################################################################################################      

tabItem("mol31",
        fluidPage(
          titlePanel("Phenotypic file"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol31', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol31', 'Header', FALSE),
              radioButtons('sepmol31', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol31', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol31",strong(code("HELP")), FALSE),
              checkboxInput("exmol31","Example", FALSE)),
            verbatimTextOutput("contentsmol31")))),

# #####################################################################################################################################################################################      

tabItem("mol32",
        fluidPage(
          titlePanel("HapMap | t(Z)"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol32', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol32', 'Header', FALSE),
              radioButtons('sepmol32', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol32', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol32",strong(code("HELP")), FALSE),
              checkboxInput("exmol32","Example", FALSE)),
            verbatimTextOutput("contentsmol32")))),

# #####################################################################################################################################################################################      

tabItem("mol33", 
        fluidPage(
          titlePanel("GWAS - Scores"),
          sidebarLayout(
            sidebarPanel(
              textInput("name.txt.mol33","Type the File Name"),
              downloadButton('downloadDatamol33', 'Download')),
            verbatimTextOutput("contentsmol33.1")))),

# #####################################################################################################################################################################################      

tabItem("mol34", 
        fluidPage(
          titlePanel("GWAS - Manhattan plot"),
          sidebarLayout(
            sidebarPanel(
              numericInput("LOD","LOD threshold",value = 4,min = 0)),
              mainPanel(plotOutput("contentsmol33.2"))))),

# #####################################################################################################################################################################################      

tabItem("mol41",
        fluidPage(
          titlePanel("Dataset"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.mol41", "Choose type of file:",
                          choices = c("Structure",
                                      "Genpop")),
              fileInput('filemol41', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput("info.mol41",strong(code("HELP")), FALSE),
              checkboxInput("exmol41","Example", FALSE)),
            verbatimTextOutput("contentsmol41")))),

# #####################################################################################################################################################################################      

tabItem("mol42",
        fluidPage(
          titlePanel("Diversity Summary"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.mol42", "Choose Results:",
                          choices = c("Summary",
                                      "Expected heterozygosity")),
                          textInput("name.txt.mol42","Type the File Name"),
                          downloadButton('downloadDatamol42', 'Download')), 
            verbatimTextOutput("contentsmol42")))),

# #####################################################################################################################################################################################      

tabItem("mol43",
        fluidPage(
          titlePanel("Graphs"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.mol43", "Choose Results:",
                          choices = c("Number of alleles per loci",
                                      "Expected Heterozygosity per loci",
                                      "Observed Heterozygosity per loci",
                                      "Heterozygosity: expected-observed",
                                      "F index:  F = (He - Ho)/He",
                                      "Expected Heterozygosity per population",
                                      "Population sample size",
                                      "Principal component analysis by population 1",
                                      "Principal component analysis by population 2",
                                      "Principal component analysis by individual",
                                      "Dendrogram by individual",
                                      "Dendrogram by population",
                                      "Population Structure",
                                      "Neighbor joining dendrogram")),
              uiOutput("ui.mol43"),
              uiOutput("ui.mol43.2"),
              uiOutput("ui.mol43.3"),
              uiOutput("ui.mol43.4")), 
            mainPanel(plotOutput("contentsmol43"))))),

# #####################################################################################################################################################################################      

tabItem("mol51",
        fluidPage(
          titlePanel("Subpopulation Groups"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol51', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol51', 'Header', FALSE),
              radioButtons('sepmol51', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol51', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol51",strong(code("HELP")), FALSE),
              checkboxInput("exmol51","Example", FALSE)),
            verbatimTextOutput("contentsmol51")))),

# #####################################################################################################################################################################################      

tabItem("mol52",
        fluidPage(
          titlePanel("Z Matrix"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol52', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol52', 'Header', FALSE),
              radioButtons('sepmol52', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol52', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol52",strong(code("HELP")), FALSE),
              checkboxInput("exmol52","Example", FALSE)),
            verbatimTextOutput("contentsmol52")))),

# #####################################################################################################################################################################################      

tabItem("mol53",
        fluidPage(
          titlePanel("PopGen Analysis"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.mol53", "Choose Results:",
                          choices = c("Whole",
                                      "By Group")),
              textInput("name.txt.mol53","Type the File Name"),
              downloadButton('downloadDatamol53', 'Download'),
              checkboxInput("info.mol53",strong(code("HELP")), FALSE)),
            verbatimTextOutput("contentsmol53")))),

# ######################################################################################################################################################################################

tabItem("mol61",
        fluidPage(
          titlePanel("Dataset"),
          sidebarLayout(
            sidebarPanel(
              fileInput('filemol61', 'Choose File.txt',
                        accept=c('text/csv', 
                                 'text/comma-separated-values,text/plain', 
                                 '.csv',".txt")),
              tags$hr(),
              checkboxInput('headermol61', 'Header', FALSE),
              radioButtons('sepmol61', 'Separator',
                           c(Comma=',',
                             Semicolon=';',
                             Tab='\t'),
                           ','),
              radioButtons('quotemol61', 'Quote',
                           c(None='',
                             'Double Quote'='"',
                             'Single Quote'="'"),
                           ''),
              checkboxInput("info.mol61",strong(code("HELP")), FALSE),
              checkboxInput("exmol61","Example", FALSE)),
            verbatimTextOutput("contentsmol61")))),

# #####################################################################################################################################################################################

tabItem("mol62",
        fluidPage(
          titlePanel("Graphs"),
          sidebarLayout(
            sidebarPanel(
              selectInput("res.mol62", "Choose Results:",
                          choices = c("Population Size",
                                      "Find Clusters",
                                      "Scatterplot",
                                      "Compoplot")),
              uiOutput("res.mol63")),
            mainPanel(plotOutput("contentsmol62"))))),

# #####################################################################################################################################################################################      
      
# Information

# #####################################################################################################################################################################################      

      tabItem("info1",         
              fluidPage(
                titlePanel("About Be-Breeder 2.0"),
                br(),
                # fluidRow(
"Be-Breeder is a R/Shiny application for phenotypic and molecular analyses in the plant breeding context, developed by Prof. Dr. Roberto Fritsche Neto, 
Dr. Filipe Inacio Matias, and Ph.D. candidate Julia Silva Morosini, who integrate the Laboratory of Allogamous Plant Breeding group - Department of Genetics,
'Luiz de Queiroz' College of Agriculture, University of Sao Paulo, Brazil.",
h1(),
"In its newest version, Be-Breeder 2.0, the platform is designed to consistently address a great variety of phenotypic, genetic, and statistical subjects,
supported by numeric and graphic outputs. Be-Breeder 2.0 is freely accessible and the source code is hosted at",
a(href="https://github.com/filipematias23/Be-Breeder/", "https://github.com/filipematias23/Be-Breeder/", target="_blank"), "."
                # pre(includeText("BeBreederAbout.txt")))
)),

      tabItem("info2",
              fluidPage(
                titlePanel("Citation"),
                br(),
"1. Matias, Filipe Inacio;  Granato, Italo Stefanine Correa; Fritsche-Neto, Roberto. (2018). Be-Breeder: an R/Shiny application for phenotypic data analyses in plant breeding.
Crop Breeding and Applied Biotechnology, 18:241-243.",
a(href="http://dx.doi.org/10.1590/1984-70332018v18n2s36", "http://dx.doi.org/10.1590/1984-70332018v18n2s36", target="_blank"),
h1(),
"2. Matias, Filipe Inácio, Granato, Italo Stefanine Correa, Dequigiovanni, Gabriel, & Fritsche-Neto, Roberto. (2017). Be-Breeder - an application for analysis of genomic data in plant breeding. 
Crop Breeding and Applied Biotechnology, 17(1):54-58",
a(href="https://dx.doi.org/10.1590/1984-70332017v17n1n8", "https://dx.doi.org/10.1590/1984-70332017v17n1n8", target="_blank"),
h1(),
"3. Fritsche-Neto, Roberto, & Matias, Filipe Inácio. (2016). Be-Breeder - Learning: a new tool for teaching and learning plant breeding principles. 
Crop Breeding and Applied Biotechnology, 16(3):240-245",
a(href="https://dx.doi.org/10.1590/1984-70332016v16n3n36", "https://dx.doi.org/10.1590/1984-70332016v16n3n36", target="_blank"),
h1()
    )),
  
      tabItem("info3",
              fluidPage(
                titlePanel("Contact"),
                h1(),
              #img(src="Allogamous.jpg", height = 100, width = 180, align = "left"),
              a(id = "web_button", class = "btn action_button",
                href= "http://www.genetica.esalq.usp.br/alogamas/index2.html",
                img(src = "Logolab.png", height = 130),
                target="_blank"),
      
              #infoBox(hr(),"Allogamous Plant Breeding Laboratory",icon=icon("close"), fill = F, color = "olive", width = 5,
                      #href = "http://www.genetica.esalq.usp.br/alogamas/index2.html"),
              br(),
              #strong(a(href="http://www.genetica.esalq.usp.br/alogamas/index2.html", "Allogamous Plant Breeding Laboratory", width=20)),
              
              h1(),
                strong("Roberto Fritsche Neto"), "- roberto.neto@usp.br",
              a(id = "web_button", class = "btn action_button",
                href= "https://orcid.org/0000-0003-4310-0047",
                img(src = "orcid.png", height = 28),
                target="_blank"), 
              a(id = "web_button", class = "btn action_button",
                href= "http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4715585D3",
                img(src = "lattes.png", height = 22),
                target="_blank"),
                
              br(),
                strong("Filipe Inácio Matias"), "- filipematias23@usp.br",
              a(id = "web_button", class = "btn action_button",
                href= "https://orcid.org/0000-0002-4414-2866",
                img(src = "orcid.png", height = 28),
                target="_blank"), 
              a(id = "web_button", class = "btn action_button",
                href= "http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4439660Y8",
                img(src = "lattes.png", height = 22),
                target="_blank"),    
              
              br(),
                strong("Júlia Silva Morosini"), "- julia.morosini@usp.br",
              a(id = "web_button", class = "btn action_button",
                href= "https://orcid.org/0000-0002-4106-6941",
                img(src = "orcid.png", height = 28),
                target="_blank"), 
              a(id = "web_button", class = "btn action_button",
                href= "http://buscatextual.cnpq.br/buscatextual/visualizacv.do?id=K4405222U1",
                img(src = "lattes.png", height = 22),
                target="_blank")
              ))
               # fluidRow(
                #  br(),
              #pre(includeText("BeBreederContact.txt")))))
# #####################################################################################################################################################################################

    )
  )
  
# #####################################################################################################################################################################################
  
  shinyApp(
    ui = dashboardPage(header, sidebar, body),
    
######################
### Shiny - Server ###
######################

    server = function(input, output) { 
      options(shiny.maxRequestSize=50*1024^2)
      
# #####################################################################################################################################################################################

# Learning      
      
# #####################################################################################################################################################################################

      # Lean.1
      
      output$cont.Lear.1.0 <- renderTable(rownames = T, bordered = T, colnames = T, digits = 4,{
        f.1<-input$Lear1.2
        a1<-2*f.1 # Va. Among
        a2<-(1-f.1) # Va. Within
        a3<-(f.1)*(1-f.1) # Vd. Among
        a4<-(1-f.1) # Vd. Within
        M<-t(cbind(round(a1,3),round(a2,3),round(a3,3),round(a4,3)))
        sigmaSq = paste('\u03c3\u00B2', rep(c("a", "d"), each=2))
        nat = rep(c("among", "wihtin:"),2)
        rownames(M)<- paste(sigmaSq, nat, sep=" ")
        colnames(M)<-c("Variance:")
        return(M)})
      
      output$cont.Lear.1 <- renderPlot({
        X<-matrix(c(0,100,0,
                    25,50,25,
                    37.5,25,37.5,
                    43.75,12.5,43.75,
                    46.875,6.25,46.875,
                    48.4375,3.125,48.4375,
                    49.21875,1.5625,49.21875,
                    49.60938,0.78125,49.60938,
                    49.80469,0.390625,49.80469,
                    50,0,50),3,10)
        out.Lear.1<-reactive({
          switch(input$Lear1.1,
                 "F1"=X[,1],
                 "F2"=X[,2],
                 "F3"=X[,3],
                 "F4"=X[,4],
                 "F5"=X[,5],
                 "F6"=X[,6],
                 "F7"=X[,7],
                 "F8"=X[,8],
                 "F9"=X[,9],
                 "F.inf"=X[,10])})
        barplot(as.vector(out.Lear.1()), names.arg=c("aa","Aa","AA"), 
                ylab=c("Genotype Frequency (%)"), 
                xlab="Genotype",  
                cex.names=1.2, 
                ylim = c(0,100),
                col=c("brown1","royalblue2","green3"))
        add_legend <- function(...) {
          opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                      mar=c(0, 0, 0, 0), new=TRUE)
          on.exit(par(opar))
          plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
          legend(...)}
        add_legend("top",
                   legend=out.Lear.1(),
                   pch=19,
                   col=c("brown1","royalblue2","green3"),
                   border = "black",
                   box.col="white",
                   bty="n",
                   bg="n",
                   text.width=c(0.3,0.3,0.3),
                   horiz=TRUE)})
      
# #####################################################################################################################################################################################

      # Lean.3
      
      classes<-function(ngd, ngp, nga){
        vgd<-matrix(rep(sample(c(3, 3, 3, 1), 1000, replace = T), ngd), ngd, 1000)
        vgp<-matrix(rep(sample(c(3, 2.5, 2.5, 1), 1000, replace = T), ngp), ngp, 1000)
        vga<-matrix(rep(sample(c(3, 2, 2, 1), 1000, replace = T), nga), nga, 1000)
        vgt<-apply(rbind(vgd, vgp, vga), 2, sum)
        nclasses<-length(unique(vgt))
        MAIN<-paste("Distribution of ",nclasses," Classes",sep = "")
        hist(vgt, xlab = "Genotypic value", main = MAIN , col = "palegreen3")
        return(paste("number of classes:", nclasses))}
      output$cont.Lear.3 <- renderPlot({ 
        classes(input$Lear3.1,input$Lear3.2,input$Lear3.3)})
      
# #####################################################################################################################################################################################
      
      # Lean.4
      
      output$cont.Lear.4 <- renderTable(rownames = T, colnames = F, digits = 0,{
        IH<-((2^((input$Lear4.2)-1)-1)/2^(input$Lear4.2))^(input$Lear4.1)
        NP<-cbind(round((log(1-input$Lear4.3)/log(1-IH)),0))
        rownames(NP)<-c("Number of plants that should be evaluated:")
        return(NP)})
      
# #####################################################################################################################################################################################
      
      # Lean.5
      
      cont.Lear.5<-reactive({
        AA<-input$Lear5.2
        Aa<-input$Lear5.3
        aa<-input$Lear5.4
        N<-sum(AA,Aa,aa)
        Total<-N
        if(N==0){return(NULL)}
        VG.AA<-input$Lear5.5
        VG.Aa<-input$Lear5.6
        VG.aa<-input$Lear5.7
        Mean.0<-round(sum((AA*VG.AA),(Aa*VG.Aa),(aa*VG.aa))/N,3)
        p<-round((AA+Aa/2)/N,3)
        q<-round(1-p,3)
        AA.Esp<-round((p^2+p*q*input$Lear5.41)*N,0)
        Aa.Esp<-round((2*p*q-2*p*q*input$Lear5.41)*N,0)
        aa.Esp<-round((q^2+p*q*input$Lear5.41)*N,0)
        a <- (VG.AA - VG.aa)/2
        d <- VG.Aa - (VG.AA + VG.aa)/2
        alpha <- a + (p - q)*d 
        Va = 2*p*q*alpha^2
        Vd = 2*p*q*d^2
        A.effect <- -q*alpha 
        a.effect <- p*alpha
        if(input$Lear5.41==0){f<-(2*p*q - Aa/N)/(2*p*q)}
        if(input$Lear5.41!=0){f<-input$Lear5.41}
        qui.cal<-round(sum((AA-AA.Esp)^2/AA.Esp,(Aa-Aa.Esp)^2/Aa.Esp,(aa-aa.Esp)^2/aa.Esp),3)
        qui.tab<-6.63
        Sel.AA<-input$Lear5.8
        Sel.Aa<-input$Lear5.9
        Sel.aa<-input$Lear5.10
        Ns<-sum(Sel.AA,Sel.Aa,Sel.aa)
        h2<-input$Lear5.11
        ps<-round((2*Sel.AA+Sel.Aa)/(2*Ns),3)
        qs<-round((2*Sel.aa+Sel.Aa)/(2*Ns),3)
        AA.Sel<-round((ps^2+ps*qs*input$Lear5.41)*Ns,0)
        Aa.Sel<-round((2*ps*qs-2*ps*qs*input$Lear5.41)*Ns,0)
        aa.Sel<-round((qs^2+ps*qs*input$Lear5.41)*Ns,0)
        new.alpha <- a + (ps - qs)*d 
        new.Va <- 2*ps*qs*new.alpha^2 
        new.Vd <- 2*ps*qs*d^2 
        new.A.effect <- -qs*new.alpha 
        new.a.effect <- ps*new.alpha
        inic.1<-c(alpha,Va,Vd,A.effect,a.effect)
        inic.2<-c(new.alpha,new.Va,new.Vd,new.A.effect,new.a.effect)
        if(input$Lear5.41==0){fs<-(2*ps*qs - Aa.Sel/Ns)/(2*ps*qs)}
        if(input$Lear5.41!=0){fs<-input$Lear5.41}
        Mean.Sel<-round(sum((AA.Sel*VG.AA),(Aa.Sel*VG.Aa),(aa.Sel*VG.aa))/sum(AA.Sel,Aa.Sel,aa.Sel),3)
        DS<-round(Mean.Sel-Mean.0,3)
        GS<-round(DS*h2,3)
        LIS.0 = as.data.frame(cbind(p,q,a,d,round((AA.Esp/N),3),round((Aa.Esp/N),3),round((aa.Esp/N),3),round(f,3)))
        names(LIS.0)<-c(" p"," q"," a"," d"," D"," H"," R"," F")
        LIS.1<-as.data.frame(cbind(rbind(AA,Aa,aa,Total),rbind(AA.Esp,Aa.Esp,aa.Esp,Total)))
        names(LIS.1)<-c("Observed","Expected")
        LIS.2<-cbind(qui.cal,qui.tab)
        colnames(LIS.2)<-c("X-squared Cal","X-squared Tab")
        LIS.3<-as.data.frame(cbind(Mean.0,Mean.Sel,DS,GS))
        colnames(LIS.3)<-c("Original Mean","Selected Mean","Selection Differential","Response to Selection")
        LIS.4<-as.data.frame(cbind(ps,qs,round((AA.Sel/Ns),3),round((Aa.Sel/Ns),3),round((aa.Sel/Ns),3),round(fs,3)))
        colnames(LIS.4)<-c("New p"," New q"," New D"," New H"," New R"," F")
        LIS.4.5<-as.data.frame(cbind(round(inic.1,3),round(inic.2,3)))
        colnames(LIS.4.5)<-c("Original", "Selected")
        sigmaSq = '\u03c3\u00B2'
        alpha = '\u03b1'
        leg = paste(sigmaSq,c("a", "d"), sep = "")
        rownames(LIS.4.5)<-c(alpha,leg, "A effect","a effect")
        LIST.Lear5<-list(LIS.0,LIS.1,LIS.2,LIS.4,LIS.4.5,LIS.3)
        return(LIST.Lear5)
      })
        
      output$cont.Lear.5.0 <- renderTable(rownames = T, colnames = T,{
        LIS.0<-as.data.frame(cont.Lear.5()[[1]], row.names = "")
        return(LIS.0)})
        
      output$cont.Lear.5.1 <- renderTable(rownames = T, colnames = T,{
        LIS.1<-as.data.frame(cont.Lear.5()[[2]], row.names = c("AA", "Aa", "aa", "Total"))
        return(LIS.1)})
      
      output$cont.Lear.5.2 <- renderTable(rownames = T, colnames = T,{
        LIS.2<-as.data.frame(cont.Lear.5()[[3]], row.names = "")
        return(LIS.2)})
        
      output$cont.Lear.5.3 <- renderTable(rownames = T, colnames = T,{
        LIS.3<-as.data.frame(cont.Lear.5()[[6]], row.names = "")
        return(LIS.3)})
        
      output$cont.Lear.5.4 <- renderTable(rownames = T, colnames = T,{
        LIS.4<-as.data.frame(cont.Lear.5()[[4]], row.names = "")
        return(LIS.4)})
      
      output$cont.Lear.5.5 <- renderTable(rownames = T, colnames = T,{
        LIS.4.5<-as.data.frame(cont.Lear.5()[[5]])
        return(LIS.4.5)})
      
      # name of each table
      output$phr1 <- renderUI({
        test1 <- strong("Frequency")
        HTML(paste(test1, sep = '<br/>'))})
      
      output$phr2 <- renderUI({
        test1 <- strong("Table Information")
        HTML(paste(test1, sep = '<br/>'))})
      
      output$phr3 <- renderUI({
        test1 <- strong("X-square Test")
        HTML(paste(test1, sep = '<br/>'))})
      
      output$phr4 <- renderUI({
        test1 <- strong("Selection Information")
        HTML(paste(test1, sep = '<br/>'))})
      
      output$phr5 <- renderUI({
        test1 <- strong("Selection Frequency")
        HTML(paste(test1, sep = '<br/>'))})
      
      output$phr6 <- renderUI({
        test1 <- strong("Allelic Effect Substitution")
        HTML(paste(test1, sep = '<br/>'))})
      
      
# #####################################################################################################################################################################################
      
      # Lean.6
      
      var1<-function(p, a, d){
        p<-as.numeric(p)  
        a<-as.numeric(a)
        d<-as.numeric(d)
        Var.g<-function(p, a, d){
          Vg<-2*p*(1-p)*(a+((1-p)-p)*d)^2+(2*p*(1-p)*d)^2
          return(Vg)}
        vg<-Var.g(p, a, d)  
        Var.a<-function(p, a, d){
          Va<-2*p*(1-p)*(a+((1-p)-p)*d)^2
          return(Va)}
        va<-Var.a(p, a, d)
        Var.d<-function(p, d){
          Vd<-(2*p*(1-p)*d)^2
          return(Vd)}
        vd<-Var.d(p, d)
        gmd<-d/a
        alpha <- a + (p - (1-p))*d
        alpha <- as.numeric(alpha)
        A.effect <- -(1-p)*alpha
        a.effect <- p*alpha
        curve(Var.g(x, a, d), from = 0, to = 1, col="black", xlab = "p", ylab = "Variances")
        curve(Var.a(x, a, d), from = 0, to = 1, add = T, col="blue")
        curve(Var.d(x, d), from = 0, to = 1, add = T, col="red")
        legend("topright", c("Vg","Va","Vd"), lty=c(1,1), lwd=c(2.5,2.5, 2.5),col=c("black", "blue","red"))}
      var2<-function(p, a, d){
        p<-as.numeric(p)  
        a<-as.numeric(a)
        d<-as.numeric(d)
        Var.g<-function(p, a, d){
          Vg<-2*p*(1-p)*(a+((1-p)-p)*d)^2+(2*p*(1-p)*d)^2
          return(Vg)}
        vg<-Var.g(p, a, d)  
        Var.a<-function(p, a, d){
          Va<-2*p*(1-p)*(a+((1-p)-p)*d)^2
          return(Va)}
        va<-Var.a(p, a, d)
        Var.d<-function(p, d){
          Vd<-(2*p*(1-p)*d)^2
          return(Vd)}
        vd<-Var.d(p, d)
        gmd<-d/a
        alpha <- a + (p - (1-p))*d
        alpha <- as.numeric(alpha)
        A.effect <- -(1-p)*alpha
        a.effect <- p*alpha
        output<-as.data.frame(rbind("Vg"=vg, "Va"=va, "Vd"=vd, "add"=gmd, "A.effect"=A.effect, "a.effect"=a.effect))
        colnames(output)<-c("Results")
        sigmaSq = '\u03c3\u00B2'
        legi = paste(sigmaSq,c("g", "a", "d"), sep = "")
        rownames(output) = c(legi, "Average degree of dominance", "A effect", "a effect")
        return(output)}
      
      output$cont.Lear.6.1 <- renderTable(rownames = T, colnames = T, digits = 4,{
        d<-input$Lear6.1
        p<-input$Lear6.2
        a<-input$Lear6.3
        #sigmaSq = '\u03c3\u00B2'
        #legi = paste(sigmaSq,c("g", "a", "d"), sep = "")
        #rownames(output) = c(legi, "average degree of dominance", "A effect", "a effect")
        return(var2(p, a, d))})
      
      output$cont.Lear.6 <- renderPlot({
        d<-input$Lear6.1
        p<-input$Lear6.2
        a<-input$Lear6.3
        var1(p, a, d)})
      
# #####################################################################################################################################################################################
      
      # Lean.7
      
      output$cont.Lear.7 <- renderTable(rownames = F, colnames = T, digits = 4,{
        prog<-input$Lear7.1
        allel.0<-as.numeric(unlist(strsplit(input$Lear7.2,",")))
        if(length(allel.0)>prog|length(allel.0)<prog){
          return(paste("The number of frequencies should be equal to the number of alleles =",prog,sep = " "))}
        if(sum(allel.0)!=1){
          return(paste("The sum of frequencies should be equal to 1.0"))}
        allel<-allel.0
        div<-(sum(allel^2))
        RH<-1-div
        Var.G<-div-div^2
        PB<-as.data.frame(cbind(RH,Var.G))
        names(PB)<-c("Residual Heterosigosy (RH)","Variance of Heterosigosy")
        return(PB)})
      
# #####################################################################################################################################################################################
      
      # Lean.8
      
      cont.Lear.8a <- reactive({
        if(input$Lear8.1=='MI/MI'){D1=0;c=(1/4);Ne=4;}
        if(input$Lear8.1=='MI/S1'){D1=0;c=(1/2);Ne=1;}
        if(input$Lear8.1=='IC/IC'){D1=0;c=(1/2);Ne=2;}
        if(input$Lear8.1=='IC/S1'){D1=0;c=(1/2);Ne=1;}
        if(input$Lear8.1=='S1/S1'){D1=(0.5);c=1;Ne=1;}
        if(input$Lear8.1=='S2/S2'){D1=(5/4);c=(3/2);Ne=(2/3);}
        I<-function(persel){
          estimate<-dnorm(qnorm(1-(persel/100), 0, 1), 0, 1)/(persel/100)
          return(estimate)}
        i<-I(input$Lear8.5)
        SP<-(input$Lear8.5)/100
        total<-input$Lear8.2
        NE<-Ne*SP*total
        f<-1/(2*NE)
        DE<-input$Lear8.3
        Va<-input$Lear8.4
        GS<-(i*c*(Va+f*D1))-(DE*6*i)/(2*NE)
        GS.T<-GS/input$Lear8.7
        LISTA.8<-list(GS,GS.T,NE,f)
        names(LISTA.8)<-c("Response to Selection per cycle","Response to Selection per year","Effective Population Size (Ne)","F Wright")
        return(LISTA.8)})
      
     # output$cont.Lear.8 <- renderPrint({
     #    return(cont.Lear.8a())})
     
     output$cont.Lear.8a1 <- renderTable(rownames = T, colnames = F, digits = 4,{
       GS<-as.matrix(cont.Lear.8a())
       return(GS)})
     
     output$cont.Lear.8g <- renderPlot({
       source("missVG.R")
       if(cont.Lear.8a()[[3]]==0){return(NULL)}
       missVg(cont.Lear.8a()[[3]], input$Lear8.6)})
      
# #####################################################################################################################################################################################
      
     # Lean.8.1
      
      cont.Lear.8.1a <- reactive({
        if(input$Lear81.1=='MI/MI'){c=(1/8);Ne=4;}
        if(input$Lear81.1=='MI/S1'){c=(1/4);Ne=1;}
        if(input$Lear81.1=='TC/MI'){c=(1/16);Ne=4;}
        if(input$Lear81.1=='IC/S1'){c=(1/4);Ne=1;}
        I<-function(persel){
          estimate<-dnorm(qnorm(1-(persel/100), 0, 1), 0, 1)/(persel/100)
          return(estimate)}
        i1<-I(input$Lear81.51)
        SP1<-(input$Lear81.51)/100
        i2<-I(input$Lear81.52)
        SP2<-(input$Lear81.52)/100
        total.1<-input$Lear81.21
        total.2<-input$Lear81.22
        NE.1<-Ne*SP1*total.1
        NE.2<-Ne*SP2*total.2
        NE.T<-as.data.frame(cbind(NE.1,NE.2))
        names(NE.T)<-c("Group 1","Group 2")
        f.1<-1/(2*NE.1)
        f.2<-1/(2*NE.2)
        f.T<-as.data.frame(cbind(f.1,f.2))
        names(f.T)<-c("Group 1","Group 2")
        Va.1<-input$Lear81.41
        Va.2<-input$Lear81.42
        GS<-(i1*c*Va.1)+(i2*c*Va.2)
        GS.T<-GS/input$Lear81.7
        LISTA.81<-list(GS,GS.T,NE.T,f.T)
        names(LISTA.81)<-c("Response to Selection per cycle","Response to Selection per Year","Effective Population Size (Ne)","Wright's inbreeding coefficient (F)")
        return(LISTA.81)})
     
     # output$cont.Lear.8.1 <- renderPrint({
     #   return(cont.Lear.8.1a())})

     output$cont.Lear.8.1 <- renderTable(rownames = T, colnames = F, digits = 4,{
       list8 = t(as.data.frame(cont.Lear.8.1a()))
       rownames(list8) = c("Response to Selection per cycle",
                           "Response to Selection per Year",
                           "Effective Population Size (Ne) - Group 1",
                           "Effective Population Size (Ne) - Group 2",
                           "F Wright - Group 1",
                           "F Wright - Group 2")
       return(list8)})
     
     output$cont.Lear.81g <- renderPlot({
       source("missVGrec.R")
       if(cont.Lear.8.1a()[[3]]==0){return(NULL)}
       if(input$Lear81.21==0|input$Lear81.22==0){return(NULL)}
       missVg(cont.Lear.8.1a()[[3]][,1],cont.Lear.8.1a()[[3]][,2],input$Lear81.6)})
      
# #####################################################################################################################################################################################
      
     # Lean.9.1
      
      out.Lear.91.1 <- reactive({
        if(input$ex.L.9!=FALSE){
          TABE.1<-matrix(c(8,15,14,11,15,4,18,12,14,18,6,10,11,12,10,3),4,4)
          colnames(TABE.1)<-c("A","B","C","D")
          rownames(TABE.1)<-c("A","B","C","D")
          return(TABE.1)}
        inFile <- input$file.L9
        if (is.null(inFile))
          return(NULL)
        TABE<-read.table(inFile$datapath, header=input$header.L9, sep=input$sep.L9, 
                         quote=input$quote.L9)
        rownames(TABE)<-TABE[,1]
        TABE.2<-TABE[,-1]
        return(TABE.2)})
      
     out.Lear.91.2 <- reactive({
        X<-out.Lear.91.1()
        if (is.null(X))
          return(NULL)
        Hib.T<-NULL
        for(i in 1:dim(X)[1]){
          for(u in 1:dim(X)[1]){
            for(t in 1:dim(X)[1]){
              if(colnames(X)[i]!=colnames(X)[u]&
                 colnames(X)[i]!=colnames(X)[t]&
                 colnames(X)[u]!=colnames(X)[t]){
                HT<-paste("TH(",colnames(X)[i],colnames(X)[u],")",colnames(X)[t],sep="")
                x.1<-X[row.names(X)==colnames(X)[i],colnames(X)==colnames(X)[t]]
                x.2<-X[row.names(X)==colnames(X)[u],colnames(X)==colnames(X)[t]]
                HT.Predict<-(1/2)*(x.1+x.2)
                Hib.T<-rbind(Hib.T,cbind(HT,HT.Predict))
                colnames(Hib.T) = c("TH", 'Value')
              }}}}
        return(Hib.T)})  
      
     out.Lear.91.3 <- reactive({
        X<-out.Lear.91.1()
        if (is.null(X))
          return(NULL) 
        Hib.D<-NULL
        for(i in 1:dim(X)[1]){
          for(u in 1:dim(X)[1]){
            for(t in 1:dim(X)[1]){
              for(w in 1:dim(X)[1]){
                if(colnames(X)[i]!=colnames(X)[u]&
                   colnames(X)[i]!=colnames(X)[t]&
                   colnames(X)[i]!=colnames(X)[w]&
                   colnames(X)[u]!=colnames(X)[t]&
                   colnames(X)[u]!=colnames(X)[w]&
                   colnames(X)[t]!=colnames(X)[w]){
                  HD<-paste("DH(",colnames(X)[i],colnames(X)[u],")(",colnames(X)[t],colnames(X)[w],")",sep="")
                  x.1<-X[row.names(X)==colnames(X)[i],colnames(X)==colnames(X)[t]]
                  x.2<-X[row.names(X)==colnames(X)[i],colnames(X)==colnames(X)[w]]
                  x.3<-X[row.names(X)==colnames(X)[u],colnames(X)==colnames(X)[t]]
                  x.4<-X[row.names(X)==colnames(X)[u],colnames(X)==colnames(X)[w]]
                  HD.Predict<-(1/4)*(x.1+x.2+x.3+x.4)
                  Hib.D<-rbind(Hib.D,cbind(HD,HD.Predict))
                  colnames(Hib.D) = c("DH", 'Value')
                }}}}}
        return(Hib.D)})  
     
     output$cont.Lear.9.1 <- renderTable(colnames = T, rownames = T,{
        switch(input$Lear91.out,
               "Table"=out.Lear.91.1(),
               "Three-way Cross Hybrids Prediction"=out.Lear.91.2(),
               "Double Cross Hybrids Prediction"=out.Lear.91.3())})
      
# #####################################################################################################################################################################################
  
      # Lean.9.2
     
     out.Lear.92.One<-reactive({
        n<-input$Lear9.2.2
        HS<-(n*(n-1))/2
        HT<-(n*(n-1)*(n-2))/2
        HD<-(n*(n-1)*(n-2)*(n-3))/8
        HIB<-as.data.frame(rbind(n,HS,HT,HD))
        rownames(HIB)<-c("Number of Lines","SH","TH","DH")
        return(HIB)})
      out.Lear.92.Two<-reactive({
        a<-input$Lear9.2.3
        b<-input$Lear9.2.4
        HS.1<-a*b
        HT.1<-a*(a-1)*b+b*(b-1)*a
        HD.1<-a*(a-1)*b*(b-1)
        HIB<-as.data.frame(rbind(a,b,HS.1,HT.1,HD.1))
        rownames(HIB)<-c("Lines Group 1","Lines Group 2","SH","TH","DH")
        return(HIB)})
      
      output$cont.Lear.9.2 <- renderTable(rownames = T, colnames = F, digits = 0,{
        switch(input$Lear9.2.1,
               "One Heterotic Group" = out.Lear.92.One(),
               "Two Heterotic Groups" = out.Lear.92.Two())
               })
      
# #####################################################################################################################################################################################
      
      # Lean.10
      
      gxa1 <- function(ga11, ga12, ga21, ga22, ga31, ga32, CV, R){
        fen<-t(matrix(c(ga11, ga12, ga21, ga22, ga31, ga32),2,3))
        rownames(fen)<-c("G1", "G2", "G3")
        colnames(fen)<-c("Env1", "Env2")
        Mg<-apply(fen, 1, mean)
        Menv<-apply(fen, 2, mean)
        Tg<-apply(fen, 1, sum)
        Tenv<-apply(fen, 2, sum)
        Yield<-c(fen[1,1], fen[1,2], fen[2,1], fen[2,2], fen[3,1], fen[3,2])
        Genotype<-as.factor(c("G1", "G1", "G2","G2", "G3", "G3"))
        Environment<-as.factor(c("Env1", "Env2","Env1", "Env2","Env1", "Env2"))
        mgxa<-data.frame(Genotype,Environment,Yield)
        X<-mean(Yield)
        QMR<-(CV*X/100)^2
        C<-sum(Yield)^2/6
        SQTotal<-(sum(Yield^2)-C)*R
        SQG<-(sum(Tg^2)/2-C)*R
        SQE<-(sum(Tenv^2)/3-C)*R
        SQGE<-SQTotal-SQG-SQE
        QMG<-SQG/2
        QME<-SQE
        QMGE<-SQGE/2
        VG<-round((QMG-QMR)/(2*R),2)
        VE<-round((QME-QMR)/(3*R),2)
        VGE<-round((QMGE-QMR)/R,2)  
        if(VG<0){VG<-0}
        if(VE<0){VE<-0}
        if(VGE<0){VGE<-0}
        Vr<-round(QMR,2)
        h2g<-round(VG/(VG+VGE/2+Vr/(R*2)),2)
        varcomp<-c(VG, VE, VGE, Vr, h2g)
        names(varcomp)<-c('Genetic var.', "Env. var.", "GxE var.", "Residual var.", "h2")
        # sigmaSq = '\u03c3\u00B2'
        # leg = paste0(sigmaSq,c("g", "e", "gxe", "res"))
        # names(varcomp) = c(leg, parse(text='h^2'))
        Lista1<-list(Mg, Menv,varcomp)
        # names(Lista1)<-c("Genotype Mean:", "Environment Mean:","Variances:")
        return(Lista1)}
      
      gxa2<-function(ga11, ga12, ga21, ga22, ga31, ga32, CV, R){
        fen<-t(matrix(c(ga11, ga12, ga21, ga22, ga31, ga32),2,3))
        rownames(fen)<-c("G1", "G2", "G3")
        colnames(fen)<-c("Env1", "Env2")
        Mg<-apply(fen, 1, mean)
        Menv<-apply(fen, 2, mean)
        Tg<-apply(fen, 1, sum)
        Tenv<-apply(fen, 2, sum)
        Yield<-c(fen[1,1], fen[1,2], fen[2,1], fen[2,2], fen[3,1], fen[3,2])
        Genotype<-as.factor(c("G1", "G1", "G2","G2", "G3", "G3"))
        Environment<-as.factor(c("Env1", "Env2","Env1", "Env2","Env1", "Env2"))
        mgxa<-data.frame(Genotype,Environment,Yield)
        Fig<-interaction.plot(mgxa$Environment, mgxa$Genotype, mgxa$Yield, trace.label= "Genotype", xlab = "Environments", ylab = "Yield", col = 2)
        return(Fig)}
  
      output$cont.Lear.10a.1 <- renderTable(colnames = T, rownames = F, align = "c", digits = 4,{
        tst10 = gxa1(input$Lear10.1, input$Lear10.2, 
                     input$Lear10.3, input$Lear10.4,
                     input$Lear10.5, input$Lear10.6,
                     input$Lear10.7, input$Lear10.8)
        return(list(t(as.data.frame(tst10[[1]]))))})
      
      output$cont.Lear.10a.2 <- renderTable(align = 'c',{
        tst10 = gxa1(input$Lear10.1, input$Lear10.2, 
                     input$Lear10.3, input$Lear10.4,
                     input$Lear10.5, input$Lear10.6,
                     input$Lear10.7, input$Lear10.8)
        return(list(t(as.data.frame(tst10[[2]]))))})  
      
      output$cont.Lear.10a.3 <- renderTable(colnames = T, align = "c",{
        tst10 = gxa1(input$Lear10.1, input$Lear10.2, 
                     input$Lear10.3, input$Lear10.4,
                     input$Lear10.5, input$Lear10.6,
                     input$Lear10.7, input$Lear10.8)
        names(tst10[[3]]) = c('Genetic', "Environmental", "GxE", "Residual", "Broad Heritability")
        return(list(t(as.data.frame(tst10[[3]]))))})
        # add.to.row = list(pos = list(0),
        #                   command = c("\u03c3\u00B2", '\u03c3\u00B2e','\u03c3\u00B2e',
        #                   '\u03c3\u00B2e', '\u03c3\u00B2e')))
    
    # Table headers  
    
      output$phr10.1 <- renderUI({
        table10.1 <- strong("Genotype Mean")
        HTML(paste(table10.1, sep = '<br/>'))})
      
      output$phr10.2 <- renderUI({
        table10.2 <- strong("Environment Mean")
        HTML(paste(table10.2, sep = '<br/>'))})
      
      output$phr10.3 <- renderUI({
        table10.3 <- strong("Variances")
        HTML(paste(table10.3, sep = '<br/>'))})

      output$cont.Lear.10.2 <- renderPlot({
        return(gxa2(input$Lear10.1, input$Lear10.2, 
                   input$Lear10.3, input$Lear10.4,
                   input$Lear10.5, input$Lear10.6,
                   input$Lear10.7, input$Lear10.8))})
      
# #####################################################################################################################################################################################
      
      # Lean.11
      
      output$cont.Lear.11 <- renderTable(colnames = F, rownames = T, digits = 4,{
        heterose<-function(Ng, gmd, dg){
          L1<-matrix(sample(c("A", "a"), Ng, replace = T), Ng, 1)
          if(dg!=0&dg!=1){
            n<-Ng-(dim(L1)[1]*dg)
            L1.1<-L1[1:n,]
            L1.2<-L1[(n+1):dim(L1)[1],]
            L1.3<-NULL
            for(u in 1:length(L1.2)){
              if(L1.2[u]=="a"){L1.3[u]<-"A"}
              if(L1.2[u]=="A"){L1.3[u]<-"a"}} 
            L2.1<-c(L1.1,L1.3)
            L2<-L2.1}
          if(dg==0){L2<-L1}
          if(dg==1){
            L2<-NULL
            for(i in 1:length(L1)){
              if(L1[i]=="a"){L2[i]<-"A"}
              if(L1[i]=="A"){L2[i]<-"a"}}}
          coinc<-L1==L2
          dist<-round(length(coinc[coinc==FALSE])/Ng,1)
          HS<-paste(L1, L2, sep = "")
          L1[L1=="A"]<-2
          L1[L1=="a"]<-0
          L2[L2=="A"]<-2
          L2[L2=="a"]<-0
          HS[HS=="AA"]<-2
          HS[HS=="aA"]<-1+1*gmd
          HS[HS=="Aa"]<-1+1*gmd
          HS[HS=="aa"]<-0
          Vg.L1<-sum(as.numeric(L1))
          Vg.L2<-sum(as.numeric(L2))
          Vg.F1<-sum(as.numeric(HS))
          PS<-max(Vg.L2,Vg.L1)
          PM<-(Vg.L2+Vg.L1)/2
          H<-(Vg.F1-PM)/PM*100
          Hb<-(Vg.F1-PS)/PS*100
          LIST.11<-data.frame(Vg.L1,Vg.L2,Vg.F1,H,Hb)
          names(LIST.11)<-c("Genetic variance Line 1","Genetic variance Line 2","Genetic variance F1","Heterosis (%)","Heterobeltiosis (%)")
          return(t(LIST.11))}
        if(input$Lear11.1==0) return(print("WARNING: 'Number of genes' must be different from zero"))
        return(heterose(input$Lear11.1, input$Lear11.2, input$Lear11.3))})

# #####################################################################################################################################################################################
      
      # Lean.12
      
      tester1 <- function(p, r, a, d, f){
        vt <- function(p, r, a, d, f){  
          estimation <- (p*(1-p)/2)*(1+f)*(a+(1-2*r)*d)^2  
          return(estimation)}
        par(mfrow=c(2,2))
        curve(vt(p = x, r = r, a = a, d = d, f = f), from = 0, to = 1, xlab = "p", ylab = "Genetic variablity", col = "red")
        curve(vt(p = p, r = x, a = a, d = d, f = f), from = 0, to = 1, xlab = "r", ylab = "Genetic variablity", col = "red")
        curve(vt(p = p, r = r, a = a, d = x, f = f), from = 0, to = 2, xlab = "d", ylab = "Genetic variablity", col = "red")
        curve(vt(p = p, r = r, a = a, d = d, f = x), from = 0, to = 1, xlab = "F", ylab = "Genetic variablity", col = "red")}
      tester2 <- function(p, r, a, d, f){
      vt <- function(p, r, a, d, f){
          estimation <- (p*(1-p)/2)*(1+f)*(a+(1-2*r)*d)^2  
          return(estimation)}
      vt.est <- vt(p = p, r = r, a = a, d = d, f = f)
        names(vt.est)<-"Genetic variability"
       return(vt.est)}
      output$cont.Lear.12 <- renderPlot({
        p<-input$Lear12.1
        r<-input$Lear12.2
        a<-input$Lear12.3
        d<-input$Lear12.4
        f<-input$Lear12.5
        tester1(p, r, a, d, f)})
      output$cont.Lear.12.1 <- renderPrint({
        p<-input$Lear12.1
        r<-input$Lear12.2
        a<-input$Lear12.3
        d<-input$Lear12.4
        f<-input$Lear12.5
        return(tester2(p, r, a, d, f))})
      
# #####################################################################################################################################################################################
      
      # Lean.13
      
      n.gen1 <- function(q, s){
        ng <- function(q, s){(q - s)/(q*s)}
        curve(ng(q = q, s = x), from = 1, to = 0, xlab = "Selection intensity", ylab = "Number of generations", col = "red")
        generations <- ng(q, s)
        return(generations)}
      n.gen2 <- function(q, s){ng <- function(q, s){(q - s)/(q*s)}
        generations <- ng(q, s)
        names(generations)<-"Generations to eliminate 'a':"
        return(generations)}
      out.Lear.13.1<-reactive({
        q<-input$Lear13.1
        s<-input$Lear13.2
        treze<-c(q,s)
        return(treze)})
      output$cont.Lear.13 <- renderPlot({
        treze<-out.Lear.13.1()
        if(treze[1]==0){plot.new()}
        else(n.gen1(treze[1],treze[2]))})
      output$cont.Lear.13.1 <- renderPrint({
        treze<-out.Lear.13.1()
        if(treze[1]==0){return(cat("WARNING: 'q' shoud be different from zero"))}
        else(n.gen2(treze[1],treze[2]))})
         
# #####################################################################################################################################################################################
      
      # Lean.14
      
      RSi1 <- function(hx, rxy, sel, size){
        i <- NULL
        p <- sel/size
        i.std <- dnorm(qnorm(p))/p
        Ac <- sqrt(hx)
        p2 <- (sel + 1/2)/(size + sel / (2 * size))
        i.s <- dnorm(qnorm(p2))/p2
        if (size >= 1000) {i <- i.std}
        if (size < 1000) {i <- i.s}
        response <- function(Ac = Ac, i = i, rxy = rxy){
          gs <- i*Ac*rxy}
        curve(response(Ac = x, rxy = rxy, i = i), from = 0, to = 1, xlab = "parameter", ylab = "Indirect response (additive deviation)", col = "red")
        curve(response(Ac = Ac, rxy = x, i = i), from = 0, to = 1, xlab = "parameter", ylab = "Indirect response (additive deviation)", col = "blue", add = TRUE)
        legend("bottom", xpd = TRUE, bty = "n", horiz = TRUE, c("Accuracy", "Corr.xy"), pch=c(1,1), text.col = c("red", "blue"), cex =1)
        GS <- response(Ac, i, rxy)
        return(GS)}
      RSi2 <- function(hx, rxy, sel, size){
        i <- NULL
        p <- sel/size
        i.std <- dnorm(qnorm(p))/p
        Ac <- sqrt(hx)
        p2 <- (sel + 1/2)/(size + sel / (2 * size))
        i.s <- dnorm(qnorm(p2))/p2
        if (size >= 1000) {i <- i.std}
        if (size < 1000) {i <- i.s}
        response <- function(Ac = Ac, i = i, rxy = rxy){
        gs <- i*Ac*rxy}
        GS <- response(Ac, i, rxy)
        names(GS) <- "Response to Indirect Selection: "
        return(GS)}
      output$cont.Lear.14 <- renderPlot({
        hx <- input$Lear14.1
        rxy <- input$Lear14.2
        size <- input$Lear14.3
        sel <- input$Lear14.4
        if(size==0){plot.new()}
        else(RSi1(hx, rxy, sel, size))})
      output$cont.Lear.14.1 <- renderPrint({
        hx <- input$Lear14.1
        rxy <- input$Lear14.2
        sel <- input$Lear14.4
        size <- input$Lear14.3
        if(size==0){return("Number of individuals evaluated should be different from zero")}
        else(return(RSi2(hx, rxy, sel, size)))})
      
# #####################################################################################################################################################################################      
     
      # Lean.15
        
      output$cont.Lear.15 <- renderTable(rownames = T, colnames = F,digits = 4,{
        n <- input$Lear15.1
        i <- input$Lear15.2
        h <- input$Lear15.3
      res.effect <- function(n, i, h){
        i <- i/100*n  
         se <- sqrt(1/h - 1)
         ee <- function(n, i, h){
         g <- as.data.frame(rnorm(n, 0, 1))
         set.seed(NULL)
         rownames(g) <- paste("T", 1:n, sep = "")
         colnames(g)[1] <- "genotype" 
            e <- as.data.frame(rnorm(n, 0, se))
            set.seed(NULL)
            rownames(e) <- paste("T", 1:n, sep = "")
            colnames(e)[1] <- "environment" 
            p <- g + e
            colnames(p) <- "phenotype"
            jointdata <- cbind(p, g, e)
            jointdata$gid <- rownames(g)
            selected <- jointdata[order(jointdata[,1], decreasing = TRUE)[1:i],"gid"]
            best <- jointdata[order(jointdata[,2], decreasing = TRUE)[1:i], "gid"]
            left <- best[!jointdata[order(jointdata[,2], decreasing = TRUE)[1:i], "gid"] %in% selected]
            pull <- selected[!selected %in% best]
            output <- t(data.frame(
              "Real accuracy" = sum(selected %in% jointdata[order(jointdata[,2], decreasing = TRUE)[1:i], "gid"])/length(selected),
              "Estimated accuracy" = round(sqrt(h),2),
              "Residual on wrong selected" = round(mean(jointdata[pull,3]),2), #foram beneficados pelo ambiente em 3.10
              "Residual on wrong discarted" = round(mean(jointdata[left, 3]),2),
              "Residual on right discarted" = round(mean(jointdata[selected,3]),2),
              "GV of wrong selected" = round(mean(jointdata[pull,2]),2),
              "GV of wrong discarted" = round(mean(jointdata[left, 2]),2), # esta dando maior que o de baixo, o que indica que estamos eleinando coisas boas!
              "GV of right selected" = round(mean(jointdata[selected,2]),2)))
            rownames(output) = c("Real accuracy", "Estimated accuracy", "Residual on wrong selected", "Residual on wrong discarted", "Residual on right discarted",
                                 "GV of wrong selected", "GV of wrong discarted",  "GV of right selected")
            return(output)}
          out <- ee(n, i, h)
          return(out)}
      if(input$Lear15.1==0){return(paste("WARNING: 'Number of individuals evaluated' should be different from zero"))}
      else(return(res.effect(n, i, h)))})

# #####################################################################################################################################################################################            
      
      # Lean.16
      
      replicates1 <- function(plots, sel, rep, h){
        Cj <- plots/rep
        add.rep <- function(plots, sel, rep, h){
          v <- sel/plots
          Cj <- plots/rep
          vj <- sel/Cj
          zo <- qnorm(v, lower.tail = F)
          zj <- qnorm(vj, lower.tail = F)
          RS <- (exp(-.5*(zj^2 - zo^2)) / rep * sqrt(rep/(1 + h*(rep - 1))) - 1)*100
          return(RS)}
        plot(y = add.rep(plots = plots, sel = sel, rep = 1:5, h = h), x = 1:5, xlab = "Number of replicates", ylab = "Effect on RS (%)", col = "red", type = "l")
        RSres <- add.rep(plots, sel, rep, h)
        return(cat( rep, "replicates increase", round(RSres, 2), "% the RS"))}
      replicates2 <- function(plots, sel, rep, h){
        Cj <- plots/rep
        add.rep <- function(plots, sel, rep, h){
          v <- sel/plots
          Cj <- plots/rep
          vj <- sel/Cj
          zo <- qnorm(v, lower.tail = F)
          zj <- qnorm(vj, lower.tail = F)
          RS <- (exp(-.5*(zj^2 - zo^2)) / rep * sqrt(rep/(1 + h*(rep - 1))) - 1)*100
          return(RS)}
       RSres <- add.rep(plots, sel, rep, h)
        return(cat( rep, "replicates increase", round(RSres, 2), "% the RS"))}
      out.Lear.16.1<-reactive({
      rep<-input$Lear16.1
      plots<-input$Lear16.2
      sel<-input$Lear16.3
      h<-input$Lear16.4
      a1<-c(rep,plots,sel,h)
      return(a1)})
      output$cont.Lear.16 <- renderPlot({
        a1<-out.Lear.16.1()
        if(0%in%a1|a1[2]<=a1[3]){plot.new()}
        else(replicates1(a1[2],a1[3],a1[1],a1[4]))})
      output$cont.Lear.16.1 <- renderPrint({
        a1<-out.Lear.16.1()
        if(a1[2]<=a1[3]){return(cat("WARNING: Values should be different from zero and the number of plots should be greater than the number of genotypes selected"))}
        else(return(replicates2(a1[2],a1[3],a1[1],a1[4])))})
      
# #####################################################################################################################################################################################                  
      
      # Lean.17
      
      budget1 <- function(size, sel, h, rep, time, cost){
        internal <- function(size, sel, h, rep, time, cost){
          i <- NULL
          p <- sel/size
          i.std <- dnorm(qnorm(p))/p
          p2 <- (sel + 1/2)/(size + sel / (2 * size))
          i.s <- dnorm(qnorm(p2))/p2
          if (size >= 1000) {i <- i.std}
          if (size < 1000) {i <- i.s}
          Vg <- h
          Ve <- (1-h)*rep 
          h2 <- Vg/(Vg+Ve/rep)
          Ac <- sqrt(h2)
          RS <- i*Ac*sqrt(Vg) 
          RS.t <- RS/time
          Total.cost <- size*cost*rep
          How.fair <- RS.t/Total.cost
          output <- How.fair
          return(output)}
        par(mfrow = c(1,3))  
        curve(internal(size = size, sel = sel, h = h, rep = x, time = time, cost = cost), from = 1, to = 5, col = "red", xlab = "Number of replicates", ylab = "Gain/Cost")
        curve(internal(size = size, sel = sel, h = h, rep = rep, time = x, cost = cost), from = .2, to = 7, col = "red", xlab = "Years per cycle", ylab = "Gain/Cost")  
        curve(internal(size = size, sel = sel, h = x, rep = rep, time = time, cost = cost), from = 0, to = 1, col = "red", xlab = "Heritability", ylab = "Gain/Cost")
        balance <- internal(size, sel, h, rep, time, cost)
        return(cat("The ratio between genetic gain and the total cost per cycle of selection is", balance))}
      budget2 <- function(size, sel, h, rep, time, cost){
        internal <- function(size, sel, h, rep, time, cost){
          i <- NULL
          p <- sel/size
          i.std <- dnorm(qnorm(p))/p
          p2 <- (sel + 1/2)/(size + sel / (2 * size))
          i.s <- dnorm(qnorm(p2))/p2
          if (size >= 1000) {i <- i.std}
          if (size < 1000) {i <- i.s}
          Vg <- h
          Ve <- (1-h)*rep 
          h2 <- Vg/(Vg+Ve/rep)
          Ac <- sqrt(h2)
          RS <- i*Ac*sqrt(Vg) 
          RS.t <- RS/time
          Total.cost <- size*cost*rep
          How.fair <- RS.t/Total.cost
          output <- How.fair
          return(output)}
        balance <- internal(size, sel, h, rep, time, cost)
        return(cat("Genetic gain / total cost =", balance))}
      output$cont.Lear.17 <- renderPlot({
        size<-input$Lear17.1
        sel<-input$Lear17.2
        h<-input$Lear17.3
        rep<-input$Lear17.4
        time<-input$Lear17.5
        cost<-input$Lear17.6
        if(input$Lear17.1==0|input$Lear17.6==0){return(NULL)}
        budget1(input$Lear17.1, input$Lear17.2, 
                input$Lear17.3, input$Lear17.4,
                input$Lear17.5, input$Lear17.6)})
      output$cont.Lear.17.1 <- renderPrint({
        size<-input$Lear17.1
        sel<-input$Lear17.2
        h<-input$Lear17.3
        rep<-input$Lear17.4
        time<-input$Lear17.5
        cost<-input$Lear17.6
        if(input$Lear17.1==0|input$Lear17.6==0){return("'Number of genotypes evaluated' and 'Cost to evaluate one plot' should be different from zero")}
        return(budget2(input$Lear17.1, input$Lear17.2, 
                       input$Lear17.3, input$Lear17.4,
                       input$Lear17.5, input$Lear17.6))})
      
# #####################################################################################################################################################################################                        
      
      # Lean.18
      
      output$cont.Lear.18 <- renderTable({
      output$ui.lear18.1 <- renderUI({
          sliderInput("Lear18.4", "Frequency of haplotype ab", min=0, max=1, value=(1-(input$Lear18.1+input$Lear18.2+input$Lear18.3)))})
        output$ui.lear18.2 <- renderUI({  
        sliderInput("Lear18.6", "Frequency of a", min=0, max=1, value=(1-input$Lear18.5))})
      output$ui.lear18.3 <- renderUI({  
          sliderInput("Lear18.8", "Frequency of b", min=0, max=1, value=(1-input$Lear18.7))})
        AB<-input$Lear18.1
        Ab<-input$Lear18.2
        aB<-input$Lear18.3
        ab<-input$Lear18.4
        A<-input$Lear18.5
        a<-input$Lear18.6
        B<-input$Lear18.7
        b<-input$Lear18.8
      LD <- function(AB, Ab, aB, ab, A, a, B, b){
        D <- AB*ab - Ab*aB
        phase <- NULL
        if (D > 0) {phase <- paste("The loci tend to be in cis because D is positive =", D)}
        if (D < 0) {phase <- paste("The loci tend to be in trans because D is negative =", D)}
        if (D == 0) {phase <- paste("The loci are independent because D is", D)}
        r2 <- round(D^2/(A*a*B*b),2)
        output <- list(phase,r2)
        names(output) <- c("Phase", "LD")
        return(output)}
      if((AB+Ab+aB+ab)!=1){return("The sum of haplotypes must be equal to 1")}
      if((A+a)!=1){return("The sum of alleles from locus A must be equal to 1")}
      if((B+b)!=1){return("The sum of alleles from locus B must be equal to 1")}
      return(LD(AB, Ab, aB, ab, A, a, B, b))})

      
# #####################################################################################################################################################################################                        
      
      # Lean.19
      
      generations1 <- function(LDlose, c){
        overtime <- function(LDlose, c){round(log(1 - LDlose/100) / log(1 - c), 1)}
        t <- overtime(LDlose, c) 
        output <- c(paste(t, "cycles to lose", LDlose,"% of the LD"))
        return(output)}
      generations2 <- function(LDlose, c){
        overtime <- function(LDlose, c){round(log(1 - LDlose/100) / log(1 - c), 1)}
        t <- overtime(LDlose, c) 
        par(mfrow=c(1,2))
        curve(overtime(LDlose = LDlose, c = x), from = 0, to = .5, xlab = "recombination fraction", ylab = paste("Generations to lose", LDlose, "% of the LD"), col = "red")
        curve(overtime(LDlose = x, c = c), from = 0, to = 100, xlab = "LD lost (%)", ylab = paste("Generations to lose", LDlose, "% of the LD"), col = "red")}
      output$cont.Lear.19 <- renderPlot({
        if(input$Lear19.1==0|input$Lear19.2==0){(plot.new())}
        else(generations2(input$Lear19.1, input$Lear19.2))})
      output$cont.Lear.19.1 <- renderPrint({
        return(generations1(input$Lear19.1, input$Lear19.2))})
      
# #####################################################################################################################################################################################      
      
# Phenotipic Breeding
      
# #####################################################################################################################################################################################
      
      # SubItem 11
      
      out.11.1 <- reactive({
        if(input$ex.11!=FALSE){
          TABE.1<-read.table("ex11.dat",header = TRUE, sep = "\t")
          return(TABE.1)}
        inFile <- input$file1
        if (is.null(inFile))
          return(NULL)
        TABE<-read.table(inFile$datapath, header=input$header, sep=input$sep, 
                         quote=input$quote)
        return(TABE)})
      out.11.2 <- reactive({
        switch(input$table.1,
               "Table" = out.11.1(),
               "Data Structure" = str(out.11.1()),
               "Names"=names(out.11.1()))})
      output$contents11 <- renderPrint({
        if(input$info.11!=TRUE)
          return (out.11.2())
        else
          includeScript("BeBreederDataFileInformation.txt")})
      
# #####################################################################################################################################################################################
      
      # SubItem 12
      
      out.12.1 <- reactive({
        if(is.null(out.11.1()))
          return(NULL)
        data12<-out.11.1()
        for(i in 1:dim(data12)[2]){
          if(names(data12)[i]=="y")
            data12$y<-as.numeric(as.character(data12$y))
          else
            data12[,i]<-as.factor(data12[,i])}
        if(input$res.12.name!=FALSE)
          return(names(out.11.1()))
        out.mod<-reactive({
          if(input$res.12.run){
            formula(isolate(input$model))}
          else
            NULL})
        if(is.null(out.mod()))
          return(NULL)
        inMod <- out.mod()
        if(input$ml != TRUE){
          mod.12<-lmer(inMod, REML=FALSE, data=data12, na.action =  na.omit)}
        if(input$ml != FALSE){
          mod.12<-lmer(inMod, REML=TRUE, data=data12, na.action =  na.omit)}
        PM<-as.matrix(cbind(order(ranef(mod.12)$Genotype,decreasing =T),mean(data12$y)+ranef(mod.12)$Genotype[order(ranef(mod.12)$Genotype,decreasing =T),]))
        colnames(PM)<-c("Genotype","Predicted")
        if(input$res.12.fix != TRUE){
          PM.1<-PM[1:round((dim(PM)[1])*input$res.12.2/100,0),]
          BLUP.12<-ranef(mod.12)
          BLUE.12<-NULL
          AM.12<-NULL}
        else {
          PM.1<-NULL
          BLUP.12<-ranef(mod.12)
          BLUP.12<-BLUP.12[names(BLUP.12)!="Genotype"]
          BLUE.12<-as.matrix(tapply(data12$y,data12$Genotype,mean)-mean(data12$y))
          colnames(BLUE.12)<-c("BLUE")
          AM.12<-as.matrix(tapply(data12$y,data12$Genotype,mean))
          colnames(AM.12)<-c("Adj. Means")}
        MEAN.12<- as.data.frame(as.numeric(fixef(mod.12)[1]))
        colnames(MEAN.12)<-c("Mean")
        ANOVA.12<-anova(mod.12,test="F")
        switch(input$res.12.1,
               "Summary" = summary(mod.12),
               "ANOVA"= ANOVA.12,
               "Analysis of Deviance"= rand(mod.12),
               "BLUP" = BLUP.12,
               "BLUE" = BLUE.12,
               "Adjusted Means"= AM.12,
               "Mean" = MEAN.12,
               "Predicted Means"= PM.1)})
      output$contents12 <- renderPrint({
        if(input$info.12!=TRUE)
          return(out.12.1())
        else
          includeScript("BeBreederEstatisticalInformation.txt")})
      output$downloadData12 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.12,'.txt', sep = "")},
        content = function(file) {
          write.table(out.12.1(), file, sep="\t")})
      
# ##############################################################################################################################################################################
      
      # SubItem 13
      
      out.1311.0 <- reactive({
        if(input$ex.1311!=FALSE){
          TABE.1<-read.table("griffing.txt", header = T)
          return(TABE.1)}
        inFile1311 <- input$file1311
        if (is.null(inFile1311))
          return(NULL)
        TABE<-read.table(inFile1311$datapath, header=input$header1311, sep=input$sep1311, 
                         quote=input$quote1311,stringsAsFactors=FALSE)
        return(TABE)})
      output$contents1311 <- renderPrint({
        if(input$info.1311!=TRUE)
          return (out.1311.0())
        else 
          includeScript("BeBreederDiallel.txt")})
      out.1321.0 <- reactive({
        if(input$ex.1321!=FALSE){
          TABE.1<-read.table("gardner.txt", header = T)
          return(TABE.1)}
        inFile1321 <- input$file1321
        if (is.null(inFile1321))
          return(NULL)
        TABE<-read.table(inFile1321$datapath, header=input$header1321, sep=input$sep1321, 
                         quote=input$quote1321,stringsAsFactors=FALSE)
        return(TABE)})
      output$contents1321 <- renderPrint({
        if(input$info.1321!=TRUE)
          return (out.1321.0())
        else 
          includeScript("BeBreederDiallel.txt")})
      out.1331.0 <- reactive({
        if(input$ex.1331!=FALSE){
          TABE.1<-read.table("NCII.txt", header = T)
          return(TABE.1)}
        inFile1331 <- input$file1331
        if (is.null(inFile1331))
          return(NULL)
        TABE<-read.table(inFile1331$datapath, header=input$header1331, sep=input$sep1331, 
                         quote=input$quote1331,stringsAsFactors=FALSE)
        return(TABE)})
      output$contents1331 <- renderPrint({
        if(input$info.1331!=TRUE)
          return (out.1331.0())
        else 
          includeScript("BeBreederDiallel.txt")})
      out.1311.1<-reactive({
        data<-out.1311.0()
        if (is.null(data))
          return(NULL)
        a<-griffing(data)
        switch(input$res.1312,
               "Output"= a)})
      out.1321.1<-reactive({
        data<-out.1321.0()
        if (is.null(data))
          return(NULL)
        a<-gardner(data)
        switch(input$res.1322,
               "Output"= a)})
      out.1331.1<-reactive({
        data<-out.1331.0()
        if (is.null(data))
          return(NULL)
        a<-factorial(data)
        switch(input$res.1332,
               "Output"= a)})
      output$contents1312 <- renderPrint({
        if(input$info.1312!=TRUE)
          return (out.1311.1())
        else
          includeScript("BeBreederDiallel.txt")})
      output$downloadData1312 <- downloadHandler(
        filename = function() {
          paste(input$name.txt.1312,'.txt', sep = "")},
        content = function(file) {
          write.table(out.1311.1(), file, sep="\t")})
      output$contents1322 <- renderPrint({
        if(input$info.1322!=TRUE)
          return (out.1321.1())
        else
          includeScript("BeBreederDiallel.txt")})
      output$downloadData1322 <- downloadHandler(
        filename = function() {
          paste(input$name.txt.1322,'.txt', sep = "")},
        content = function(file) {
          write.table(out.1321.1(), file, sep="\t")})
      output$contents1332 <- renderPrint({
        if(input$info.1332!=TRUE)
          return (out.1331.1())
        else
          includeScript("BeBreederDiallel.txt")})
      output$downloadData1332 <- downloadHandler(
        filename = function() {
          paste(input$name.txt.1332,'.txt', sep = "")},
        content = function(file) {
          write.table(out.1331.1(), file, sep="\t")})
     
# ##############################################################################################################################################################################
      
      # SubItem 14
      
      out.14<-reactive({ 
        if(input$ex.14!=FALSE){
          IndSel.1<-read.table("ex14.dat",header = TRUE, sep = "\t")
          return(IndSel.1)}
        inFile14 <- input$file14
        if (is.null(inFile14))
          return(NULL)
        IndSel<-read.table(inFile14$datapath, header=input$header.14, sep=input$sep.14, 
                           quote=input$quote.14)
        return(IndSel)})
      output$contents14.1<- renderPrint({
        if(input$info.141!=TRUE)
          return(out.14())
        else 
          includeScript("BeBreederDataIndexFile.txt")})
      out.14.1<-reactive({
        if(is.null(out.14()))
          return(NULL)
        IndSel<-out.14()
        if(input$res.14.name!=FALSE)
          return(names(IndSel))
        y<-as.numeric(unlist(strsplit(input$vec.14,",")))
        if(length(y)<(dim(IndSel)[2]-1)|length(y)>(dim(IndSel)[2]-1))
          return(NULL)
        A<-scale(as.matrix(IndSel[,2:dim(IndSel)[2]]))
        Genotype<-as.matrix(IndSel[,1])
        IndSelec<-function(y,A,Genotype){
          IS<-y%*%t(A)
          Sel<-cbind(Genotype[order(IS, decreasing = T)],IS[order(IS, decreasing = T)])
          return(Sel)}
        INDSEL<-IndSelec(y,A,Genotype)
        colnames(INDSEL)<-c("Genotype","Index")
        INDSEL.1<-INDSEL[1:(dim(INDSEL)[1]*input$res.14/100),]
        return(INDSEL.1)})
      output$contents14.2<- renderPrint({
        if(input$info.142!=TRUE)
          return(out.14.1())
        else 
          includeScript("BeBreederDataIndexAnalyse.txt")})
      output$downloadData14 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.14,'.txt',sep="")},
        content = function(file) {
          write.table(out.14.1(), file, sep="\t")})
      
# ##############################################################################################################################################################################
      
      # SubItem 15
      
      out.15.0 <- reactive({
        if(input$ex.15.0!=FALSE){
          TABE15.1<-read.table("ex15.dat",header = TRUE, sep = "\t")
          return(as.data.frame(TABE15.1))}
        inFile15 <- input$file15
        if (is.null(inFile15))
          return(NULL)
        TABE15<-read.table(inFile15$datapath, header=input$header15, sep=input$sep15, 
                           quote=input$quote15)
        return(as.data.frame(TABE15))})
      output$contents151 <- renderPrint({
        if(input$info.151!=TRUE)
          return (out.15.0())
        else 
          includeScript("BeBreederDataCorrelation.txt")})
      out.15.1<-reactive({
        if(is.null(out.15.0())){
          return(NULL)}
        data15<-as.matrix(out.15.0())
        name15<-colnames(out.15.0())
        data15.1<-data15[,-c(which(name15=="Geno"))]
        data15.2<-apply(data15.1,2,function(x){
          as.numeric(as.character(x))})
        COR.T<-correlation(data15.2)
        LISTA<-COR.T
        names(LISTA)<-c("Correlation","pValue","Num.Obs")
        return(LISTA)})
      output$contents15.2<- renderPrint({
        if(is.null(out.15.0())){
          return(NULL)}
        if(input$info.152!=TRUE){
          switch(input$res.15.2,
                 "Genotype Correlation"= round(out.15.1()[[1]],3),
                 "pValue"=round(out.15.1()[[2]],3))}
        else 
          includeScript("BeBreederDataCorrelation.txt")})
      output$contents15.1<- renderPlot({
        if(is.null(out.15.0())){
          return(NULL)}
        if(input$info.152!=TRUE){
          data15<-as.matrix(out.15.0())
          name15<-colnames(out.15.0())
          data15.1<-data15[,-c((name15=="Geno"))]
          data15.2<-apply(data15.1,2,function(x){
            as.numeric(as.character(x))})
          data15.2 %>% correlate() %>% network_plot(min_cor = input$res.15.3)
          }
        else 
          includeScript("BeBreederDataCorrelation.txt")})
      output$contents15.3<- renderPlot({
        if(is.null(out.15.0())){
          return(NULL)}
        if(input$info.152!=TRUE){
          corrgram(out.15.1()[[1]], type = "cor", lower.panel = panel.shade, upper.panel = panel.pie)}
        else 
          includeScript("BeBreederDataCorrelation.txt")})
      
      output$downloadData15 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.15,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.15.1()), file, sep="\t")})
      
# ##############################################################################################################################################################################
      
      # SubItem 16
      
      out.16.0 <- reactive({
        if(input$ex.16!=FALSE){
          TABE16.1<-read.table("ex14.dat",header = TRUE, sep = "\t")
          return(as.data.frame(TABE16.1))}
        inFile16 <- input$file16
        if (is.null(inFile16))
          return(NULL)
        TABE16<-read.table(inFile16$datapath, header=input$header16, sep=input$sep16, 
                           quote=input$quote16)
        return(as.data.frame(TABE16))})
      output$contents161 <- renderPrint({
        if(input$info.161!=TRUE)
          return (out.16.0())
        else 
          includeScript("BeBreederDataPathAnalyse.txt")})
      out.16.1<-reactive({
        if(is.null(out.16.0()))
          return(NULL)
        data16<-as.matrix(out.16.0())
        name16<-names(out.16.0())
        if(input$res.16.name!=FALSE)
          return(list(name16))
        y.test<-input$name.16
        if(!(y.test%in%name16))
          return(NULL)
        Main.Trait<-data16[,which(name16==input$name.16)]
        x16<-data16[,-c(which(name16=="Genotype"),which(name16==input$name.16))]
        cor.y<-correlation(Main.Trait,x16)$correlation
        cor.x<-correlation(x16)$correlation
        PATH16<-path.analysis(cor.x,cor.y)
        COR.T<-correlation(data16[,-c(which(name16=="Genotype"))])
        N16<-paste("Path Analysis: ",input$name.16)
        LISTA<-list(PATH16,COR.T)
        names(LISTA)<-c(N16,"Traits Correlation")
        return(LISTA)})
      output$contents16.2<- renderPrint({
        if(input$info.162!=TRUE){
          switch(input$res.16.2,
                 "Path Analysis"= out.16.1()[[1]],
                 "Traits Correlation"= out.16.1()[[2]])}
        else 
          includeScript("BeBreederDataPathAnalyse.txt")})
      output$downloadData16 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.16,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.16.1()), file, sep="\t")})
      
# ##############################################################################################################################################################################
      
      # SubItem 17
      
      out.17.0 <- reactive({
        if(input$ex.17!=FALSE){
          TABE17<-read.table("ex17.dat",header = TRUE, sep = "\t")}
        else{
          inFile17 <- input$file17
          if (is.null(inFile17))
            return(NULL)
          TABE17<-read.table(inFile17$datapath, header=input$header17, sep=input$sep17,
                             quote=input$quote17)}
        rownames(TABE17)<-TABE17$Genotype
        name.17<-names(TABE17)
        TABE17.1<-TABE17[,-(which(name.17=="Genotype"))]
        return(TABE17.1)})
      output$contents171 <- renderPrint({
        if(input$info.171!=TRUE)
          return (out.17.0())
        else 
          includeScript("BeBreederbiplot.txt")})
     out.17.1<-reactive({
      if(is.null(out.17.0()))
        return(NULL)
      data17<-as.matrix(out.17.0())
      SUMMARY<- summary(prcomp(data17,scale = TRUE))
      SCORES.GEN<- summary(prcomp(data17,scale = TRUE))$x
      SCORES.ENV<- summary(prcomp(data17,scale = TRUE))$rotation
      LIST.17<-list(SUMMARY,SCORES.GEN,SCORES.ENV)
      return(LIST.17)})
    output$contents173<- renderPlot({
      if(is.null(out.17.0()))
        return(NULL)
      data17<-as.matrix(out.17.0())
      name17<-colnames(data17)
      CP<-(summary(prcomp(data17,scale = TRUE))[[1]])^2
      CP1<-paste("PC 1 (",(round(CP[1]/sum(CP),3))*100,"%",")", sep="")
      CP2<-paste("PC 2 (",(round(CP[2]/sum(CP),3))*100,"%",")", sep="")
      biplot(prcomp(data17,scale = TRUE),
             ylab=CP2,xlab=CP1)
      abline(v=0,h=0)})
    out.17.3<-reactive({
      switch(input$res.172,
             "Summary"= out.17.1()[[1]],
             "Scores of G"= out.17.1()[[2]],
             "Scores of E"= out.17.1()[[3]])})
     output$contents172<- renderPrint({
        if(input$info.172!=TRUE)
          return(out.17.3())
        else 
          includeScript("BeBreederbiplot.txt")})
      output$downloadData172 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.172,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.17.3()), file, sep="\t")})
      
# #####################################################################################################################################################################################
      
      # Sub.174

      out.174.1<-reactive({
        if(is.null(out.17.0()))
          return(NULL)
        mydata <- na.omit(out.17.0()) 
        mydata <- scale(mydata)
        return(mydata)})
      out.174.00<-reactive({
      if(is.null(out.17.0()))
          return(NULL)
        mydata<-out.174.1()
        set.seed(0001)
        fit <- kmeans(mydata, input$num.174)
        set.seed(0002)
        d <- dist(mydata, method = "euclidean")
        set.seed(0003)
        fit2 <- hclust(d, method="ward")
        set.seed(0004)
        groups <- cutree(fit2, k=input$num.174)
        set.seed(0005)
        mydata2 <- data.frame(mydata, fit$cluster)
        if(input$summary.174=="Cluster information"){
          return(fit)}
        if(input$summary.174=="Genotype groups"){
          return(mydata2)}
        if(input$summary.174=="Cluster mean"){
          return(aggregate(mydata,by=list(fit$cluster),FUN=mean))}})
      output$contents174<-renderPrint({
        if(is.null(out.17.0()))
          return(NULL)
        out.174.00()})
      output$downloadData174 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.174,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.174.00()), file, sep="\t")})
      output$contents175<-renderPlot({
        if(is.null(out.17.0()))
          return(NULL)
        mydata<-out.174.1()
        set.seed(0001)
        fit <- kmeans(mydata, input$num.175)
        set.seed(0002)
        d <- dist(mydata, method = "euclidean")
        set.seed(0003)
        fit2 <- hclust(d, method="ward")
        set.seed(0004)
        groups <- cutree(fit2, k=input$num.175)
        if(input$sel.175!=FALSE){
        plot(fit2,xlab="Distance")
        rect.hclust(fit2, k=input$num.175, border="red")}
        if(input$sel.175==FALSE){
          plot(fit2,xlab="Distance")}})

# #####################################################################################################################################################################################
      
      # Sub.18
      
      output$ui.18.1 <- renderUI({
      switch(input$res.18,
             "Completely randomized design"={
               textInput("name.txt.18.1","Treatments",'A,B,C')},
             "Randomized complete block design"={
               textInput("name.txt.18.1","Treatments",'A,B,C,D,E')},
             "Latin square design"={
               textInput("name.txt.18.1","Treatments",'A,B,C,D')},
             "Balanced Incomplete Block Designs"={
               textInput("name.txt.18.1","Treatments",'A,B,C,D,E')},
             "Lattice designs"={
               textInput("name.txt.18.1","Treatments",'a,b,c,d,e,f,g,h,i')},
             "Alpha designs"={
               textInput("name.txt.18.1","Treatments",'a,b,c,d,e,f,g,h,i,j,k,l,m,n,o')},
             "Augmented block designs"={
               textInput("name.txt.18.1","Treatments-01",'A,B,C,D')},
             "Split plot designs"={
               textInput("name.txt.18.1","Treatments-01",'A,B,C,D')},
             "Factorial"={
               textInput("name.txt.18.1","Treatments",'3,2')})})
      output$ui.18.2 <- renderUI({
        switch(input$res.18,
               "Completely randomized design"={
                 textInput("name.txt.18.2","Replications","4, 3, 4")},
               "Randomized complete block design"={
                 numericInput("name.txt.18.2","Blocks",min = 1,value = 4)},
               "Latin square design"={},
               "Balanced Incomplete Block Designs"={},
               "Lattice designs"={
                 numericInput("name.txt.18.2","Replications",min = 2,max = 3,value = 3)},
               "Alpha designs"={
                 numericInput("name.txt.18.2","Replications",min = 1,value = 3)},
               "Augmented block designs"={
                 numericInput("name.txt.18.2","Blocks",min = 1,value = 3)},
               "Split plot designs"={
                 numericInput("name.txt.18.2","Replications",min = 1,value = 3)},
               "Factorial"={
                 numericInput("name.txt.18.2","Replications",min = 1,value = 3)})})
      output$ui.18.3 <- renderUI({
        switch(input$res.18,
               "Completely randomized design"={selectInput("name.txt.18.3", "Method for to randomize:",
                                                           choices = c("Wichmann-Hill", "Marsaglia-Multicarry", 
                                                                       "Super-Duper", "Mersenne-Twister", "Knuth-TAOCP", 
                                                                       "user-supplied", "Knuth-TAOCP-2002", "default"),width = 400)},
               "Randomized complete block design"={selectInput("name.txt.18.3", "Method for to randomize:",
                                                               choices = c("Wichmann-Hill", "Marsaglia-Multicarry", 
                                                                           "Super-Duper", "Mersenne-Twister", "Knuth-TAOCP", 
                                                                           "user-supplied", "Knuth-TAOCP-2002", "default"),width = 400)},
               "Latin square design"={selectInput("name.txt.18.3", "Method for to randomize:",
                                                  choices = c("Wichmann-Hill", "Marsaglia-Multicarry", 
                                                              "Super-Duper", "Mersenne-Twister", "Knuth-TAOCP", 
                                                              "user-supplied", "Knuth-TAOCP-2002", "default"),width = 400)},
               "Balanced Incomplete Block Designs"={
                 numericInput("name.txt.18.3","size block (k)",min = 1,value = 4)},
               "Lattice designs"={
                 selectInput("name.txt.18.3", "Method for to randomize:",
                             choices = c("Wichmann-Hill", "Marsaglia-Multicarry", 
                                         "Super-Duper", "Mersenne-Twister", "Knuth-TAOCP", 
                                         "user-supplied", "Knuth-TAOCP-2002", "default"),width = 400)},
               "Alpha designs"={
                 numericInput("name.txt.18.3","size block (k)",min = 1,value = 3)},
               "Augmented block designs"={
                 textInput("name.txt.18.3","Treatments-02",'t,u,v,w,x,y,z')},
               "Split plot designs"={textInput("name.txt.18.3","Treatments-02",'a,b,c')},
               "Factorial"={selectInput("name.txt.18.3", "Designs type:",
                                        choices = c("rcbd","crd","lsd"),width = 400)})})
      out.18.0 <- reactive({
        switch(input$res.18,
               "Completely randomized design"={
                 trt<-do.call(c,strsplit(input$name.txt.18.1,","))
                 rep<-do.call(c,strsplit(input$name.txt.18.2,","))
                 outdesign <- design.crd(trt=trt,
                                         r=as.numeric(as.character(rep)),
                                         seed=input$num.18.4,
                                         serie=0,kinds =input$name.txt.18.3)
                 book1 <- outdesign$book
                 return(book1)},
               "Randomized complete block design"={
                 trt<-do.call(c,strsplit(input$name.txt.18.1,","))
                 rep<-input$name.txt.18.2
                 outdesign <- design.rcbd(trt=trt,
                                         r=as.numeric(as.character(rep)),
                                         seed=input$num.18.4,
                                         serie=2,kinds =input$name.txt.18.3)
                book2<- zigzag(outdesign) 
                 return(book2)},
               "Latin square design"={
                 trt<-do.call(c,strsplit(input$name.txt.18.1,","))
                 outdesign <- design.lsd(trt, seed=input$num.18.4, serie=2,kinds =input$name.txt.18.3)
                 book <- zigzag(outdesign)
                 return(book)},
               "Balanced Incomplete Block Designs"={
                 trt<-do.call(c,strsplit(input$name.txt.18.1,","))
                 k<-input$name.txt.18.3
                 outdesign <- design.bib(trt,k=as.numeric(k), seed=input$num.18.4, serie=2)
                 book <- outdesign$book
                 return(book)},
               "Lattice designs"={
                 trt<-do.call(c,strsplit(input$name.txt.18.1,","))
                 rep<-input$name.txt.18.2
                 outdesign <-design.lattice(trt, r = rep, serie = 2, seed = input$num.18.4, kinds =input$name.txt.18.3)
                 book <- outdesign$book
                 return(book)},
               "Alpha designs"={
                 trt<-do.call(c,strsplit(input$name.txt.18.1,","))
                 k<-input$name.txt.18.3
                 rep<-input$name.txt.18.2
                 outdesign <- design.alpha(trt,r=as.numeric(rep),k=as.numeric(k),seed=input$num.18.4)
                 book <- outdesign$book
                 return(book)},
               "Augmented block designs"={
                 trt1 <- do.call(c,strsplit(input$name.txt.18.1,","))
                 trt2 <- do.call(c,strsplit(input$name.txt.18.3,","))
                 rep<-input$name.txt.18.2
                 outdesign <- design.dau(trt1, trt2, r=as.numeric(rep), seed=input$num.18.4, serie=2)
                 book <- zigzag(outdesign)
                 return(book)},
               "Split plot designs"={
                 trt1 <- do.call(c,strsplit(input$name.txt.18.1,","))
                 trt2 <- do.call(c,strsplit(input$name.txt.18.3,","))
                 rep<-input$name.txt.18.2
                 outdesign <-design.split(trt1,trt2,r=as.numeric(rep),serie=2,seed=input$num.18.4)
                 book <- outdesign$book
                 return(book)},
               "Factorial"={
                 trt <- do.call(c,strsplit(input$name.txt.18.1,","))
                 rep<-input$name.txt.18.2
                 outdesign <-design.ab(trt, r=rep, seed=input$num.18.4,serie=2)
                 book <- outdesign$book
                 return(book)})})
      output$contents18<-renderPrint({
        if(input$info.18!=TRUE)
          return(out.18.0())
        else 
          includeScript("BeBreederDataDesing.txt")})
      output$downloadData18 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.18,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.18.0()), file, sep="\t")})
      
# #####################################################################################################################################################################################
      
# Molecular Breeding
      
# #####################################################################################################################################################################################
      
      # Mol121
      
      out.mol121.0 <- reactive({
        if(input$exmol121!=FALSE){
          TABEmol121.1 <- read.table("dupont.txt", header = TRUE, na.strings = "-", sep = "")
          return(TABEmol121.1)}
        inFilemol121 <- input$filemol121
        if (is.null(inFilemol121))
          return(NULL)
        TABEmol121<-read.table(inFilemol121$datapath, header=input$headermol121, sep=input$sepmol121,
                              quote=input$quotemol121, na.strings =input$na.stringsmol121)
        return(TABEmol121)})
      
      output$contentsmol121 <- renderPrint({
        if(input$info.mol121!=TRUE)
          return (out.mol121.0())
        else 
          includeScript("BeBreederMol12.txt")})
      
# #####################################################################################################################################################################################
      
      # Mol1212
      
      out.mol1212.0 <- reactive({
        if(input$exmol1212!=FALSE){
          TABEmol1212.1 <- read.table("hapmap.txt", header = TRUE, na.strings = "-", sep = "")
          return(TABEmol1212.1)}
        inFilemol1212 <- input$filemol1212
        if (is.null(inFilemol1212))
          return(NULL)
        TABEmol1212<-read.table(inFilemol1212$datapath, header=input$headermol1212, sep=input$sepmol1212,
                               quote=input$quotemol1212, na.strings =input$na.stringsmol1212)
        return(TABEmol1212)})
      output$contentsmol1212 <- renderPrint({
        if(input$info.mol1212!=TRUE)
          return (out.mol1212.0())
        else 
          includeScript("BeBreederMol12.txt")})
      
# #####################################################################################################################################################################################
      
      # Mol122
      
      out.mol122<-reactive({
        if (is.null(out.mol121.0()))
          return(NULL)
        Z<-out.mol121.0()
        hap<- out.mol1212.0()
        withProgress(message = 'Making analysis', value = 0, {
          incProgress(1/2, detail = paste("Doing part", 1))
        if(is.null(hap)){
        out.1.mol11<-reactive({
          if(input$res7.mol122){
            x<-raw.data(as.matrix(Z), frame=input$res5.mol122, imput=input$res4.mol122, imput.type=input$res4.5.mol122, sweep.sample =input$res3.mol122 , call.rate =input$res2.mol122 , maf=input$res1.mol122, outfile =input$res6.mol122)
            return(x)}})}
        if(!is.null(hap)){
        out.1.mol11<-reactive({
          if(input$res7.mol122){
            x<-raw.data(as.matrix(Z), frame=input$res5.mol122, imput=input$res4.mol122, imput.type=input$res4.5.mol122, sweep.sample =input$res3.mol122 , call.rate =input$res2.mol122 , maf=input$res1.mol122, outfile =input$res6.mol122, hapmap = hap )
            return(x)}})}
        X <- out.1.mol11()
        incProgress(2/2, detail = paste("Doing part", 2))})
        return(X)})
      out.Mol122.End <- reactive({
        if (is.null(out.mol122()))
          return(NULL)
        x<-out.mol122()
          switch(input$res8.mol122,
                 "M.clean"=x$M.clean, 
                 "Hapmap"=x$Hapmap, 
                 "Report"=x$report)})
      output$contentsmol122 <- renderPrint({
        if(input$info.mol122!=TRUE)
          return (out.Mol122.End())
        else 
          includeScript("BeBreederMol12.txt")})
      output$downloadDatamol122 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.mol122,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.Mol122.End()), file, sep="\t")})
      
# #####################################################################################################################################################################################
      
      # Mol241
      
      out.mol241.0 <- reactive({
        if(input$exmol241!=FALSE){
          TABEmol241.1<-read.table("exMol241.txt",header = TRUE, sep = "\t")
          return(TABEmol241.1)}
        inFilemol241 <- input$filemol241
        if (is.null(inFilemol241))
          return(NULL)
        TABEmol241<-read.table(inFilemol241$datapath, header=input$headermol241, sep=input$sepmol241,
                              quote=input$quotemol241)
        return(as.data.frame(TABEmol241))})
      output$contentsmol241 <- renderPrint({
        if(input$info.mol241!=TRUE)
          return (out.mol241.0())
        else
          includeScript("BeBreederMol24.txt")})   
      
# #####################################################################################################################################################################################
      
      # Mol242
      
      out.mol242.0 <- reactive({
        if(input$exmol242!=FALSE){
          TABEmol242.1<-read.table("exMol242.txt",header = TRUE, sep = "\t")
          return(TABEmol242.1)}
        inFilemol242 <- input$filemol242
        if (is.null(inFilemol242))
          return(NULL)
        TABEmol242<-read.table(inFilemol242$datapath, header=input$headermol242, sep=input$sepmol242,
                               quote=input$quotemol242)
        return(as.data.frame(TABEmol242))})
      output$contentsmol242 <- renderPrint({
        if(input$info.mol242!=TRUE)
          return (out.mol242.0())
        else
          includeScript("BeBreederMol24.txt")})
      
# #####################################################################################################################################################################################
      
      # Mol1243
      
      out.mol243<-reactive({
        if (is.null(out.mol241.0()))
          return(NULL)
        if (is.null(out.mol242.0()))
          return(NULL)
        Z<-out.mol242.0()
        B<-out.mol241.0()
        Gen<-Z$Genotype
        Mar<-Z[,-(which(colnames(Z)=="Genotype"))]
        rownames(Mar)<-Gen
        efec<-B$Effect
        Mar2<-Mar[,colnames(Mar)%in%as.character(B$Marker)]
        Mar3<-Mar2[,match(colnames(Mar2),as.character(B$Marker))]
        Predit<-as.matrix(Mar3)%*%as.matrix(efec)
        SEL<-as.matrix(Predit[order(Predit,decreasing = T),])
        GS.SEL<-as.matrix(SEL[1:round((dim(SEL)[1])*input$res.mol243/100,0),])
        colnames(GS.SEL)<-"Genotype"
        return(GS.SEL)})
      output$contentsmol243 <- renderPrint({
        out.mol243()})
      output$downloadDatamol243 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.mol243,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.mol243()), file, sep="\t",col.names = F,row.names = F)})
      
# #####################################################################################################################################################################################
      
      # Mol11
      
      out.mol11.0 <- reactive({
        if(input$exmol11!=FALSE){
          TABEmol11.1 <- fread("exMol11.txt", data.table=FALSE)
          rownames(TABEmol11.1) <- TABEmol11.1[,1]
          TABEmol11.1 <- as.matrix(TABEmol11.1[,-1])
          return(as.data.frame(TABEmol11.1))}
        inFilemol11 <- input$filemol11
        if (is.null(inFilemol11))
          return(NULL)
        TABEmol11<-read.table(inFilemol11$datapath, header=input$headermol11, sep=input$sepmol11,
                              quote=input$quotemol11)
        rownames(TABEmol11) <- TABEmol11[,1]
        TABEmol11 <- as.matrix(TABEmol11[,-1])
        return(as.data.frame(TABEmol11))})
      output$contentsmol11 <- renderPrint({
        if(input$info.mol11!=TRUE)
          return (out.mol11.0()[1:15,1:15])
        else 
          includeScript("BeBreederMol13.txt")})
      
# #####################################################################################################################################################################################
      
      # Mol13 
      
      output$res3.mol13 <- renderUI({
        if(input$res1.mol13=="VanRaden"){
          selectInput("res4.mol13","Choose Results:",
                      choices = c("G.additive",
                                  "G.dominant"))}})
      out.mol13<-reactive({
        if (is.null(out.mol11.0()))
          return(NULL)
        withProgress(message = 'Making analysis', value = 0, {
          incProgress(1/2, detail = paste("Doing part", 1))
        Z<-as.matrix(out.mol11.0())
        out.1.mol11<-reactive({
          if(input$res5.mol13){
            x<-G.matrix(Z, method=input$res1.mol13, format =input$res2.mol13)
          return(x)}})
        X <- out.1.mol11()
        incProgress(2/2, detail = paste("Doing part", 2))})
        return(X)})
      out.Mol13.End <- reactive({
        x<-out.mol13()
        if(input$res1.mol13=="VanRaden"){
        switch(input$res4.mol13,
               "G.additive"={a<-x$Ga},
               "G.dominant"={a<-x$Gd})
          return(a)}
        if(input$res1.mol13!="VanRaden"){
          return (x)}})
      output$contentsmol13 <- renderPrint({
        if(input$info.mol13!=TRUE)
          return (out.Mol13.End())
        else 
          includeScript("BeBreederMol13.txt")})
      output$downloadDatamol13 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.mol13,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.Mol13.End()), file, sep="\t")})
      
# #####################################################################################################################################################################################
      
      # Mol21
      
      out.mol21.0 <- reactive({
        if(input$exmol21!=FALSE){
          TABEmol21.1<-read.table("R.dat",header = TRUE, sep = "\t")
          return(TABEmol21.1)}
        inFilemol21 <- input$filemol21
        if (is.null(inFilemol21))
          return(NULL)
        TABEmol21<-read.table(inFilemol21$datapath, header=input$headermol21, sep=input$sepmol21,
                              quote=input$quotemol21)
        return(as.data.frame(TABEmol21))
      })
      output$contentsmol21 <- renderPrint({
        if(input$info.mol21!=TRUE)
          return (out.mol21.0())
        else
          includeScript("BeBreederMol21.txt")})
      
# #####################################################################################################################################################################################
      
      # Mol22
      
      out.mol22.0 <- reactive({
        if(input$exmol22!=FALSE){
          TABEmol22.1<-read.table("M1.dat",header = TRUE, sep = "\t")
          return(TABEmol22.1)}
        inFilemol22 <- input$filemol22
        if (is.null(inFilemol22))
          return(NULL)
        TABEmol22<-read.table(inFilemol22$datapath, header=input$headermol22, sep=input$sepmol22,
                              quote=input$quotemol22)
        return(as.data.frame(TABEmol22))})
      output$contentsmol22 <- renderPrint({
        if(input$info.mol22!=TRUE)
          return (out.mol22.0())
        else
          includeScript("BeBreederMol21.txt")})

######################################################################################################################################################################################
      
      # Mol23 
      
      out.mol23<-reactive({
        if(is.null(out.mol21.0()))
         return(NULL)
        if(is.null(out.mol22.0()))
          return(NULL)
        B<-out.mol21.0()
        Z<-out.mol22.0()
        withProgress(message = 'Making analysis', value = 0, {
          incProgress(1/2, detail = paste("Doing part", 1))
        Gen<-Z[,colnames(Z)=="Genotype"]
        Mar<-Z[,-(which(colnames(Z)=="Genotype"))]
        Gen2<-B[,colnames(B)=="Genotype"]
        y<-as.matrix(B[,-(which(colnames(B)=="Genotype"))])
        Mar2<-Mar[Gen%in%as.character(Gen2),]
        M<-Mar2[match(Gen,as.character(Gen2)),]
        y<-as.matrix(as.numeric(as.character(y)))
        rownames(y)<-Gen2
        M<-apply(M, 2,function(x) as.numeric(as.character(x)))
        rownames(Mar)<-Gen
        ans1 <- mixed.solve(y,Z=M)
        ans2 <- mixed.solve(y,K=A.mat(M))
        accuracy <- cor(y,ans2$u)
        Y.est<-as.matrix(ans2$u)
        rownames(Y.est)<-Gen2
        incProgress(2/2, detail = paste("Doing part", 1))})
        switch(input$res.mol23,
               "Marker effects"= ans1$u,
               "Breeding values"= Y.est,
               "Accuracy"= accuracy)})
     output$contentsmol23<-renderPrint({
       return (out.mol23())})
      output$downloadDatamol23 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.mol23,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.mol23()), file, sep="\t")})
        
# #####################################################################################################################################################################################
      
      # Mol31
      
      out.mol31.0 <- reactive({
        if(input$exmol31!=FALSE){
          load("exMol31.dat")
          return(as.data.frame(TABEmol31.1))}
        inFilemol31 <- input$filemol31
        if (is.null(inFilemol31))
          return(NULL)
        TABEmol31<-read.table(inFilemol31$datapath, header=input$headermol31, sep=input$sepmol31, 
                           quote=input$quotemol31)
        return(as.data.frame(TABEmol31))})
      output$contentsmol31 <- renderPrint({
        if(input$info.mol31!=TRUE)
          return (out.mol31.0())
        else 
          includeScript("BeBreederMol3.txt")})
      
# #####################################################################################################################################################################################
      
      # Mol32
      
      out.mol32.0 <- reactive({
        if(input$exmol32!=FALSE){
          load("exMol32.dat")
          return(as.data.frame(TABEmol32.1))}
        inFilemol32 <- input$filemol32
        if (is.null(inFilemol32))
          return(NULL)
        TABEmol32<-read.table(inFilemol32$datapath, header=input$headermol32, sep=input$sepmol32, 
                              quote=input$quotemol32)
        return(as.data.frame(TABEmol32))})
      output$contentsmol32 <- renderPrint({
        if(input$info.mol32!=TRUE)
          return (out.mol32.0())
        else 
          includeScript("BeBreederMol3.txt")})
      
# #####################################################################################################################################################################################
      
      # Mol33  
      
      out.mol33<-reactive({
        if(is.null(out.mol31.0()))
          return(NULL)
        if(is.null(out.mol32.0()))
          return(NULL)
        pheno<-out.mol31.0()
        geno<-out.mol32.0()
        withProgress(message = 'Making analysis', value = 0, {
          incProgress(1/2, detail = paste("Doing part", 1))
        scores <- GWAS(pheno,geno,plot=F)
        incProgress(2/2, detail = paste("Doing part", 2))})
        return(scores)})
      output$contentsmol33.1<-renderPrint({
        if(is.null(out.mol33()))
          return(NULL)
        out.mol33()})
      output$contentsmol33.2<-renderPlot({
        if(is.null(out.mol33()))
          return(NULL)
        scores<-out.mol33()
        Xpos<-list()
        for(i in 1:length(unique(scores$chrom))){
          cr<-subset(scores,chrom==i)
          pcr<-dim(cr)[1]/2
          Xpos[[i]]<-cr$pos==pcr}
        plot(abs(scores$y),col=scores$chrom, xlab ="Chromosomes",ylab = "-log10(p)",xaxt='n',pch=19)
        axis(side=1, at=which(unlist(Xpos)==TRUE), labels=c(labels=unique(scores$chrom)))
        abline(h=input$LOD, col = 2)})
      output$downloadDatamol33 <- downloadHandler(
        filename = function() { 
          paste(input$name.txt.mol33,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.mol33()), file, sep="\t")})
      
# #####################################################################################################################################################################################
      
      # Mol41
      
      out.mol41.0 <- reactive({
        if(input$exmol41!=FALSE){
          if(input$res.mol41=="Structure"){
            a= read.table("SSR.str") 
            n=nrow(a)
            popn=length(levels(a[,2]))
            c=(ncol(a)-2)/2
            TABEmol41.1=read.structure("SSR.str", n.ind = n, n.loc = c, col.lab = 1, col.pop = 2, onerowperind = T, row.marknames = 0, col.others = 0)}
          if(input$res.mol41=="Genpop"){
            TABEmol41.1=read.genepop("SNP.gen") }
          return(TABEmol41.1)}
        inFilemol41 <- input$filemol41
        if (is.null(inFilemol41))
          return(NULL)
        if(input$res.mol41=="Structure"){
          a<-read.table(inFilemol41$datapath)
          n=nrow(a)
          c=(ncol(a)-2)/2
          ale.name<-runif(1)
          write.table(a,paste(ale.name,".str",sep = ""),col.names = FALSE,row.names = FALSE,quote = FALSE)
          TABEmol41=read.structure(paste(ale.name,".str",sep = ""), n.ind = n, n.loc = c, col.lab = 1, col.pop = 2, onerowperind = T, row.marknames = 0, col.others = 0)
          file.remove(paste(ale.name,".str",sep = ""))
          rm(a,n, c, ale.name)}
        if(input$res.mol41=="Genpop"){
          b<-read.table(inFilemol41$datapath,sep="?")
          ale.name2<-runif(1)
          write.table(b,paste(ale.name2,".gen",sep = ""),col.names = FALSE,row.names = FALSE,quote = FALSE)
          TABEmol41=read.genepop(paste(ale.name2,".gen",sep = ""))
        file.remove(paste(ale.name2,".gen",sep = ""))
        rm(b,ale.name2)}
        return(TABEmol41)})
      out.mol.end.41<-reactive({
        if (is.null(out.mol41.0()))
          return(NULL)
        as.data.frame(out.mol41.0())})
      output$contentsmol41 <- renderPrint({
        if(input$info.mol41!=TRUE)
          return (out.mol.end.41())
        else 
          includeScript("BeBreederMol4.txt")})
      
# #####################################################################################################################################################################################      
      
      # Mol42
     
      out.mol42<-reactive({
        if (is.null(out.mol41.0()))
          return(NULL)
        obj<-out.mol41.0()
        switch(input$res.mol42,
               "Summary"= {
                 sum<-summary(obj)
                 S<-list("Number of individuals:"= sum$n,
                         "Group sizes:"= sum$n.by.pop,
                         "Number of alleles per loci:"= sum$loc.n.all,
                         "Number of alleles per group:"= sum$pop.n.all,
                         "Percentage of missing data:"= sum$NA.perc,
                         "Observed heterozygosity:"= sum$Hobs,
                         "Expected heterozygosity:"= sum$Hexp)
                 return(S)},
               "Expected heterozygosity"= Hs(obj))})
      output$contentsmol42 <- renderPrint({
        out.mol42()})
      output$downloadDatamol42 <- downloadHandler(
        filename = function() {
          paste(input$name.txt.mol42,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.mol42()), file, sep="\t")})
      
# #####################################################################################################################################################################################      
      
      # Mol43
      
      output$ui.mol43 <- renderUI({
      if(input$res.mol43=="Principal component analysis by population 1"){
        textInput("title.mol43","Main Title:")}})
      output$ui.mol43.2 <- renderUI({
        if(input$res.mol43=="Population Structure"){
          textInput("name.txt.mol44","Type the File Name")}})
      output$ui.mol43.3 <- renderUI({
        if(input$res.mol43=="Population Structure"){
          downloadButton('downloadDatamol44', 'Download')}})
      output$ui.mol43.4 <- renderUI({
        if(input$res.mol43=="Principal component analysis by population 1"){
          numericInput("res4.mol43", "Label size:",min = 0,max = 5,value = 0.8,step = 0.1)}})
      output$contentsmol43 <- renderPlot({
        if (is.null(out.mol41.0()))
          return(NULL)
        obj<-out.mol41.0()
        sum=summary(obj)
        He=Hs(obj)
        X <- scaleGen(obj, NA.method="mean")
        pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=4)
        popn<-length(unique(obj$pop))
        add_legend <- function(...) {
          opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                      mar=c(0, 0, 0, 0), new=TRUE)
          on.exit(par(opar))
          plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
          legend(...)}
        switch(input$res.mol43,
               "Number of alleles per loci"={barplot(obj$loc.n.all, ylab="Number of alleles per loci")},
               "Expected Heterozygosity per loci"={barplot(sum$Hexp, main="Expected Heterozygosity per loci", ylab="He", xlab="Loci")},
               "Expected Heterozygosity per population"={barplot(He, main="Expected Heterozygosity per population", ylab="He", xlab="Population")},
               "Observed Heterozygosity per loci"={barplot(sum$Hobs, main="Observed Heterozygosity per loci", ylab="Ho", xlab="Loci")},
               "Heterozygosity: expected-observed"={barplot(sum$Hexp-sum$Hobs, main="Heterozygosity: expected-observed", ylab="Hexp - Hobs")},
               "F index:  F = (He - Ho)/He"={barplot(sum$Hexp-sum$Hobs/sum$Hexp, main="F index:  F = (He - Ho)/He", ylab= "F")},
               "Population sample size"={
                 barplot(sum$n.by.pop, main="Population sample size", ylab="Number of genotypes",las=3)},
               "Principal component analysis by population 1"={
               s.class(pca1$li, obj$pop,xax=1,yax=2, col=rainbow(popn, alpha = 1), axesell=FALSE, pch=19, grid=FALSE, cellipse=0, clabel=input$res4.mol43, cstar=0, addaxes=1)
               title(input$title.mol43)
               add.scatter.eig(pca1$eig, 3,1,2)},
               "Principal component analysis by individual"={
                 s.label(pca1$li, sub="CA 1-2",cpoint="", boxes=F, clabel=0.7)},
               "Dendrogram by individual"={dist <- aboot(obj, tree="upgma", dist = nei.dist, missing="mean", sample = 100, cutoff = 0)},
               "Dendrogram by population"={
                 obj_pop=genind2genpop(obj)
                 dist_pop <- aboot(obj_pop, tree="upgma", dist = rogers.dist, sample = 10, cutoff = 50)},
               "Principal component analysis by population 2"={
                 obj_pop=genind2genpop(obj)
                 X <- scaleGen(obj_pop, NA.method="mean")
                 pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=4)
                 s.label(pca1$li, sub="CA 1-2",cpoint="", boxes=T, clabel=1)
               },
               "Population Structure"={
                 dapc1 <- dapc(obj, n.da=100, n.pca=nrow(obj@tab)/3)
                 temp <- optim.a.score(dapc1)
                 dapc2 <- dapc(obj, n.da=100, n.pca=temp$best)
                 compoplot(dapc2, lab="", posi=list(x=12,y=-.01), cleg=.7)
                 output$downloadDatamol44 <- downloadHandler(
                   filename = function() {
                     paste(input$name.txt.mol44,'.txt',sep="")},
                   content = function(file) {
                     write.table(capture.output(dapc2$posterior), file, sep="\t",row.names = FALSE,col.names = FALSE)})},
               "Neighbor joining dendrogram"={
                 X <- scaleGen(obj, NA.method="mean")
                 D <- dist(X)
                 tre <- nj(D)
                 myCol=rainbow(popn)[as.integer(pop(obj))]
                 plot(tre, type = "unr", show.tip.lab = FALSE)
                 tiplabels(col = myCol, pch = 20)
                 add_legend("topright", legend=unique(obj$pop), pch=20, 
                            col=unique(myCol),
                            horiz=FALSE, bty='n', cex=0.8)})}) 

# #####################################################################################################################################################################################      
      
      # Mol51
      
      out.mol51<-reactive({ 
        if(input$exmol51!=FALSE){
          TAB.mol51.1<-read.table("PS.txt",header = T)
          return(TAB.mol51.1)}
        inFilemol51 <- input$filemol51
        if (is.null(inFilemol51))
          return(NULL)
        TAB.mol51.0<-read.table(inFilemol51$datapath, header=input$headermol51, sep=input$sepmol51, 
                           quote=input$quotemol51)
        return(TAB.mol51.0)})
      output$contentsmol51 <- renderPrint({
        if(input$info.mol51!=TRUE){
          if (is.null(out.mol51())){
            return(NULL)}
          return (as.data.frame(out.mol51()))}
          else 
          includeScript("BeBreederMol5.txt")})
      
# #####################################################################################################################################################################################      
      
      # Mol52
      
      out.mol52<-reactive({ 
        if(input$exmol52!=FALSE){
          TAB.mol52.1<-read.table("Z.txt")
          return(TAB.mol52.1)}
        inFilemol52 <- input$filemol52
        if (is.null(inFilemol52))
          return(NULL)
        TAB.mol52<-read.table(inFilemol52$datapath, header=input$headermol52, sep=input$sepmol52, 
                                quote=input$quotemol52)
        return(TAB.mol52)})
      output$contentsmol52 <- renderPrint({
        if(input$info.mol52!=TRUE)
          return (out.mol52())
        else 
          includeScript("BeBreederMol5.txt")})
      
# #####################################################################################################################################################################################      
      
      # Mol53
      
      out.mol53.0<-reactive({
        if(is.null(out.mol51()))
          return(NULL)
        if(is.null(out.mol52()))
          return(NULL)
        PS1<-out.mol51()
        row.names(PS1)<-PS1[,1]
        PS<-as.data.frame(as.character(PS1[,-1]))
        colnames(PS)<-"PS"
        Z<-as.matrix(out.mol52())
        a.53<-popgen(Z,subgroups=PS)
        return(a.53)})
      out.mol53<-reactive({
        if (is.null(out.mol53.0()))
          return(NULL)
        obj.53<-out.mol53.0()
        switch(input$res.mol53,
               "Whole"=obj.53$whole,
               "By Group"=obj.53$bygroup)})
      output$contentsmol53 <- renderPrint({
        if(input$info.mol53!=TRUE)
          return (out.mol53())
        else 
          includeScript("BeBreederMol5.txt")})
      output$downloadDatamol53 <- downloadHandler(
        filename = function() {
          paste(input$name.txt.mol53,'.txt',sep="")},
        content = function(file) {
          write.table(capture.output(out.mol53()), file, sep="\t")})
      
# #####################################################################################################################################################################################      
      
      # Mol61
      
      out.mol61<-reactive({ 
        if(input$exmol61!=FALSE){
          TAB.mol61.1<-read.table("TABEmol61.txt")
          return(TAB.mol61.1)}
        inFilemol61 <- input$filemol61
        if (is.null(inFilemol61))
          return(NULL)
        TAB.mol61<-read.table(inFilemol61$datapath, header=input$headermol61, sep=input$sepmol61, 
                              quote=input$quotemol61)
        return(TAB.mol61)})
      output$contentsmol61 <- renderPrint({
        if(input$info.mol61!=TRUE)
          return (out.mol61()[1:20,1:10])
        else 
          includeScript("BeBreederMol6.txt")})
      
# #####################################################################################################################################################################################      
      
      # Mol62
      
           output$res.mol63<-renderUI({
             switch(input$res.mol62,
                    "Population Size"={},
                    "Find Clusters"={},
                    "Scatterplot"={},
                    "Compoplot"={checkboxInput("res.mol64","Identification", FALSE)})})
            output$contentsmol62 <- renderPlot({
        if (is.null(out.mol61()))
          return(NULL)
        tab1<-out.mol61()
        pop<-tab1[,1]
        X<-as.data.frame(t(tab1[,-1]))
        XX<-apply(X,2,function(x){
          as.numeric(as.character(x))})
        rownames(XX)<-rownames(X)
        tab<-t(XX)
        loc.fac<-colnames(tab)
        loc.n.all<-apply(XX,1,function(x1){
          length(unique(x1))})
        names(loc.n.all)<-colnames(tab)
        names(pop)<-rownames(tab)
        x<-as.genind(tab)
        x$loc.n.all<-loc.n.all
        x$pop<-as.factor(pop)
        grp <- find.clusters(x,  n.pca = length(unique(pop)), n.clust = length(unique(pop)))
        dapc1 <- dapc(x, grp$grp,  n.pca=length(unique(pop)), n.da=length(unique(pop)))
        switch(input$res.mol62,
               "Population Size"={barplot(table(pop(x)), col=funky(17), las=3,xlab="", ylab="Sample size")},
               "Find Clusters"={table.value(table(pop(x), grp$grp), col.lab=paste("Cluster", 1:length(unique(pop))),
                                                         row.lab=rownames(table(pop(x), grp$grp)))},
                 "Scatterplot"={
                   scatter(dapc1, posi.da="bottomright", bg="white", pch=17:(17+length(unique(pop))))},
                 "Compoplot"={
                   if(input$res.mol64!=FALSE){
                     set.seed(1)
                     compoplot(dapc1, posi="bottomright",
                                                        txt.leg=paste("Cluster", 1:length(unique(pop))), lab=rownames(dapc1$tab),
                                                        ncol=1, xlab="", col=funky(length(unique(pop))),legend = T)}
                   if(input$res.mol64==FALSE){
                     set.seed(1)
                     compoplot(dapc1, posi="bottomright",
                                                        txt.leg=paste("Cluster", 1:length(unique(pop))), lab="",
                                                        ncol=1, xlab="", col=funky(length(unique(pop))),legend = T)}})}) 
      
# #####################################################################################################################################################################################      
    })
  
###############
### The End ###
###############
  