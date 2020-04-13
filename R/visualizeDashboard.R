#' visualizeDashboard
#' @param SNPs
#' @param SNPSummary
#' @keywords internal
#' @importFrom grDevices colorRampPalette
#' @importFrom reshape2 dcast
#' @importFrom shinythemes shinytheme
#' @import rio
#' @importFrom plotly plotlyOutput renderPlotly plot_ly layout
#' @import conquer.d3js
#' @import htmlwidgets
#' @import shiny
#' @importFrom DT DTOutput renderDT
#' @importFrom BioCircos BioCircosOutput renderBioCircos
#' @importFrom BiocGenerics toTable
#' @return [[NULL]]
visualizeDashboard <- function(SNPs,SNPSummary){

  KEGG_DATA <- prepare_KEGG(species = "hsa",
                                              KEGG_Type = "KEGG",
                                              keyType = "kegg")
  KEGG_DATA$ENSG <- lapply(X = KEGG_DATA$PATHID2EXTID,
                           FUN = entrezToENSEMBL)

  ensg_symb <- merge(BiocGenerics::toTable(org.Hs.egENSEMBL),BiocGenerics::toTable(org.Hs.egSYMBOL),by="gene_id")


  shiny::addResourcePath("logo", directoryPath = system.file("logo", package = "CONQUER"))
  #load QTLs
  pQTLs <- conquer.db::pQTLs
  meQTLs <- conquer.db::meQTLs
  miQTLex <- conquer.db::miQTLexperiment
  miQTLpred <- conquer.db::miQTLpredict

  ui <- shiny::navbarPage(title = shiny::div(shiny::img(src = "logo/CONQUER.png", style="margin-top:-10px;")),shiny::tags$head(HTML("<title>test</title>")),
                          selected = "Tissue Specific",
                          theme = shinythemes::shinytheme("flatly"),
                          shiny::tabPanel("Tissue Specific",
                                          shiny::fluidPage(
                                            shiny::sidebarLayout(
                                              shiny::sidebarPanel(
                                                width = 3,
                                                shiny::selectInput(inputId = "tissueSel",
                                                                   label = "Select Tissue:",
                                                                   choices = colnames(SNPSummary) %>% sort()),
                                                shiny::uiOutput("responsiveUI_C")
                                              ),
                                              shiny::mainPanel(
                                                shiny::tabsetPanel(id = "tis_spec",
                                                  shiny::tabPanel("Overview",
                                                           shiny::br(),
                                                           shiny::br(),
                                                           shiny::fluidRow(
                                                             shiny::column(12,align='center',
                                                                           conquer.d3js::ConquerRingOutput("RingPlot",
                                                                                                           width = 900,
                                                                                                           height = 900)
                                                                           )
                                                           )

                                                           ),
                                                  shiny::tabPanel("Modules",value = "mod",
                                                                  shiny::br(),
                                                                  shiny::br(),
                                                                  shiny::fluidRow(shiny::column(6,
                                                                                                plotly::plotlyOutput("moduleHeat")),
                                                                                  shiny::column(6,
                                                                                                DT::DTOutput("moduleTable"))
                                                                  ),
                                                                  shiny::br(),
                                                                  shiny::br(),
                                                                  shiny::fluidRow(shiny::column(12,
                                                                                         DT::DTOutput("eQTL_SNP_Table")

                                                                    )
                                                                  ),
                                                                  shiny::fluidRow(shiny::tags$h3("Module-SNP association:")),
                                                                  shiny::br(),
                                                                  shiny::fluidRow(
                                                                    shiny::column(7,DT::DTOutput("eQTL_check")),
                                                                    shiny::column(5,conquer.d3js::ConquerViolinOutput("module_Violin"))
                                                                    )
                                                                    )
                                                )
                                              )
                                            )

                                          )),
                          shiny::tabPanel("All Tissues",
                                          shiny::fluidPage(
                                            shiny::sidebarLayout(
                                              shiny::sidebarPanel(
                                                width = 3
                                                ),
                                              shiny::mainPanel(
                                                shiny::tabsetPanel(
                                                  shiny::tabPanel("Pathway Overview",
                                                           shiny::fluidRow(
                                                             if(!is.null(SNPSummary))
                                                             {
                                                               shiny::column(12,align = 'center',
                                                                             conquer.d3js::ConquerEdgeOutput("EdgePlot",
                                                                                                             width = 1100,
                                                                                                             height = 1100))
                                                             }else{
                                                               shiny::h3("No summary file provided")
                                                             }
                                                             )),
                                                  shiny::tabPanel("pQTL Overview",
                                                                  shiny::br(),
                                                                  DT::DTOutput("pQTLOverview"),
                                                                  shiny::br(),
                                                                  shiny::br(),
                                                                  DT::DTOutput("pQTLOverview_LD")),
                                                  shiny::tabPanel("meQTL Overview",
                                                                  shiny::br(),
                                                                  DT::DTOutput("meQTLOverview"),
                                                                  shiny::br(),
                                                                  shiny::br(),
                                                                  DT::DTOutput("meQTLOverview_LD")),
                                                  shiny::tabPanel("miQTL Overview",
                                                                  shiny::tabsetPanel(
                                                                    shiny::tabPanel("Experimental",
                                                                                    shiny::br(),
                                                                                    shiny::br(),
                                                                                    DT::DTOutput("mi_overview_exp"),
                                                                                    shiny::br(),
                                                                                    shiny::br(),
                                                                                    DT::DTOutput("mi_overview_exp_LD")),
                                                                    shiny::tabPanel("Predicted",
                                                                                    shiny::br(),
                                                                                    shiny::br(),
                                                                                    DT::DTOutput("mi_overview_pred"),
                                                                                    shiny::br(),
                                                                                    shiny::br(),
                                                                                    DT::DTOutput("mi_overview_pred_LD"))
                                                                  )
                                                                  )

                                                )
                                              )

                                            )

                                          )),
                          shiny::tabPanel("Single SNP",
                                          shiny::fluidPage(
                                            shiny::sidebarLayout(
                                              shiny::sidebarPanel(
                                                width = 3,
                                                shiny::selectInput(inputId = "snpSel",
                                                                   label = "Select SNP:",
                                                                   choices = names(SNPs)),
                                                shiny::uiOutput("responsiveUI_A"),
                                                shiny::uiOutput("responsiveUI_B")
                                              ),
                                              shiny::mainPanel(
                                                shiny::tabsetPanel(id = "single",
                                                  shiny::tabPanel(
                                                    title = "Linkage Disequilibrium", value="LD",
                                                    shiny::br(),
                                                    shiny::downloadButton("downloadLocus", "Download Locus Zoom"),
                                                    shiny::br(),
                                                    shiny::fluidRow(
                                                      shiny::column(12, align = 'center',
                                                                    shiny::uiOutput("LocusHeader"),
                                                                    conquer.d3js::ConquerLocusZoomOutput("LocusZoom",
                                                                                                         width = 500,
                                                                                                         height = 600)
                                                                    )
                                                    ),
                                                    shiny::downloadButton("downloadLD", "Download Linkage Disequilibrium"),
                                                    shiny::br(),
                                                    shiny::br(),
                                                    shiny::fluidRow(
                                                      shiny::column(12, align = 'center',
                                                                    DT::DTOutput("LDTable")
                                                                    )
                                                    ),
                                                    shiny::fluidRow(
                                                      shiny::column(5, align = 'left',
                                                                    shiny::uiOutput("PMIDHeader"),
                                                                    shiny::br(),
                                                                    shiny::br(),
                                                                    DT::DTOutput("PMIDTable")
                                                      )
                                                    )
                                                  ),
                                                  shiny::tabPanel(title="Chromosomal interactions",value = "CI",
                                                                  shiny::br(),
                                                                  shiny::uiOutput("CIMessage"),
                                                                  BioCircos::BioCircosOutput(outputId = "Circos",width = 1000,height = 1000)),
                                                  shiny::tabPanel(title = "Chromatin States", value = "CS",
                                                                  shiny::br(),
                                                                  shiny::fluidRow(
                                                                    shiny::column(1,shiny::uiOutput("placeHolder_downloadData")),
                                                                    shiny::column(1,shiny::uiOutput("placeHolder_downloadPlot"))
                                                                    ),
                                                                  shiny::br(),
                                                                  plotly::plotlyOutput("chromatinStates",height = 1100)),
                                                  shiny::tabPanel(title="QTLs", value = "QTLs",
                                                                  shiny::tabsetPanel(id="QTLs_tab",
                                                                                     shiny::tabPanel(title="eQTLs",value = "eQTLs",
                                                                                                     shiny::fluidRow(shiny::tags$h1("Hive plot"),
                                                                                                                     shiny::downloadButton("downloadHive", "Download Hive plot"),
                                                                                                                     shiny::br(),
                                                                                                                     shiny::br(),
                                                                                                                     conquer.d3js::ConquerHiveOutput("hive_plot",
                                                                                                                                                     width = 900,
                                                                                                                                                     height = 700)
                                                                                                     )
                                                                                                     ,
                                                                                                     shiny::fluidRow(
                                                                                                       shiny::tags$h1("eQTLs"),
                                                                                                       shiny::downloadButton("downloadeQTLs", "Download eQTLs"),
                                                                                                       shiny::br(),
                                                                                                       shiny::br(),
                                                                                                       shiny::column(7,
                                                                                                                     shiny::br(),
                                                                                                                     DT::DTOutput("eQTLsTable")
                                                                                                       ),
                                                                                                       shiny::column(5,
                                                                                                                     shiny::br(),
                                                                                                                     shiny::uiOutput("DownloadViolinButton"),
                                                                                                                     shiny::br(),
                                                                                                                     conquer.d3js::ConquerViolinOutput("cis_Violin")
                                                                                                       ))
                                                                                                     ),
                                                                                     shiny::tabPanel(title="pQTLs",
                                                                                     shiny::column(12,
                                                                                                   shiny::br(),
                                                                                                   DT::DTOutput("pQTLsTable")
                                                                                                   )
                                                                                     ),
                                                                                     shiny::tabPanel(title="miQTLs",
                                                                                                     shiny::fluidRow(
                                                                                                       shiny::column(6,
                                                                                                                     shiny::h4("Predicted miQTLs"),
                                                                                                                     shiny::br(),
                                                                                                                     DT::DTOutput("miQTLsPred_table")),
                                                                                                       shiny::column(6,
                                                                                                                     shiny::h4("Experimental determined miQTLs"),
                                                                                                                     shiny::br(),
                                                                                                                     DT::DTOutput("miQTLsExp_table"))
                                                                                                     )),
                                                                                     shiny::tabPanel(title="meQTLs",
                                                                                                     shiny::column(12,
                                                                                                                   shiny::br(),
                                                                                                                   DT::DTOutput("meQTLsTable"))))),
                                                  shiny::tabPanel(title="Gene expression",value = "GeX",
                                                                  shiny::fluidRow(
                                                                    shiny::tags$h1("Gene expression"),
                                                                    shiny::br(),
                                                                    shiny::fluidRow(
                                                                      shiny::column(1,shiny::uiOutput("geneExpr_downloadData")),
                                                                      shiny::column(1,shiny::uiOutput("geneExpr_downloadPlot"))
                                                                    ),
                                                                    shiny::br(),
                                                                    plotly::plotlyOutput("expressionHeatmap",
                                                                                         width = 1350,
                                                                                         height = 1000)
                                                                  )
                                                                  )
                                                )

                                              )
                                            )

                                          ))
                          )




  # Define server logic required to draw a histogram
  server <- function(input, output) {
    ####Tissue Specific####
    output$RingPlot <- conquer.d3js::renderConquerRing({
      shiny::req(input$tissueSel)
      conquer.d3js::ConquerRing(SNPSummary,KEGG_DATA=KEGG_DATA,tissue = input$tissueSel,hoverID = "#sss")
    })

    output$responsiveUI_C <- shiny::renderUI({
      tissue <- req(input$tissueSel)
      if(input$tis_spec == "mod" | input$tis_spec == "check"){
        shiny::selectInput(inputId = "ModuleSel",
                           label = "Select Module:",
                           choices = names(SNPSummary[["canOR",tissue]]))
      }
    })

    output$moduleHeat <- plotly::renderPlotly({
      tissue <- shiny::req(input$tissueSel)
      module <- shiny::req(input$ModuleSel)
      genes <- SNPSummary[["Module_Genes",tissue]][module] %>% unlist()
      geneExpression <- SNPSummary[["Expression",tissue]]
      corMatrix <- cor(geneExpression[,genes],method = "spearman")

      colnames(corMatrix) <- ensg_symb[match(colnames(corMatrix),ensg_symb$ensembl_id),"symbol"]
      rownames(corMatrix) <- ensg_symb[match(rownames(corMatrix),ensg_symb$ensembl_id),"symbol"]
      plotly::plot_ly(
        x = rownames(corMatrix),y = colnames(corMatrix),zmin = -1, zmax = 1,colors = colorRampPalette(c("blue","white","red"))(100),
        z = t(as.matrix(corMatrix)), type = "heatmap") %>% plotly::layout(title =sprintf("Correlation matrix of module %s (%s)",module,tissue))
    })


    output$moduleTable <-  DT::renderDT({
      tissue <- shiny::req(input$tissueSel)
      module <- shiny::req(input$ModuleSel)
      data <- SNPSummary[["canOR",tissue]][[module]]
      disp_data <- cbind(signif(data$OR,2),
                         signif(data$low_CI,2),
                         signif(data$up_CI,2),
                         data$bg,
                         signif(data$P.Val,2),
                         data$intersect)
      rownames(disp_data) <- data$PathwayName
      colnames(disp_data) <- c("OR","low_CI","up_CI","background","P.Val","intersect")
      disp_data
    })

    output$eQTL_SNP_Table <- DT::renderDT({
      tissue <- shiny::req(input$tissueSel)
      module <- shiny::req(input$ModuleSel)
      cbind(SNPSummary[["Module_SNPs_eQTLs",tissue]][[module]][,c("gene","SNP")],
            signif(SNPSummary[["Module_SNPs_eQTLs",tissue]][[module]][,c("pValue","pValueThreshold","Pval.ratio")],2))

    },selection = 'single')



    moduleQTLdata <- shiny::reactive({
      tissue <- shiny::req(input$tissueSel)
      module <- shiny::req(input$ModuleSel)
      rowSel <- shiny::req(input$eQTL_SNP_Table_rows_selected)
      data <- cbind(SNPSummary[["Module_SNPs_eQTLs",tissue]][[module]][,c("gene","SNP")],
                    signif(SNPSummary[["Module_SNPs_eQTLs",tissue]][[module]][,c("pValue","pValueThreshold","Pval.ratio")],2))
      SelectedSNP <- data[rowSel,"SNP"]
      genes <- SNPSummary[["Module_Genes",tissue]][[module]]
      genes_symb <- ensg_symb[match(genes,ensg_symb$ensembl_id),"symbol"]
      genes_versioned <- GeneNameToVersionedID(genes_symb)
      bulk <- get_eQTL_bulk(genesx = genes_versioned ,lead = SelectedSNP ,tissues = tissue)
      bulk
    })

    output$eQTL_check <- DT::renderDT({
      bulk <- moduleQTLdata()
      bulk <- cbind(bulk[,c("snpId","geneSymbol")],signif(bulk[,c("pValue","pValueThreshold")],2))
      bulk
    },selection = 'single')

    output$module_Violin <- conquer.d3js::renderConquerViolin({
      row_selected <- shiny::req(input$eQTL_check_rows_selected)
      data <- moduleQTLdata()
      SNP <- data[row_selected,"snpId"]
      gene <- strsplit(data[row_selected,"gencodeId"],"[.]")[[1]][1]
      tissue <- shiny::req(input$tissueSel)
      conquer.d3js::ConquerViolin(SNP = SNP, gene = gene,tissue = tissue)
    })



    ####ALL TISSUES####
    output$EdgePlot <- conquer.d3js::renderConquerEdge({
      conquer.d3js::ConquerEdge(SNPSummary)
    })

    LDSNPs <- shiny::reactive({
      LDSNPs <- lapply(names(SNPs),function(SNP){
        LDtable <- SNPs[[SNP]]$LD
        LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
        as.character(LDSNPs)
      })
      names(LDSNPs) <- names(SNPs)
      LDSNPs <- stack(LDSNPs)
      colnames(LDSNPs) <- c("LDSNP", "leadingSNP")
      return(LDSNPs)
      })

    ###pQTLs
    output$pQTLOverview <- DT::renderDT({
      AllTissuespQTLsData()
    },options=list(scrollX=T),selection = "single")

    AllTissuespQTLsData <- shiny::reactive({
      LDSNPs <- LDSNPs()
      all <- unique(c(LDSNPs$LDSNP, LDSNPs$leadingSNP))
      output <- pQTLs[pQTLs$rsID %in% all,]
      return(output)
    })

    output$pQTLOverview_LD <- DT::renderDT({
      row_selected <- shiny::req(input$pQTLOverview_rows_selected)
      LDSNPs <- LDSNPs()
      pQTLs <- AllTissuespQTLsData()
      selectedSNP <- pQTLs[row_selected,"rsID"]
      leadingSNP <- LDSNPs[LDSNPs$LDSNP == selectedSNP,"leadingSNP"]
      fullLD <- SNPs[[leadingSNP]]$LD
      fullLD <- fullLD[fullLD$variation == selectedSNP, ]
      fullLD <- cbind(leadingSNP, fullLD)
      return(fullLD)
    })

    ###meQTLs
    AllTissuesMeQTLsData <- shiny::reactive({
      LDSNPs <- LDSNPs()
      meQTLs <- conquer.db::meQTLs
      all <- unique(c(LDSNPs$LDSNP, LDSNPs$leadingSNP))
      output <- meQTLs[meQTLs$rsID %in% all, ]
      return(output)
    })

    output$meQTLOverview <- DT::renderDT({
      AllTissuesMeQTLsData()
    },options=list(scrollX=T),selection = "single")

    output$meQTLOverview_LD <- DT::renderDT({
      row_selected <- shiny::req(input$meQTLOverview_rows_selected)
      LDSNPs <- LDSNPs()
      meQTLs <- AllTissuesMeQTLsData()
      selectedSNP <- meQTLs[row_selected,"rsID"]
      leadingSNP <- LDSNPs[LDSNPs$LDSNP == selectedSNP,"leadingSNP"]
      fullLD <- SNPs[[leadingSNP]]$LD
      fullLD <- fullLD[fullLD$variation == selectedSNP, ]
      fullLD <- cbind(leadingSNP, fullLD)
      return(fullLD)
    })

    ###miQTLS
    #Experimental
    AllTissuesMiQTLexperiment <-  shiny::reactive({
      LDSNPs <- LDSNPs()
      all <- unique(c(LDSNPs$LDSNP, LDSNPs$leadingSNP))
      miQTLs <- conquer.db::miQTLexperiment
      miQTLs <- miQTLs[miQTLs$SNP %in% all,]
      return(miQTLs)
    })

    output$mi_overview_exp <- DT::renderDT({
      AllTissuesMiQTLexperiment()
    },options=list(scrollX=T),selection = "single")

    output$mi_overview_exp_LD <- DT::renderDT({
      row_selected <- shiny::req(input$mi_overview_exp_rows_selected)
      LDSNPs <- LDSNPs()
      miQTLs <- AllTissuesMiQTLexperiment()
      selectedSNP <- miQTLs[row_selected,"SNP"]
      selectedSNP <- c(selectedSNP) %>% unlist()
      leadingSNP <- LDSNPs[LDSNPs$LDSNP == selectedSNP | LDSNPs$leadingSNP == selectedSNP ,"leadingSNP"]
      fullLD <- SNPs[[leadingSNP]]$LD
      fullLD <- fullLD[fullLD$variation == selectedSNP, ]
      fullLD <- cbind(leadingSNP, fullLD)
      return(fullLD)
    },options=list(scrollX=T),selection = "single")

    #Predicted
    AllTissuesMiQTLpredict <- shiny::reactive({
      LDSNPs <- LDSNPs()
      all <- unique(c(LDSNPs$LDSNP, LDSNPs$leadingSNP))
      miQTLs <- conquer.db::miQTLpredict
      miQTLs <- miQTLs[miQTLs$SNP %in% all,]
      return(miQTLs)
    })

    output$mi_overview_pred <- DT::renderDT({
      AllTissuesMiQTLpredict()[,c("SNP","Gene", "miRNA", "Celltype","Change","Effect","RefAllele")]
    },options=list(scrollX=T),selection = "single")


    output$mi_overview_pred_LD <- DT::renderDT({
      row_selected <- shiny::req(input$mi_overview_pred_rows_selected)
      print(row_selected)
      LDSNPs <- LDSNPs()
      miQTLs <- AllTissuesMiQTLpredict()
      selectedSNP <- miQTLs[row_selected,"SNP"]
      selectedSNP <- c(selectedSNP) %>% unlist() %>% unname() %>% as.character()
      print(selectedSNP)
      leadingSNP <- LDSNPs[LDSNPs$LDSNP == selectedSNP | LDSNPs$leadingSNP == selectedSNP ,"leadingSNP"]
      fullLD <- SNPs[[leadingSNP]]$LD
      fullLD <- fullLD[fullLD$variation == selectedSNP, ]
      fullLD <- cbind(leadingSNP, fullLD)
      return(fullLD)
    })

    ####SINGLE SNP####
    output$responsiveUI_A <- shiny::renderUI({

      if(input$single == "QTLs"){
        checkboxGroupInput("check_eQTL", label = "Type of eQTLs: ",
                           choices = list("cis-eQTLs" = 1, "trans-eQTLs" = 2),
                           selected = c(1,2))
      }else if(input$single == "GeX"){
        checkboxGroupInput("check_eQTL", label = "Type of eQTLs: ",
                           choices = list("cis-eQTLs" = 1, "trans-eQTLs" = 2),
                           selected = c(1,2))
      }
    })

    output$responsiveUI_B <- shiny::renderUI({
      if(input$single == "GeX"){
        checkboxInput("log10Trans", "Transform Log10", value = FALSE, width = NULL)
      }else if (input$single == "CI"){
        groups <- unique(conquer.db::ChromatinGroups$States$group)
        shiny::selectInput(inputId = "chromTissue",
                           label = "Select Tissue:",
                           choices = sort(groups))
      }else if(input$single == "CS" ){
        groups <- unique(conquer.db::ChromatinGroups$States$group)
        shiny::selectInput(inputId = "chromTissue",
                           label = "Select Tissue:",
                           choices = c("All", sort(groups)))
      }
    })


    # Load selected SNP into the environment
    output$LocusZoom <- conquer.d3js::renderConquerLocusZoom({
      conquer.d3js::ConquerLocusZoom(SNPs[[input$snpSel]])
    })
    output$LocusHeader <- shiny::renderUI({
      shiny::tags$h2(sprintf("Linkage Disequilibrium for: %s",input$snpSel))
    })

    output$downloadLocus <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_LocusZoom", ".html", sep = "")
      },
      content = function(file) {
        htmlwidgets::saveWidget(widget = conquer.d3js::ConquerLocusZoom(SNPs[[input$snpSel]], width = 500), file = file)
      })




    ########################LD Table########################
    #Data
    LDTable <- reactive({
      SelectedSNP <- shiny::req(input$snpSel)
      data <- SNPs[[input$snpSel]]$LD
      data <- data[order(data$r2,decreasing = T),]
      data
    })

    #Table
    output$LDTable <- DT::renderDT({
      data <- LDTable()
      data <- data[order(data$r2,decreasing = T),c("variation","start","consequence_type","r2","clinical_significance")]
      data
    },selection = "single")

    #Download
    output$downloadLD <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_Linkage_Disequilibrium", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(LDTable(), file)
      }
    )
    ########################END########################



    output$PMIDTable <- DT::renderDT({
      shiny::req(input$LDTable_rows_selected)
      data <- SNPs[[input$snpSel]]$LD
      data <- data[order(data$r2,decreasing = T),c("variation","start")]
      sel <- data[input$LDTable_rows_selected,"variation"]
      json <- jsonlite::fromJSON(sprintf("https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/%s",gsub("rs","",sel)))
      pmidData <- data.frame("PMID"= paste0("<a target='_blank' href='",
                                     sprintf("https://www.ncbi.nlm.nih.gov/pubmed/%s",
                                             json$citations),
                                     "'>",
                                     json$citations,"</a>"))

      pmidData
    },escape = FALSE,selection = "single")


    output$PMIDHeader <- shiny::renderUI({
      shiny::req(input$LDTable_rows_selected)
      data <- SNPs[[input$snpSel]]$LD
      data <- data[order(data$r2,decreasing = T),c("variation","start")]
      sel <- data[input$LDTable_rows_selected,"variation"]
      shiny::tags$h3(sprintf("Publications in which %s is mentioned: ", sel))
    })



    ########################eQTL Table########################
    #Data
    eQTLsData <- shiny::reactive({
      if(nrow(SNPs[[input$snpSel]]$eQTLs) == 0){
        cis <- SNPs[[input$snpSel]]$eQTLs
      }else{
        cis <- SNPs[[input$snpSel]]$eQTLs
        cis$type <- "cis"
      }
      if(nrow(SNPs[[input$snpSel]]$eQTLsTrans) == 0){
        trans <- SNPs[[input$snpSel]]$eQTLsTrans
      }else{
        trans <- SNPs[[input$snpSel]]$eQTLsTrans
        trans$type <- "trans"
      }
      shiny::req(input$check_eQTL)

      if(length(input$check_eQTL) == 2){
        data <- rbind(cis,trans)
      } else if(input$check_eQTL == 1 ){
        data <- cis
      } else if (input$check_eQTL == 2) {
        data <- trans
      }
      data <- data[order(data$pValue,decreasing = F),]
      data
    })

    #Table
    output$eQTLsTable <- DT::renderDT({
      data <- eQTLsData()
      data <- cbind(data[,c("SNP","gene","tissue")],signif(data[,c("pValue","Pval.ratio")],2))
    },options=list(scrollX=T),selection = "single")

    #Download
    output$downloadeQTLs <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_eQTLs", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(eQTLsData(), file)
      }
    )
    ########################END########################
    output$pQTLsTable <- DT::renderDT({
      LDtable <- SNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewpQTLs <- pQTLs[pQTLs$rsID %in% LDSNPs,]
      if(nrow(ViewpQTLs) == 0){

      }else{
        ViewpQTLs
      }
    },options=list(scrollX=T),selection = "single")



    output$meQTLsTable <- DT::renderDT({
      LDtable <- SNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewmeQTLs <- meQTLs[meQTLs$rsID %in% LDSNPs,]
      if(nrow(ViewmeQTLs) == 0){

      }else{
        ViewmeQTLs
      }
    },options=list(scrollX=T),selection = "single")



    output$miQTLsPred_table <- DT::renderDT({
      LDtable <- SNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewmiQTLsPred <- miQTLpred[miQTLpred$SNP %in% LDSNPs, c("SNP","Gene","miRNA","Celltype","Change","Effect")]
      if(is.null(ViewmiQTLsPred)){

      }else{
        ViewmiQTLsPred
      }
    })

    output$miQTLsExp_table <- DT::renderDT({
      LDtable <- SNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewmiQTLsExp <- miQTLex[miQTLex$SNP %in% LDSNPs,]

      if(nrow(ViewmiQTLsExp) == 0){

      }else{
        ViewmiQTLsExp
      }
    },options=list(scrollX=T),selection = "single")


    output$CIMessage <- shiny::renderUI({
      SNP <- shiny::req(input$snpSel)
      tissue <- shiny::req(input$chromTissue)
      if(is.null(SNPs[[SNP]]$chromInt)){
        shiny::h3(sprintf("There are no known interactions around %s", SNP))
      }else{
        NULL
      }
    })

    output$Circos <- BioCircos::renderBioCircos({
      SNP <- shiny::req(input$snpSel)
      tissue <- shiny::req(input$chromTissue)
      if(is.null(SNPs[[SNP]]$chromInt)){
        return(NULL)
      }else{
        ConquerCircos(SNPs[[SNP]],tissue = tissue)
      }
    })


    ########################ChromatinStates########################
    chromatinStatesData <- shiny::reactive({
      ChromatinGroups <- conquer.db::ChromatinGroups$States
      chromSelect <- req(input$chromTissue)
      if(chromSelect == "All"){
        statesSel <- NULL
      }else{
        statesSel <- ChromatinGroups[ChromatinGroups$group == chromSelect, "name"]
      }
      output <- PrepareChromatinStates(SNPs[[input$snpSel]],statesSel)
      return(output)
    })
    output$chromatinStates <- plotly::renderPlotly({
      data <- chromatinStatesData()
      PlotChromatinStates(StatesPlotData = data[[1]], fillScale = data[[2]],jit = data[[3]])
    })
    output$placeHolder_downloadData <- shiny::renderUI({
      data <- chromatinStatesData()
      if(!is.null(data)){
        shiny::downloadButton("chromatinStates_downloadData", "Data")
      }
    })
    output$placeHolder_downloadPlot <- shiny::renderUI({
      data <- chromatinStatesData()
      if(!is.null(data)){
        shiny::downloadButton("chromatinStates_downloadPlot", "Plot")
      }
    })

    #Download
    output$chromatinStates_downloadData <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_chromatinStates", ".xlsx", sep = "")
      },
      content = function(file) {
        data <- chromatinStatesData()[[1]]
        data <- data[,c("seqnames","start","end","width","strand","target","sample","group")]
        rio::export(data, file)
      }
    )

    output$chromatinStates_downloadPlot <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_chromatinStates", ".html", sep = "")
      },
      content = function(file) {
        data <- chromatinStatesData()
        htmlwidgets::saveWidget(PlotChromatinStates(StatesPlotData = data[[1]], fillScale = data[[2]],jit = data[[3]]), file = file)
      }
    )
    ########################END########################
    ########################VIOLIN########################
    #Data
    violinPlotData <- reactive({
      shiny::req(input$eQTLsTable_rows_selected)
      SNP <- input$snpSel
      tmp_eQTL <- eQTLsData()
      tissue <-   tmp_eQTL[input$eQTLsTable_rows_selected,"tissue"]
      gene <- tmp_eQTL[input$eQTLsTable_rows_selected,"gencodeId"]
      gene <- strsplit(gene,"[.]")[[1]][1]
      list("SNP" = SNP, "gene" = gene, "tissue" = tissue)
    })
    #Plot
    output$cis_Violin <- conquer.d3js::renderConquerViolin({
      shiny::req(input$eQTLsTable_rows_selected)
      data <- violinPlotData()
      conquer.d3js::ConquerViolin(SNP = data[["SNP"]], gene = data[["gene"]],tissue = data[["tissue"]])
    })
    #Download button
    output$DownloadViolinButton <- shiny::renderUI({
      shiny::req(input$eQTLsTable_rows_selected)
      shiny::downloadButton("downloadVioloinPlot", "Download Violin")
    })
    #Actual download
    output$downloadVioloinPlot <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_",violinPlotData()[["gene"]],"_",violinPlotData()[["tissue"]],"_violin", ".html", sep = "")
      },
      content = function(file) {
        htmlwidgets::saveWidget(conquer.d3js::ConquerViolin(SNP = violinPlotData()[["SNP"]], gene = violinPlotData()[["gene"]],tissue = violinPlotData()[["tissue"]]), file = file)
      })
    ########################END########################



    ########################Hive Plot########################
    output$hive_plot <- conquer.d3js::renderConquerHive({
      selec <- shiny::req(input$snpSel)
      cis_trans <- shiny::req(input$check_eQTL)
      if(length(cis_trans) == 2){
        cis <- T
        trans <- T
      } else if(cis_trans == 1 ){
        cis <- T
        trans <- F
      }else if (cis_trans == 2) {
        cis <- F
        trans <- T
      }else{
        cis <- F
        trans <- F
      }
      data <- SNPs[[selec]]

      if(cis & nrow(data$eQTLs) == 0){
        cis <- F
      }
      if(trans & nrow(data$eQTLsTrans) == 0){
        trans <- F
      }

      if(!cis & !trans){
        return(NULL)
      }else{
        conquer.d3js::ConquerHive(data,cis,trans)
      }


    })

    output$downloadHive <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_Hiveplot", ".html", sep = "")
      },
      content = function(file) {
        htmlwidgets::saveWidget(
          widget = conquer.d3js::ConquerHive(SNPData = SNPs[[shiny::req(input$snpSel)]],
                                             cis = T, trans = T, width = 900, height = 700),
          file = file)
      })
    ########################END########################

    #"geneExpr_downloadData"
    #"geneExpr_downloadPlot"

    geneExpressionData <- shiny::reactive({
      check_eQTL <- shiny::req(input$check_eQTL)
      SNP <- shiny::req(input$snpSel)
      transform <- FALSE
      transform <- input$log10Trans
      cisGenes <- SNPs[[SNP]]$eQTLs$gencodeId %>% unique()
      transGenes <- SNPs[[SNP]]$eQTLsTrans$gencodeId %>% unique()
      if(length(check_eQTL) == 2){
        genes <- c(cisGenes, transGenes)
      } else if(check_eQTL == 1 ){
        genes <- cisGenes
      }else if (check_eQTL == 2) {
        genes <- transGenes
      }
      geneExpr<- SNPs[[SNP]]$geneExpr
      geneExpr <- geneExpr[geneExpr$gencodeId %in% genes,]
      geneExpr<- reshape2::dcast(geneExpr, formula = geneSymbol ~ tissueSiteDetailId,value.var = "median")
      rownames(geneExpr) <- geneExpr$geneSymbol
      geneExpr$geneSymbol <- NULL
      if(transform){geneExpr <- log10(geneExpr+1)}
      return(geneExpr)
    })

    output$expressionHeatmap <- plotly::renderPlotly({
      geneExpr <- geneExpressionData()
      plotly::plot_ly(
        x = rownames(geneExpr),y = colnames(geneExpr), colors = colorRampPalette(c("blue","white","red"))(100),
        z = t(as.matrix(geneExpr)), type = "heatmap") %>% plotly::layout(title ="Median Expression (TPM)")
    })

    output$geneExpr_downloadData <- shiny::renderUI({
      shiny::downloadButton("geneExpr_Data", "Data")
    })
    output$geneExpr_downloadPlot <- shiny::renderUI({
      shiny::downloadButton("geneExpr_Plot", "Plot")
    })

    output$geneExpr_Data <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_geneExpression", ".xlsx", sep = "")
      },
      content = function(file) {
        data <- geneExpressionData()
        genes <- rownames(data)
        data <-  cbind(genes,data)
        rio::export(data, file)
      }
    )

    output$geneExpr_Plot <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_geneExpression", ".html", sep = "")
      },
      content = function(file) {
        data <- geneExpressionData()
        plot <- plotly::plot_ly(
          x = rownames(data),y = colnames(data), colors = colorRampPalette(c("blue","white","red"))(100),
          z = t(as.matrix(data)), type = "heatmap") %>% plotly::layout(title ="Median Expression (TPM)")
        htmlwidgets::saveWidget(plot, file = file)
      }
    )


  }
  return(shiny::shinyApp(ui = ui, server = server))

}



