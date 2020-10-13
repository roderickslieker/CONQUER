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
#' @importFrom shinyjs useShinyjs toggle
#' @import shinycssloaders
#' @return [[NULL]]
visualizeDashboard <- function(loadedSNPs, SNPSummary, ColocSummary, tissues=NULL){
  currversion.db <- rio::import("https://raw.githubusercontent.com/roderickslieker/CONQUER.db/master/DESCRIPTION", nrow=1, skip=3, format='\t', header=F)[1,2]
  currversion.d3 <- rio::import("https://raw.githubusercontent.com/roderickslieker/CONQUER.d3/master/DESCRIPTION", nrow=1, skip=3, format='\t', header=F)[1,2]
  currversion <- rio::import("https://raw.githubusercontent.com/roderickslieker/CONQUER/master/DESCRIPTION", nrow=1, skip=3, format='\t', header=F)[1,2]

  if(packageVersion("conquer.db") == currversion.db)
  {
  }else{
    stop("You are not using the latest version of CONQUER.db, please update from GitHub to prevent broken links (devtools::install.github('roderickslieker/conquer.db')")
  }

  if(packageVersion("CONQUER") == currversion)
  {
  }else{
    stop("You are not using the latest version of CONQUER, please update from GitHub to prevent broken links (devtools::install.github('roderickslieker/CONQUER')")
  }

  if(packageVersion("conquer.d3js") == currversion.d3)
  {
  }else{
    stop("You are not using the latest version of CONQUER.d3js, please update from GitHub to prevent broken links (devtools::install.github('roderickslieker/conquer.d3')")
  }
  #KEGG_DATA <- CONQUER:::prepare_KEGG(species = "hsa",
  #                                KEGG_Type = "KEGG",
  #                                keyType = "kegg")
  #KEGG_DATA$ENSG <- lapply(X = KEGG_DATA$PATHID2EXTID,
  #                         FUN = entrezToENSEMBL)
  ensg_symb <- merge(BiocGenerics::toTable(org.Hs.egENSEMBL),BiocGenerics::toTable(org.Hs.egSYMBOL),by="gene_id")

  #Logo
  shiny::addResourcePath("logo", directoryPath = system.file("logo", package = "CONQUER"))

  #LD
  ld.snps <- getAllSNPLD(SNPs = loadedSNPs)

  cat("Loading data.....","\n")
  qtls <- c("pqtls","meQTLs","miQTLexperiment","miQTLpredict","mqtls_LC",
            "mqtls_NG","sqtls1","sqtls2","sqtls3","sqtls4","lqtls")


  for(qtl in qtls)
  {
    cat(sprintf("..Loading %s.....",qtl),"\n")
    temp <- get(qtl)
    if(length(grep("sqtl", qtl)) != 0)
    {
      temp <- temp[temp$variant_id %in% ld.snps$id,]
      temp$LeadSNP <- ld.snps[match(temp$variant_id, ld.snps$id),"leadingSNP"]
    }else if(length(grep("miQTL", qtl)) !=0){
      temp <- temp[temp$SNP %in% ld.snps$LDSNP,]
      temp$LeadSNP <- ld.snps[match(temp$SNP, ld.snps$LDSNP),"leadingSNP"]
    }else{
      temp <- temp[temp$rsID %in% ld.snps$LDSNP,]
      temp$LeadSNP <- ld.snps[match(temp$rsID, ld.snps$LDSNP),"leadingSNP"]
    }
    qtl <- paste0(qtl, "_internal")
    assign(qtl, temp, envir = baseenv())
    rm(temp)
  }

  sqtls_internal <- rbind(sqtls1_internal,sqtls2_internal,sqtls3_internal,sqtls4_internal)

  cat(sprintf("Starting dashboard.....",qtl))


  buttonStyle <-  "display:block;
      height: 32px;
      width: 32px;
      padding-top:0px;
      padding-left:0px;
      padding-bottom:0px;
      background-color: #A4A4A4;
      padding-right:0px;
      color: white;
      font-size: 20px;
      vertical-align: middle;
      border-radius: 50%;
      border: -px solid #A4A4A4"

  ui <- shiny::navbarPage(title = shiny::div(shiny::img(src = "logo/CONQUER.png", style="margin-top:-10px;")),
                          shiny::tags$head(shiny::HTML("<title>test</title>")),
                          windowTitle = shiny::HTML("CONQUER"),
                          selected = "Modules",
                          theme = shinythemes::shinytheme("flatly"),
                          shiny::tabPanel("Modules",
                                          shiny::fluidPage(
                                            shiny::sidebarLayout(
                                              shiny::sidebarPanel(
                                                width = 3,
                                                shiny::selectInput(inputId = "tissueSel",
                                                                   label = "Select Tissue:",
                                                                   choices = colnames(SNPSummary) %>% sort()),
                                                shiny::uiOutput("responsiveUI_C"),
                                                shinyjs::useShinyjs(),
                                                shiny::tags$style(shiny::HTML("
                                                          .btn {
                                                            display:block;
                                                            height: 50px;
                                                            width: 50px;
                                                            fill:red;
                                                            padding-top:0px;
                                                            padding-left:0px;
                                                            padding-bottom:0px;
                                                            padding-right:0px;
                                                            color:#A4A4A4;
                                                            border-radius: 50%;
                                                            border: 1px solid #ECF0F1;
                                                            font-size: 30px;


                                                            }

                                                            ")),
                                                shiny::actionButton("button", "", icon=shiny::icon("info-circle"), style = "background-color: #ECF0F1"),
                                                shinyjs::hidden(
                                                  shiny::div(id='text_div',
                                                             shiny::htmlOutput("text")
                                                  )
                                                ),
                                              ),
                                              shiny::mainPanel(
                                                shiny::tabsetPanel(id = "tis_spec",
                                                                   shiny::tabPanel("Overview", value = "overview",
                                                                                   shiny::br(),
                                                                                   shiny::br(),
                                                                                   shiny::fluidRow(
                                                                                     shiny::column(12,align='center',
                                                                                                   conquer.d3js::ConquerRingOutput("RingPlot",
                                                                                                                                   width = 900,
                                                                                                                                   height = 900) %>% withSpinner(type=7)
                                                                                     )
                                                                                   )

                                                                   ),
                                                                   shiny::tabPanel("Modules",value = "mod",
                                                                                   shiny::br(),
                                                                                   shiny::br(),
                                                                                   shiny::fluidRow(shiny::column(6,
                                                                                                                 plotly::plotlyOutput("moduleHeat")%>% withSpinner(type=7)
                                                                                   ),
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
                                                                   ),
                                                                   shiny::tabPanel("Advanced Enrichment",value = "advenr",
                                                                                   shiny::br(),
                                                                                   shiny::br(),
                                                                                   shiny::fluidRow(shiny::column(12,
                                                                                                                 plotly::plotlyOutput("moduleEnrichment")%>% withSpinner(type=7)
                                                                                   )
                                                                                   )
                                                                   )
                                                )
                                              )
                                            )

                                          )),
                          shiny::tabPanel("All SNPs",
                                          shiny::fluidPage(
                                            shiny::sidebarLayout(
                                              shiny::sidebarPanel(
                                                width = 3,

                                                shiny::actionButton("button2", "", icon=shiny::icon("info-circle"), style = "background-color: #ECF0F1"),
                                                shinyjs::hidden(
                                                  shiny::div(id='text_div2',
                                                             shiny::htmlOutput("text2")
                                                  )
                                                )
                                              ),
                                              shiny::mainPanel(
                                                shiny::tabsetPanel(id = "pws",
                                                                   shiny::tabPanel("KEGG Pathways", value = "pway",
                                                                                   shiny::fluidRow(
                                                                                     if(!is.null(SNPSummary))
                                                                                     {
                                                                                       shiny::column(12,align = 'center',
                                                                                                     conquer.d3js::ConquerEdgeOutput("EdgePlot",
                                                                                                                                     width = 1100,
                                                                                                                                     height = 1100)%>% withSpinner(type=7))
                                                                                     }else{
                                                                                       shiny::h3("No summary file provided")
                                                                                     }
                                                                                   )),
                                                                   shiny::tabPanel("KEGG Disease", value = "pway",
                                                                                   shiny::fluidRow(
                                                                                     if(!is.null(SNPSummary))
                                                                                     {
                                                                                       shiny::column(12,align = 'center',
                                                                                                     conquer.d3js::ConquerEdgeOutput("EdgePlotDisease",
                                                                                                                                     width = 1100,
                                                                                                                                     height = 1100))
                                                                                     }else{
                                                                                       shiny::h3("No summary file provided")
                                                                                     }
                                                                                   )),
                                                                   shiny::tabPanel("DNAm QTL", value = "meqtl",
                                                                                   shiny::br(),
                                                                                   customDownloadbutton("downloadallmeqtls", "", icon="file-excel", style = buttonStyle,
                                                                                                        class="btn btn-default htmlwidgets-download-link"),
                                                                                   DT::DTOutput("meQTLOverview"),
                                                                                   shiny::br(),
                                                                                   shiny::br(),
                                                                                   DT::DTOutput("meQTLOverview_LD")),
                                                                   shiny::tabPanel("miRNA QTL", value = "miqtl",
                                                                                   shiny::tabsetPanel(
                                                                                     shiny::tabPanel("Experimental",
                                                                                                     shiny::br(),
				                                                                                     customDownloadbutton("downloadallmiqtlsexp", "", icon="file-excel", style = buttonStyle,
				                                                                                                        class="btn btn-default shiny-download-link"),
				                                                                                                     shiny::br(),
                                                                                                     DT::DTOutput("mi_overview_exp"),
                                                                                                     shiny::br(),
                                                                                                     shiny::br(),
                                                                                                     DT::DTOutput("mi_overview_exp_LD")),
                                                                                     shiny::tabPanel("Predicted",
                                                                                                     shiny::br(),
                                                                                                     customDownloadbutton("downloadallmiqtlspred", "", icon="file-excel", style = buttonStyle,
				                                                                                                        class="btn btn-default shiny-download-link"),

                                                                                                     shiny::br(),
                                                                                                     DT::DTOutput("mi_overview_pred"),
                                                                                                     shiny::br(),
                                                                                                     shiny::br(),
                                                                                                     DT::DTOutput("mi_overview_pred_LD"))
                                                                                   )
                                                                   ),
                                                                   shiny::tabPanel("Protein QTL",value = "pqtl",
                                                                                   shiny::br(),
                                                                                   customDownloadbutton("downloadallpqtls", "", icon="file-excel", style = buttonStyle,
                                        				                                                        class="btn btn-default shiny-download-link"),
                                                                                   shiny::br(),
                                                                                   DT::DTOutput("pQTLOverview"),
                                                                                   shiny::br(),
                                                                                   shiny::br(),
                                                                                   DT::DTOutput("pQTLOverview_LD")),
                                                                   shiny::tabPanel("Splicing QTL",value = "sqtl",
                                                                                   shiny::br(),
																				   customDownloadbutton("downloadallsqtls", "", icon="file-excel", style = buttonStyle,
																											class="btn btn-default shiny-download-link"),
                                                                                   shiny::br(),
                                                                                   DT::DTOutput("sQTLOverview")),
                                                                   shiny::tabPanel("Lipid QTL",value = "lqtl",
                                                                                   shiny::br(),
                                                                                   customDownloadbutton("downloadalllqtls", "", icon="file-excel", style = buttonStyle,
                                                                                                        class="btn btn-default shiny-download-link"),
                                                                                   shiny::br(),
                                                                                   DT::DTOutput("lQTLOverview")),
                                                                   shiny::tabPanel("Metabolite QTL",value = "mqtl",
                                                                                   shiny::tabsetPanel(
                                                                                     shiny::tabPanel("Nightingale",
                                                                                                     shiny::br(),
                                                                                                     customDownloadbutton("downloadallmqtlsnigh", "", icon="file-excel", style = buttonStyle,
                                                                                                        class="btn btn-default shiny-download-link"),
                                                                                                     shiny::br(),
                                                                                                     DT::DTOutput("mqtls_overview_ng")),
                                                                                     shiny::tabPanel("Multiplatform",
                                                                                                     shiny::br(),
                                                                                                     customDownloadbutton("downloadallmqtlsmulti", "", icon="file-excel", style = buttonStyle,
                                                                                                        class="btn btn-default shiny-download-link"),
                                                                                                     shiny::br(),
                                                                                                     DT::DTOutput("mqtls_overview_lc"))
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
                                                                   choices = names(loadedSNPs)),
                                                shiny::uiOutput("responsiveUI_A"),
                                                shiny::uiOutput("responsiveUI_B"),
                                                shiny::uiOutput("responsiveUI_D"),
                                                shiny::uiOutput("responsiveUI_E"),
                                                shiny::actionButton("button3", "", icon=shiny::icon("info-circle"), style = "background-color: #ECF0F1",),
                                                shinyjs::hidden(
                                                  shiny::div(id='text_div3',
                                                             shiny::htmlOutput("text3")
                                                  )
                                                )),
                                              shiny::mainPanel(
                                                shiny::tabsetPanel(id = "single",
                                                                   shiny::tabPanel(
                                                                     title = "Linkage Disequilibrium", value="LD",
                                                                     shiny::br(),
                                                                     customDownloadbutton("downloadLocus", "", icon="cloud-download",
                                                                                          style = buttonStyle,
                                                                                          class="btn btn-default shiny-download-link"),
                                                                     shiny::fluidRow(
                                                                       shiny::column(12, align = 'center',
                                                                                     shiny::uiOutput("LocusHeader"),
                                                                                     conquer.d3js::ConquerLocusZoomOutput("LocusZoom",
                                                                                                                          width = 500,
                                                                                                                          height = 600)
                                                                       )
                                                                     ),
                                                                     shiny::br(),
                                                                     customDownloadbutton("downloadLD", "", icon="file-excel", style = buttonStyle,
                                                                                          class="btn btn-default shiny-download-link"),
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
                                                                                   shiny::uiOutput("downloadCPlot"),
                                                                                   BioCircos::BioCircosOutput(outputId = "Circos",width = 1000,height = 1000)),
                                                                   shiny::tabPanel(title = "Chromatin States", value = "CS",
                                                                                   shiny::br(),
                                                                                   shiny::fluidRow(
                                                                                     shiny::column(1,shiny::uiOutput("placeHolder_downloadData")),
                                                                                     shiny::column(1,shiny::uiOutput("placeHolder_downloadPlot"))
                                                                                   ),
                                                                                   shiny::br(),
                                                                                   plotly::plotlyOutput("chromatinStates",height = 1100) %>% withSpinner(type=7)),
                                                                   shiny::tabPanel(title="QTLs", value = "QTLs",
                                                                                   shiny::tabsetPanel(id="QTLs_tab",
                                                                                                      shiny::tabPanel(title="eQTLs",value = "eqtls_sub",
                                                                                                                      shiny::fluidRow(shiny::tags$h1("Hive plot"),
                                                                                                                                      customDownloadbutton("downloadHive", "", icon="cloud-download",
                                                                                                                                                           style = buttonStyle,
                                                                                                                                                           class="btn btn-default shiny-download-link"),
                                                                                                                                      shiny::br(),
                                                                                                                                      shiny::br(),
                                                                                                                                      conquer.d3js::ConquerHiveOutput("hive_plot",
                                                                                                                                                                      width = 900,
                                                                                                                                                                      height = 700) %>% withSpinner(type=7)
                                                                                                                      ),
                                                                                                                      shiny::fluidRow(
                                                                                                                        shiny::tags$h1("eQTLs"),
                                                                                                                        customDownloadbutton("downloadeQTLs", "", icon="cloud-download",
                                                                                                                                             style = buttonStyle,
                                                                                                                                             class="btn btn-default shiny-download-link"),
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
                                                                                                                                      conquer.d3js::ConquerViolinOutput("cis_Violin") %>% withSpinner(type=7)
                                                                                                                        ))),
                                                                                                      shiny::tabPanel(title="Colocalization", value = "coloc_sub",
                                                                                                                      shiny::tags$h3("Bayes Factor colocalization analysis"),
                                                                                                                      shiny::br(),
                                                                                                                      shiny::br(),
                                                                                                                      shiny::column(1,shiny::uiOutput("placeHolder_downloadPlot_coloc")),
                                                                                                                      shiny::br(),
                                                                                                                      shiny::br(),
                                                                                                                      shiny::fluidRow(shiny::column(width = 12,
                                                                                                                                                    withSpinner(plotly::plotlyOutput("moduleColoc", width = "650px", height="1000px"), type=7))
                                                                                                                      )),
                                                                                                      shiny::tabPanel(title="sQTLs", value = 'sQTLs_sub',
                                                                                                                      shiny::column(12,
                                                                                                                                    shiny::br(),
                                                                                                                                    DT::DTOutput("sQTLsExp_table"))),
                                                                                                      shiny::tabPanel(title="miQTLs", value = "miqtl_sub",
                                                                                                                      shiny::fluidRow(
                                                                                                                        shiny::column(6,
                                                                                                                                      shiny::h4("Predicted miQTLs"),
                                                                                                                                      shiny::br(),
                                                                                                                                      DT::DTOutput("miQTLsPred_table")),
                                                                                                                        shiny::column(6,
                                                                                                                                      shiny::h4("Experimental determined miQTLs"),
                                                                                                                                      shiny::br(),
                                                                                                                                      DT::DTOutput("miQTLsExp_table")))),
                                                                                                      shiny::tabPanel(title="meQTLs", value = 'meqtl_sub',
                                                                                                                      shiny::column(12,
                                                                                                                                    shiny::br(),
                                                                                                                                    DT::DTOutput("meQTLsTable"))),
                                                                                                      shiny::tabPanel(title="pQTLs", value = "pqtls_sub",
                                                                                                                      shiny::column(12,
                                                                                                                                    shiny::br(),
                                                                                                                                    DT::DTOutput("pqtlsTable"))),
                                                                                                      shiny::tabPanel(title="lQTLs", value = 'lqtl_sub',
                                                                                                                      shiny::column(12,
                                                                                                                                    shiny::br(),
                                                                                                                                    DT::DTOutput("lQTLsExp_table"))),
                                                                                                      shiny::tabPanel(title="mQTLs (NG)", value = 'mqtlng_sub',
                                                                                                                      shiny::column(12,
                                                                                                                                    shiny::br(),
                                                                                                                                    DT::DTOutput("mQTLsNGExp_table"))),
                                                                                                      shiny::tabPanel(title="mQTLs (multi)", value = 'mqtllc_sub',
                                                                                                                      shiny::column(12,
                                                                                                                                    shiny::br(),
                                                                                                                                    DT::DTOutput("mQTLsLCExp_table")))
                                                                                   )),
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

                                          )
                          ),
                          shiny::tabPanel("Colocalization",
                                          shiny::fluidPage(
                                            shiny::sidebarLayout(
                                              shiny::sidebarPanel(
                                                width = 3,
                                                shiny::selectInput(inputId = "snpc",
                                                                   label = "Select SNP:",
                                                                   choices = names(loadedSNPs) %>% sort()),
                                                shiny::textInput(inputId = "gencodec",
                                                                 label = "Provide versioned gencodeId (gene.1)"),
                                                shiny::selectInput(inputId = "tissueSelc",
                                                                   label = "Select Tissue:",
                                                                   choices = gtexTissuesV8 %>% sort()),
                                                shiny::actionButton("run", "Go"),
                                                shiny::actionButton("button", "", icon=shiny::icon("info-circle"), style = "background-color: #ECF0F1"),
                                                shinyjs::hidden(
                                                  shiny::div(id='text_div',
                                                             shiny::htmlOutput("text6")
                                                  )
                                                ),
                                              ),
                                              shiny::mainPanel(
                                                shiny::tabsetPanel(id = "coloc_main",
                                                                   shiny::tabPanel("Overview", value = "coloc_overview",
                                                                                   shiny::br(),
                                                                                   shiny::br(),
                                                                                   shiny::fluidRow(
                                                                                     shiny::column(12,align='center',
                                                                                                   withSpinner(
                                                                                                     plotly::plotlyOutput("single_coloc_plot",
                                                                                                                          width = 500,
                                                                                                                          height = 500), type=7)

                                                                                     )
                                                                                   )

                                                                   )
                                                )
                                              )
                                            )

                                          )
                          ),
                          shiny::tabPanel("About",
                                          shiny::fluidPage(
                                            shiny::sidebarLayout(
                                              shiny::sidebarPanel(
                                                width = 3,

                                              ),
                                              shiny::mainPanel(
                                                shiny::tabsetPanel(id = "about_tab",
                                                                   shiny::tabPanel("About", value = "about_overview",
                                                                                   shiny::br(),
                                                                                   shiny::br(),
                                                                                   shiny::fluidRow(
                                                                                     shiny::column(12,
                                                                                                   shiny::includeMarkdown(paste0(path.package("CONQUER"), "/About.md"))
                                                                                     )
                                                                                   )
                                                                   )
                                                )
                                              )
                                            )
                                          )
                          )
  )

  # Define server logic required to draw a histogram
  server <- function(input, output) {

    #### Info labels ####
    observeEvent(input$button, {
      shinyjs::toggle('text_div')
      if(input$tis_spec == "mod")
      {
        output$text <- renderUI({shiny::HTML("Select a tissue of interest and a module from the dropdown menu. <br><br>
          Top left: Heatmap of the correlation of the genes. <br><br>
          Top right: Enriched pathways in this module <br><br>
          Bottom: Table of eQTLs in this module")})
      }else if(input$tis_spec == "overview"){
        output$text <- renderUI({shiny::HTML("Select a tissue of interest from the dropdown menu. Click on the dots to navigate across to see more information about the modules, enriched pathways, genes and SNPs")})
      }
    })

    observeEvent(input$button2, {
      shinyjs::toggle('text_div2')
      if(input$pws == "pway")
      {
        output$text2 <- shiny::renderUI({shiny::HTML("This figure shows the pathways are tissue-specific and  tissue-shared. Click a tissue to see the enriched pathways. Click a pathway to see in which tissues the pathway is enriched.")})
      }else if(input$pws == "meqtl"){
        output$text2 <- shiny::renderUI({shiny::HTML("This table shows the DNA methylation QTLs in whole blood based on the BIOS Consortium data. For more information see the about tab.")})
      }else if(input$pws == "miqtl"){
        output$text2 <- shiny::renderUI({shiny::HTML("These tables contain the miRNA QTLs, both the experimentally determined and predicted QTLs. For more information see the about tab")})
      }else if(input$pws == "pqtl"){
        output$text2 <- shiny::renderUI({shiny::HTML("This table shows the pQTLs in plasma. For more information see the about tab.")})
      }
    })

    observeEvent(input$button3, {
      shinyjs::toggle('text_div3')
      if(input$single == "LD"){
        output$text3 <- shiny::renderUI({shiny::HTML("<h4>Linkage Disequilibrium</h4>
        Plot shows correlation of nearby SNPs, recombination rate and genes witin the region. Hovering over SNP shows position, R<sup>2</sup> with leading SNP and consequence type.<br><br>
        The table shows the linkage disequilibrium data of selected SNP. Contains IDs, position, consequence type, R<sup>2</sup> and clinincal significance of nearby SNPs.")})
      }else if(input$single == "CI"){
        output$text3 <- shiny::renderUI({shiny::HTML("<h4>Chromosomal interaction</h4>
         Circos plot of chromosomal interaction near selected SNP and in selected tissue. Also shows chromatin state segmentations of selected tissue. Outer ring contains LD-Block (R<sup>2</sup> > 0.8).<br><br>
         Hovering over link shows tissue and the linked locations. For linked genes the gene symbols and ensembl gene IDs are given. Hovering over chromatin state segmentations shows location and chromatin state.")})
      }else if(input$single == "CS"){
        output$text3 <- shiny::renderUI({shiny::HTML("<h4>Chromatin state segmentations</h4>
        Plot of the chromatin state segmentations near selected SNP and in selected tissue. <br> By default all tissues are selected. <br><br> Data and plot can separately be downloaded")})
      }else if(input$single == "QTLs"){
        output$text3 <- shiny::renderUI({shiny::HTML("<h3>Quantitative trait loci (QTLs)</h3>
                         <h4> Expression QTLs (eQTLs, GTEx v8) </h4>
                         eQTLs tab shows a hive plot with three axis (SNP - top, tissues - left, genes - right), hovering over links shows connection between nodes. Table shows all calculated eQTLs. Clicking on a row
                         in the table generates a violin plot of the normalized expression of the respective gene and genotypes of selected SNP.<br><br>
                         <h4> Colocalization </h4>
                         The y-axis of the colocalization figure shows the posterior probability (PP) that an SNP is the causal variant for the observed eQTL. The PP is based on the normalized effect size and standard error from precalculated eQTLs from GTEx. Note that new eQTLs may be found in the eQTL tab as CONQUER additionally tests SNPs. The black dot indicates the lead SNP.
                         <h4> protein QTLs (pQTLs) </h4>
                         Table of plasma pQTLs for selected SNP. <br><br>
                         <h4> microRNA QTLs (miQTLs) </h4>
                         Table of predicted miQTLs for selected SNP (left). Table of experimental determined miQTLs for selected SNP (left<br><br>
                         <h4> methylation QTLs (meQTLs). </h4>
                         Table of plasma meQTLs.")})
      }else if(input$single == "GeX"){
        output$text3 <- shiny::renderUI({shiny::HTML("<h4>Gene expression</h4>
        Plot shows normalized gene expression in all available tissues in GTEx v8 for genes near the selected SNP.")})
      }
    })






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

    output$responsiveUI_C2 <- shiny::renderUI({
      AllGenes <- getGenes(shiny::req(input$snpSel), loadedSNPs)
      shiny::selectInput(inputId = "genecoloc",
                         label = "Select gene:",
                         choices = AllGenes)
    })


    single_coloc <- eventReactive(input$run, {
      snp.in <- shiny::req(input$snpc)
      gencode.in <- shiny::req(input$gencodec)
      tissue.in <- shiny::req(input$tissueSelc)
      out <- list(gencode.in, snp.in, tissue.in)
      out
    })



    output$single_coloc_plot <- plotly::renderPlotly({
      out <- single_coloc()
      out <- getColocalizationSingle(gencodeId = out[[1]], leadSNP=out[[2]],
                                     tissue=out[[3]], loadedSNPs=loadedSNPs)
      return(out)
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

    output$moduleEnrichment <- plotly::renderPlotly({
      tissue <- shiny::req(input$tissueSel)
      plotAdvancedEnrichment(tissue = tissue, interactive=TRUE,SNPSummary=SNPSummary)
    })


    output$moduleColoc <- plotly::renderPlotly({
      if(!is.null(ColocSummary[[input$snpSel]]))
      {
        if(is.null(tissues))
        {
          plotColoc(shiny::req(input$snpSel), all.coloc=ColocSummary, loadedSNPs=loadedSNPs, filter=TRUE)
        }else{
          plotColoc(shiny::req(input$snpSel), all.coloc=ColocSummary, loadedSNPs=loadedSNPs, filter=FALSE, tissues=tissues)
        }
      }
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
      conquer.d3js::ConquerEdge(SNPSummary, disease="no")
    })

    output$EdgePlotDisease <- conquer.d3js::renderConquerEdge({
      conquer.d3js::ConquerEdge(SNPSummary, disease="yes")
    })

    LDSNPs <- shiny::reactive({
      LDSNPs <- getAllSNPLD(loadedSNPs)
    })

    ###pQTLs
    output$pQTLOverview <- DT::renderDT({
      AllTissuespQTLsData()
    },options=list(scrollX=T),selection = "single")


    AllTissuespQTLsData <- shiny::reactive({
      output <- pqtls_internal
      return(output)
    })

    ## sQTLs
    output$sQTLOverview <- DT::renderDT({
      AllTissuessQTLsData()
    },options=list(scrollX=T),selection = "single")

    AllTissuessQTLsData <- shiny::reactive({
      output <- sqtls_internal
      return(output)
    })

    ## mQTLs NG
    output$mqtls_overview_ng <- DT::renderDT({
      AllTissuesmQTLNGsData()
    },options=list(scrollX=T),selection = "single")

    AllTissuesmQTLNGsData <- shiny::reactive({
      output <- mqtls_NG_internal
      return(output)
    })

    ## mQTLs Multi
    output$mqtls_overview_lc <- DT::renderDT({
      AllTissuesmQTLLCsData()
    },options=list(scrollX=T),selection = "single")

    AllTissuesmQTLLCsData <- shiny::reactive({
      output <- mqtls_LC_internal
      return(output)
    })

    ## lQTLs
    output$lQTLOverview <- DT::renderDT({
      AllTissueslQTLsData()
    },options=list(scrollX=T),selection = "single")

    AllTissueslQTLsData <- shiny::reactive({
      return(lqtls_internal)
    })


###########################################################
#Download all QTLs
###########################################################

	#MeQTls
    output$downloadallmeqtls <- shiny::downloadHandler(
      filename = function() {
        paste("Methylation QTLs", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(meQTLs_internal, file)
      }
    )
    #miQTLs exp
    output$downloadallmiqtlsexp <- shiny::downloadHandler(
      filename = function() {
        paste("MiRNA_Experimentally_QTLs", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(miQTLexperiment_internal, file)
      }
    )
    #miQTls pred
    output$downloadallmiqtlspred <- shiny::downloadHandler(
      filename = function() {
        paste("MiRNA_Predicted_QTLs", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(miQTLpredict_internal, file)
      }
    )
    #pqtls
    output$downloadallpqtls <- shiny::downloadHandler(
      filename = function() {
        paste("Protein_QTLs", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(pqtls_internal, file)
      }
    )
    #lipidQTls
    output$downloadalllqtls <- shiny::downloadHandler(
      filename = function() {
        paste("Lipid_QTLs", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(lqtls_internal, file)
      }
    )
    #sQTLs
    output$downloadallsqtls <- shiny::downloadHandler(
      filename = function() {
        paste("Splicing_QTLs", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(sqtls_internal, file)
      }
    )
    #mQTLs nightingale
    output$downloadallmqtlsnigh <- shiny::downloadHandler(
      filename = function() {
        paste("Metabolite_QTLs_Nightingale", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(mqtls_NG_internal, file)
      }
    )
    #mQTLs multiplatform
    output$downloadallmqtlsmulti <- shiny::downloadHandler(
      filename = function() {
        paste("Metabolite_QTLs_MultiPlatform", ".xlsx", sep = "")
      },
      content = function(file) {
        rio::export(mqtls_LC_internal, file)
      }
    )

    output$pQTLOverview_LD <- DT::renderDT({
      row_selected <- shiny::req(input$pQTLOverview_rows_selected)
      LDSNPs <- LDSNPs()
      pQTLs.sub <- AllTissuespQTLsData()
      selectedSNP <- pQTLs.sub[row_selected,"rsID"]
      leadingSNP <- LDSNPs[LDSNPs$LDSNP == selectedSNP,"leadingSNP"]
      fullLD <- loadedSNPs[[leadingSNP]]$LD
      fullLD <- fullLD[fullLD$variation == selectedSNP, ]
      fullLD <- cbind(leadingSNP, fullLD)
      return(fullLD)
    })

    ###meQTLs
    AllTissuesMeQTLsData <- shiny::reactive({
      output <- meQTLs_internal
      return(output)
    })

    output$meQTLOverview <- DT::renderDT({
      AllTissuesMeQTLsData()
    },options=list(scrollX=T),selection = "single")

    output$meQTLOverview_LD <- DT::renderDT({
      row_selected <- shiny::req(input$meQTLOverview_rows_selected)
      LDSNPs <- LDSNPs()
      meQTLs.sub <- AllTissuesMeQTLsData()
      selectedSNP <- meQTLs.sub[row_selected,"rsID"]
      leadingSNP <- LDSNPs[LDSNPs$LDSNP == selectedSNP,"leadingSNP"]
      fullLD <- loadedSNPs[[leadingSNP]]$LD
      fullLD <- fullLD[fullLD$variation == selectedSNP, ]
      fullLD <- cbind(leadingSNP, fullLD)
      return(fullLD)
    })

    ###miQTLS
    #Experimental
    AllTissuesMiQTLexperiment <-  shiny::reactive({
      output <- miQTLexperiment_internal
      return(output)
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
      fullLD <- loadedSNPs[[leadingSNP]]$LD
      fullLD <- fullLD[fullLD$variation == selectedSNP, ]
      fullLD <- cbind(leadingSNP, fullLD)
      return(fullLD)
    },options=list(scrollX=T),selection = "single")

    #Predicted
    AllTissuesMiQTLpredict <- shiny::reactive({
      return(miQTLpredict_internal)
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
      fullLD <- loadedSNPs[[leadingSNP]]$LD
      fullLD <- fullLD[fullLD$variation == selectedSNP, ]
      fullLD <- cbind(leadingSNP, fullLD)
      return(fullLD)
    })

    ####SINGLE SNP####
    output$responsiveUI_A <- shiny::renderUI({

      if(input$single == "QTLs" & input$QTLs_tab == "eqtls_sub"){
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

    output$responsiveUI_D <- shiny::renderUI({
      if(input$single == "QTLs" & input$QTLs_tab == "coloc_sub"){

        shiny::sliderInput(inputId = "Colocwidth", label = "Export width", min = 5, max = 15, value = 10, step = 1)
      }
    })


    output$responsiveUI_E <- shiny::renderUI({
      if(input$single == "QTLs" & input$QTLs_tab == "coloc_sub"){

        shiny::sliderInput(inputId = "Colocheight", label = "Export height", min = 10, max = 25, value = 12, step = 1)
      }
    })


    # Load selected SNP into the environment
    output$LocusZoom <- conquer.d3js::renderConquerLocusZoom({
      conquer.d3js::ConquerLocusZoom(loadedSNPs[[input$snpSel]])
    })
    output$LocusHeader <- shiny::renderUI({
      shiny::tags$h2(sprintf("Linkage Disequilibrium for: %s",input$snpSel))
    })

    output$downloadLocus <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_LocusZoom", ".html", sep = "")
      },
      content = function(file) {
        htmlwidgets::saveWidget(widget = conquer.d3js::ConquerLocusZoom(loadedSNPs[[input$snpSel]], width = 500), file = file)
      })




    ########################LD Table########################
    #Data
    LDTable <- reactive({
      SelectedSNP <- shiny::req(input$snpSel)
      data <- loadedSNPs[[input$snpSel]]$LD
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
      data <- loadedSNPs[[input$snpSel]]$LD
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
      data <- loadedSNPs[[input$snpSel]]$LD
      data <- data[order(data$r2,decreasing = T),c("variation","start")]
      sel <- data[input$LDTable_rows_selected,"variation"]
      shiny::tags$h3(sprintf("Publications in which %s is mentioned: ", sel))
    })



    ########################eQTL Table########################
    #Data
    eQTLsData <- shiny::reactive({
      if(nrow(loadedSNPs[[input$snpSel]]$eQTLs) == 0){
        cis <- loadedSNPs[[input$snpSel]]$eQTLs
      }else{
        cis <- loadedSNPs[[input$snpSel]]$eQTLs
        cis$type <- "cis"
      }
      if(nrow(loadedSNPs[[input$snpSel]]$eQTLsTrans) == 0){
        trans <- loadedSNPs[[input$snpSel]]$eQTLsTrans
      }else{
        trans <- loadedSNPs[[input$snpSel]]$eQTLsTrans
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



    ########################pQTL Table########################
    output$pQTLsTable <- DT::renderDT({
      LDtable <- loadedSNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewpQTLs <- pqtls_internal[pqtls_internal$rsID %in% LDSNPs,]
      return(ViewpQTLs)
    },options=list(scrollX=T),selection = "single")
    ########################END########################


    ########################meQTL Table########################
    output$meQTLsTable <- DT::renderDT({
      LDtable <- loadedSNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewmeQTLs <- meQTLs_internal[meQTLs_internal$rsID %in% LDSNPs,]
      return(ViewmeQTLs)
    },options=list(scrollX=T),selection = "single")

    ########################END########################

    ########################miQTL Table########################

    output$miQTLsPred_table <- DT::renderDT({
      LDtable <- loadedSNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewmiQTLsPred <- miQTLpredict_internal[, c("SNP","Gene","miRNA","Celltype","Change","Effect")]
      return(ViewmiQTLsPred)
    })
    ########################END########################

    ########################miQTL Table########################

    output$miQTLsExp_table <- DT::renderDT({
      LDtable <- loadedSNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewmiQTLsExp <- miQTLexperiment_internal[miQTLexperiment_internal$SNP %in% LDSNPs,]
      return(ViewmiQTLsExp)
    },options=list(scrollX=T),selection = "single")
    ########################END########################

    ########################sQTL Table########################
    output$sQTLsExp_table <- DT::renderDT({
      LDtable <- loadedSNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,]
      IDs <- paste0("chr",LDSNPs$chr, "_", LDSNPs$start)
      ViewsQTLsExp <- sqtls_internal[sqtls_internal$variant_id %in% IDs,]
      return(ViewsQTLsExp)
    },options=list(scrollX=T),selection = "single")
    ########################END########################

    ########################lQTL Table########################
    output$lQTLsExp_table <- DT::renderDT({
      LDtable <- loadedSNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewlQTLsExp <- lqtls_internal[lqtls_internal$rsID %in% LDSNPs,]
      return(ViewlQTLsExp)
    },options=list(scrollX=T),selection = "single")

    ########################END########################

    ########################mQTL NG########################
    output$mQTLsNGExp_table <- DT::renderDT({
      LDtable <- loadedSNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewmngQTLsExp <- mqtls_NG_internal[mqtls_NG_internal$rsID %in% LDSNPs,]
      return(ViewmngQTLsExp)
    },options=list(scrollX=T),selection = "single")
    ########################END########################

    ########################mQTL LC########################
    output$mQTLsLCExp_table <- DT::renderDT({
      LDtable <- loadedSNPs[[input$snpSel]]$LD
      LDSNPs <- LDtable[LDtable$r2 >= 0.8,"variation"]
      ViewmlcQTLsExp <- mqtls_LC_internal[mqtls_LC_internal$rsID %in% LDSNPs,]
      return(ViewmlcQTLsExp)
    },options=list(scrollX=T),selection = "single")
    ########################END########################


    ########################Circos########################


    #Download button
    output$downloadCPlot <- shiny::renderUI({
      #shiny::downloadButton(outputId = "circosdown", "Download CC")

      customDownloadbutton("circosdown", "", icon="cloud-download",
                           style = buttonStyle,
                           class="btn btn-default shiny-download-link")

    })



    #Actual download
    output$circosdown <- shiny::downloadHandler(
      filename = function() {
        paste(shiny::req(input$snpSel),"_",shiny::req(input$chromTissue),"_CircosPlot", ".html", sep = "")
        #return("Test.html")
      },
      content = function(file) {
        htmlwidgets::saveWidget(widget =ConquerCircos(SNPData = loadedSNPs[[shiny::req(input$snpSel)]], tissue = shiny::req(input$chromTissue)), file=file)
      })
    ########################END########################






    ########################    ########################
    output$CIMessage <- shiny::renderUI({
      SNP <- shiny::req(input$snpSel)
      tissue <- shiny::req(input$chromTissue)
      if(is.null(loadedSNPs[[SNP]]$chromInt)){
        shiny::h3(sprintf("There are no known interactions around %s", SNP))
      }else{
        NULL
      }
    })

    output$Circos <- BioCircos::renderBioCircos({
      SNP <- shiny::req(input$snpSel)
      tissue <- shiny::req(input$chromTissue)
      if(is.null(loadedSNPs[[SNP]]$chromInt)){
        return(NULL)
      }else{
        ConquerCircos(loadedSNPs[[SNP]],tissue = tissue)
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
      output <- PrepareChromatinStates(loadedSNPs[[input$snpSel]],statesSel)
      return(output)
    })
    output$chromatinStates <- plotly::renderPlotly({
      data <- chromatinStatesData()
      PlotChromatinStates(StatesPlotData = data[[1]], fillScale = data[[2]],jit = data[[3]])
    })
    output$placeHolder_downloadData <- shiny::renderUI({
      data <- chromatinStatesData()
      if(!is.null(data)){
        customDownloadbutton("chromatinStates_downloadData", "", icon="file-excel",
                             style = buttonStyle,
                             class="btn btn-default shiny-download-link")
      }
    })
    output$placeHolder_downloadPlot <- shiny::renderUI({
      data <- chromatinStatesData()
      if(!is.null(data)){
        customDownloadbutton("chromatinStates_downloadPlot", "", icon="cloud-download",
                             style = buttonStyle,
                             class="btn btn-default shiny-download-link")
      }
    })

    output$placeHolder_downloadPlot_coloc <- shiny::renderUI({
      #data <- chromatinStatesData()
      if(!is.null(ColocSummary[[input$snpSel]])){
        customDownloadbutton("Colocalization_downloadPlot", "", icon="cloud-download",
                             style = buttonStyle,
                             class="btn btn-default shiny-download-link")
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
        paste(input$snpSel,"_chromatinStates", ".pdf", sep = "")
      },
      content = function(file) {
        data <- chromatinStatesData()
        pdf(file, width=8, height=12)
        PlotChromatinStates(StatesPlotData = data[[1]], fillScale = data[[2]],jit = data[[3]], interactive=FALSE) %>% print()
        dev.off()
      }
    )
    ########################END########################

    ########################VIOLIN########################
    output$Colocalization_downloadPlot <- shiny::downloadHandler(
      filename = function() {
        paste(input$snpSel,"_colocalization", ".pdf", sep = "")
      },
      content = function(file) {
        out.width.coloc <- req(input$Colocwidth) %>% as.numeric()
        out.height.coloc <- req(input$Colocheight) %>% as.numeric()

        pdf(file, width = out.width.coloc, height = out.height.coloc)
        plotColoc(shiny::req(input$snpSel), all.coloc=ColocSummary, loadedSNPs=loadedSNPs, filter=FALSE, tissues=tissues, interactive=F) %>% print()
        dev.off()
      }
    )
    ######################################################

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
      customDownloadbutton("downloadViolinPlot", "", icon="cloud-download",
                             style = buttonStyle,
                             class="btn btn-default shiny-download-link")
      })

    #Actual download
    output$downloadViolinPlot <- shiny::downloadHandler(
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
      data <- loadedSNPs[[selec]]

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
          widget = conquer.d3js::ConquerHive(SNPData = loadedSNPs[[shiny::req(input$snpSel)]],
                                             cis = T, trans = T, width = 900, height = 700),
          file = file)
      })
    ########################END########################

    geneExpressionData <- shiny::reactive({
      check_eQTL <- shiny::req(input$check_eQTL)
      SNP <- shiny::req(input$snpSel)
      transform <- FALSE
      transform <- input$log10Trans
      cisGenes <- loadedSNPs[[SNP]]$eQTLs$gencodeId %>% unique()
      transGenes <- loadedSNPs[[SNP]]$eQTLsTrans$gencodeId %>% unique()
      if(length(check_eQTL) == 2){
        genes <- c(cisGenes, transGenes)
      } else if(check_eQTL == 1 ){
        genes <- cisGenes
      }else if (check_eQTL == 2) {
        genes <- transGenes
      }
      geneExpr<- loadedSNPs[[SNP]]$geneExpr
      geneExpr <- geneExpr[geneExpr$gencodeId %in% genes,drop=F,]
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
      #shiny::downloadButton("geneExpr_Data", "Data")
      customDownloadbutton("geneExpr_Data", "", icon="file-excel",
                           style = buttonStyle,
                           class="btn btn-default shiny-download-link")
    })

    output$geneExpr_downloadPlot <- shiny::renderUI({
      #shiny::downloadButton("geneExpr_Plot", "Plot")
      customDownloadbutton("geneExpr_Plot", "", icon="cloud-download",
                           style = buttonStyle,
                           class="btn btn-default shiny-download-link")
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


customDownloadbutton <- function (outputId, label, icon = NULL, width = NULL, class = "btn btn-default shiny-download-link", ...)
{
  aTag <- shiny::tags$a(id = outputId, class = paste(class,
                                                     class), href = "", target = "_blank", download = NA,
                        shiny::icon(icon), label, ...)
}



validateIcon <- function (icon)
{
  if (is.null(icon) || identical(icon, character(0))) {
    return(icon)
  }
  else if (inherits(icon, "shiny.tag") && icon$name ==
           "i") {
    return(icon)
  }
  else {
    stop("Invalid icon. Use Shiny's 'icon()' function to generate a valid icon")
  }
}

