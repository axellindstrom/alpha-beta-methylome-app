library(shiny)
library(shinythemes)
library(ggplot2)
library(qs)
library(zip)
library(plotly)



library(showtext)
# Load a font that supports β
font_add_google("Noto Sans", "noto")
showtext_auto()
source('./Scripts/functions.R')

ui <- suppressWarnings(navbarPage(
  
  
  # Title of the app
  title = "alpha-beta-methylome",
  
  # Theme used
  theme = shinytheme("flatly"),
  
  
  # Initiate a sidebar 
  sidebarLayout(position = "left",
                sidebarLayout(
                  
                  # Add a side menu
                  sidebarPanel(
                    
                    # Add drop-down menu to choose coverage cut off for WGBS
                    selectInput('Coverage', 'Min coverage', choices = c('5', '10'), width = '200px', selected = '10'),
                    
                    # Add drop-down menu to choose comparison
                    selectInput('Comparison', 'Comparison', choices = c('Age (α-cells)',
                                                                        'Age (β-cells)',
                                                                        'α- vs β-cells',
                                                                        'ND vs pre-T2D/T2D (α-cells)', 
                                                                        'ND vs pre-T2D/T2D (β-cells)',
                                                                        'Sex (α-cells)',
                                                                        'Sex (β-cells)'), width = '200px', selected = "Alpha"),
                    
                    selectizeInput(inputId = 'Gene',
                                   label = 'Gene',
                                   choices = NULL,
                                   selected = NULL,
                                   options = list(placeholder = 'Enter gene symbol'),
                                   width = '200px'),
                    
                    selectInput('p_value', 'P-value cutoff', choices = c(1, 0.05, 0.01, 0.001), width = '200px', selected = "Alpha"),
                    
                    selectInput('Extend_search', 'Extended chr search (bp)', choices = c(0, 200, 1000, 3000, 5000), width = '200px', selected = "Alpha"),
                    
                    uiOutput("downloadButtonContainer"),
                    
                    p(HTML('<p style="text-align:justify; font-size:13px;">
                              You are free to use and adapt the outputs for any purpose, as long as you provide appropriate credit. Please cite the accompanying scientific article: 
                              <a href="https://www.nature.com/articles/s42255-026-01498-9" target="_blank">Cell-specific DNA methylation in human alpha and beta cells regulates gene expression in type 2 diabetes</a>
                           </p>')),
                    
                    width = '2'),
                    
                  
                  mainPanel(
                    tabsetPanel(
                      tabPanel('DNA Methylation',
                               column(
                                 12,
                                 tags$head(
                                   tags$style(HTML("
                                   /* 1. Hide the labels (min, max, and current value) */
                                   .irs-single { display: none !important; }
                                  
                                  /* 2. Hide the grid/tick marks */
                                  .irs-grid { display: none !important; }
                                  
                                  /* 3. Remove the color from the selection bar (make it transparent or match background) */
                                  .irs-bar {
                                    background: transparent !important;
                                    border-top: 1px solid transparent !important;
                                    border-bottom: 1px solid transparent !important;
                                  }
                                  
                                  /* 4. Optional: Make the background line uniform */
                                  .irs-line { background: #eee !important; }
                                    "))
                                 ),
                               plotOutput('Meth_plot'),
                               uiOutput("scroll_slider")),
                               tableOutput('Meth_table')),
                      
                      tabPanel('Gene Expression',
                               plotOutput("expression_plot"),
                               tableOutput('expression_table')),
                      tabPanel('Info',
                               fluidRow(
                                 column(12,
                                        p(HTML('
                                        <h3> About </h3>
                                        <p style="text-align:justify;">
                                        The alpha-beta-methylome web application is a comprehensive open resource based on the whole genome bisulfite sequencing and RNA-seq data from sorted human pancreatic islet α- and β-cells included in 
                                        <a href="https://www.nature.com/articles/s42255-026-01498-9" target="_blank">Cell-specific DNA methylation in human alpha and beta cells regulates gene expression in type 2 diabetes</a>. alpha-beta-methylome lets the user explore cell type-, T2D-, age-, and sex-associations in DNA methylation and gene expression.
                                        </p>
                                        <h3> How To </h3>
                                        <p style="text-align:justify;">
                                        <b> Min coverage: </b> Choose the minimum sequencing coverage per CpG site in the DNA methylation analysis. <b> Note! </b> Choosing a higher cutoff will result in higher resolution but may also reduce the number of CpG sites and/or samples included in the analysis.
                                        </p>
                                        <br>
                                        <p style="text-align:justify;">
                                        <b> Comparison: </b> Choose a comparison of interest.
                                        </p>
                                          <ul>
                                            <li>Age in α-cells</li>
                                            <li>Age in β-cells</li>
                                            <li>α- vs β-cells</li>
                                            <li>ND vs pre-T2D/T2D in α-cells</li>
                                            <li>ND vs pre-T2D/T2D in β-cells</li>
                                            <li>sex in α-cells</li>
                                            <li>sex in β-cells</li>
                                          </ul>
                                        <p style="text-align:justify;">
                                          In the age comparison a linear model adjusted for sex, BMI and number of days cultured is fitted to the data from individual CpG sites and gene expression. The plots and tables display the β-coefficient (slope), standard error for each CpG site and gene together with the t-statistics, and <i>p</i>-value.
                                          <br>
                                          <br>
                                          For the α- vs β-cell comparison, Wilcoxon matched pairs signed rank test is performed on individual CpG sites. Only samples with a coverage ≥<b> Min coverage </b> are included in the analysis and only CpG sites with data from at least 3 sample pairs are included in the analysis. Differential expression analysis was performed using DESeq2 (Wald test). The plots and tables display the mean β-value/expression and standard error of the mean for each CpG site and gene together with the <i>p</i>-value.
                                          <br>
                                          <br>
                                          In the ND vs T2D comparison a Wilcoxon Rank-Sum test is performed on individual CpG sites. Only samples with a coverage ≥<b> Min coverage </b>are included in the analysis and only CpG sites with data from at least 3 samples in each group are included in the analysis. Differential expression analysis was performed using DESeq2 (Wald test). The plots and tables display the mean β-value/expression and standard error of the mean for each CpG site and gene together with the <i>p</i>-value.
                                          <br>
                                          <br>
                                          In the sex comparison a Wilcoxon Rank-Sum test is performed on individual CpG sites. Only samples with a coverage ≥<b> Min coverage </b> are included in the analysis and only CpG sites with data from at least 3 samples in each group are included in the analysis. Differential expression analysis was performed using DESeq2 (Likelihood Ratio Test). The plots and tables display the mean β-value/expression and standard error of the mean for each CpG site and gene together with the <i>p</i>-value.
                                        </p>
                                        <br>
                                        <p style="text-align:justify;">
                                        <b> Gene: </b> Enter a gene symbol.
                                        </p>
                                        <br>
                                        <p style="text-align:justify;">
                                          <b> P-value cutoff: </b> <i>p</i>-value cutoff to be used to only show results with a <i>p</i>-value lower than:
                                        </p>
                                          <ul>
                                            <li>1</li>
                                            <li>0.05</li>
                                            <li>0.01</li>
                                            <li>0.001</li>
                                          </ul>
                                        <p style="text-align:justify;">
                                          <b> Extended chr search (bp): </b> Extend the region around the selected gene, to include CpG sites, upstream of the transcription start site and downstream of the transcription end site by:
                                        </p>
                                          <ul> 
                                            <li>0 bp</li>
                                            <li>200 bp</li>
                                            <li>1000 bp</li>
                                            <li>3000 bp</li>
                                            <li>5000 bp</li>
                                          </ul>
                                        <p style="text-align:justify;">
                                        For more information regarding sample characteristics and data availability see accompanying 
                                        <a href="https://www.nature.com/articles/s42255-026-01498-9" target="_blank">scientific article</a>.
                                        </p>
                                        <h3> Output Licensing </h3>
                                          <p style="text-align:justify;">
                                            All tables and plots generated by this application are licensed under a
                                            <a href="https://creativecommons.org/licenses/by/4.0/" target="_blank" rel="license">Creative Commons Attribution 4.0 International License</a> 
                                            <img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">
                                            <img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">
                                          </p>
                                        <br>
                                        <br>
                                        <p style="text-align:justify;">
                                        You are free to use and adapt the outputs for any purpose, as long as you provide appropriate credit. Please cite the accompanying scientific article: 
                                        <a href="https://www.nature.com/articles/s42255-026-01498-9" target="_blank">Cell-specific DNA methylation in human alpha and beta cells regulates gene expression in type 2 diabetes</a>
                                        </p>
                                               '))
                                 )
                               ),
                      ),
                    ))
                ),
                mainPanel(hr(),
                          div(class = "footer",
                              p(HTML('All tables and plots generated by this application are licensed under a
                              <a href="https://creativecommons.org/licenses/by/4.0/" target="_blank" rel="license">Creative Commons Attribution 4.0 International License</a> <img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">
                              '), style = "font-size:11px;")
                          ))
  )
  
))

