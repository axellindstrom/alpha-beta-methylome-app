

server <- function(input, output, session) {
  annotation <- qread('./Data/Annotations/annotations.qs')
  Genes_symbols <- annotation$Symbol
  Genes_symbols <- Genes_symbols[order(Genes_symbols)]
  updateSelectizeInput(session, 'Gene',
                       choices = Genes_symbols,
                       selected = character(0), # Clears previous selection
                       server = TRUE)
  
  
  
  # Read in stats for selected chromosome
  chromosome_data <- reactive({
    req(input$Comparison) # Ensure is not NULL 
    selected_Coverage <- input$Coverage
    selected_Comparison <- input$Comparison
    selected_Gene <- input$Gene
    selected_p_value <- as.numeric(input$p_value)
    
    # Get chromosome number
    chr <- gsub('chr', '', annotation[annotation$Symbol == selected_Gene, 'Chr'])
    
    # Load in selected stats data 
    path_to_chromosome_data <- paste0(coverage_map[[selected_Coverage]], chromosome_files[[selected_Comparison]], 'chromosome_', chr, '.qs') #chromosome_', chr, '.qs')
    chromosome_data <- qread(as.character(path_to_chromosome_data))
    
  })
  
  expression <- reactive({
    req(input$Comparison) # Ensure is not NULL 
    selected_Comparison <- input$Comparison
    selected_Gene <- input$Gene
    
    
    # Load in selected stats data 
    path_to_expression_data <- paste0('./Data/Expression/', chromosome_files[[selected_Comparison]], 'expression.qs') #chromosome_', chr, '.qs')
    expression_data <- qread(path_to_expression_data)
    expression_data <- expression_data[expression_data$gene_name == selected_Gene, ]
    
    
  })
  
  
  gene_methylation <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene
    if (!is.null(selected_Gene) && selected_Gene != "") { # Check for NULL and empty string
      # If a gene is selected, create a data frame to display it
      chr <- annotation[annotation$Symbol == selected_Gene, 'Chr']
      Start <- as.numeric(annotation[annotation$Symbol == selected_Gene, 'Start'])
      End <- as.numeric(annotation[annotation$Symbol == selected_Gene, 'End'])
      extened <- as.numeric(input$Extend_search)
      if(input$Comparison == 'Age (α-cells)' | input$Comparison == 'Age (β-cells)'){
        age_data <- chromosome_data()
        age_data <- age_data[!is.nan(age_data$Age_p_value), ]
        gene_methylation <- age_data[as.numeric(age_data$CpG) <= (End + extened) & as.numeric(age_data$CpG) >= (Start - extened), ,drop = FALSE]
        
      }else{
        gene_methylation <- chromosome_data()[,as.numeric(colnames(chromosome_data())) <= (End + extened) & as.numeric(colnames(chromosome_data())) >= (Start - extened), drop = FALSE]
        as.data.frame(gene_methylation)
      }
    }
  })
  
  
  chr_num <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene
    if (!is.null(selected_Gene) && selected_Gene != "") { # Check for NULL and empty string
      chr <- annotation[annotation$Symbol == selected_Gene, 'Chr']
    }
  })
  
  
  start_pos <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene
    if (!is.null(selected_Gene) && selected_Gene != "") { # Check for NULL and empty string
      Start <- as.numeric(annotation[annotation$Symbol == selected_Gene, 'Start'])
    }
  })
  
  
  end_pos <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene
    if (!is.null(selected_Gene) && selected_Gene != "") { # Check for NULL and empty string
      End <- as.numeric(annotation[annotation$Symbol == selected_Gene, 'End'])
    }
  })
  
  
  extended_pos <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene
    if (!is.null(selected_Gene) && selected_Gene != "") { # Check for NULL and empty string
      extened <- as.numeric(input$Extend_search)
    }
  })
  
  
  meth_plot <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene
    
        if (!is.null(selected_Gene) && selected_Gene != "") { # Check for NULL and empty string
          
          if(input$Comparison == 'Age (α-cells)' | input$Comparison == 'Age (β-cells)'){
            
            
            filtered_data_stats <- gene_methylation()
            filtered_data_stats <- filtered_data_stats[!is.nan(filtered_data_stats$Age_p_value),, drop = FALSE]
            gene_methylation <- filtered_data_stats[as.numeric(filtered_data_stats$Age_p_value) < as.numeric(input$p_value), ,drop = FALSE]
            total_sites <- nrow(gene_methylation)
            values <- c(gene_methylation$Age_Estimate + gene_methylation$Age_Std.Error, gene_methylation$Age_Estimate - gene_methylation$Age_Std.Error)
            minval <- min(values)
            maxval <- max(values)
            
            
            if(total_sites> 30){
              if(is.null(input$pos)){
                xmin <- which(gene_methylation$CpG == gene_methylation$CpG[1])
              }else{
                if(input$pos %in% gene_methylation$CpG){
                  xmin <- which(gene_methylation$CpG == input$pos)
                }else{
                  xmin <- which(gene_methylation$CpG == gene_methylation$CpG[1])
                }
              }
            if((xmin + 30) > nrow(gene_methylation)){
              xmax <- nrow(gene_methylation)
            }else{
              xmax <- xmin + 29
            }
              if(xmax - xmin < 29){
                gene_methylation <- gene_methylation[(xmax - 29):xmax, , drop = FALSE]
              }else{
                gene_methylation <- gene_methylation[xmin:xmax, , drop = FALSE]
              }
            }
            
            if(nrow(gene_methylation) > 0){
              chr <- chr_num()
              start <- (start_pos() - extended_pos())
              end <- (end_pos() + extended_pos())
              n_sites <- nrow(gene_methylation)
              p <- age_methylation_plot(gene_methylation, chr, start, end, selected_Gene, total_sites, minval, maxval)
            }else{
              if(as.numeric(input$p_value) < 1){
                p <- create_empty_plot(paste("No significant results for:", selected_Gene))
              }else{
                p <- create_empty_plot(paste("No available data for:", selected_Gene))
              }
            }
          }else{
            filtered_data_stats <- gene_methylation()
            gene_methylation <- filtered_data_stats[, (as.numeric(filtered_data_stats['p_value', ]) < as.numeric(input$p_value)) & (!is.nan(unlist(filtered_data_stats['p_value', ]))), drop = FALSE]
            total_sites <- ncol(gene_methylation)
            
            if(total_sites > 30){
              if(is.null(input$pos)){
                xmin <- which(colnames(gene_methylation) == colnames(gene_methylation)[1])
              }else{
                if(input$pos %in% colnames(gene_methylation)){
                  xmin <- which(colnames(gene_methylation) == input$pos)
                }else{
                  xmin <- which(colnames(gene_methylation) == colnames(gene_methylation)[1])
                }
              }
              

              if((xmin + 30) > ncol(gene_methylation)){
                xmax <- ncol(gene_methylation)
              }else{
                xmax <- xmin + 29
              }
              
              if(xmax - xmin < 29){
                gene_methylation <- gene_methylation[,(xmax - 29):xmax , drop = FALSE]
              }else{
                gene_methylation <- gene_methylation[,xmin:xmax, drop = FALSE]
              }
            }
            
            if(ncol(gene_methylation) > 0){
              chr <- chr_num()
              start <- (start_pos() - extended_pos())
              end <- (end_pos() + extended_pos())
              n_sites <- ncol(gene_methylation)
              p <- meth_bar_plot(gene_methylation, chr, start, end, selected_Gene, total_sites)
            }else{
              if(as.numeric(input$p_value) < 1){
                p <- create_empty_plot(paste("No significant result for:", selected_Gene))
              }else{
                p <- create_empty_plot(paste("No available data for:", selected_Gene))
              }
            }
          }
        }
         
        return(p)
  })
  
  
  meth_table <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene
    if (!is.null(selected_Gene) && selected_Gene != "") { # Check for NULL and empty string
      filtered_data_stats <- gene_methylation()
      
      
      
      if(input$Comparison == 'Age (α-cells)' | input$Comparison == 'Age (β-cells)'){
        filtered_data_stats <- as.data.frame(filtered_data_stats[as.numeric(filtered_data_stats$Age_p_value) < as.numeric(input$p_value), , drop = FALSE])
        if(nrow(filtered_data_stats) > 0){
          rownames(filtered_data_stats) <- filtered_data_stats$CpG
          filtered_data_stats <- filtered_data_stats[,-1]
          colnames(filtered_data_stats) <- c('β-coefficient', 'Std Error', 't-value', 'p-value')
          filtered_data_stats
        }
      }else{
        filtered_data_stats <- as.data.frame(t(filtered_data_stats[, (as.numeric(filtered_data_stats['p_value', ]) < as.numeric(input$p_value)) & (!is.nan(unlist(filtered_data_stats['p_value', ]))), drop = FALSE]))
        if(nrow(filtered_data_stats) > 0){
          
          
          filtered_data_stats[,4] <- filtered_data_stats[,4]
          filtered_data_stats[,2] <- filtered_data_stats[,2]
          
          filtered_data_stats[,4]  <-  gsub('-', '±', filtered_data_stats[,4] )
          filtered_data_stats[,2] <-  gsub('-', '±', filtered_data_stats[,2])
          colnames(filtered_data_stats) <- c(paste0('Mean (', gsub('.*_', '', colnames(filtered_data_stats)[1]),')'), 
                                             paste0('SEM (', gsub('.*_', '', colnames(filtered_data_stats)[2]), ')'),
                                             paste0('Mean (', gsub('.*_', '', colnames(filtered_data_stats)[3]), ')'),
                                             paste0('SEM (', gsub('.*_', '', colnames(filtered_data_stats)[4]), ')'),
                                             'p-value')
          filtered_data_stats
        }
      }
    }
  })
  
  
  epression_plot <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene
    
    
        if(selected_Gene %in% expression()$gene_name){
          if(input$Comparison == 'Age (α-cells)' | input$Comparison == 'Age (β-cells)'){
            if(nrow(expression()[expression()$gene_name == selected_Gene, ]) > 0 && as.numeric(expression()[expression()$gene_name == selected_Gene, 'pvalue']) < as.numeric(input$p_value)){
              p <- age_expression_plot(expression()[expression()$gene_name == selected_Gene, ], selected_Gene)
            }else{
              if(as.numeric(input$p_value) < 1){
                p <- create_empty_plotly(paste("No significant result for:", selected_Gene))
              }else{
                p <- create_empty_plotly(paste("No available data for:", selected_Gene))
              }
            }
          }else{
            if(as.numeric(expression()[expression()$gene_name == selected_Gene, 'pvalue']) < as.numeric(input$p_value)){
              p <- expression_bar_plot(expression(), selected_Gene)
            }else{
              if(as.numeric(input$p_value) < 1){
                p <- create_empty_plot(paste("No significant result for:", selected_Gene))
              }else{
                p <- create_empty_plot(paste("No available data for:", selected_Gene))
              }
            }
          }
          
          
          
          # If there is no available data for selected gene, render emtpy plot
        }else if(selected_Gene %in% Genes_symbols & !(selected_Gene %in% expression()$gene_name)){
          p <- create_empty_plot(paste('No data available for:', selected_Gene))
        }
        
        
        return(p)
  })
  
  
  epression_table <- reactive({
    req(input$Gene) # Ensure Comparison is not NULL before proceeding
    selected_Gene <- input$Gene

    if(selected_Gene %in% expression()$gene_name){
      if(input$Comparison == 'Age (α-cells)' | input$Comparison == 'Age (β-cells)'){
        if(nrow(expression()[expression()$gene_name == selected_Gene, ]) > 0 && as.numeric(expression()[expression()$gene_name == selected_Gene, 'pvalue']) < as.numeric(input$p_value)){

          table_data <- suppressWarnings(expression()[expression()$gene_name == selected_Gene, drop = FALSE])
          table_data <- table_data[,c('beta_coefficient', 'SE', 'pvalue')]
          colnames(table_data) <- c('β-coefficient', 'Std Error', 'p-value')
          table_data
        }
        
      }else{
        table_data <- suppressWarnings(expression()[expression()$gene_name == selected_Gene, drop = FALSE])
        table_data <- table_data[, c(2,4,3,5,1), drop = FALSE]
        
        if(input$Comparison == 'Sex (α-cells)' | input$Comparison == 'Sex (β-cells)'){
          colnames(table_data) <- c(paste0('Mean (', gsub('Mean_', '', colnames(table_data)[1]), ')'),
                                    paste0('SEM (', gsub('Mean_', '', colnames(table_data)[1]), ')'),
                                    paste0('Mean (', gsub('Mean_', '', colnames(table_data)[3]), ')'),
                                    paste0('SEM (', gsub('Mean_', '', colnames(table_data)[3]), ')'),
                                    'p-value')
        }else{
          colnames(table_data) <- c(paste0('Mean (', gsub('Mean_', '', colnames(table_data)[1]), ')'),
                                    paste0('SEM (', gsub('Mean_', '', colnames(table_data)[1]), ')'),
                                    paste0('Mean (', gsub('Mean_', '', colnames(table_data)[3]), ')'),
                                    paste0('SEM (', gsub('Mean_', '', colnames(table_data)[3]), ')'),
                                    'p-value')
        }
        table_data
      }
      if(as.numeric(expression()[expression()$gene_name == selected_Gene, 'pvalue']) < as.numeric(input$p_value)){
        table_data
      }
    }
  })
  
  
  # Observe changes in the 'Gene' dropdown and render it as a table
  observeEvent(input$Gene, {
    req(input$Gene)
    
    # Render plot of methylation level in selected gene
    output$Meth_plot <- renderPlot({
      meth_plot()
      
    },outputArgs = list(dev = "svg"))
    
    # Render table of methylation level in selected gene
    output$Meth_table <- renderTable({
      round_values(meth_table())
      
    }, ignoreNULL = TRUE, rownames = TRUE) # Set ignoreNULL to FALSE to react to initial NULL/empty state
    
    # Render plot of expression level in selected gene
    output$expression_plot <- renderPlot({
        epression_plot()
      
    }, outputArgs = list(dev = "svg"))
    
    # Render table of expression level in selected gene
    output$expression_table <- renderTable({
      round_values(epression_table())
      
      
    }, ignoreNULL = TRUE) 
    
    
    
  }, ignoreNULL = TRUE) # Set ignoreNULL to FALSE to react to initial NULL/empty state
  
  
  output$download <- downloadHandler(
    filename = function() {
      paste0('alpha-beta-methylome', '_', Sys.Date(), ".zip")
    },
    content = function(file) {
      temp_dir <- tempdir()
      
      folder_name <- 'Data' # Name of the folder inside the zip
      project_path <- file.path(temp_dir, folder_name)
      dir.create(project_path, recursive = TRUE)
      on.exit(unlink(paste0(temp_dir, '/Data'), recursive = TRUE, force = TRUE))
      
      # Save methylation data
      if(input_is_valid()[[2]]){
          
          # File path to the methylation plot
          file1_path <- file.path(project_path, paste0(input$Gene, "_methylation_plot.pdf"))
          # Save methylation plot
          ggsave(plot = meth_plot(), filename = file1_path, width = 14, height = 8, dpi = 150, units = "in")
          
          # File path to the methylation table
          file2_path <- file.path(project_path, paste0(input$Gene, "_methylation_stats.txt"))
          
          # Save methylation table
          write.table(meth_table(), file2_path, sep = '\t', col.names = TRUE, quote = FALSE, row.names = TRUE)
        # }
      }
      
      
      # Save expression data
      if(input_is_valid()[[3]]){
        
        file3_path <- file.path(project_path, paste0(input$Gene, "_expression_plot.pdf"))
        # File path to the expression table
        file4_path <- file.path(project_path, paste0(input$Gene, "_expression_stats.txt"))
        
        # Save methylation plot
        ggsave(plot = epression_plot(), filename = file3_path, width = 14, height = 8, dpi = 150, units = "in")
        
        # Save expression table
        write.table(epression_table(), file4_path, sep = '\t', col.names = TRUE, quote = FALSE, row.names = FALSE)
      }
      
      
      # Zip the folder
      zipr(zipfile = file, files = project_path, recurse = TRUE) # recurse = TRUE is important for folders
      
    },
    contentType = "application/zip"
  )
  
  
  input_is_valid <- reactive({
    
    data_output <- list('data_available' = FALSE, 'meth_data' =FALSE, 'expression_data' = FALSE)
    
    E <- !is.null(epression_table())
    M <- !is.null(meth_table())
    
    if((E && M) || (E || M)){
      data_output[['data_available']] <- TRUE
      if(M){
        data_output[['meth_data']] <- TRUE
      }
      
      if(E){
        data_output[['expression_data']] <- TRUE
      }
    }
    return(data_output)
  })
  
  #############################################################
  output$downloadButtonContainer <- renderUI({
    if (input_is_valid()[[2]] || input_is_valid()[[3]]) {
      # If there is data available, show the download button
      downloadButton("download", "Download")
    }
  })
  
  observeEvent(event_data("plotly_relayout", source = "meth_plot"), {
    relayout_data <- event_data("plotly_relayout", source = "meth_plot")
    
    
    # Check if the X-axis range was modified
    if (!is.null(relayout_data[["xaxis.range"]])) {
      xmin <- relayout_data[["xaxis.range"]][[1]]
      xmax <- relayout_data[["xaxis.range"]][[2]]
      current_width <- xmax - xmin
      
      # If the user tried to zoom (width != 30), force it back to 30
      if (abs(current_width - 30) > 0.1) {
        new_xmax <- xmin + 30

        
        #Use Proxy to snap the UI back to a 30-point window
        plotlyProxy("meth_plot_output", session) %>%
          plotlyProxyInvoke("relayout", list(
            "xaxis.range" = c(xmin, new_xmax)
          ))
      }
    }
  })
  
  data_reactive <- reactive({
    req(input$Gene)
    if(input$Comparison == 'Age (α-cells)' | input$Comparison == 'Age (β-cells)'){
      filtered_data_stats <- gene_methylation()
      filtered_data_stats <- filtered_data_stats[!is.nan(filtered_data_stats$Age_p_value),, drop = FALSE]
      gene_methylation <- filtered_data_stats[as.numeric(filtered_data_stats$Age_p_value) < as.numeric(input$p_value), ,drop = FALSE]
    }else{
      filtered_data_stats <- gene_methylation()
      gene_methylation <- filtered_data_stats[, (as.numeric(filtered_data_stats['p_value', ]) < as.numeric(input$p_value)) & (!is.nan(unlist(filtered_data_stats['p_value', ]))), drop = FALSE]
      }
    })
  
  observeEvent(data_reactive(), {
    
    if(input$Comparison == 'Age (α-cells)' | input$Comparison == 'Age (β-cells)'){
      n <- nrow(data_reactive()) - 31
    }else{
      n <- ncol(data_reactive()) - 31
    }
    
    updateSliderInput(
      session,
      inputId = "pos",
      min   = 1,
      max   = n,
      value = 1
    )
  })
  
  output$scroll_slider <- renderUI({
    if(input$Comparison == 'Age (α-cells)' | input$Comparison == 'Age (β-cells)'){
      if(nrow(data_reactive()) > 30){
        cpg_range <- c(data_reactive()$CpG[1:(length(data_reactive()$CpG) - 30)], data_reactive()$CpG[length(data_reactive()$CpG)])
        shinyWidgets::sliderTextInput(
          inputId = "pos",
          label = "<-- Slide to view more -->",
          choices = cpg_range,  # custom allowed values
          selected = cpg_range[1],
          grid = TRUE,     # optional grid for steps
          width = "100%"
        )
      }
    }else{
      if(ncol(data_reactive()) > 30){
        cpg_range <- c(colnames(data_reactive())[1:(length(colnames(data_reactive())) - 30)],colnames(data_reactive())[length(colnames(data_reactive()))])
        shinyWidgets::sliderTextInput(
          inputId = "pos",
          label = "<-- Slide to view more -->",
          choices = cpg_range,  # custom allowed values
          selected = cpg_range[1],
          grid = TRUE,     # optional grid for steps
          width = "100%"
        )
      }
    }
  })
}






