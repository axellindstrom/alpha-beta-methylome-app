

# Map to correct path based on input options
coverage_map <- list('5' = './Data/Results_5/', '10' = './Data/Results_10/') # change to ../Data/Results_5/
chromosome_files <- list('ND vs pre-T2D/T2D (α-cells)' = 'alpha/ctrl_vs_T2D/',
                         'ND vs pre-T2D/T2D (β-cells)' = 'beta/ctrl_vs_T2D/',
                         'α- vs β-cells' = 'alpha_vs_beta/',
                         'Sex (α-cells)' = 'alpha/Sex/',
                         'Sex (β-cells)' = 'beta/Sex/',
                         'Age (α-cells)' = 'alpha/Age/',
                         'Age (β-cells)' = 'beta/Age/')

# Plot methylation level at CpG sites
meth_bar_plot <- function(chr_stats, chr, start, end, selected_Gene, total_sites){
  chr_stats <- as.data.frame(t(chr_stats))
  
  chr_stats$CpG <- gsub('X', '', rownames(chr_stats))
  
  control_group <- chr_stats[,c(1,2,5,6)]
  comp_group <- chr_stats[,c(3,4,5,6)]
  colnames(control_group) <- c('mean', 'SEM', 'p_value', 'CpG')
  colnames(comp_group) <- c('mean', 'SEM', 'p_value', 'CpG')
  
  control_group$Group <- rep(gsub('.*_', '', colnames(chr_stats)[1]), nrow(control_group))
  comp_group$Group <- rep(gsub('.*_', '', colnames(chr_stats)[3]), nrow(control_group))
  
  chr_stats <- rbind(control_group, comp_group)
  if('Control' %in% chr_stats$Group){
    chr_stats$Group <- gsub('Control', 'ND', chr_stats$Group)
  }
  
  if('T2D' %in% chr_stats$Group){
    chr_stats$Group <- gsub('T2D', 'pre-T2D/T2D', chr_stats$Group)
  }
  
  if('Males' %in% chr_stats$Group){
    chr_stats$Group <- factor(chr_stats$Group, levels = c('Males', 'Females'))
  }

  chr_stats$mean <- round(chr_stats$mean*100,0)
  chr_stats$SEM <- round(chr_stats$SEM*100,2)
  p <- ggplot(chr_stats, aes(fill = Group, y = mean, x = CpG)) +
    geom_bar(position = position_dodge(), stat = 'identity') +
    geom_errorbar(aes(x=CpG, ymin=mean-SEM, ymax=mean+SEM), position = position_dodge(0.9), width=0.4, colour="black", alpha=0.7, linewidth=0.3) +
    guides(fill=guide_legend(title="")) +
    scale_fill_manual(values=c("mediumblue", "deeppink1")) +
    theme_classic(base_size = 20) + 
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 12))+ #8.5)) +
    coord_cartesian(ylim=c(0,max(c(100, max(chr_stats$mean + chr_stats$SEM))))) + 
    scale_y_continuous(expand=c(0,0)) +
    ylab('Methylation (%)') +
    xlab('Position') +
    ggtitle(paste0(chr, ': ', start, ' - ', end, ' (', selected_Gene, ')'))

  
  return(p)
  
}



# ß-coefficients of age in regards to methylation level at CpG sites
age_methylation_plot <- function(chr_stats, chr, start, end, selected_Gene, total_sites, minval, maxval){
  

  lims <- c(minval - 0.0005, maxval + 0.0005)
  

  
  chr_stats$Group <- sign(chr_stats$Age_Estimate)
  chr_stats$Group[chr_stats$Group == 1] <- 'Pos'
  chr_stats$Group[chr_stats$Group == -1] <- 'Neg'
  chr_stats$Group <- factor(chr_stats$Group , levels = c('Pos', 'Neg'))
  chr_stats$CpG <- as.factor(chr_stats$CpG)
  
  p <- ggplot(chr_stats, aes(fill = Group, y = Age_Estimate, x = CpG)) +
      geom_bar(stat = "identity", aes(fill = Group)) +
      geom_errorbar(aes(x=CpG, ymin=Age_Estimate-Age_Std.Error, ymax=Age_Estimate+Age_Std.Error), position = position_dodge(0.9), width=0.4, colour="black", alpha=0.7, linewidth=0.3) +
      scale_fill_manual(values=c('Neg' = "brown3", 'Pos' = "steelblue2")) +
      guides(fill=guide_legend(title="")) +
      theme_classic(base_size = 20) + 
      theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = 15)) + #, legend.position="none") + #8.5)) +
      coord_cartesian(ylim=lims) + 
      scale_y_continuous(expand=c(0,0)) +
      ylab('β-coefficient') +
      xlab('Position') +
      ggtitle(paste0(chr, ': ', start, ' - ', end, ' (', selected_Gene, ')'))
      ggtitle(selected_Gene) +
      geom_hline(yintercept = 0, color = 'grey54')
  
  return(p)
}


age_expression_plot <- function(expression_data, selected_Gene){
  
  
  if((expression_data$beta_coefficient - expression_data$SE - 0.0005) <= 0 && (expression_data$beta_coefficient + expression_data$SE + 0.0005) >0){
    lims <- c((expression_data$beta_coefficient - expression_data$SE - 0.0005), (expression_data$beta_coefficient + expression_data$SE + 0.0005))
  }else if((expression_data$beta_coefficient - expression_data$SE - 0.0005) <= 0 && (expression_data$beta_coefficient + expression_data$SE + 0.0005) <= 0){
    lims <- c((expression_data$beta_coefficient - expression_data$SE - 0.0005), 0)
  }else{
    lims <- c(0, (expression_data$beta_coefficient + expression_data$SE + 0.0005))
  }
  

  
  
  
  expression_data$Group <- as.factor(sign(expression_data$beta_coefficient))
  expression_data$gene_name <- as.factor(expression_data$gene_name)
  
  
  p <- ggplot(expression_data, aes(fill = Group, y = beta_coefficient, x = gene_name)) +
    geom_bar(stat = "identity", aes(fill = Group)) +
    geom_errorbar(aes(x=gene_name, ymin=beta_coefficient-SE, ymax=beta_coefficient+SE), position = position_dodge(0.9), width=0.4, colour="black", alpha=0.7, linewidth=0.3) +
    guides(fill=guide_legend(title="")) +
    scale_fill_manual(values=c('-1' = "brown3", '1' = "steelblue2")) +
    theme_classic(base_size = 20) + 
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), legend.position="none") + #8.5)) +
    coord_cartesian(ylim=lims) + 
    scale_y_continuous(expand=c(0,0)) +
    ylab('β-coefficient') +
    xlab('Gene') +
    ggtitle(selected_Gene) +
    geom_hline(yintercept = 0, color = 'grey54')
  
  return(p)
}

# Plot expression of genes
expression_bar_plot <- function(expression_data, gene_name){

  exp_data <- expression_data[expression_data$gene_name == gene_name, , drop = FALSE]
  exp_data <- exp_data[1, , drop = FALSE]

  

  
  control_group <- exp_data[,c(1,2,4)]
  comp_group <- exp_data[,c(1,3,5)]

  
  
  control_group$Group <- gsub('Mean_', '', colnames(control_group)[2])
  comp_group$Group <- gsub('Mean_', '', colnames(comp_group)[2])
  
  colnames(control_group) <- c('pvalue', 'mean', 'SEM', 'Group')
  colnames(comp_group) <- c('pvalue', 'mean', 'SEM', 'Group')
  
  
  expr_stats <- rbind(control_group, comp_group)
  
  if('Control' %in% expr_stats$Group){
    expr_stats$Group <- gsub('Control', 'ND', expr_stats$Group)
  }
  
  if('Males' %in% expr_stats$Group){
    expr_stats$Group <- factor(expr_stats$Group, levels = c('Males', 'Females'))
  }

  p <- ggplot(expr_stats, aes(fill = Group, y = mean, x = Group)) +
    geom_bar(position = position_dodge(), stat = 'identity') +
    geom_errorbar(aes(ymin=mean-SEM, ymax=mean+SEM), position = position_dodge(0.9), width=0.4, colour="black", alpha=0.7, linewidth=0.3) +
    guides(fill=guide_legend(title="")) +
    scale_fill_manual(values=c("mediumblue", "deeppink1")) +
    theme_classic(base_size = 20) + 
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))+ #8.5)) +
    coord_cartesian(ylim=c(0,max(expr_stats$mean + expr_stats$SEM))) + 
    scale_y_continuous(expand=c(0,0)) +
    ylab('Count per million') +
    xlab('Group') +
    ggtitle(gene_name)
  
  return(p)
  
}


round_values <- function(input_table, digits = 3){
  r_names <- rownames(input_table)
  table_data <- data.frame(lapply(input_table, function(x) {
    x <- as.numeric(x)
    if (is.numeric(x)) formatC(x, format = "f", digits = digits) else x
  }), check.names = FALSE)
  
  rownames(table_data) <- r_names
  return(table_data)
  
}






create_empty_plot <- function(message) {
  df <- data.frame()
  p <- ggplot(df) +
    geom_point() + xlim(0, 10) + ylim(0, 10) +
    annotate("text", x=5, y=5, label= paste(message), size = 10) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(p)
}




