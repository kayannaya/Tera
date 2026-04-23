library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library("biomaRt")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("stringr")
library('corrr')
library('ggpubr')
library('pheatmap')
library(stringr)
library(zFPKM)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(immunedeconv)
library(Rtsne)
library(tidyr)
library(stringr)
library(xlsx)

library(AnnotationDbi)
library(org.Hs.eg.db)

## EXTRACTING DATA FROM TCGA
## VARIANCE
varFilt <- function(matrix, thres) {
  CVFILTER <- thres
  mean_nolym <- apply(matrix, 1, mean)
  var_nolym <- apply(matrix, 1, var)
  cv_nolym <- abs(var_nolym / mean_nolym)
  filt_genes <- subset(cv_nolym, cv_nolym > CVFILTER)
  
  exprZero <- matrix
  expr.dat <- matrix[rownames(matrix) %in% names(filt_genes), ]
  exprZero <-
    subset(exprZero,!(rownames(exprZero) %in% names(filt_genes)))
  
  return(expr.dat)
}

## BATCH EFFECT REMOVAL
removeBatchEffect = function(expr.dat, thres = 0, md) {
  expr.dat <- log2(expr.dat + 1)
  
  CVFILTER <- thres
  mean_nolym <- apply(expr.dat, 1, mean)
  var_nolym <- apply(expr.dat, 1, var)
  cv_nolym <- abs(var_nolym / mean_nolym)
  filt_genes <- subset(cv_nolym, cv_nolym > CVFILTER)
  
  exprZero <- expr.dat
  expr.dat <- expr.dat[rownames(expr.dat) %in% names(filt_genes), ]
  exprZero <-
    subset(exprZero,!(rownames(exprZero) %in% names(filt_genes)))
  
  expr.limma = limma::removeBatchEffect(as.matrix(expr.dat), meta$batch)
  expr.limma = rbind(expr.limma, exprZero)
  return(expr.limma)
}

## MIR31HG DATA
get_gene = function(x, data) {
  row = grep(x, rownames(data))
  data = as.data.frame(data)
  gene = data[row, ]
  gene = t(gene)
  colnames(gene) = x
  
  return(gene)
}


## zFPKM

activeG <- function(se) {
  activeGenes <- assay(se, "zFPKM") %>%
    mutate(gene = rownames(assay(se, "zFPKM"))) %>%
    gather(barcode, zfpkm,-gene) %>%
    left_join(dplyr::select(as.data.frame(colData(se)), barcode, sample_type), by =
                "barcode") %>%
    group_by(gene, sample_type) %>%
    summarize(median_zfpkm = median(zfpkm)) %>%
    ungroup() %>%
    mutate(active = (median_zfpkm > -3)) %>%
    filter(active) %>%
    dplyr::select(gene) %>%
    distinct()
  seActive <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(assay(se, "tpm")[activeGenes$gene,])),
                                   colData = colData(se))
  active = assay(seActive, "counts")
  return(active)
}


## CONVERT ENSEMBL TO SYMBOL
convert <- function(x) {
  x = x
  rownames(x) <- str_replace(rownames(x),
                             pattern = ".[0-9]+$",
                             replacement = "") # replaces everything after . with nothing
  rownames(x) <- mapIds(
    org.Hs.eg.db,
    keys = rownames(x),
    # column input
    column = 'SYMBOL',
    # Desired output
    keytype = 'ENSEMBL',
    # column input type
    multiVals = "first" # selects first to match
  )
  return(x)
}

## IMMUNE DECONVOLUTION
TME <- function(expr, method) {
  imm = immunedeconv::deconvolute(gene_expression =  expr, method)
  imm = t(imm)
  colnames(imm) = imm[1, ]
  imm = imm[-1, ]
  
  imm = as.data.frame(imm)
  samples = rownames(imm)
  
  imm = sapply(imm, as.numeric)
  rownames(imm) = samples
  
  return(imm)
}

remove_symbols <- function(mat) {
  colnames(mat) <-
    gsub("[[:punct:]]", "", colnames(mat))  # Removing symbols using gsub
  return(mat)
}

frac_plot = function(data_imm) {
  plot = data_imm %>%
    gather(sample, fraction,-cell_type) %>%
    # plot as stacked bar chart
    ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(limits = rev(levels(data_imm)))
  return(plot)
}

# Normalize fractions
frac_norm = function(data_TME) {
  test = data_TME[, 1:ncol(data_TME)-1]
  rowsums = rowSums(test)
  to_multiply = 1 / rowsums
  out = test * as.vector(t(to_multiply))
  out = cbind(out, data_TME[, ncol(data_TME)])
  colnames(out)[ncol(out)] = "uncharacterized cells"
  out = na.omit(out)
  return(out)
}

normfrac_plot = function(data_frac_norm) {
  data_frac_norm = as.data.frame(t(data_frac_norm))
  data_frac_norm$cell_type = rownames(data_frac_norm)
  plot = data_frac_norm %>%
    gather(sample, fraction, -cell_type) %>%
    # plot as stacked bar chart
    ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_brewer(palette = "Paired") +
    scale_x_discrete(limits = rev(levels(out)))
  return(plot)
}

get_plot_data = function(geneTarget,
                         input_list,
                         gene_col_name,
                         sample_subset = NULL) {
  out_list = list()
  if (is.null(sample_subset)) {
    for (i in names(input_list)) {
      print(all(rownames(geneTarget) == rownames(input_list[[i]])))
      o = as.data.frame(cbind(input_list[[i]], setNames(geneTarget, gene_col_name)))
      names(o)[length(names(o))]<-gene_col_name
      out_list[[i]] = o
    }
  } else{
    for (i in names(input_list)) {
      dat = input_list[[i]][sample_subset,]
      gene = as.data.frame(geneTarget[sample_subset,])
      colnames(gene) = gene_col_name
      
      print(all(rownames(gene) == rownames(dat)))
      o = as.data.frame(cbind(dat, gene))
      out_list[[i]] = o
    }
  }
  return(out_list)
}

get_by_cells = function(plot_data,
                        geneTarget,
                        method_list,
                        cell_list,
                        sample_subset=NULL) {
  by_cells = list()
  if (is.null(sample_subset)) {
    for (i in cell_list) {
      by_cells[[i]] = as.data.frame(geneTarget)
      for (deconv in method_list) {
        print(all(rownames(by_cells[[i]]) == rownames(plot_data[[deconv]])))
        data = plot_data[[deconv]]
        cells = data %>% select(contains(i))
        if (ncol(cells) != 0) {
          colnames(cells) = paste(colnames(cells), suffix = deconv, sep = "-")
          by_cells[[i]] = as.data.frame(cbind(by_cells[[i]], cells))
        } else{
          print(i)
        }
      }
    }
  } else{
    for (i in cell_list) {
      by_cells[[i]] = as.data.frame(geneTarget[sample_subset,])
      for (deconv in method_list) {
        print(all(rownames(by_cells[[i]]) == rownames(plot_data[[deconv]][sample_subset, ])))
        data = plot_data[[deconv]][sample_subset,]
        cells = data %>% select(contains(i))
        if (ncol(cells) != 0) {
          colnames(cells) = paste(colnames(cells), suffix = deconv, sep = "-")
          by_cells[[i]] = as.data.frame(cbind(by_cells[[i]], cells))
        } else{
          print(i)
        }
      }
    }
  }
  return(by_cells)
}

plot_multiple_y <-
  function(data,
           x_column,
           y_columns,
           xlab_text = "X Axis Label",
           ylab_text = "Y Axis Label",
           method = "loess") {
    library(ggplot2)
    
    # Initialize an empty list to store the plots
    plots_list <- list()
    # Loop through each y column and generate a plot with LOESS, correlation, and p-value
    if (method == "loess") {
      for (y_col in y_columns) {
        loess_fit <- loess(data[[y_col]] ~ data[[x_column]])
        fitted_values <- predict(loess_fit)
        correlation <- cor(data[[y_col]], fitted_values)
        correlation_test <- cor.test(data[[y_col]], fitted_values)
        p_value <- correlation_test$p.value
        data_plot <-
          data.frame(x = data[[x_column]], y = data[[y_col]])
        
        plot = ggplot(data_plot, aes(x = x, y = y)) +
          geom_point() +
          geom_smooth(span = 0.7, lwd = 1) +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = paste("Correlation (LOESS):", round(correlation, 2)),
            hjust = 1.1,
            vjust = 2,
            size = 4
          ) +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = paste("p-value:", format(p_value, scientific = TRUE)),
            hjust = 1.1,
            vjust = 3.7,
            size = 4,
            color = "black"
          ) +
          theme_bw() +
          xlab(xlab_text) +
          ylab(y_col)
        
        plots_list[[y_col]] <- plot
      }
    }
    else if (method == "lm") {
      for (y_col in y_columns) {
        data_plot <- data.frame(x = data[[x_column]], y = data[[y_col]])
        
        plot = ggplot(data_plot, aes(x = x, y = y)) +
          geom_point() +
          geom_smooth(method = "lm",
                      lwd = 1,
                      color = "red") +
          stat_cor() +
          theme_bw() +
          xlab(xlab_text) +
          ylab(y_col)
        plots_list[[y_col]] <- plot
        
      }
    }
    else {
      "Method required, options: lm, loess"
    }
    return(plots_list)
    
  }

plot_multiple_x <-
  function(data,
           x_columns,
           y_column,
           xlab_text = "X Axis Label",
           ylab_text = "Y Axis Label",
           method = "loess") {
    library(ggplot2)
    
    # Initialize an empty list to store the plots
    plots_list <- list()
    # Loop through each y column and generate a plot with LOESS, correlation, and p-value
    if (method == "loess") {
      for (x_col in x_columns) {
        loess_fit <- loess(data[[y_column]] ~ data[[x_col]])
        fitted_values <- predict(loess_fit)
        correlation <- cor(data[[y_column]], fitted_values)
        correlation_test <- cor.test(data[[y_column]], fitted_values)
        p_value <- correlation_test$p.value
        data_plot <-
          data.frame(x = data[[x_col]], y = data[[y_column]])
        
        plot = ggplot(data_plot, aes(x = x, y = y)) +
          geom_point() +
          geom_smooth(span = 0.7, lwd = 1) +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = paste("Correlation (LOESS):", round(correlation, 2)),
            hjust = 1.1,
            vjust = 2,
            size = 4
          ) +
          annotate(
            "text",
            x = Inf,
            y = Inf,
            label = paste("p-value:", format(p_value, scientific = TRUE)),
            hjust = 1.1,
            vjust = 3.7,
            size = 4,
            color = "black"
          ) +
          theme_bw() +
          xlab(x_col) +
          ylab(ylab_text)
        
        plots_list[[x_col]] <- plot
      }
    }
    else if (method == "lm") {
      for (x_col in x_columns) {
        data_plot <- data.frame(x = data[[x_col]], y = data[[y_column]])
        
        plot = ggplot(data_plot, aes(x = x, y = y)) +
          geom_point() +
          geom_smooth(method = "lm",
                      lwd = 1,
                      color = "red") +
          stat_cor() +
          theme_bw() +
          xlab(x_col) +
          ylab(ylab_text)
        plots_list[[x_col]] <- plot
        
      }
    }
    else {
      "Method required, options: lm, loess"
    }
    return(plots_list)
    
  }

#
# plot_lm_multiple_y_with_pvalue <- function(data, x_column, y_columns, xlab_text = "X Axis Label", ylab_text = "Y Axis Label") {
#   library(ggplot2)
#
#   # Initialize an empty list to store the plots
#   plots_list <- list()
#
#   # Loop through each y column and generate a plot with LOESS, correlation, and p-value
#   for (y_col in y_columns) {
#
#     data_plot <- data.frame(x = data[[x_column]], y = data[[y_col]])
#
#     plot = ggplot(data_plot, aes(x = x, y = y)) +
#       geom_point() +
#       geom_smooth(method = "lm", lwd = 1, color = "red") +
#       stat_cor() +
#       theme_bw() +
#       xlab(xlab_text) +
#       ylab(ylab_text)
#     #+ggtitle(paste(title_text, "(", y_col, ")", sep = ""))
#
#     # Store each plot in the list
#     plots_list[[y_col]] <- plot
#   }
#
#   # Return the list of plots
#   return(plots_list)
# }
#
# plot_lm_multiple_x_with_pvalue <- function(data, x_columns, y_column, xlab_text = "X Axis Label", ylab_text = "Y Axis Label") {
#   library(ggplot2)
#
#   # Initialize an empty list to store the plots
#   plots_list <- list()
#
#   # Loop through each x column and generate a plot with linear regression line, correlation, and p-value
#   for (x_col in x_columns) {
#     data_plot <- data.frame(x = data[[x_col]], y = data[[y_column]])
#
#     plot <- ggplot(data_plot, aes(x = x, y = y)) +
#       geom_point() +
#       geom_smooth(method = "lm", lwd = 1, color = "red") +
#       stat_cor() +
#       theme_bw() +
#       xlab(paste(xlab_text, "(", x_col, ")", sep = "")) +
#       ylab(ylab_text) +
#       ggtitle(paste("Relationship between", x_col, "and", y_column))
#
#     # Store each plot in the list
#     plots_list[[paste(x_col, y_column, sep = "_")]] <- plot
#   }
#
#   # Return the list of plots
#   return(plots_list)
# }
#
# plot_loess_multiple_y_with_pvalue <- function(data, x_column, y_columns, xlab_text = "X Axis Label", ylab_text = "Y Axis Label", title_text = "Your Plot Title") {
#   library(ggplot2)
#
#   # Initialize an empty list to store the plots
#   plots_list <- list()
#
#   # Loop through each y column and generate a plot with LOESS, correlation, and p-value
#   for (y_col in y_columns) {
#     # Perform LOESS regression
#     loess_fit <- loess(data[[y_col]] ~ data[[x_column]])
#
#     # Predict values using LOESS fit
#     fitted_values <- predict(loess_fit)
#
#     # Calculate correlation between 'y' and fitted values
#     correlation <- cor(data[[y_col]], fitted_values)
#
#     # Perform correlation test to get p-value
#     correlation_test <- cor.test(data[[y_col]], fitted_values)
#     p_value <- correlation_test$p.value
#
#     # Create data frame for plotting
#     data_plot <- data.frame(x = data[[x_column]], y = data[[y_col]])
#
#     # Create the scatterplot with LOESS line, annotate with the correlation value and p-value
#     plot <- ggplot(data_plot, aes(x = x, y = y)) +
#       geom_point() + # Scatterplot of your data points
#       geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + # LOESS line
#       # Annotate with the correlation value and p-value
#       annotate("text", x = Inf, y = Inf, label = paste("Correlation (LOESS):", round(correlation, 2)),
#                hjust = 1, vjust = 1, size = 4) +
#       annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, scientific = TRUE)),
#                hjust = 1, vjust = 0.8, size = 4, color = "red") +
#       xlab(xlab_text) +
#       ylab(ylab_text) +
#       ggtitle(paste(title_text, "(", y_col, ")", sep = ""))
#
#     # Store each plot in the list
#     plots_list[[y_col]] <- plot
#   }
#
#   # Return the list of plots
#   return(plots_list)
# }
#
# plot_loess_multiple_x_with_pvalue <- function(data, x_columns, y_column, xlab_text = "X Axis Label", ylab_text = "Y Axis Label", title_text = "Your Plot Title") {
#   library(ggplot2)
#
#   # Initialize an empty list to store the plots
#   plots_list <- list()
#
#   # Loop through each x column and generate a plot with LOESS, correlation, and p-value
#   for (x_col in x_columns) {
#     # Perform LOESS regression
#     loess_fit <- loess(data[[y_column]] ~ data[[x_col]])
#
#     # Predict values using LOESS fit
#     fitted_values <- predict(loess_fit)
#
#     # Calculate correlation between 'y' and fitted values
#     correlation <- cor(data[[y_column]], fitted_values)
#
#     # Perform correlation test to get p-value
#     correlation_test <- cor.test(data[[y_column]], fitted_values)
#     p_value <- correlation_test$p.value
#
#     # Create data frame for plotting
#     data_plot <- data.frame(x = data[[x_col]], y = data[[y_column]])
#
#     # Create the scatterplot with LOESS line, annotate with the correlation value and p-value
#     plot <- ggplot(data_plot, aes(x = x, y = y)) +
#       geom_point() + # Scatterplot of your data points
#       geom_smooth(method = "loess", formula = y ~ x, se = FALSE) + # LOESS line
#       # Annotate with the correlation value and p-value
#       annotate("text", x = Inf, y = Inf, label = paste("Correlation (LOESS):", round(correlation, 2)),
#                hjust = 1, vjust = 1, size = 4) +
#       annotate("text", x = Inf, y = Inf, label = paste("p-value:", format(p_value, scientific = TRUE)),
#                hjust = 1, vjust = 0.8, size = 4, color = "red") +
#       xlab(paste(xlab_text, "(", x_col, ")", sep = "")) +
#       ylab(ylab_text) +
#       ggtitle(paste(title_text, "(", y_column, " vs ", x_col, ")", sep = ""))
#
#     # Store each plot in the list
#     plots_list[[x_col]] <- plot
#   }
#
#   # Return the list of plots
#   return(plots_list)
# }


# melted_data <- df %>%
#   pivot_longer(cols = -c(samples, ENSG00000171889),
#                names_to = "Category",
#                values_to = "Value") %>%
#   filter(!is.na(Value))  # Filter out NA values if any
#
# # Create a separate dataframe for Gene_Expression
# gene_expression_data <- df %>%
#   select(samples, ENSG00000171889) %>%
#   mutate(Category = "ENSG00000171889",
#          Value = ENSG00000171889) %>%
#   select(Category, samples, Value)
# final_data <- bind_rows(melted_data, gene_expression_data)
#
#
# ggplot(final_data, aes(x = ENSG00000171889, y = Value, color = Category)) +
#   geom_point() +
#   labs(x = "Cell Proportion", y = "mir31hg Gene Expression", color = "Cell Category") +
#   theme_minimal()
#
# # Conversion of matrix to dataframe
# tsne_plot <- data.frame(x = tsne_out$Y[,1],
#                         y = tsne_out$Y[,2])
#
# # Plotting the plot using ggplot() function
# ggplot2::ggplot(tsne_plot,label=Species) + geom_point(aes(x=x,y=y))
#
# ## LIMMA PIPELINE
# limma_pipeline = function(
    #     tcga_data,
#     condition_variable,
#     reference_group=NULL){
#
#   design_factor = colData(tcga_data)[, condition_variable, drop=T]
#
#   group = factor(design_factor)
#   if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
#
#   design = model.matrix(~ group)
#
#   dge = DGEList(counts=assay(tcga_data),
#                 samples=colData(tcga_data),
#                 genes=as.data.frame(rowData(tcga_data)))
#
#   # filtering
#   keep = filterByExpr(dge,design)
#   dge = dge[keep,,keep.lib.sizes=FALSE]
#   rm(keep)
#
#   # Normalization (TMM followed by voom)
#   dge = calcNormFactors(dge)
#   v = voom(dge, design, plot=TRUE)
#
#   # Fit model to data given design
#   fit = lmFit(v, design)
#   fit = eBayes(fit)
#
#   # Show top genes
#   topGenes = topTable(fit, coef=ncol(design), number=100, sort.by="p")
#
#   return(
#     list(
#       voomObj=v, # normalized data
#       fit=fit, # fitted model and statistics
#       topGenes=topGenes # the 100 most differentially expressed genes
#     )
#   )
# }
#
# limma_res = limma_pipeline(
#   tcga_data=tcga_data,
#   condition_variable="definition",
#   reference_group="Solid Tissue Normal"
# )
#
# plot_PCA = function(voomObj, condition_variable){
#   group = factor(voomObj$targets[, condition_variable])
#   pca = prcomp(t(voomObj$E))
#   # Take PC1 and PC2 for the plot
#   plot(pca$x[,1:2],col=group, pch=19)
#   # include a legend for points
#   legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
#   return(pca)
# }
# res_pca = plot_PCA(limma_res$voomObj, "definition")
