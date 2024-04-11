rm(list = ls())
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(funkyheatmap)
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/benchmark.RData')
Scenario = c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3))
Sample_size = rep(c(200, 300, 400), 3)

####FDR
bama_FDR_matrix = cbind(t(bama_FDR_matrix[complete.cases(bama_FDR_matrix),]),Scenario,Sample_size)
bama_FDR_matrix = cbind(bama_FDR_matrix,Model = rep('bama',9))
# HIMA_FDR_matrix = cbind(t(HIMA_FDR_matrix[complete.cases(HIMA_FDR_matrix),]),Scenario,Sample_size)
mednet_FDR_matrix = cbind(t(mednet_FDR_matrix[complete.cases(mednet_FDR_matrix),]),Scenario,Sample_size)
mednet_FDR_matrix = cbind(mednet_FDR_matrix,Model = rep('medNet',9))
medfix_FDR_matrix = cbind(t(medfix_FDR_matrix[complete.cases(medfix_FDR_matrix),]),Scenario,Sample_size)
medfix_FDR_matrix = cbind(medfix_FDR_matrix,Model = rep('MedFix',9))
medfix_FDR_matrix_long = pivot_longer(as.data.frame(medfix_FDR_matrix), cols = starts_with("V"), 
                                  names_to = "Variable", values_to = "Value")
mednet_FDR_matrix_long = pivot_longer(as.data.frame(mednet_FDR_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
bama_FDR_matrix_long = pivot_longer(as.data.frame(bama_FDR_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
m_FDR_long = rbind(medfix_FDR_matrix_long,mednet_FDR_matrix_long,bama_FDR_matrix_long)

p_1 = ggplot(as.data.frame(m_FDR_long), aes(x = Sample_size, y = as.numeric(Value), fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#B0B0B0', '#E29998', '#67A2A3')) +
  theme_minimal() +
  theme_bw()+
  coord_flip() +
  theme(text = element_text(size = 12), 
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mediator FDR Comparison Across Different Methods", 
       x = "Sample Size", 
       y = "FDR Value",
       fill = "Model")


####FPR
bama_FPR_matrix = cbind(t(bama_FPR_matrix[complete.cases(bama_FPR_matrix),]),Scenario,Sample_size)
bama_FPR_matrix = cbind(bama_FPR_matrix,Model = rep('bama',9))
# HIMA_FPR_matrix = cbind(t(HIMA_FPR_matrix[complete.cases(HIMA_FPR_matrix),]),Scenario,Sample_size)
mednet_FPR_matrix = cbind(t(mednet_FPR_matrix[complete.cases(mednet_FPR_matrix),]),Scenario,Sample_size)
mednet_FPR_matrix = cbind(mednet_FPR_matrix,Model = rep('medNet',9))
medfix_FPR_matrix = cbind(t(medfix_FPR_matrix[complete.cases(medfix_FPR_matrix),]),Scenario,Sample_size)
medfix_FPR_matrix = cbind(medfix_FPR_matrix,Model = rep('MedFix',9))
medfix_FPR_matrix_long = pivot_longer(as.data.frame(medfix_FPR_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
mednet_FPR_matrix_long = pivot_longer(as.data.frame(mednet_FPR_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
bama_FPR_matrix_long = pivot_longer(as.data.frame(bama_FPR_matrix), cols = starts_with("V"), 
                               names_to = "Variable", values_to = "Value")
m_FPR_long = rbind(medfix_FPR_matrix_long,mednet_FPR_matrix_long,bama_FPR_matrix_long)

p_2 = ggplot(as.data.frame(m_FPR_long), aes(x = Sample_size, y = as.numeric(Value), fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#B0B0B0', '#E29998', '#67A2A3')) +
  theme_minimal() +
  theme_bw()+
  coord_flip() +
  theme(text = element_text(size = 12), 
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mediator FPR Comparison Across Different Methods", 
       x = "Sample Size", 
       y = "FPR Value",
       fill = "Model")

####TPR
bama_TPR_matrix = cbind(t(bama_TPR_matrix[complete.cases(bama_TPR_matrix),]),Scenario,Sample_size)
bama_TPR_matrix = cbind(bama_TPR_matrix,Model = rep('bama',9))
# HIMA_TPR_matrix = cbind(t(HIMA_TPR_matrix[complete.cases(HIMA_TPR_matrix),]),Scenario,Sample_size)
mednet_TPR_matrix = cbind(t(mednet_TPR_matrix[complete.cases(mednet_TPR_matrix),]),Scenario,Sample_size)
mednet_TPR_matrix = cbind(mednet_TPR_matrix,Model = rep('medNet',9))
medfix_TPR_matrix = cbind(t(medfix_TPR_matrix[complete.cases(medfix_TPR_matrix),]),Scenario,Sample_size)
medfix_TPR_matrix = cbind(medfix_TPR_matrix,Model = rep('MedFix',9))
medfix_TPR_matrix_long = pivot_longer(as.data.frame(medfix_TPR_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
mednet_TPR_matrix_long = pivot_longer(as.data.frame(mednet_TPR_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
bama_TPR_matrix_long = pivot_longer(as.data.frame(bama_TPR_matrix), cols = starts_with("V"), 
                               names_to = "Variable", values_to = "Value")
m_TPR_long = rbind(medfix_TPR_matrix_long,mednet_TPR_matrix_long,bama_TPR_matrix_long)

p_3 = ggplot(as.data.frame(m_TPR_long), aes(x = Sample_size, y = as.numeric(Value), fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values = c('#B0B0B0', '#E29998', '#67A2A3')) +
  theme_minimal() +
  theme_bw()+
  coord_flip() +
  theme(text = element_text(size = 12), 
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Mediator TPR Comparison Across Different Methods", 
       x = "Sample Size", 
       y = "TPR Value",
       fill = "Model")

cowplot::plot_grid(p_2,p_3,p_1,ncol = 1,nrow = 3)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/benchmark_comparison.png', width = 8, height = 15)

##TPR
mednet_TPR_matrix = as.data.frame(mednet_TPR_matrix)
mednet_TPR_matrix[,1:(ncol(mednet_TPR_matrix)-3)] <- t(as.data.frame(apply(mednet_TPR_matrix[,1:(ncol(mednet_TPR_matrix)-3)],1, function(x) as.numeric(x))))
medfix_TPR_matrix = as.data.frame(medfix_TPR_matrix)
medfix_TPR_matrix[,1:(ncol(medfix_TPR_matrix)-3)] <- t(as.data.frame(apply(medfix_TPR_matrix[,1:(ncol(medfix_TPR_matrix)-3)],1, function(x) as.numeric(x))))
bama_TPR_matrix = as.data.frame(bama_TPR_matrix)
bama_TPR_matrix[,1:(ncol(bama_TPR_matrix)-3)] <- t(as.data.frame(apply(bama_TPR_matrix[,1:(ncol(bama_TPR_matrix)-3)],1, function(x) as.numeric(x))))
mean_TPR_bama = cbind(Value = rowMeans(bama_TPR_matrix[,1:(ncol(bama_TPR_matrix)-3)]),bama_TPR_matrix[,(ncol(bama_TPR_matrix)-2):ncol(bama_TPR_matrix)])
mean_TPR_medfix = cbind(Value = rowMeans(medfix_TPR_matrix[,1:(ncol(medfix_TPR_matrix)-3)]),medfix_TPR_matrix[,(ncol(medfix_TPR_matrix)-2):ncol(medfix_TPR_matrix)])
mean_TPR_mednet = cbind(Value = rowMeans(mednet_TPR_matrix[,1:(ncol(mednet_TPR_matrix)-3)]),mednet_TPR_matrix[,(ncol(mednet_TPR_matrix)-2):ncol(mednet_TPR_matrix)])
mean_TPR = rbind(mean_TPR_bama,mean_TPR_medfix,mean_TPR_mednet)
mean_TPR = cbind(mean_TPR, Type = rep('TPR',nrow(mean_TPR)))
##FPR
mednet_FPR_matrix = as.data.frame(mednet_FPR_matrix)
mednet_FPR_matrix[,1:(ncol(mednet_FPR_matrix)-3)] <- t(as.data.frame(apply(mednet_FPR_matrix[,1:(ncol(mednet_FPR_matrix)-3)],1, function(x) as.numeric(x))))
medfix_FPR_matrix = as.data.frame(medfix_FPR_matrix)
medfix_FPR_matrix[,1:(ncol(medfix_FPR_matrix)-3)] <- t(as.data.frame(apply(medfix_FPR_matrix[,1:(ncol(medfix_FPR_matrix)-3)],1, function(x) as.numeric(x))))
bama_FPR_matrix = as.data.frame(bama_FPR_matrix)
bama_FPR_matrix[,1:(ncol(bama_FPR_matrix)-3)] <- t(as.data.frame(apply(bama_FPR_matrix[,1:(ncol(bama_FPR_matrix)-3)],1, function(x) as.numeric(x))))
mean_FPR_bama = cbind(Value = rowMeans(bama_FPR_matrix[,1:(ncol(bama_FPR_matrix)-3)]),bama_FPR_matrix[,(ncol(bama_FPR_matrix)-2):ncol(bama_FPR_matrix)])
mean_FPR_medfix = cbind(Value = rowMeans(medfix_FPR_matrix[,1:(ncol(medfix_FPR_matrix)-3)]),medfix_FPR_matrix[,(ncol(medfix_FPR_matrix)-2):ncol(medfix_FPR_matrix)])
mean_FPR_mednet = cbind(Value = rowMeans(mednet_FPR_matrix[,1:(ncol(mednet_FPR_matrix)-3)]),mednet_FPR_matrix[,(ncol(mednet_FPR_matrix)-2):ncol(mednet_FPR_matrix)])
mean_FPR = rbind(mean_FPR_bama,mean_FPR_medfix,mean_FPR_mednet)
mean_FPR = cbind(mean_FPR, Type = rep('FPR',nrow(mean_FPR)))
#FDR
mednet_FDR_matrix = as.data.frame(mednet_FDR_matrix)
mednet_FDR_matrix[,1:(ncol(mednet_FDR_matrix)-3)] <- t(as.data.frame(apply(mednet_FDR_matrix[,1:(ncol(mednet_FDR_matrix)-3)],1, function(x) as.numeric(x))))
medfix_FDR_matrix = as.data.frame(medfix_FDR_matrix)
medfix_FDR_matrix[,1:(ncol(medfix_FDR_matrix)-3)] <- t(as.data.frame(apply(medfix_FDR_matrix[,1:(ncol(medfix_FDR_matrix)-3)],1, function(x) as.numeric(x))))
bama_FDR_matrix = as.data.frame(bama_FDR_matrix)
bama_FDR_matrix[,1:(ncol(bama_FDR_matrix)-3)] <- t(as.data.frame(apply(bama_FDR_matrix[,1:(ncol(bama_FDR_matrix)-3)],1, function(x) as.numeric(x))))
mean_FDR_bama = cbind(Value = rowMeans(bama_FDR_matrix[,1:(ncol(bama_FDR_matrix)-3)]),bama_FDR_matrix[,(ncol(bama_FDR_matrix)-2):ncol(bama_FDR_matrix)])
mean_FDR_medfix = cbind(Value = rowMeans(medfix_FDR_matrix[,1:(ncol(medfix_FDR_matrix)-3)]),medfix_FDR_matrix[,(ncol(medfix_FDR_matrix)-2):ncol(medfix_FDR_matrix)])
mean_FDR_mednet = cbind(Value = rowMeans(mednet_FDR_matrix[,1:(ncol(mednet_FDR_matrix)-3)]),mednet_FDR_matrix[,(ncol(mednet_FDR_matrix)-2):ncol(mednet_FDR_matrix)])
mean_FDR = rbind(mean_FDR_bama,mean_FDR_medfix,mean_FDR_mednet)
mean_FDR = cbind(mean_FDR, Type = rep('FDR',nrow(mean_FDR)))

mean_value = rbind(mean_TPR,mean_FPR,mean_FDR)

column_info <- tribble(
  ~id,     ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "Model",    "Method",             "Method",                         "text",       NA,          list(hjust = 0, width = 3),
  "Scenario",    "Experimental Setting",             "Scenario",            "text",       NA,          list(hjust = 0, width = 3),
  "Sample_size",   "Experimental Setting",      "Sample Size",           "bar",        "palette1",  lst(),
  "TPR",   "Mediator",      "TPR",           "circle",        "Score",  list(hjust = 0, width = 1.8,scale = T),
  "FPR",   "Mediator",      "FPR",           "circle",        "Score",  list(hjust = 0, width = 1.8,scale = T),
  "FDR",   "Mediator",      "FDR",           "circle",     "Score",  list(hjust = 0, width = 1.8,scale = T),
)
column_groups <- tribble(
  ~Experiment,    ~Category,                                      ~group,                   ~palette,
  "Experiment",       "\n",                                           "Method",  "overall",
  "Experiment",      "\n",             "Experimental Setting",          "overall",
  "Scaled mean value",     "Mediator",                           "Mediator",       "benchmark2",
)
palettes <- tribble(
  ~palette,             ~colours,
  "palette1",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "Score",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-8:-9]))(101),
  "overall",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "benchmark",          grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "benchmark2", grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-8:-9]))(101),
)




# Assuming mean_value is your data frame
wide_df <- pivot_wider(mean_value, 
                       names_from = Type, 
                       values_from = Value, 
                       values_fill = list(Value = NA))
wide_df$Sample_size <- as.numeric(wide_df$Sample_size)
wide_df = wide_df%>% arrange(Scenario,Model,Sample_size)
p = funky_heatmap(wide_df, column_info = column_info,palettes = palettes,column_groups = column_groups)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/benchmark_heatmap.pdf',p, device = cairo_pdf, width = 6, height = 15)











