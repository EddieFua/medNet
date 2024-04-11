rm(list = ls())
library(funkyheatmap)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)

load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/all_method.RData')
result_matrix_2 = matrix(NA, ncol = 5, nrow = 27)
result_matrix_2[, 1] = rowMeans(e_TP_matrix, na.rm = T)
result_matrix_2[, 2] = rowMeans(e_FP_matrix, na.rm = T)
result_matrix_2[, 3] = rowMeans(m_TP_matrix, na.rm = T)
result_matrix_2[, 4] = rowMeans(m_FP_matrix, na.rm = T)
result_matrix_2[, 5] = rowMeans(m_FDR_matrix, na.rm = T)

new_dat = data.frame(t(result_matrix_2))

column_info <- tribble(
  ~id,     ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "id",    "Method",             "\n",                         "text",       NA,          list(hjust = 0, width = 3),
  "Scenario",    "Experimental Setting",             "Scenario",            "text",       NA,          list(hjust = 0, width = 3),
  "Sample_size",   "Experimental Setting",      "Sample Size",           "bar",        "palette1",  lst(),
  "e_TP",   "Exposure",      "TPR",           "bar",        "palette2",  list(hjust = 0, width = 2,scale = F),
  "e_FP",   "Exposure",      "FPR",           "bar",        "palette2",  list(hjust = 0, width = 2,scale = T),
  "m_TP",   "Mediator",      "TPR",           "circle",        "Score",  list(hjust = 0, width = 1.8,scale = T),
  "m_FP",   "Mediator",      "FPR",           "circle",        "Score",  list(hjust = 0, width = 1.8,scale = T),
  "m_FDR",   "Mediator",      "FDR",           "circle",        "Score",  list(hjust = 0, width = 1.8,scale = T),
)
column_groups <- tribble(
  ~Experiment,    ~Category,                                      ~group,                   ~palette,
  "Experiment",       "\n",                                           "Method",  "overall",
  "Experiment",      "\n",             "Experimental Setting",          "overall",
  "Mean value",     "Exposure",                                   "Exposure",       "benchmark",
  "Scaled mean value",     "Mediator",                           "Mediator",       "benchmark2",
)
new_dat = result_matrix_2
colnames(new_dat) = c('e_TP', 'e_FP', 'm_TP', 'm_FP', 'm_FDR')

Scenario = rep(c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3)), 3)
Sample_size = rep(c(200, 300, 400), 9)
Sample_size = as.numeric(Sample_size) # Ensure Sample_size is numeric

new_dat = as.data.frame(cbind(new_dat, Scenario, Sample_size))
id = c(c('LR_200_1', 'LR_300_1','LR_400_1','LR_200_2', 'LR_300_2','LR_400_2','LR_200_3', 'LR_300_3','LR_400_3'),
       c('SVM_200_1', 'SVM_300_1','SVM_400_1','SVM_200_2', 'SVM_300_2','SVM_400_2','SVM_200_3', 'SVM_300_3','SVM_400_3'),
       c('RF_200_1', 'RF_300_1','RF_400_1','RF_200_2', 'RF_300_2','RF_400_2','RF_200_3', 'RF_300_3','RF_400_3'))
# id = c(c('LR_1_200_0.6', 'LR_1_300_0.6','LR_1_400_0.6','LR_2_200_0.6', 'LR_2_300_0.6','LR_2_400_0.6','LR_3_200_0.6', 'LR_3_300_0.6','LR_3_400_0.6'),
#        c('LR_1_200_0.65', 'LR_1_300_0.65','LR_1_400_0.65','LR_2_200_0.65', 'LR_2_300_0.65','LR_2_400_0.65','LR_3_200_0.65', 'LR_3_300_0.65','LR_3_400_0.65'),
#        c('LR_1_200_0.7', 'LR_1_300_0.7','LR_1_400_0.7','LR_2_200_0.7', 'LR_2_300_0.7','LR_2_400_0.7','LR_3_200_0.7', 'LR_3_300_0.7','LR_3_400_0.7'))
new_dat = cbind(id, new_dat)
new_dat = new_dat %>% mutate(Sample_size = as.numeric(Sample_size)) # Convert Sample_size to numeric
library(dplyr)
row_info = new_dat%>%dplyr::select(Scenario,id)
colnames(row_info) = c('group','id')
row_groups = tribble(~group,~Group,
                     'Scenario 1','Scenario 1',
                     'Scenario 2','Scenario 2',
                     'Scenario 3','Scenario 3')
palettes <- tribble(
  ~palette,             ~colours,
  "palette1",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "palette2",          grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "Score",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-8:-9]))(101),
  "overall",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "benchmark",          grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "benchmark2", grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-8:-9]))(101),
)
new_dat <- new_dat %>%
  mutate(
    e_TP = as.numeric(e_TP), 
    e_FP = as.numeric(e_FP), 
    m_TP = as.numeric(m_TP), 
    m_FP = as.numeric(m_FP), 
    m_FDR = as.numeric(m_FDR), 
    Scenario = as.factor(Scenario),
    Sample_size = as.numeric(Sample_size), 
  )
new_dat = new_dat%>% arrange(Scenario,id)
row_info = row_info%>% arrange(group,id)

p = funky_heatmap(new_dat, column_info = column_info, expand = list(xmax = 4),row_info = row_info,row_groups = row_groups,palettes = palettes,column_groups = column_groups,
                  col_annot_offset = 3.2)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/method_heatmap.pdf',p, device = cairo_pdf, width = 8, height = 15)


library(tidyr)
library(ggplot2)
library(dplyr)
model = c(rep('LR',9),rep('SVM',9),rep('RF',9))
# model = c(rep('T = 0.6',9),rep('T = 0.65',9),rep('T = 0.7',9))
m_FP_matrix = cbind(m_FP_matrix,Sample_size,model)
m_FP_matrix_long <- pivot_longer(as.data.frame(m_FP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")


p_1 = ggplot(as.data.frame(m_FP_matrix_long), aes(x = Sample_size, y = as.numeric(Value), fill = model)) +
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
  labs(title = "Mediator FPR Comparison Across T", 
       x = "Sample Size", 
       y = "FPR Value",
       fill = "Model")


m_TP_matrix = cbind(m_TP_matrix,Sample_size,model)
m_TP_matrix_long <- pivot_longer(as.data.frame(m_TP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")


p_2 = ggplot(as.data.frame(m_TP_matrix_long), aes(x = Sample_size, y = as.numeric(Value), fill = model)) +
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
  labs(title = "Mediator TPR Comparison Across T", 
       x = "Sample Size", 
       y = "TPR Value",
       fill = "Model")


m_FDR_matrix = cbind(m_FDR_matrix,Sample_size,model)
m_FDR_matrix_long <- pivot_longer(as.data.frame(m_FDR_matrix), cols = starts_with("V"), 
                                  names_to = "Variable", values_to = "Value")


p_3 = ggplot(as.data.frame(m_FDR_matrix_long), aes(x = Sample_size, y = as.numeric(Value), fill = model)) +
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
  labs(title = "Mediator FDR Comparison Across T", 
       x = "Sample Size", 
       y = "FDR Value",
       fill = "Model")

cowplot::plot_grid(p_1,p_2,p_3,ncol = 1,nrow = 3)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/method_comparison.png', width = 8, height = 15)







