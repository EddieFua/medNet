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
  "id",    "",             "",                         "text",       NA,          list(hjust = 0, width = 3),
  "Scenario",    "",             "Scenario",            "text",       NA,          list(hjust = 0, width = 2),
  "Sample_size",   "",      "Sample Size",           "bar",        "palette1",  lst(),
  "e_TP",   "Exposure",      "TPR",           "bar",        "palette2",  lst(),
  "e_FP",   "Exposure",      "FPR",           "bar",        "palette2",  lst(),
  "m_TP",   "Mediator",      "TPR",           "circle",        "palette3",  lst(),
  "m_FP",   "Mediator",      "FPR",           "circle",        "palette3",  lst(),
  "m_FDR",   "Mediator",      "FDR",           "circle",        "palette3",  lst(),
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
new_dat = cbind(id, new_dat)
new_dat = new_dat %>% mutate(Sample_size = as.numeric(Sample_size)) # Convert Sample_size to numeric

row_info = new_dat%>%select(Scenario,id)
colnames(row_info) = c('group','id')
row_groups = tribble(~group,~Group,
                        'Scenario 1','Scenario 1',
                        'Scenario 2','Scenario 2',
                        'Scenario 3','Scenario 3')
palettes <- tribble(
  ~palette,             ~colours,
  "palette1",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101),
  "palette2",          grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636")))(101),
  "palette3",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Reds")[-8:-9]))(101),
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
row_groups
funky_heatmap(new_dat, column_info = column_info, expand = list(xmax = 4),row_info = row_info,row_groups = row_groups,palettes = palettes)

funky_heatmap(new_dat, column_info = column_info,palettes = palettes)


library(tidyr)
library(ggplot2)
library(dplyr)
model = c(rep('LR',9),rep('SVM',9),rep('RF',9))
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
  labs(title = "Mediator FPR Comparison Across Models", 
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
  labs(title = "Mediator TPR Comparison Across Models", 
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
  labs(title = "Mediator FDR Comparison Across Models", 
       x = "Sample Size", 
       y = "FDR Value",
       fill = "Model")

cowplot::plot_grid(p_1,p_2,p_3,ncol = 1,nrow = 3)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/model_comparison.png', width = 8, height = 15)







