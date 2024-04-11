rm(list = ls())
cbind_column_16  = cbind(method = c(rep('SVM',9),rep('RF',9),rep('SVM',9),rep('RF',9)),H = c(rep(12,18),rep(16,18)),
                         Scenario = rep(c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3)), 4),
                         Sample_size = rep(c(200, 300, 400), 12))
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/SVM_RF_12_16.RData')
e_FP_matrix_16 = e_FP_matrix
e_TP_matrix_16 = e_TP_matrix
m_FDR_matrix_16 = m_FDR_matrix
m_FP_matrix_16 = m_FP_matrix
m_TP_matrix_16 = m_TP_matrix
e_FP_matrix_16 = cbind(e_FP_matrix_16,cbind_column_16)
e_TP_matrix_16 = cbind(e_TP_matrix_16,cbind_column_16)
m_FDR_matrix_16 = cbind(m_FDR_matrix_16,cbind_column_16)
m_FP_matrix_16 = cbind(m_FP_matrix_16,cbind_column_16)
m_TP_matrix_16 = cbind(m_TP_matrix_16,cbind_column_16)


cbind_column_18  = cbind(method = c(rep('LR',9),rep('SVM',9),rep('RF',9)),H = rep(18,27), Scenario = rep(c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3)), 3),
                         Sample_size = rep(c(200, 300, 400), 9))
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/LR_SVM_RF_18.RData')
e_FP_matrix_18 = e_FP_matrix
e_TP_matrix_18 = e_TP_matrix
m_FDR_matrix_18 = m_FDR_matrix
m_FP_matrix_18 = m_FP_matrix
m_TP_matrix_18 = m_TP_matrix
e_FP_matrix_18 = cbind(e_FP_matrix_18,cbind_column_18)
e_TP_matrix_18 = cbind(e_TP_matrix_18,cbind_column_18)
m_FDR_matrix_18 = cbind(m_FDR_matrix_18,cbind_column_18)
m_FP_matrix_18 = cbind(m_FP_matrix_18,cbind_column_18)
m_TP_matrix_18 = cbind(m_TP_matrix_18,cbind_column_18)




cbind_column_14  = cbind(method = c(rep('LR',9),rep('SVM',9),rep('RF',9)),H = rep(14,27), Scenario = rep(c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3)), 3),
                         Sample_size = rep(c(200, 300, 400), 9))
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/all_method.RData')
e_FP_matrix_14 = e_FP_matrix
e_TP_matrix_14 = e_TP_matrix
m_FDR_matrix_14 = m_FDR_matrix
m_FP_matrix_14 = m_FP_matrix
m_TP_matrix_14 = m_TP_matrix
e_FP_matrix_14 = cbind(e_FP_matrix_14,cbind_column_14)
e_TP_matrix_14 = cbind(e_TP_matrix_14,cbind_column_14)
m_FDR_matrix_14 = cbind(m_FDR_matrix_14,cbind_column_14)
m_FP_matrix_14 = cbind(m_FP_matrix_14,cbind_column_14)
m_TP_matrix_14 = cbind(m_TP_matrix_14,cbind_column_14)


cbind_column_LR  = cbind(method = rep('LR',27),H = c(rep(12,9),rep(14,9),rep(16,9)), Scenario = rep(c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3)), 3),
                         Sample_size = rep(c(200, 300, 400), 9))
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/all_H.RData')
e_FP_matrix_LR = e_FP_matrix
e_TP_matrix_LR = e_TP_matrix
m_FDR_matrix_LR = m_FDR_matrix
m_FP_matrix_LR = m_FP_matrix
m_TP_matrix_LR = m_TP_matrix
e_FP_matrix_LR = cbind(e_FP_matrix_LR,cbind_column_LR)
e_TP_matrix_LR = cbind(e_TP_matrix_LR,cbind_column_LR)
m_FDR_matrix_LR = cbind(m_FDR_matrix_LR,cbind_column_LR)
m_FP_matrix_LR = cbind(m_FP_matrix_LR,cbind_column_LR)
m_TP_matrix_LR = cbind(m_TP_matrix_LR,cbind_column_LR)
e_FP_matrix_LR = e_FP_matrix_LR[as.numeric(e_FP_matrix_LR[,22])!=14,]
e_TP_matrix_LR = e_TP_matrix_LR[as.numeric(e_TP_matrix_LR[,22])!=14,]
m_FDR_matrix_LR = m_FDR_matrix_LR[as.numeric(m_FDR_matrix_LR[,22])!=14,]
m_FP_matrix_LR = m_FP_matrix_LR[as.numeric(m_FP_matrix_LR[,22])!=14,]
m_TP_matrix_LR = m_TP_matrix_LR[as.numeric(m_TP_matrix_LR[,22])!=14,]

e_FP_matrix = rbind(e_FP_matrix_14,e_FP_matrix_16,e_FP_matrix_LR,e_FP_matrix_18)
e_TP_matrix = rbind(e_TP_matrix_14,e_TP_matrix_16,e_TP_matrix_LR,e_TP_matrix_18)
m_FDR_matrix = rbind(m_FDR_matrix_14,m_FDR_matrix_16,m_FDR_matrix_LR,m_FDR_matrix_18)
m_FP_matrix = rbind(m_FP_matrix_14,m_FP_matrix_16,m_FP_matrix_LR,m_FP_matrix_18)
m_TP_matrix = rbind(m_TP_matrix_14,m_TP_matrix_16,m_TP_matrix_LR,m_TP_matrix_18)




library(tidyr)
library(ggplot2)
library(dplyr)

m_FP_matrix_long <- pivot_longer(as.data.frame(m_FP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
m_FP_matrix_long = m_FP_matrix_long[m_FP_matrix_long$Sample_size==400,]
p_1 = ggplot(as.data.frame(m_FP_matrix_long), aes(x = H, y = as.numeric(Value), fill = method)) +
  geom_boxplot() +
  facet_wrap(~Scenario, strip.position = "bottom") +
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


m_TP_matrix_long <- pivot_longer(as.data.frame(m_TP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
m_TP_matrix_long = m_TP_matrix_long[m_TP_matrix_long$Sample_size==400,]

p_2 = ggplot(as.data.frame(m_TP_matrix_long), aes(x = H, y = as.numeric(Value), fill = method)) +
  geom_boxplot() +
  facet_wrap(~Scenario, strip.position = "bottom") +
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
  labs(#title = "Mediator TPR Comparison Across Models", 
       x = "Sample Size", 
       y = "TPR Value",
       fill = "Model")


m_FDR_matrix_long <- pivot_longer(as.data.frame(m_FDR_matrix), cols = starts_with("V"), 
                                  names_to = "Variable", values_to = "Value")
m_FDR_matrix_long = m_FDR_matrix_long[m_FDR_matrix_long$Sample_size==400,]

p_3 = ggplot(as.data.frame(m_FDR_matrix_long), aes(x = H, y = as.numeric(Value), fill = method)) +
  geom_boxplot() +
  facet_wrap(~Scenario, strip.position = "bottom") +
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
  labs(#title = "Mediator FDR Comparison Across Models", 
       x = "Sample Size", 
       y = "FDR Value",
       fill = "Model")

cowplot::plot_grid(p_1,p_2,p_3,ncol = 1,nrow = 3)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/model_comparison_all.png', width = 15, height = 15)













