rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/bama_nonlinear.RData')
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/nonlinear.RData')

mednet_TPR_matrix = matrix(ncol = 20, nrow = 18)
mednet_FPR_matrix = matrix(ncol = 20, nrow = 18)
mednet_FDR_matrix = matrix(ncol = 20, nrow = 18)
for (i in 1:20) {
  mednet_TPR_matrix[,i] = res[[i]][[3]]
  mednet_FPR_matrix[,i] = res[[i]][[4]]
  mednet_FDR_matrix[,i] = res[[i]][[5]]
}

cbind_col = cbind(method = c(rep('SVM',9),rep('RF',9)),
                  Scenario = rep(c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3)), 2),
                  Sample_size = rep(c(200, 300, 400), 6))
mednet_TPR_matrix = cbind(mednet_TPR_matrix,cbind_col)
mednet_FPR_matrix = cbind(mednet_FPR_matrix,cbind_col)
mednet_FDR_matrix = cbind(mednet_FDR_matrix,cbind_col)
cbind_col = cbind(method = rep('bama',9),
                  Scenario = rep(c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3)), 1),Sample_size = rep(c(200, 300, 400), 3))
bama_TPR_matrix = cbind(t(bama_TPR_matrix),cbind_col)
bama_FPR_matrix = cbind(t(bama_FPR_matrix),cbind_col)
bama_FDR_matrix = cbind(t(bama_FDR_matrix),cbind_col)

m_FP_matrix = rbind(mednet_FPR_matrix,bama_FPR_matrix)
m_TP_matrix = rbind(mednet_TPR_matrix,bama_TPR_matrix)
m_FDR_matrix = rbind(mednet_FDR_matrix,bama_FDR_matrix)
library(tidyr)
library(ggplot2)
library(dplyr)
base_size= 14
m_FP_matrix_long <- pivot_longer(as.data.frame(m_FP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
p_1 = ggplot(as.data.frame(m_FP_matrix_long), aes(x =Sample_size , y = as.numeric(Value), fill = method)) +
  geom_boxplot() +
  facet_wrap(~Scenario, strip.position = "bottom") +
  scale_fill_manual(values = c('#B0B0B0', '#E29998', '#67A2A3')) +
  theme_minimal() +
  theme_bw()+
  coord_flip() +
  theme(plot.title = element_text(size = base_size * 1.2, face = "bold", hjust = 0.5),
        axis.title = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size),
        # legend.position = "",
        strip.text.x = element_text(size = base_size, face = "bold"),
        axis.text.x = element_text(vjust = 0.5, size = base_size))+
  labs(#title = "Performance comparison between medNet and bama in nonlinear environment", 
       x = "Sample Size", 
       y = "FPR Value",
       fill = "Model")+labs(fill = NULL)+
  theme(strip.text.x = element_text(size = base_size * 1, face = "bold"))


m_TP_matrix_long <- pivot_longer(as.data.frame(m_TP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")


p_2 = ggplot(as.data.frame(m_TP_matrix_long), aes(x =Sample_size , y = as.numeric(Value), fill = method)) +
  geom_boxplot() +
  facet_wrap(~Scenario, strip.position = "bottom") +
  scale_fill_manual(values = c('#B0B0B0', '#E29998', '#67A2A3')) +
  theme_minimal() +
  theme_bw()+
  coord_flip() +
  theme(plot.title = element_text(size = base_size * 1.2, face = "bold", hjust = 0.5),
        axis.title = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size),
        # legend.position = "",
        strip.text.x = element_text(size = base_size, face = "bold"),
        axis.text.x = element_text(vjust = 0.5, size = base_size))+
  labs(#title = "Mediator TPR Comparison Across Models", 
    x = "Sample Size", 
    y = "TPR Value",
    fill = "Model")+labs(fill = NULL)+
  theme(strip.text.x = element_text(size = base_size * 1, face = "bold"))


m_FDR_matrix_long <- pivot_longer(as.data.frame(m_FDR_matrix), cols = starts_with("V"), 
                                  names_to = "Variable", values_to = "Value")


p_3 = ggplot(as.data.frame(m_FDR_matrix_long), aes(x =Sample_size , y = as.numeric(Value), fill = method)) +
  geom_boxplot() +
  facet_wrap(~Scenario, strip.position = "bottom") +
  scale_fill_manual(values = c('#B0B0B0', '#E29998', '#67A2A3')) +
  theme_minimal() +
  theme_bw()+
  coord_flip() +
  theme(plot.title = element_text(size = base_size * 1.2, face = "bold", hjust = 0.5),
        axis.title = element_text(size = base_size, face = "bold"),
        axis.text = element_text(size = base_size),
        legend.title = element_blank(),
        legend.text = element_text(size = base_size),
        # legend.position = "",
        strip.text.x = element_text(size = base_size, face = "bold"),
        axis.text.x = element_text(vjust = 0.5, size = base_size))+
  labs(#title = "Mediator FDR Comparison Across Models", 
    x = "Sample Size", 
    y = "FDR Value",
    fill = "Model")+labs(fill = NULL)+
  theme(strip.text.x = element_text(size = base_size * 1, face = "bold"))

cowplot::plot_grid(p_1,p_2,p_3,ncol = 1,nrow = 3)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/nonlinear_comparison_all.png', width = 15, height = 15)













