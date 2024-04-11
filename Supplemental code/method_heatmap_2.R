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


e_TP_matrix_long <- pivot_longer(as.data.frame(e_TP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
e_TP_matrix_long = e_TP_matrix_long[e_TP_matrix_long$H==18,]
e_FP_matrix_long <- pivot_longer(as.data.frame(e_FP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
e_FP_matrix_long = e_FP_matrix_long[e_FP_matrix_long$H==18,]

m_TP_matrix_long <- pivot_longer(as.data.frame(m_TP_matrix), cols = starts_with("V"), 
                                 names_to = "Variable", values_to = "Value")
m_TP_matrix_long = m_TP_matrix_long[m_TP_matrix_long$H==18,]
m_FP_matrix_long <- pivot_longer(as.data.frame(m_FP_matrix), cols = starts_with("V"), 
                                  names_to = "Variable", values_to = "Value")
m_FP_matrix_long = m_FP_matrix_long[m_FP_matrix_long$H==18,]
m_FDR_matrix_long <- pivot_longer(as.data.frame(m_FDR_matrix), cols = starts_with("V"), 
                                  names_to = "Variable", values_to = "Value")
m_FDR_matrix_long = m_FDR_matrix_long[m_FDR_matrix_long$H==18,]

p_1 = ggplot(as.data.frame(m_FP_matrix_long), aes(x = Sample_size, y = as.numeric(Value), fill = method)) +
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
  labs(#title = "Mediator FPR Comparison Across Models", 
       x = "Sample Size", 
       y = "FPR Value",
       fill = "Model")


p_2 = ggplot(as.data.frame(m_TP_matrix_long), aes(x = Sample_size, y = as.numeric(Value), fill = method)) +
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
  labs(#title = "Mediator TPR Comparison Across Models", 
       x = "Sample Size", 
       y = "TPR Value",
       fill = "Model")



p_3 = ggplot(as.data.frame(m_FDR_matrix_long), aes(x = Sample_size, y = as.numeric(Value), fill = method)) +
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
  labs(#title = "Mediator FDR Comparison Across Models", 
       x = "Sample Size", 
       y = "FDR Value",
       fill = "Model")

cowplot::plot_grid(p_1,p_2,p_3,ncol = 1,nrow = 3)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/model_comparison.png', width = 8, height = 15)


m_TP_matrix_long = m_TP_matrix_long[m_TP_matrix_long$H==18,]
m_FP_matrix_long = m_FP_matrix_long[m_FP_matrix_long$H==18,]
m_FDR_matrix_long = m_FDR_matrix_long[m_FDR_matrix_long$H==18,]

library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
Scenario = rep(c(rep('Scenario 1', 3), rep('Scenario 2', 3), rep('Scenario 3', 3)), 3)
Sample_size = rep(c(200, 300, 400), 9)
Sample_size = as.numeric(Sample_size)
id = c(c('LR_200_1', 'LR_300_1','LR_400_1','LR_200_2', 'LR_300_2','LR_400_2','LR_200_3', 'LR_300_3','LR_400_3'),
       c('SVM_200_1', 'SVM_300_1','SVM_400_1','SVM_200_2', 'SVM_300_2','SVM_400_2','SVM_200_3', 'SVM_300_3','SVM_400_3'),
       c('RF_200_1', 'RF_300_1','RF_400_1','RF_200_2', 'RF_300_2','RF_400_2','RF_200_3', 'RF_300_3','RF_400_3'))
new_dat = as.data.frame(cbind(id, Scenario,Sample_size))
method = sub("^(.*?)_.*", "\\1", id)

e_TP = c()
for (i in 1:nrow(new_dat)) {
  idx = e_TP_matrix_long$method==method[i]&e_TP_matrix_long$Scenario==new_dat$Scenario[i]&e_TP_matrix_long$Sample_size==new_dat$Sample_size[i]
  e_TP = c(e_TP,mean(as.numeric(e_TP_matrix_long$Value[idx])))
}
e_FP = c()
for (i in 1:nrow(new_dat)) {
  idx = e_FP_matrix_long$method==method[i]&e_FP_matrix_long$Scenario==new_dat$Scenario[i]&e_FP_matrix_long$Sample_size==new_dat$Sample_size[i]
  e_FP = c(e_FP,mean(as.numeric(e_FP_matrix_long$Value[idx])))
}
m_TP = c()
for (i in 1:nrow(new_dat)) {
  idx = m_TP_matrix_long$method==method[i]&m_TP_matrix_long$Scenario==new_dat$Scenario[i]&m_TP_matrix_long$Sample_size==new_dat$Sample_size[i]
  m_TP = c(m_TP,mean(as.numeric(m_TP_matrix_long$Value[idx])))
}
m_FP = c()
for (i in 1:nrow(new_dat)) {
  idx = m_FP_matrix_long$method==method[i]&m_FP_matrix_long$Scenario==new_dat$Scenario[i]&m_FP_matrix_long$Sample_size==new_dat$Sample_size[i]
  m_FP = c(m_FP,mean(as.numeric(m_FP_matrix_long$Value[idx])))
}
m_FDR = c()
for (i in 1:nrow(new_dat)) {
  idx = m_FDR_matrix_long$method==method[i]&m_FDR_matrix_long$Scenario==new_dat$Scenario[i]&m_FDR_matrix_long$Sample_size==new_dat$Sample_size[i]
  m_FDR = c(m_FDR,mean(as.numeric(m_FDR_matrix_long$Value[idx])))
}

new_dat = cbind(new_dat,e_TP,e_FP,m_TP,m_FP,m_FDR)

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
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/main_paper_heatmap.pdf',p, device = cairo_pdf, width = 8, height = 15)









