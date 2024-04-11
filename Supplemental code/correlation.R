rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/breast cancer res/cleandata.RData')
corr_real1 = c()
for (i in 1:ncol(dat$mediator)){
  corr_real1 = c(corr_real1,cor(dat$response,dat$mediator[,i]))
}

load('/Users/fuyinghao/Documents/Mediator analysis/Second data/clean data.RData')
corr_real2 = c()
for (i in 1:ncol(mediator)){
  corr_real2 = c(corr_real2,cor(response,mediator[,i]))
}


corr_simu_1 = c()
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/simu_1.RData')
idx1 = which(data[[2]]$exposure_mediator==1,arr.ind = T)[,2]
simu_data = data[[2]]$mediator[,idx1]
corr_simu_1_idx = c()
for (i in 1:ncol(simu_data)){
  corr_simu_1_idx = c(corr_simu_1_idx,lm(data[[2]]$response~simu_data[,i])$coefficients[2])
}

simu_data = data[[2]]$mediator
corr_simu_1 = c()
for (i in 1:ncol(simu_data)){
  corr_simu_1 = c(corr_simu_1,lm(data[[2]]$response~simu_data[,i])$coefficients[2])
}
corr_simu_2 = c()
idx2 = which(data[[5]]$exposure_mediator==1,arr.ind = T)[,2]
simu_data = data[[5]]$mediator[,idx2]
corr_simu_2_idx = c()
for (i in 1:ncol(simu_data)){
  corr_simu_2_idx = c(corr_simu_2_idx,lm(data[[5]]$response~simu_data[,i])$coefficients[2])
}
simu_data = data[[5]]$mediator
corr_simu_2 = c()
for (i in 1:ncol(simu_data)){
  corr_simu_2 = c(corr_simu_2,lm(data[[5]]$response~simu_data[,i])$coefficients[2])
}
corr_simu_3 = c()
idx3 = which(data[[8]]$exposure_mediator==1,arr.ind = T)[,2]
simu_data = data[[8]]$mediator[,idx3]
corr_simu_3_idx = c()
for (i in 1:ncol(simu_data)){
  corr_simu_3_idx = c(corr_simu_3_idx,lm(data[[8]]$response~simu_data[,i])$coefficients[2])
}
simu_data = data[[8]]$mediator
corr_simu_3 = c()
for (i in 1:ncol(simu_data)){
  corr_simu_3 = c(corr_simu_3,lm(data[[8]]$response~simu_data[,i])$coefficients[2])
}


max_length <- max(length(corr_real1), length(corr_real2), length(corr_simu_3))
dat <- data.frame(
  METABRIC = c(corr_real1, rep(NA, max_length - length(corr_real1))),
  Metabolomics = c(corr_real2, rep(NA, max_length - length(corr_real2))),
  Scenario1 = c(corr_simu_1, rep(NA, max_length - length(corr_simu_1))),
  Scenario2 = c(corr_simu_2, rep(NA, max_length - length(corr_simu_2))),
  Scenario3 = c(corr_simu_3, rep(NA, max_length - length(corr_simu_3)))
)

max_length <- max(length(corr_simu_1_idx), length(corr_simu_2_idx), length(corr_simu_3_idx))
dat_point = data.frame(
  Scenario1 = c(corr_simu_1_idx, rep(NA, max_length - length(corr_simu_1_idx))),
  Scenario2 = c(corr_simu_2_idx, rep(NA, max_length - length(corr_simu_2_idx))),
  Scenario3 = c(corr_simu_3_idx, rep(NA, max_length - length(corr_simu_3_idx)))
)
library(reshape2)
long_dat <- melt(dat, variable.name = "Dataset", value.name = "Correlation", na.rm = TRUE)
long_dat_point <- melt(dat_point, variable.name = "Dataset", value.name = "Correlation", na.rm = TRUE)


p <- ggplot(long_dat, aes(x = Dataset, y = Correlation, fill = Dataset)) +
  geom_boxplot(outlier.shape = 20, outlier.size = 2, outlier.stroke = 0.5) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.5) +
  # geom_point(data = long_dat_point, aes(x = Dataset, y = Correlation), color = "red", size = 3) +
  scale_fill_manual(values = c("#5D9AD3", "#C4121A", "#DE722A",'#4B3D7C','7F7F7F')) +
  labs(title = "Correlation Between Mediators and Response",
       subtitle = "Comparative Analysis Across Datasets",
       x = "",
       y = "Correlation Coefficient") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.line = element_line(size = 0.5, color = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p = ggplot(long_dat, aes(x = Dataset, y = Correlation, fill = Dataset)) +
  geom_boxplot(outlier.shape = NA, outlier.size = 2, outlier.stroke = 0.5) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("#5D9AD3", "#C4121A", "#DE722A",'#4B3D7C','7F7F7F')) +
  labs(title = "Correlation Between Mediators and Response",
       subtitle = "Comparative Analysis Across Datasets",
       x = "",
       y = "Correlation Coefficient") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.line = element_line(size = 0.5, color = "black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )



# Save the plot
ggsave("/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/correlation_comparison.pdf", p, width = 10, height = 10)


