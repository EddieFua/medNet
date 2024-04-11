rm(list = ls())
load('/Users/fuyinghao/Documents/Mediator analysis/Revision/res/Metablomic/ten_times.RData')
intersect_pathway_num = matrix(NA,nrow = 3,ncol = 4)
intersect_meta_num = matrix(NA,nrow = 3,ncol = 4)
intersect_pathway_num[1,] = colMeans(idx1_matrix)
intersect_pathway_num[2,] = colMeans(idx2_matrix)
intersect_pathway_num[3,] = colMeans(idx3_matrix)
intersect_meta_num[1,] = colMeans(idx1_matrix_meta)
intersect_meta_num[2,] = colMeans(idx2_matrix_meta)
intersect_meta_num[3,] = colMeans(idx3_matrix_meta)


disturbed_rate = seq(0.2,0.8,0.2)

data = data.frame(disturbed_rate, intersect_meta_num = intersect_meta_num[1,])

# metarate the line plot for metas
meta_plot_1 <- ggplot(data, aes(x = disturbed_rate, y = intersect_meta_num)) +
  geom_line(color = 'blue') + # Line color
  geom_point(color = 'red') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection meta Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Metabolites")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


# Assuming 'intersect_pathway_num' is another vector corresponding to 'disturbed_rate'
# Create a data frame for pathway data
pathway_data <- data.frame(disturbed_rate, intersect_pathway_num = intersect_pathway_num[1,])

# metarate the line plot for pathways
pathway_plot_1 <- ggplot(pathway_data, aes(x = disturbed_rate, y = intersect_pathway_num)) +
  geom_line(color = 'green') + # Line color
  geom_point(color = 'orange') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection Pathway Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Pathways")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


data = data.frame(disturbed_rate, intersect_meta_num = intersect_meta_num[2,])

# metarate the line plot for metas
meta_plot_2 <- ggplot(data, aes(x = disturbed_rate, y = intersect_meta_num)) +
  geom_line(color = 'blue') + # Line color
  geom_point(color = 'red') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection meta Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Metabolites")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


# Assuming 'intersect_pathway_num' is another vector corresponding to 'disturbed_rate'
# Create a data frame for pathway data
pathway_data <- data.frame(disturbed_rate, intersect_pathway_num = intersect_pathway_num[2,])

# metarate the line plot for pathways
pathway_plot_2 <- ggplot(pathway_data, aes(x = disturbed_rate, y = intersect_pathway_num)) +
  geom_line(color = 'green') + # Line color
  geom_point(color = 'orange') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection Pathway Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Pathways")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


data = data.frame(disturbed_rate, intersect_meta_num = intersect_meta_num[3,])

# metarate the line plot for metas
meta_plot_3 <- ggplot(data, aes(x = disturbed_rate, y = intersect_meta_num)) +
  geom_line(color = 'blue') + # Line color
  geom_point(color = 'red') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection meta Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Metabolites")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )


# Assuming 'intersect_pathway_num' is another vector corresponding to 'disturbed_rate'
# Create a data frame for pathway data
pathway_data <- data.frame(disturbed_rate, intersect_pathway_num = intersect_pathway_num[3,])

# metarate the line plot for pathways
pathway_plot_3 <- ggplot(pathway_data, aes(x = disturbed_rate, y = intersect_pathway_num)) +
  geom_line(color = 'green') + # Line color
  geom_point(color = 'orange') + # Point color
  theme_minimal() + # Minimal theme similar to Nature style
  theme_bw()+
  labs(
    # title = "Relationship Between Disturbed Rate and Intersection Pathway Count",
    x = "Disturbed Rate",
    y = paste0("Number of Intersecting Pathways")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    axis.text = element_text(size = 12), # Adjust text size
    axis.title = element_text(size = 14), # Adjust title size
    legend.position = "none" # Hide legend if not needed
  )

cowplot::plot_grid(pathway_plot_1,meta_plot_1,pathway_plot_2,meta_plot_2,pathway_plot_3,meta_plot_3,ncol = 2,nrow = 3)
ggsave('/Users/fuyinghao/Documents/Mediator analysis/Revision/figure/disturbed_network_meta.png', width = 12, height = 14)
