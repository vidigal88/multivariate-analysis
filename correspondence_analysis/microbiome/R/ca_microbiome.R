
# Load libraries ----------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(ca)
library(ggrepel)

# Load data ---------------------------------------------------------------

#' Data with the counts of the 182 species in the 26 samples

load("correspondence_analysis/microbiome/data/GP.OTU.RData")

#' Data with information about the type of sample

load("correspondence_analysis/microbiome/data/GP.SampleData.RData")

#' Data with the Phylum and Class labels of the 182 species. The rows of this
#' matrix correspond to the rows in the GP.OTU matrix

load("correspondence_analysis/microbiome/data/GP.TaxData.RData")


# Transform species/sample data into a long format ------------------------

specie_sample <- GP.OTU %>% 
    rownames_to_column(var = "specie") %>%
    pivot_longer(
        -specie, 
        names_to = "sample", 
        values_to = "count"
    )


# Transform sample data into a long format --------------------------------

sample_data <- GP.SampleData %>% 
    data.frame() %>% 
    rownames_to_column(var = "sample") %>% 
    rename(sample_type = SampleType)


# Transform taxonomy data into a long format ------------------------------

tax_data <- GP.TaxData %>% 
    rownames_to_column(var = "specie") %>% 
    rename(phylum = Phylum,
           class = Class)


# Calculate total by phylum/sample type -----------------------------------

phylum_sample_type_data <- specie_sample %>%
    left_join(sample_data, by = "sample") %>% 
    left_join(tax_data, by = "specie") %>% 
    mutate(sample_type = sample_type %>% as.character()) %>% 
    select(phylum, sample_type, count) %>% 
    group_by(phylum, sample_type) %>% 
    summarise(total_count = sum(count)) %>%
    ungroup() %>% 
    pivot_wider(names_from = sample_type, 
                values_from = total_count) %>% 
    column_to_rownames(var = 'phylum')


# Calculate total by class/sample type ------------------------------------

class_sample_type_data <- specie_sample %>%
    left_join(sample_data, by = "sample") %>% 
    left_join(tax_data, by = "specie") %>% 
    mutate(sample_type = sample_type %>% as.character()) %>% 
    select(class, sample_type, count) %>% 
    group_by(class, sample_type) %>% 
    summarise(total_count = sum(count)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = sample_type, 
                values_from = total_count) %>% 
    column_to_rownames(var = 'class')


# Bind class and phylum datasets ------------------------------------------

class_philum_sample <- rbind(class_sample_type_data, 
                             phylum_sample_type_data)


# Calculate sum of human/non-human samples for each class/phylum ----------

class_philum_sample <- class_philum_sample %>% 
    mutate(
        Human = Feces + Mock  + Skin + Tongue,
        `Non-Human` = Freshwater + `Freshwater (creek)` + Ocean + `Sediment (estuary)` + Skin  + Soil
    )


# CA ----------------------------------------------------------------------

#' Phylum and Human/Non-Human as supplementary variables as they are
#' additional info.

ca_summary <- summary(ca(class_philum_sample, suprow = 17:21, supcol = 10:11))

# Default CA map ----------------------------------------------------------

ca_plot_dim <- plot(ca(class_philum_sample, suprow = 17:21, supcol = 10:11))

# CA map with ggplot2 -----------------------------------------------------

ca_dimensions <- data.frame(dim_1 = c(ca_plot_dim$rows[,1],
                                      ca_plot_dim$cols[,1]),
                            dim_2 = c(ca_plot_dim$rows[,2],
                                      ca_plot_dim$cols[,2]),
                            variable = c(rep('Class', 16),
                                         rep('Phylum', 5),
                                         rep('Sample Type', 9),
                                         rep('Human/Non-Human', 2))) %>% 
    rownames_to_column(var = 'label') %>% 
    mutate(variable = variable %>% as_factor(),
           variable = fct_relevel(variable))


# CA map with ggplot2 -----------------------------------------------------

ca_map_ggplot <- ca_dimensions %>%
    ggplot(aes(x = dim_1, y = dim_2, 
               colour = variable, 
               group = variable)) +
    geom_vline(xintercept = 0, lty = "dashed", alpha = .5) +
    geom_hline(yintercept = 0, lty = "dashed", alpha = .5) +
    geom_point(size = 2.5, shape = 16) +
    scale_color_manual(values = c("#63872c", "#244a16", "#edb62a", "#c6412a")) +
    geom_text_repel(aes(label = label, size = 24), show.legend = F, fontface = 'bold') + 
    labs(
        x = paste0("Dimension 1 - ", round(ca_summary$scree[1, 3], 1), "%"), 
        y = paste0("Dimension 2 - ", round(ca_summary$scree[2, 3], 1), "%")) +
    theme_bw() +
    theme(plot.title = element_text(size = 16),
          plot.subtitle = element_text(size = 14),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          panel.background = element_rect(fill = "#fefaf2", color = "#fefaf2"),
          plot.background = element_rect(fill = "#fefaf2", color = "#fefaf2"),
          legend.background = element_rect(fill = "#fefaf2", color = "#fefaf2")
    )

ggsave("correspondence_analysis/microbiome/plot/ca_map_ggplot.png", 
       ca_map_ggplot, 
       device = ragg::agg_png, res = 400, units = "in", h = 8, w = 12)


