library(data.table)
library(tidyverse)
library(dbscan)
library(scales)
library(ggthemes)
library(gghighlight)
library(vegan)
library(viridis)
library(cowplot)

#
# C. difficile plot
#

colours <- c("#ff7f00","#cab2d6","#6a3d9a","#fdbf6f","#ffff99","#1f78b4","#b15928","#33a02c","#e31a1c","#a6cee3","#fb9a99","#b2df8a")

# https://figshare.com/ndownloader/files/30449916
kraken <- read.table("File1_full_krakenbracken.txt", header = T, sep = "\t", quote = "", comment.char = "")

perc_cols <- seq(3, 101, 2)

kraken_abundance <- kraken[,perc_cols]
kraken_abundance[is.na(kraken_abundance)] <- 0
simpsons_diversity <- diversity(kraken_abundance, index = "simpson")

cdiff_perc <- function(x) {
  perc <- 0
  for (col in seq(2, 101, 2)) {
    if (x[col] == "Clostridioides difficile") {
      perc <- as.numeric(x[col + 1])
      break
    }
  }
  perc
}

cdiff <- apply(kraken, 1, cdiff_perc)

bac616k_embedding <- fread("data/616k_v4.embedding.txt.bz2") %>% as_tibble()
bac616k_names <- fread("data/616k_mandrake.names.txt.gz", sep='\t', header = FALSE) %>% as_tibble()
bac616k_species <- fread("data/616k_species.txt.gz", header = FALSE) %>% as_tibble()
bac616k_embedding$diversity <- simpsons_diversity[match(bac616k_names$V1, kraken$sample_id)]
bac616k_embedding$cdiff <- cdiff[match(bac616k_names$V1, kraken$sample_id)]
bac616k_embedding$species <- bac616k_species$V2[match(bac616k_names$V1, bac616k_species$V1)]

set.seed(1234)
p1 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=diversity)) + 
  geom_point(alpha=0.5, size=0.1) +
  theme_bw(base_size = 20) +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  labs(colour = "Simpson's diversity") +
  scale_colour_viridis() +
  xlim(-10, 5) +
  ylim(-5, 10) +
  xlab('SCE DIM 1') + ylab('SCE DIM 2')

p1_legend <- get_legend(p1 + 
                          theme(plot.background = element_blank(),
                                legend.background = element_blank(),
                                legend.title = element_text(size = 10),
                                legend.text = element_text(size = 10),
                                legend.position = "right"))

p2 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=diversity)) + 
  geom_point(alpha=0.5, size=0.1) +
  gghighlight(species == "Clostridioides difficile") +
  theme_bw(base_size = 20) +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  labs(colour = "Simpson's diversity") +
  scale_colour_viridis() +
  xlim(-10, 5) +
  ylim(-5, 10) +
  xlab('SCE DIM 1') + ylab('SCE DIM 2')

p3 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=diversity)) + 
  geom_point(alpha=0.5, size=0.1) +
  gghighlight(species != "Clostridioides difficile") +
  theme_bw(base_size = 20) +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  labs(colour = "Simpson's diversity") +
  scale_colour_viridis() +
  xlim(-10, 5) +
  ylim(-5, 10) +
  xlab('') + ylab('')

p4 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=cdiff)) + 
  geom_point(alpha=0.05, size=0.1) +
  theme_bw(base_size = 20) +
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  labs(colour = "C. difficile %") +
  scale_colour_viridis(direction = -1) +
  xlim(-10, 5) +
  ylim(-5, 10) +
  xlab('') + ylab('')

p4_legend <- get_legend(p4 + 
                          theme(plot.background = element_blank(),
                                legend.background = element_blank(),
                                legend.title = element_text(size = 10),
                                legend.text = element_text(size = 10),
                                legend.position = "right"))

row1 <- plot_grid(p1, p2, p3, p1_legend, rel_widths = c(1, 1, 1, .3), nrow = 1, ncol = 4, labels = c("A", "B", "C"), label_size = 12)
row2 <- plot_grid(NULL, p4, p4_legend, NULL, rel_widths = c(1, 1, .3, 1), nrow = 1, ncol = 4, labels = c("", "D"))
plot_grid(row1, row2, ncol = 1)

ggsave('./figures/bac616k_cdiff.png', width = 18, height = 10)

#
# Salmonella plot
#

salmonella_hiercc <- read.table("data/hiercc_filtered.txt", header = F, sep = "\t", col.names = c("name", "hcc900", "hcc2850"))

bac616k_embedding$hcc900 <- salmonella_hiercc$hcc900[match(bac616k_names$V1, salmonella_hiercc$name)]
bac616k_embedding$hcc2850 <- salmonella_hiercc$hcc2850[match(bac616k_names$V1, salmonella_hiercc$name)]

tb <- table(bac616k_embedding$hcc900)
bac616k_embedding$major_hcc900 <- as.character(bac616k_embedding$hcc900)
bac616k_embedding$major_hcc900[bac616k_embedding$major_hcc900 %in% names(tb)[tb<1000]] <- 'other'

salmonella_embedding <- fread("data/hiercc_mandrake.embedding.txt.bz2") %>% as_tibble()
salmonella_names <- fread("data/hiercc_mandrake.names.txt.gz", sep='\t', header = FALSE) %>% as_tibble()
salmonella_embedding$hcc900 <- salmonella_hiercc$hcc900[match(salmonella_names$V1, salmonella_hiercc$name)]
salmonella_embedding$major_hcc900 <- as.character(salmonella_embedding$hcc900)
salmonella_embedding$major_hcc900[salmonella_embedding$major_hcc900 %in% names(tb)[tb<1000]] <- 'other'
set.seed(1234)
l <- length(unique(bac616k_embedding$hcc900))
all_cols = rainbow(l, s=.6, v=.9)[sample(1:l,l)]
l <- length(unique(bac616k_embedding$major_hcc900))
major_cols = rainbow(l, s=.6, v=.9)[sample(1:l,l)]

set.seed(1234)
sp1 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=factor(hcc2850))) + 
  geom_point(alpha=0.5, size=0.1) +
  gghighlight(!is.na(hcc2850)) +
  scale_colour_manual(values=sample(colours, length(unique(bac616k_embedding$hcc2850)), replace = TRUE)) +
  theme_bw(base_size = 20) +
  labs(colour = "HC2850") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(plot.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none') +
  xlab('SCE DIM 1') + ylab('SCE DIM 2')

sp1_legend <- get_legend(sp1 + 
                           theme(plot.background = element_blank(),
                                 legend.background = element_blank(),
                                 legend.title = element_text(size = 10),
                                 legend.text = element_text(size = 10),
                                 legend.position = "right"))

sp2 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=factor(hcc900))) + 
  geom_point(alpha=0.5, size=0.1) +
  gghighlight(!is.na(hcc900)) +
  scale_colour_manual(values=cols) +
  theme_bw(base_size = 20) +
  theme(plot.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none') +
  xlab('SCE DIM 1') + ylab('SCE DIM 2')

set.seed(1234)
sp3 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=factor(major_hcc900))) + 
  geom_point(alpha=0.5, size=0.1) +
  gghighlight(!is.na(major_hcc900) & major_hcc900 != 'other') +
  scale_colour_manual(values=sample(colours, length(unique(bac616k_embedding$major_hcc900)), replace = TRUE)) +
  theme_bw(base_size = 20) +
  labs(colour = "HC900") +
  theme(plot.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none') +
  xlab('SCE DIM 1') + ylab('SCE DIM 2')

set.seed(1234)
sp4 <- ggplot(salmonella_embedding, aes(x=V1, y=V2, colour=factor(major_hcc900))) + 
  geom_point(alpha=0.5, size=0.3) +
  gghighlight(!is.na(major_hcc900) & major_hcc900 != 'other') +
  theme_bw(base_size = 20) +
  scale_colour_manual(values=sample(colours, length(unique(salmonella_embedding$major_hcc900)), replace = TRUE)) +
  labs(colour = "HC900") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(plot.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none') +
  xlab('SCE DIM 1') + ylab('SCE DIM 2')

sp4_legend <- get_legend(sp4 + 
                          theme(plot.background = element_blank(),
                                legend.background = element_blank(),
                                legend.title = element_text(size = 10),
                                legend.text = element_text(size = 10),
                                legend.position = "right"))

sp_row1 <- plot_grid(sp1, sp2, sp1_legend, rel_widths = c(1, 1, .3), nrow = 1, ncol = 3, labels = c("A", "B"), label_size = 12)
sp_row2 <- plot_grid(sp3, sp4, sp4_legend, rel_widths = c(1, 1, .3), nrow = 1, ncol = 3, labels = c("C", "D"), label_size = 12)
plot_grid(sp_row1, sp_row2, ncol = 1)

ggsave('./figures/bac616k_salmonella.png', width = 18, height = 10)

#
# Listeria plot
#

listeria_pp <- read.table("data/listeria_pp_assign.txt", header = F, sep = ",", col.names = c("name", "VLKC"))
bac616k_embedding$listeria_pp <- listeria_pp$VLKC[match(bac616k_names$V1, listeria_pp$name)]

# strain 5 is misassigned due to contamination
bac616k_embedding$listeria_pp[bac616k_embedding$listeria_pp == 5] <- NA

ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=factor(listeria_pp))) + 
  geom_point(alpha=0.3, size=2) +
  gghighlight(!is.na(listeria_pp),
              use_direct_label = FALSE,
              unhighlighted_params = list(size=0.2)) +
  theme_bw(base_size = 20) +
  scale_colour_brewer(type="qual") +
  theme_bw(base_size = 20) +
  theme(plot.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #legend.background = element_blank(),
        legend.position = 'right') +
  guides(colour = guide_legend(override.aes = list(size=7, alpha=1))) +
  labs(colour = "VLKC")

ggsave('./figures/bac616k_lm.png', width = 14, height = 10)

#
# Mtb plot
#

mtb_lin <- read.table("data/mtb_lineages.csv", header = F, sep = ",", col.names = c("name", "mtb_lineage"))
bac616k_embedding$mtb_lin <- mtb_lin$mtb_lineage[match(bac616k_names$V1, mtb_lin$name)]

tb_p1 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=mtb_lin)) + 
  geom_point(alpha=0.5, size=0.1) +
  gghighlight(!is.na(mtb_lin), use_direct_label = FALSE) +
  theme_bw(base_size = 20) +
  scale_colour_brewer(type="qual") +
  theme_bw(base_size = 20) +
  theme(plot.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.position = 'right') +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  labs(colour = "Lineage")

tb_p2 <- ggplot(bac616k_embedding, aes(x=V1, y=V2, colour=mtb_lin)) + 
  geom_point(alpha=1, size=1) +
  gghighlight(!is.na(mtb_lin), use_direct_label = FALSE) +
  theme_bw(base_size = 20) +
  scale_colour_brewer(type="qual") +
  theme_bw(base_size = 20) +
  theme(plot.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.position = 'right') +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  labs(colour = "Lineage") +
  xlim(20, 35) + ylim(-15, 2)

tb_p2_legend <- get_legend(tb_p2 + 
                           theme(plot.background = element_blank(),
                                legend.background = element_blank(),
                                legend.title = element_text(size = 14),
                                legend.text = element_text(size = 14),
                                legend.position = "right"))

plot_grid(tb_p1 + theme(legend.position = 'none'),
          tb_p2 + theme(legend.position = 'none'),
          tb_p2_legend,
          rel_widths = c(1, 1, .3), nrow = 1, ncol = 3, labels = c("A", "B"), label_size = 12)

ggsave('./figures/bac616k_tb.png', width = 18, height = 7)


