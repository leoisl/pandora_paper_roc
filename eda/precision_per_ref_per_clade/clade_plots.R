library("ggplot2")
library("gridExtra")
library("grid")
library("tidyverse")

###############################################################################################
# args
args <- commandArgs(trailingOnly=TRUE)
csv_file <- args[1]
plot_title <- args[2]
output_file <- args[3]

ref_colouring <- c(
  # A
  "H6 (CP010170.1)" = "#e735f32b",
  "H7 (CP010171.1)" = "#e735f32b",
  "S1 (CP010226.1)" = "#e735f32b",
  "W3110 (NC_007779.1)" = "#e735f32b",
  "VR50 (NZ_CP011134.1)" = "#e735f32b",
  "Y5 (NZ_CP013483.1)" = "#e735f32b",
  "MRSN346595 (NZ_CP018109.1)" = "#e735f32b",
  "EcRV308Chr (NZ_LM995446.1)" = "#e735f32b",

  # B1
  "S21 (CP010230.1)" = "#002fff2b",
  "09-00049 (NZ_CP015228.1)" = "#002fff2b",

  # B2
  "S88 (NC_011742.1)" = "#d4f7ffff",
  "NGF1 (NZ_CP016007.1)" = "#d4f7ffff",
  "LF82 (NC_011993.1)" = "#d4f7ffff",
  "JJ1886 (NC_022648.1)" = "#d4f7ffff",
  "NCTC_13441 (NZ_LT632320.1)" = "#d4f7ffff",
  "EC958 (NZ_HG941718.1)" = "#d4f7ffff",

  # D
  "C1 (CP010116.1)" = "#0080002c",
  "C4 (CP010121.1)" = "#0080002c",
  "MRSN346647 (CP018206.1)" = "#0080002c",
  "UMN026 (CU928163.2)" = "#0080002c",

  # F
  "ECONIH1 (NZ_CP009859.1)" = "#ff00002b",
  "ST648 (NZ_CP008697.1)" = "#ff00002b",
  "SMS-3-5 (NC_010498.1)" = "#ff00002b",
  "CE10 (NC_017646.1)" = "#ff00002b"
)
ref_ordering <- names(ref_colouring)

sample_colouring <- c(
                # A
                "Escherichia_coli_MSB1_9D" = "#FF00FF",
                "H131800734" = "#FF00FF",
                "Escherichia_coli_MINF_8D" = "#FF00FF",
                "Escherichia_coli_MSB1_1A" = "#FF00FF",
                "Escherichia_coli_MSB1_8G" = "#FF00FF",
  
                # B1
                "Escherichia_coli_MSB1_4E" = "#00ceffff",
                "063_STEC" = "#00ceffff",
                "CFT073" = "#00ceffff",
                "Escherichia_coli_MSB2_1A" = "#00ceffff",
                "Escherichia_coli_MINF_1D" = "#00ceffff", 
               
                # D
                "Escherichia_coli_MSB1_7C" = "#008000ff",
                "Escherichia_coli_MINF_7C" = "#008000ff",
                "Escherichia_coli_MSB1_6C" = "#008000ff",
                "ST38" = "#008000ff",
                "Escherichia_coli_MSB1_7A" = "#008000ff",
                "Escherichia_coli_MSB1_3B" = "#008000ff",
                
                # F
                "Escherichia_coli_MSB1_4I" = "#FF0000",
                "Escherichia_coli_MINF_1A" = "#FF0000",
                "Escherichia_coli_MINF_9A" = "#FF0000",
                "Escherichia_coli_MSB1_8B" = "#FF0000"
)
sample_ordering <- names(sample_colouring)
###############################################################################################


###############################################################################################
# input
precision_table <- read.csv(csv_file, header=TRUE)

# bugfix for medaka sample Escherichia_coli_MSB1_7A, reference NZ_CP011134.1
if (plot_title == "precision_medaka_pandora") {
	precision_table[nrow(precision_table) + 1,] = list("Medaka", "Escherichia_coli_MSB1_7A", 0.0, "VR50 (NZ_CP011134.1)", "Medaka / VR50 (NZ_CP011134.1)")
}

precision_table_for_pandora_no_denovo <- precision_table %>% filter(tool == "Pandora illumina no denovo")
if(nrow(precision_table_for_pandora_no_denovo) == 0){
  precision_table_for_pandora_no_denovo <- precision_table %>% filter(tool == "Pandora nanopore no denovo")
}

precision_table_for_pandora_with_denovo <- precision_table %>% filter(tool == "Pandora illumina with denovo")
if(nrow(precision_table_for_pandora_with_denovo) == 0){
  precision_table_for_pandora_with_denovo <- precision_table %>% filter(tool == "Pandora nanopore with denovo")
}
###############################################################################################

###############################################################################################
# processing
png(output_file, width = 2500, height = 2000)
index <- 1
plots <- list()

min_ylim <- 0.90
if (plot_title == "precision_medaka_pandora") {
	min_ylim <- 0.50
}

for (current_ref in ref_ordering) {
  if (!identical(current_ref, "PRG")) {
    precision_table_for_ref <- precision_table %>% filter(ref == current_ref)
    precision_table_for_ref$sample <- factor(precision_table_for_ref$sample, levels=sample_ordering)
    ref_color <- ref_colouring[[current_ref]]

    plot <- ggplot(data = precision_table_for_ref, aes(x=sample, y=precision, fill=sample, colour=tool)) +
        scale_fill_manual(values = sample_colouring) +
        geom_bar(stat = "identity") +
        coord_cartesian(ylim=c(min_ylim, 1)) +
        geom_line(data=precision_table_for_pandora_no_denovo, aes(x=sample, y=precision, group=1), colour="black", size = 2) +
        geom_line(data=precision_table_for_pandora_with_denovo, aes(x=sample, y=precision, group=1), colour = "orange", size = 2) +
        ggtitle(current_ref) +
        theme(
          panel.background = element_rect(fill = ref_color,
                                          colour = ref_color),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size=30),
          axis.text.x=element_blank(),
          legend.position = "none"
        )
    plots[[index]] <- plot
    index <- index + 1
  }
}

print(length(plots))
grid.arrange(grobs=plots, ncol = 5, nrow = 5, top = textGrob(plot_title,gp=gpar(fontsize=60,font=3)))

dev.off()
