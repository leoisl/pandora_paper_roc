library("ggplot2")
library("gridExtra")
library("tidyverse")

###############################################################################################
# configs
tool = "samtools_pandora"

plot_name = paste("precision_per_ref_per_clade_histogram", tool, "png", sep=".")


ref_colouring <- c(
  # A
  "CP010170.1" = "#e735f32b",
  "CP010171.1" = "#e735f32b",
  "CP010226.1" = "#e735f32b",
  "NC_007779.1" = "#e735f32b",
  "NZ_CP011134.1" = "#e735f32b",
  "NZ_CP013483.1" = "#e735f32b",
  "NZ_CP018109.1" = "#e735f32b",
  "NZ_LM995446.1" = "#e735f32b",
  
  # B1
  "CP010230.1" = "#002fff2b",
  "NZ_CP015228.1" = "#002fff2b",
  
  # B2
  "NC_011742.1" = "#d4f7ffff",
  "NZ_CP016007.1" = "#d4f7ffff",
  "NC_011993.1" = "#d4f7ffff",
  "NC_004431.1" = "#d4f7ffff",
  "NC_022648.1" = "#d4f7ffff",
  "NZ_LT632320.1" = "#d4f7ffff",
  "NZ_HG941718.1" = "#d4f7ffff",
  
  # D
  "CP010116.1" = "#0080002c",
  "CP010121.1" = "#0080002c",
  "CP018206.1" = "#0080002c",
  "CU928163.2" = "#0080002c",
  
  # F
  "NZ_CP009859.1" = "#ff00002b",
  "NZ_CP008697.1" = "#ff00002b",
  "NC_010498.1" = "#ff00002b",
  "NC_017646.1" = "#ff00002b"
)
ref_ordering = names(ref_colouring)

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
sample_ordering = names(sample_colouring)
###############################################################################################


###############################################################################################
# input
precision_table = read.csv(paste("precision_per_ref_per_clade", tool, "csv", sep="."), header=TRUE)
precision_table_for_pandora_no_denovo = precision_table %>% filter(tool == "Pandora illumina no denovo")
precision_table_for_pandora_with_denovo = precision_table %>% filter(tool == "Pandora illumina with denovo")
###############################################################################################

###############################################################################################
# processing
png(plot_name, width = 2000, height = 2000)
index = 1
plots <- list()
for (current_ref in ref_ordering) {
  if (!identical(current_ref, "PRG")) {
    precision_table_for_ref = precision_table %>% filter(ref == current_ref)
    precision_table_for_ref$sample = factor(precision_table_for_ref$sample, levels=sample_ordering)
    ref_color = ref_colouring[[current_ref]]
    
    plot <- ggplot(data = precision_table_for_ref, aes(x=sample, y=precision, fill=sample, colour=tool)) +
        scale_fill_manual(values = sample_colouring) +
        geom_bar(stat = "identity") +
        coord_cartesian(ylim=c(0.90, 1)) + 
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
    index = index + 1
  }
}


grid.arrange(grobs=plots, ncol = 5, nrow = 5)

dev.off()