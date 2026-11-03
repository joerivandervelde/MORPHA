library(R.utils)
library(data.table)
library(this.path)
library(ggplot2)
library(ggtext)
library(ggbeeswarm)
library(scales)

gene_name = "MEFV"

script_dir <- this.dir()
parent_dir <- dirname(script_dir)
setwd(dirname(script_dir))

#######################################################
# Read in the protein haplotypes for 1kGP individuals #
#######################################################
haplo_path <- "results/MEFV/MEFV_protein_haplotypes.fa"
fasta_lines <- readLines(haplo_path, warn = FALSE)
header_idx <- grep("^>", fasta_lines)

n_records <- length(header_idx)
gene_vec      <- character(n_records)
variant_vec   <- character(n_records)
n_vec         <- rep(NA_integer_, n_records)
samples_vec   <- rep(NA_character_, n_records)
sequence_vec  <- character(n_records)
variantID_vec <- rep(NA_character_, n_records)

for (i in seq_along(header_idx)) {
  h_line <- sub("^>", "", fasta_lines[header_idx[i]])
  parts <- strsplit(h_line, "\\|")[[1]]
  
  variant_id <- parts[1]
  gene_vec[i] <- gene_name
  variant_vec[i] <- variant_id
  
  n_str <- parts[grep("^n=", parts)]
  if (length(n_str) > 0) { n_vec[i] <- as.integer(sub("^n=", "", n_str[1])) }
  
  samples_str <- parts[grep("^samples=", parts)]
  if (length(samples_str) > 0) {
    samples_val <- sub("^samples=", "", samples_str[1])
    samples_vec[i] <- samples_val
  }
  
  start_seq <- header_idx[i] + 1
  end_seq <- if (i < length(header_idx)) {
    header_idx[i + 1] - 1
  } else {
    length(fasta_lines)
  }
  
  if (start_seq <= end_seq) {
    sequence_vec[i] <- paste(fasta_lines[start_seq:end_seq], collapse = "")
  } else {
    sequence_vec[i] <- ""
  }
}

all_haplotypes <- data.frame(
  gene = gene_vec, variant = variant_vec, n = n_vec, 
  variantID = as.character(variantID_vec), samples = samples_vec, 
  sequence = sequence_vec, stringsAsFactors = FALSE
)
cat("Total nr of haplotypes:", dim(all_haplotypes)[1], "\n")

############################################
# Map samples to Superpopulations          #
############################################
dt <- copy(as.data.table(all_haplotypes))
dt[, row_id := .I]
sample_long <- dt[, .(sample = trimws(unlist(strsplit(samples, ",", fixed = TRUE)))), by = row_id]
sample_long <- merge(sample_long, dt[, .(row_id, gene, variant)], by = "row_id", all.x = TRUE)

samples <- read.delim("viz/igsr_samples.tsv", sep = "\t", check.names = FALSE)
setnames(samples, names(samples), gsub(" ", "_", names(samples)))

sample_long[, sample_base := sub("\\.h[12]$", "", sample)]
sample_long[samples, on = .(sample_base == Sample_name), superpop := i.Superpopulation_code]
sample_long[samples, on = .(sample_base == Sample_name), pop := i.Population_code]

# Special handling for admixed individual
sample_long[sample == "HG01783.h1", superpop := "EUR"]
sample_long[sample == "HG01783.h2", superpop := "AFR"]
sample_long[sample == "HG01783.h1", pop := "IBS"]
sample_long[sample == "HG01783.h2", pop := "MSL"]
sample_long[, sample_base := NULL]

#######################################################
# Bridge: Aggregate sample_long into Superpop Counts  #
#######################################################
superpop_counts <- sample_long[!is.na(superpop), .N, by = .(gene, variant, superpop)]
dt_all <- dcast(superpop_counts, gene + variant ~ paste0("n_", superpop), value.var = "N", fill = 0)

required_superpops <- c("n_EUR", "n_AFR", "n_EAS", "n_SAS", "n_AMR")
for (col in required_superpops) {
  if (!col %in% names(dt_all)) dt_all[, (col) := 0]
}
setnames(dt_all, "variant", "natural_haplo")

####################
# Plot Generation  #
####################
natuhap_pal <- c("white", "#505050")
superpop_pal <- c(
  "EUR" = "#0072b2",
  "AFR" = "#d55e00",
  "EAS" = "#009e73",
  "SAS" = "#cc79a7",
  "AMR" = "#f0e442"
)

dt_all[, n_samples := n_EUR + n_AFR + n_EAS + n_SAS + n_AMR]
setkey(dt_all, NULL)
dt_all[, pair_id := .I]

superpop_cols <- c("n_EUR", "n_AFR", "n_EAS", "n_SAS", "n_AMR")
long_all <- melt(
  dt_all,
  id.vars = c("gene", "natural_haplo", "pair_id", "n_samples"),
  measure.vars = superpop_cols,
  variable.name = "superpop_raw",
  value.name = "count"
)
long_all[, superpop := sub("^n_", "", superpop_raw)]
long_all[, superpop_raw := NULL]
long_all <- long_all[count > 0]

# Build simplified panel labels
gene_labels <- long_all[, {
  base_label <- paste0("<b>", unique(gene), "</b>")
  curr_gene <- unique(gene)
  if (curr_gene == "LOC128092249") {
    panel_label <- paste0('<span>', base_label, '</span>')
  } else {
    panel_label <- base_label
  }
  .(panel_label = panel_label)
}, by = gene]

long_all <- merge(long_all, gene_labels, by = "gene", all.x = TRUE)

# Create pair_summary and sort by the true numerical value of the variant suffix
pair_summary <- long_all[, .(
  n_samples = max(n_samples)
), by = .(panel_label, natural_haplo, pair_id)]

# Extract numerical suffix for correct sorting (e.g., variant_2 comes before variant_10)
pair_summary[, sort_num := as.integer(sub(".*_([0-9]+)$", "\\1", natural_haplo))]
setorder(pair_summary, panel_label, sort_num)

# Assign clean sequential numbering 1, 2, 3, 4... per panel
pair_summary[, nathapnr := 1:.N, by = panel_label]

# Merge nathapnr back into long_all
long_all[pair_summary, on = .(pair_id, natural_haplo, panel_label), nathapnr := i.nathapnr]

setorder(pair_summary, panel_label, nathapnr)
long_all[, pair_id := factor(pair_id, levels = pair_summary$pair_id)]

axis_labels_map <- setNames(as.character(pair_summary$nathapnr), pair_summary$pair_id)
label_dt <- unique(long_all[, .(pair_id, panel_label, n_samples)])

p <- ggplot(long_all, aes(x = pair_id, y = count, fill = superpop)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8, linewidth = 0.5) +
  geom_text(
    data = label_dt, 
    aes(x = pair_id, y = 0, label = paste0("", n_samples)), 
    inherit.aes = FALSE, 
    angle = 45,           # Rotates text 90 degrees
    hjust = 1.1,          # Pushes the text down below the axis point
    vjust = 0.5,          # Centers the text vertically on the tick
    size = 12 / .pt, 
    fontface = "italic"
  ) +
  scale_fill_manual(values = superpop_pal, labels = c("EUR" = "European","AFR" = "African","EAS" = "East Asian","SAS" = "South Asian","AMR" = "Admixed American")) +
  scale_color_identity() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(labels = axis_labels_map, guide = guide_axis(check.overlap = TRUE)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Sample and superpopulation distribution of MEFV haplotypes in 1kGP",
    x = "Natural protein haplotypes from 1kGP with total sample counts in italics",
    y = "Relative 1kGP superpopulation breakdown per natural haplotype",
    fill = "Superpopulation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 21),
    plot.subtitle = element_text(size = 14),
    strip.text = element_blank(),       
    strip.background = element_blank(), 
    axis.text = element_text(size = 12, margin = margin(t = 0), color="black"),
    axis.ticks.length = unit(0, "pt"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid = element_blank(),
  )
ggsave(filename = "superpopulation_breakdown.png", plot = p, width = 500, height = 200, units = "mm", dpi = 300, limitsize = FALSE)


#########
# By population
###########
# Changed from superpop to pop
pop_counts <- sample_long[!is.na(pop), .N, by = .(gene, variant, pop)]
dt_all <- dcast(pop_counts, gene + variant ~ paste0("n_", pop), value.var = "N", fill = 0)

pop_cols <- grep("^n_", names(dt_all), value = TRUE)
setnames(dt_all, "variant", "natural_haplo")

dt_all[, n_samples := rowSums(.SD), .SDcols = pop_cols]
setkey(dt_all, NULL)
dt_all[, pair_id := .I]

long_all <- melt(
  dt_all,
  id.vars = c("gene", "natural_haplo", "pair_id", "n_samples"),
  measure.vars = pop_cols,
  variable.name = "pop_raw",
  value.name = "count"
)
long_all[, pop := sub("^n_", "", pop_raw)]
long_all[, pop_raw := NULL]
long_all <- long_all[count > 0]

# Attach superpop mapping for smart color grouping
pop_to_superpop <- unique(samples[, c("Population_code", "Superpopulation_code")])
setnames(pop_to_superpop, c("Population_code", "Superpopulation_code"), c("pop", "superpop"))
long_all[pop_to_superpop, on = .(pop), superpop := i.superpop]
long_all[pop %in% c("IBS", "MSL") & is.na(superpop), superpop := fifelse(pop == "IBS", "EUR", "AFR")]

# Build simplified panel labels
gene_labels <- long_all[, {
  base_label <- paste0("<b>", unique(gene), "</b>")
  curr_gene <- unique(gene)
  if (curr_gene == "LOC128092249") {
    panel_label <- paste0('<span>', base_label, '</span>')
  } else {
    panel_label <- base_label
  }
  .(panel_label = panel_label)
}, by = gene]

long_all <- merge(long_all, gene_labels, by = "gene", all.x = TRUE)

pair_summary <- long_all[, .(
  n_samples = max(n_samples)
), by = .(panel_label, natural_haplo, pair_id)]

pair_summary[, sort_num := as.integer(sub(".*_([0-9]+)$", "\\1", natural_haplo))]
setorder(pair_summary, panel_label, sort_num)

pair_summary[, nathapnr := 1:.N, by = panel_label]
long_all[pair_summary, on = .(pair_id, natural_haplo, panel_label), nathapnr := i.nathapnr]

setorder(pair_summary, panel_label, nathapnr)
long_all[, pair_id := factor(pair_id, levels = pair_summary$pair_id)]

axis_labels_map <- setNames(as.character(pair_summary$nathapnr), pair_summary$pair_id)
label_dt <- unique(long_all[, .(pair_id, panel_label, n_samples)])

# Generate coordinated color gradients per superpopulation
base_superpop_colors <- c("EUR" = "#0072b2", "AFR" = "#d55e00", "EAS" = "#009e73", "SAS" = "#cc79a7", "AMR" = "#f0e442")
set.seed(42)
pop_pal <- unlist(lapply(unique(long_all$superpop), function(sp) {
  pops_in_sp <- sort(unique(long_all[superpop == sp, pop]))
  n <- length(pops_in_sp)
  if (n == 1) return(setNames(base_superpop_colors[sp], pops_in_sp))
  color_ramp <- colorRampPalette(c(base_superpop_colors[sp], "white"))(n + 2)[2:(n + 1)]
  setNames(color_ramp, pops_in_sp)
}))

p <- ggplot(long_all, aes(x = pair_id, y = count, fill = pop)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8, linewidth = 0.5) +
  geom_text(
    data = label_dt, 
    aes(x = pair_id, y = 0, label = paste0("", n_samples)), 
    inherit.aes = FALSE, 
    angle = 45,           # Rotates text 90 degrees
    hjust = 1.1,          # Pushes the text down below the axis point
    vjust = 0.5,          # Centers the text vertically on the tick
    size = 12 / .pt, 
    fontface = "italic"
  ) +
  scale_fill_manual(values = pop_pal) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(labels = axis_labels_map, guide = guide_axis(check.overlap = TRUE)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Sample and population distribution of MEFV haplotypes in 1kGP",
    x = "Natural protein haplotypes from 1kGP with total sample counts in italics",
    y = "Relative 1kGP population breakdown per natural haplotype",
    fill = "Population"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 21),
    plot.subtitle = element_text(size = 14),
    strip.text = element_blank(),       
    strip.background = element_blank(), 
    axis.text = element_text(size = 12, margin = margin(t = 0), color="black"),
    axis.ticks.length = unit(0, "pt"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    panel.grid = element_blank(),
  )

ggsave(filename = "population_breakdown.png", plot = p, width = 500, height = 200, units = "mm", dpi = 300)
