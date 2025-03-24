# Load packages
library(edgeR)
library(Rsubread)
library(GenomicFeatures)
library(limma)

# Initialize file paths
programs_path <- "/Users/teaguemcc/Programs"
git_path <- file.path(programs_path, "git_clones/ECE759_project")
gtf_path <- file.path(programs_path, "data/genomes/Arabidopsis_thaliana.TAIR10.54.gtf") # only stored on cloudy server # nolint (this removed VSCode error detection)

# Load the data from csv file
data <- read.csv(file.path(git_path, "cellular_clarity/read_count_matrix.csv"), row.names = 1) # nolint
data_matrix <- as.matrix(t(data)) # Transpose where columns are samples and rows are genes for edgeR # nolint
sample_metadata <- read.csv(file.path(git_path, "cellular_clarity/sample_names.csv"), row.names = 1) # Used to group data for DGE selection # nolint
colnames(data_matrix) <- sample_metadata[colnames(data_matrix), "sample_name"] # matches sample ids and replaces them with sample names #nolint
sample_metadata <- sample_metadata[match(colnames(data_matrix), sample_metadata$sample_name), , drop = FALSE] # Reorder sample_names so rows match the column order of data_matrix  # nolint

# Define group up front
group <- factor(sample_metadata$experimental_condition)
levels(group) <- gsub(" -", "_minus", levels(group))
levels(group) <- gsub(" \\+", "_plus", levels(group))

# Check if row names in sample_names match column names in data_matrix
if (all(sample_metadata$sample_name %in% colnames(data_matrix))) {
  # Assign to the DGEList object
  dge <- DGEList(counts = data_matrix, group = group)
} else {
    stop("Mismatch: Some row names in sample_names do not match columns in data_matrix") # nolint
}
# Check that the sample_ids match in metadata and the data matrix 
# print(identical(sample_metadata$sample_name, colnames(dge$counts)))  # Check if sample order matches: Should return TRUE, did return true when tested # nolint

# Creating gene_lenths for RPKM calculation - only needed to run once, see the read.csv command in next section for results # nolint
#txdb <- makeTxDbFromGFF(gtf_path, format = "gtf") # Creates a transcript database from the gtf file # nolint
#exons_by_gene <- exonsBy(txdb, by = "gene") # nolint
#gene_lengths <- sapply(exons_by_gene, function(exons) sum(width(reduce(exons)))) # nolint
#gene_id_list <- rownames(data_matrix) # nolint
#gene_lengths_of_interest <- gene_lengths[intersect(names(gene_lengths), gene_id_list)] # nolint
#write.csv(gene_lengths_of_interest, file.path(git_path, "cellular_clarity/gene_lengths.csv")) # nolint

# Normalize Counts - use RPKM normalized data (Selene's disertation section 4.3.2) # nolint
gene_lengths <- read.csv(file.path(git_path, "cellular_clarity/gene_lengths.csv"), row.names = 1)  # nolint
gene_lengths_vector <- gene_lengths[[1]]  # extract first (or appropriate) column # nolint
names(gene_lengths_vector) <- rownames(gene_lengths) # Convert to named numeric vector # nolint
dge$genes$Length <- gene_lengths_vector[rownames(dge$counts)] # Provide gene lengths to DGE object # nolint
#rpkm_values <- rpkm(dge) # Only needed to do once, see read.csv command below # nolint
#write.csv(rpkm_values, file.path(git_path, "cellular_clarity/rpkm_values.csv")) # nolint
rpkm_values <- read.csv(file.path(git_path, "cellular_clarity/rpkm_values.csv"), row.names = 1, check.names = FALSE) # Saved rpkm_values to be in genes x samples  # nolint
# Rename columns of rpkm_values using sample_metadata$sample_name, converting to "A + 1" # nolint
#colnames(rpkm_values) <- sample_metadata[colnames(rpkm_values), "sample_name"]
# Rename dge sample columns to match rpkm_values, temporary
# colnames(dge$counts) <- sample_metadata[colnames(dge$counts), "sample_name"]
# print(sum(is.na(dge$genes$Length)))  # Should be 0, tested and this was true # nolint
# print(identical(colnames(dge$counts), colnames(rpkm_values))) # Should be TRUE, checked and TRUE # nolint
# print(identical(rownames(dge$samples), colnames(dge$counts)))   # Should be TRUE, checked and TRUE # nolint
# print(identical(sample_metadata$sample_name, colnames(data_matrix)))  # Should be TRUE, checked and TRUE # nolint

# Filter lowly expressed genes - need to use the logic presented in paper methods (Selene's disertation section 4.3.2) # nolint
# Get all sample names for timepoint A
samples_A <- grep("^A", colnames(rpkm_values), value = TRUE) # Get all sample names for timepoint A # nolint
samples_B <- grep("^B", colnames(rpkm_values), value = TRUE) # Get all sample names for timepoint B # nolint
samples_C <- grep("^C", colnames(rpkm_values), value = TRUE) # Get all sample names for timepoint C # nolint
samples_D <- grep("^D", colnames(rpkm_values), value = TRUE) # Get all sample names for timepoint D # nolint
samples_E <- grep("^E", colnames(rpkm_values), value = TRUE) # Get all sample names for timepoint E # nolint
samples_F <- grep("^F", colnames(rpkm_values), value = TRUE) # Get all sample names for timepoint F # nolint
samples_G <- grep("^G", colnames(rpkm_values), value = TRUE) # Get all sample names for timepoint G # nolint
# Make DGE objects for each timepoint, this is becuase filtering is timepoint specific # nolint
# Keep genes wtih RPKM >= 20 for at least 5 out of 6 samples
subset_dge <- function(dge, genes_keep, samples_keep) {
  new_dge <- list()
  new_dge$counts <- dge$counts[genes_keep, samples_keep, drop = FALSE]
  new_dge$samples <- dge$samples[samples_keep, , drop = FALSE]
  new_dge$samples$group <- dge$samples$group[samples_keep] # correctly define sample groups # nolint

  if (!is.null(dge$genes) && is.data.frame(dge$genes)) {
    new_dge$genes <- dge$genes[genes_keep, , drop = FALSE]
  }
  class(new_dge) <- "DGEList"
  calcNormFactors(new_dge)
}

rpkm_A <- rpkm_values[, samples_A] # Subset RPKM matrix to timepoint A # nolint
keep_A <- rowSums(rpkm_A >= 20) >= 2 # Keep genes with RPKM >= 20 in at least 5 of 6 samples (2 of 3 for A only) # nolint
dge_A <- subset_dge(dge, keep_A, samples_A) # Note that these are normalized read counts (not RPKM but normalized in the subseting function for lib size) # nolint

rpkm_B <- rpkm_values[, samples_B] # nolint
keep_B <- rowSums(rpkm_B >= 20) >= 5 # nolint
dge_B <- subset_dge(dge, keep_B, samples_B)

rpkm_C <- rpkm_values[, samples_C] # nolint
keep_C <- rowSums(rpkm_C >= 20) >= 5 # nolint
dge_C <- subset_dge(dge, keep_C, samples_C)

rpkm_D <- rpkm_values[, samples_D] # nolint
keep_D <- rowSums(rpkm_D >= 20) >= 5 # nolint
dge_D <- subset_dge(dge, keep_D, samples_D)

rpkm_E <- rpkm_values[, samples_E] # nolint
keep_E <- rowSums(rpkm_E >= 20) >= 5 # nolint
dge_E <- subset_dge(dge, keep_E, samples_E)

rpkm_F <- rpkm_values[, samples_F] # nolint
keep_F <- rowSums(rpkm_F >= 20) >= 5 # nolint
dge_F <- subset_dge(dge, keep_F, samples_F)

rpkm_G <- rpkm_values[, samples_G] # nolint
keep_G <- rowSums(rpkm_G >= 20) >= 5 # nolint
dge_G <- subset_dge(dge, keep_G, samples_G)

# Analysis to identify differentially expressed genes: 
# Create a design matrix consistent with comparing Fe - to Fe + samples
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Estimate the dispersion before statistical testing
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Define contrasts for differential expression analysis
contrast_matrix <- makeContrasts(
  E_minus_vs_E_plus = E_minus - E_plus,
  B_minus_vs_B_plus = B_minus - B_plus,
  D_minus_vs_D_plus = D_minus - D_plus,
  G_minus_vs_G_plus = G_minus - G_plus,
  C_minus_vs_C_plus = C_minus - C_plus,
  F_minus_vs_F_plus = F_minus - F_plus,
  G_minus_vs_G_plus = G_minus - G_plus,
  levels = design
)

# Perform statistical testing for differential expression analysis
qlf_E <- glmQLFTest(fit, contrast = contrast_matrix[, "E_minus_vs_E_plus"])
qlf_B <- glmQLFTest(fit, contrast = contrast_matrix[, "B_minus_vs_B_plus"])
qlf_D <- glmQLFTest(fit, contrast = contrast_matrix[, "D_minus_vs_D_plus"])
qlf_G <- glmQLFTest(fit, contrast = contrast_matrix[, "G_minus_vs_G_plus"])
qlf_C <- glmQLFTest(fit, contrast = contrast_matrix[, "C_minus_vs_C_plus"])
qlf_F <- glmQLFTest(fit, contrast = contrast_matrix[, "F_minus_vs_F_plus"])

# Extract the significant DEGs based on the papers stated criteria 
deg_path <- file.path(git_path, "cellular_clarity/DEGs")

results_E <- topTags(qlf_E, n = Inf)$table
degs_E <- results_E[results_E$FDR < 0.05 & abs(results_E$logFC) > 0.75, ]
write.csv(degs_E, file.path(deg_path, "DEGs_E_minus_vs_E_plus.csv"))

results_B <- topTags(qlf_B, n = Inf)$table
degs_B <- results_B[results_B$FDR < 0.05 & abs(results_B$logFC) > 0.75, ]
write.csv(degs_B, file.path(deg_path, "DEGs_B_minus_vs_B_plus.csv"))

results_D <- topTags(qlf_D, n = Inf)$table
degs_D <- results_D[results_D$FDR < 0.05 & abs(results_D$logFC) > 0.75, ]
write.csv(degs_D, file.path(deg_path, "DEGs_D_minus_vs_D_plus.csv"))

results_G <- topTags(qlf_G, n = Inf)$table
degs_G <- results_G[results_G$FDR < 0.05 & abs(results_G$logFC) > 0.75, ]
write.csv(degs_G, file.path(deg_path, "DEGs_G_minus_vs_G_plus.csv"))

results_C <- topTags(qlf_C, n = Inf)$table
degs_C <- results_C[results_C$FDR < 0.05 & abs(results_C$logFC) > 0.75, ]
write.csv(degs_C, file.path(deg_path, "DEGs_C_minus_vs_C_plus.csv"))

results_F <- topTags(qlf_F, n = Inf)$table
degs_F <- results_F[results_F$FDR < 0.05 & abs(results_F$logFC) > 0.75, ]
write.csv(degs_F, file.path(deg_path, "DEGs_F_minus_vs_F_plus.csv"))
print('script is complete')