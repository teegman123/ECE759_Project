library(edgeR)

# Load the data from csv file
data <- read.csv('/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/read_count_matrix.csv', row.names = 1)
data_matrix <- as.matrix(t(data)) # Transpose where columns are samples and rows are genes for edgeR
sample_names <- read.csv('/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/sample_names.csv', row.names = 1)

# Check if row names in sample_names match column names in data_matrix
if (all(rownames(sample_names) %in% colnames(data_matrix))) {
    # Define experimental conditions
    group <- factor(sample_names$experimental_condition)
    # Assign to the DGEList object
    dge <- DGEList(counts = data_matrix, group = group)
} else {
    stop("Mismatch: Some row names in sample_names do not match columns in data_matrix")
}
# Check that group levels match expected conditions
#table(group)
#print(identical(rownames(sample_names), colnames(dge$counts)))  # Check if sample order matches: Should return TRUE

# Filter lowly expressed genes
keep <- filterByExpr(dge, group = group)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize Counts - accounts for differences in library sizes, do not do this, instead use RPKM normalized data
#dge <- normLibSizes(dge)

# Create a design matrix consistent with comparing Fe - to Fe + samples
# Convert group names into valid R variable names while keeping Fe+ and Fe- distinct
levels(group) <- gsub(" -", "_minus", levels(group))
levels(group) <- gsub(" \\+", "_plus", levels(group))  # Escape + because it's a regex special character

# Create the design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Estimate the dispersion before statistical testing
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Define contrasts for differential expression analysis
library(limma)

# Define valid contrasts
contrast_matrix <- makeContrasts(
  E_minus_vs_E_plus = E_minus - E_plus,
  B_minus_vs_B_plus = B_minus - B_plus,
  D_minus_vs_D_plus = D_minus - D_plus,
  G_minus_vs_G_plus = G_minus - G_plus,
  C_minus_vs_C_plus = C_minus - C_plus,
  F_minus_vs_F_plus = F_minus - F_plus,
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
results_E <- topTags(qlf_E, n = Inf)$table
degs_E <- results_E[results_E$FDR < 0.05 & abs(results_E$logFC) > 0.75, ]
write.csv(degs_E, "/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/DEGs/DEGs_E_minus_vs_E_plus.csv")

results_B <- topTags(qlf_B, n = Inf)$table
degs_B <- results_B[results_B$FDR < 0.05 & abs(results_B$logFC) > 0.75, ]
write.csv(degs_B, "/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/DEGs/DEGs_B_minus_vs_B_plus.csv")

results_D <- topTags(qlf_D, n = Inf)$table
degs_D <- results_D[results_D$FDR < 0.05 & abs(results_D$logFC) > 0.75, ]
write.csv(degs_D, "/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/DEGs/DEGs_D_minus_vs_D_plus.csv")

results_G <- topTags(qlf_G, n = Inf)$table
degs_G <- results_G[results_G$FDR < 0.05 & abs(results_G$logFC) > 0.75, ]
write.csv(degs_G, "/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/DEGs/DEGs_G_minus_vs_G_plus.csv")

results_C <- topTags(qlf_C, n = Inf)$table
degs_C <- results_C[results_C$FDR < 0.05 & abs(results_C$logFC) > 0.75, ]
write.csv(degs_C, "/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/DEGs/DEGs_C_minus_vs_C_plus.csv")

results_F <- topTags(qlf_F, n = Inf)$table
degs_F <- results_F[results_F$FDR < 0.05 & abs(results_F$logFC) > 0.75, ]
write.csv(degs_F, "/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/DEGs/DEGs_F_minus_vs_F_plus.csv")
print('script is complete')