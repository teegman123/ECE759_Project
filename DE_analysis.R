library(edgeR)

# Load the data from csv file
data <- read.csv('/Users/teaguemcc/Programs/git_clones/ECE759_Project/cellular_clarity/read_count_matrix.csv', row.names = 1)
data_matrix <- as.matrix(t(data)) # Transpose where columns are samples and rows are genes for edgeR

# Create DGEList Object for edgeR
dge <- DGEList(counts = data_matrix)