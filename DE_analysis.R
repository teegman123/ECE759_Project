library(edgeR)

# Load the data from csv file
git_path = '/Users/teaguemcc/Programs/git_clones/ECE759_Project'
data <- read.csv(git_path+'/cellular_clarity/read_count_matrix.csv', row.names = 1)
data_matrix <- as.matrix(t(data)) # Transpose where columns are samples and rows are genes for edgeR
sample_names <- read.csv(git_path+'/cellular_clarity/sample_names.csv')

# Create DGEList Object for edgeR
dge <- DGEList(counts = data_matrix)