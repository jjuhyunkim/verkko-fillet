library(SVbyEye)
library(stringr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

paf.file <- args[1]
outPrefix <- args[2]
minlen <- ifelse(length(args) >= 3, as.numeric(args[3]), 500000)

print(paf.file)
print(outPrefix)
print(minlen)

# Load the PAF file
paf.table <- read.csv(paf.file, header = FALSE, sep = '\t')

# Select relevant columns
paf.table <- paf.table[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12)]
colnames(paf.table) <- c('q.name', 'q.len', 'q.start', 'q.end', 'strand', 't.name', 't.len', 't.start', 't.end', 'n.match', 'mapq')

# Filter data based on query names
paf.table_dam <- paf.table[str_starts(paf.table$q.name, "dam"), ]
paf.table_sire <- paf.table[str_starts(paf.table$q.name, "sire"), ]

# Plot sire
plotGenome(paf.table = paf.table_sire, min.query.aligned.bp = minlen)
ggsave(filename = paste0(outPrefix, "_sire.png"))

# Plot dam
plotGenome(paf.table = paf.table_dam, min.query.aligned.bp = minlen)
ggsave(filename = paste0(outPrefix, "_dam.png"))

