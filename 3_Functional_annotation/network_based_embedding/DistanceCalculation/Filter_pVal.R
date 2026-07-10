suppressMessages(library(data.table))


args = commandArgs(trailingOnly=TRUE)

tem <- fread(args[1], header = T) 

tem <- subset(tem, classic <= 0.05)

tem[,'FDR'] <- p.adjust(tem$classic, method = 'fdr')
tem <- as.data.frame(subset(tem, FDR <= 0.05))
print(tem)