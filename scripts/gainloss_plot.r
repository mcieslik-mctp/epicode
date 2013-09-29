library(ggplot2)
library(reshape2)
library(stringr)
library(rjson)

args = commandArgs(TRUE)

file.name = args[1]
plot.name = args[2]

tbl = read.table(file.name, header=TRUE, sep="\t")
tbl$code = factor(1:nrow(tbl))
tbl =melt(tbl, id.vars="code")
colnames(tbl) = c("code", "mark", "value")

tbl$dir = ""
tbl[grep(".g", tbl$mark), "dir"] = "gain" 
tbl[grep(".l", tbl$mark), "dir"] = "loss" 
tbl$mark = str_replace(tbl$mark, ".g", "")
tbl$mark = str_replace(tbl$mark, ".l", "")

plt = ggplot(tbl) +
  aes(x=mark, y=value, fill=mark) + facet_grid(code ~ dir) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=-90, hjust=0.5, vjust=0.5), axis.title = element_blank())

ggsave(plot.name, plt)
