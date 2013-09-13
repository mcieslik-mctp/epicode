X11.options(type="nbcairo")
options(stringsAsFactors = FALSE)
options(max.print=100)
library(ggplot2)
library(reshape2)
library(stringr)
library(rjson)

args = commandArgs(FALSE)
script.name = str_split(tail(str_split(args[8], "/")[[1]], n=1), "\\.")[[1]][1]
args = args[-1:-9]

if (any(args!="")) {
  file.name = args[1]
  plot.name = args[2]
}

if (is.na(plot.name)) {
  plot.name = str_replace(file.name,
    str_sub(file.name, nchar(file.name)-2, nchar(file.name)), "png")
}

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
