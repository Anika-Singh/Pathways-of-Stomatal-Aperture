rownames(GSM1272576_sample_table)[0] <- "id"


row <- rownames(GSM1272576_sample_table)

total <- merge( GSM1272576_sample_table,GSM1272577_sample_table,by="ID_REF")
total <- merge(  total, GSM1272578_sample_table,by="ID_REF")
total <- merge( total, GSM1272579_sample_table,by="ID_REF")
total <- merge( total, GSM1272580_sample_table, by="ID_REF")
total <- merge( total,GSM1272581_sample_table,by="ID_REF")


colnames(total) <- c("gene_id","root.MES.rep1","root.MES.rep2","root.MES.rep3" , "root.ATP.rep1" , "root.ATP.rep2" , "root.ATP.rep3" )

total <- merge( total,Correspondance,by="gene_id")

colnames(Correspondance) <- c("gene_id", "gene")

rownames(total) <- total$gene

total$gene <- NULL

total$gene_id <- NULL

total[,7] <- NULL
source("http://bioconductor.org/biocLite.R")

install.packages("limma")

biocLite("limma")
library("limma")


root.MES.mean = apply(total[, c("root.MES.rep1", "root.MES.rep2" , "root.MES.rep3")], 1, mean)

root.ATP.mean = apply(total[, c("root.ATP.rep1", "root.ATP.rep2", "root.ATP.rep3")], 1, mean)

df = data.frame(root.ATP.mean, root.MES.mean)

trmean = apply(df, 2, mean, trim=0.02)

mean.of.trmeans = mean(trmean)

exprSet.trmean = df / trmean * mean.of.trmeans

root.MES.mean = apply(root.MES.mean, 2, mean, trim=0.02)

plotMA(total)

total = exprSet.trmean

difference = total$root.ATP.mean - total$root.MES.mean

total = cbind(total, difference)

first.day.MIT = total

exprSet <- total

dataset.1 <- exprSet[1, c("brain.1", "brain.2")]

dataset.2 <- exprSet[1, c("brain.1", "brain.2")]

