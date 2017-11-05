source("http://bioconductor.org/biocLite.R")
biocLite("affy")
library("affy")

#Reading the data 
affy.data = ReadAffy()
eset.mas5 = mas5(affy.data)
exprSet.nologs = exprs(eset.mas5)
colnames(exprSet.nologs)

#Logarithm of the Affymetrix data
exprSet = log(exprSet.nologs, 2)
exprSet <- as.data.frame(exprSet)

colnames(exprSet) = c("Root.MES.1", "Root.MES.2", 
                             "Root.MES.3", "Root.ATP.1",
                             "Root.ATP.2", "Root.ATP.3", 
                             "liver.1", "liver.2", "delete", "s" , "wq" ,"q" )
exprSet$liver.1 <- NULL
exprSet$liver.2 <- NULL
exprSet$s <- NULL
exprSet$q <- NULL
exprSet$wq <- NULL
exprSet$delete <- NULL

#RMA analysis 
# 
# total <- merge( GSM1272576_sample_table,GSM1272577_sample_table,by="ID_REF")
# total <- merge(  total, GSM1272578_sample_table,by="ID_REF")
# total <- merge( total, GSM1272579_sample_table,by="ID_REF")
# total <- merge( total, GSM1272580_sample_table, by="ID_REF")
# total <- merge( total,GSM1272581_sample_table,by="ID_REF")

# Get the actual A/P calls
data.mas5calls = mas5calls(affy.data)
data.mas5calls.calls = exprs(data.mas5calls)


# total <- exprSet
# 
# colnames(total) <- c("gene_id","root.MES.rep1","root.MES.rep2","root.MES.rep3" , "root.ATP.rep1" , "root.ATP.rep2" , "root.ATP.rep3" )
# 
# total <- merge( total,Correspondance,by="gene_id")
# 
# colnames(Correspondance) <- c("gene_id", "gene")

# rownames(total) <- total$gene

# total$gene <- NULL
# 
# total$gene_id <- NULL
# 
# total[,7] <- NULL

biocLite("limma")
library("annotate")

#Normalization

total <- exprSet

root.MES.mean = apply(total[, c("Root.MES.1", "Root.MES.2" , "Root.MES.3")], 1, mean)
root.ATP.mean = apply(total[, c("Root.ATP.1", "Root.ATP.2", "Root.ATP.3")], 1, mean)
df = data.frame(root.ATP.mean, root.MES.mean)
trmean = apply(df, 2, mean, trim=0.02)
mean.of.trmeans = mean(trmean)
exprSet = exprSet / trmean * mean.of.trmeans
plotMA(exprSet)
exprSet = log(total, 2)
difference = df$root.ATP.mean - df$root.MES.mean
total = cbind(df, difference)
first.day.MIT = total

#Identifying differentially expressed genes
dataset.1 <- exprSet[1, c("Root.MES.1", "Root.MES.2", "Root.MES.3")]
dataset.2 <- exprSet[1, c("Root.ATP.1", "Root.ATP.2", "Root.ATP.3")]

#Get t test for 1st gene
t.test.gene.1 = t.test(dataset.1, dataset.2, "two.sided")
t.test.gene.1$p.value
#t test for the rest of the genes
ATP.p.value.all.genes = apply(exprSet, 1, function(x) { t.test(x[1:3], x[4:7]) $p.value } )
AP = apply(data.mas5calls.calls, 1, paste, collapse="")
exprSet.present = exprSet[genes.present,]
ATP.raw.pvals.present = ATP.p.value.all.genes[genes.present]
ATP.fdr.pvals.present = p.adjust(ATP.raw.pvals.present, method="fdr")
ATP.fdr.pvals.present.sorted = ATP.fdr.pvals.present[order(ATP.fdr.pvals.present)]

expression.plus.pvals = cbind(exprSet.present, ATP.raw.pvals.present, 
                              ATP.fdr.pvals.present)

ATP.DE.probesets = names(ATP.raw.pvals.present[ATP.raw.pvals.present < 0.01])

all.data = first.day.MIT
ss = rownames(all.data)
volcano = cbind(all.data[, c("difference")],ATP.raw.pvals.present)

ATP.DE.log2.ratios <- as.data.frame(ATP.DE.log2.ratios )[which(data > 0)]

rownames(ATP.DE.log2.ratios) <- ATP.DE.probesets

colnames(ATP.DE.log2.ratios) <- "data"

#Plotting the data
write.table(ATP.DE.log2.ratios, "ATP.DE.log2.ratios.txt", sep="\t", quote=F)

root = all.data[, "root.MES.mean"]
ATP.root = all.data[,  "root.ATP.mean"]
A = (root + ATP.root) / 2
M = ATP.root - root
plot(A, M)
plot(A, M, main="MA plot of ATP treatment vs Wild type", pch=19, cex=0.2, col="red")

log2.ratios = first.day.MIT$difference
log2.ratios = expression.plus.pvals[, "Root.ATP.1"] -  expression.plus.pvals[, "Root.MES.1"]
p.values = expression.plus.pvals[, "ATP.raw.pvals.present"]

sep1 = log2(log2.ratios)


plot(log2.ratios, -log(p.values, 10) ,pch=19, cex=0.2, col="red")


suppressMessages(library("ath1121501.db"))
getEG(ATP.DE.probesets, "ath1121501.db")

install.packages("ggplot2")
library("heatmap.2")
heatmap(as.matrix(all.data))
as.matrix(all.data)

write.table(all.data, "data", sep="\t", quote=F)

#Plotting the differentially expressed sets
heatmap(as.matrix(exprSet[ 
                          c("Root.MES.1","Root.MES.2","Root.MES.3","Root.ATP.1","Root.ATP.2","Root.ATP.3")]
))


cool = head(exprSet)

heatmap(as.matrix(exprSet[ATP.DE.probesets, c("Root.MES.1","Root.MES.2","Root.MES.3","Root.ATP.1","Root.ATP.2","Root.ATP.3")]))


rownames(exprSet)


source("http://www.bioconductor.org/biocLite.R") 
biocLite("limma") 
biocLite("affy") 
biocLite("hgu95av2") 
biocLite("estrogen") 
biocLite("ath1121501.db") 
biocLite("simpleaffy") 
biocLite("annotate") 
biocLite("XML")

ll <- as.matrix(getEG(ATP.DE.probesets, "ath1121501.db"))


p <- plot_ly(data = exprSet, x = root.MES.mean, y = -log10(root.MES.mean), mode = "markers")


