#!/usr/bin/Rscript

library('limma')

probe.file <- 'NBB.Sample.Probe.txt'
sample.file <- 'NBB.Samplesheet.csv'

expr <- read.ilmn(files=probe.file)

expressed.x <- rowSums(expr$other$Detection < 0.05) >= 3
expressed.y <- colSums(expr$other$Detection < 0.05) >= 3
filtered <- expr[expressed.x,expressed.y]
expr.neqc <- neqc(filtered)

samples <- read.delim(sample.file, skip=7, header=T, sep=',', stringsAsFactors=F)

expr$E[1:5,1:5]
samples[1:5,1:5]

#save(expr.neqc, samples, file='NBB.neqc.RData')

intersection <- intersect(colnames(expr.neqc), samples$Sample_Name)
expr.neqc.int <- expr.neqc[,colnames(expr.neqc) %in% intersection]
samples.int <- samples[samples$Sample_Name %in% intersection,]

samples.int <<- samples.int[match(colnames(expr.neqc.int$E), samples.int$Sample_Name) ,]

stopifnot(nrow(samples.int) == ncol(expr.neqc.int))
avg.expr.neqc.int <<- avereps(expr.neqc.int, ID=expr.neqc.int$genes$SYMBOL)

groups <- as.factor(samples.int$Sample_Group)
design <- model.matrix(~0+groups)

fit <- lmFit(expr.neqc.int, design)
contrasts <- makeContrasts(groupsAD-groupsCtrl, levels=design)

groups <- as.factor(samples.int$ApoE)
design <- model.matrix(~0+groups)


new <- expr.neqc.int[,-c(301,302,303,304,305,306,307,308)]
fit <- lmFit(new, design)

#22 23 24 33 34 44
contrasts <- makeContrasts((groups44+groups34+groups24)/3-(groups33+groups23+groups22)/3, levels=design)#,
contrasts <- makeContrasts(groups44-(groups34+groups24+groups33+groups23+groups22)/5, levels=design)#, 

contrasts

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=T)
summary(decideTests(fit2, method='global'))

topTable(fit2, coef=1, p.value=0.05, number=10)


genes = c('TSPO','BGN', 'GNA01', 'RTN1')
genes <- intersect(genes, rownames(avg.expr.neqc.int$E))

bplot <- function(gene) {
	png(file=paste(gene,'png',sep='.'))
	boxplot(avg.expr.neqc.int$E[gene,1:300] ~ as.factor(samples.int$Braak[1:300]))
	points(avg.expr.neqc.int$E[gene,1:300] ~ as.factor(samples.int$Braak[1:300]))
	dev.off()
}

for (gene in genes) {
	bplot(gene)
}

correlations <- apply(avg.expr.neqc.int$E, 1, function (g) { cor(avg.expr.neqc.int$E['ADAM17',], g) })
