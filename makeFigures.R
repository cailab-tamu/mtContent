# Size of the figures should be a single (86 mm) or double (178 mm) column width.
library(ggplot2)
library(ggjoy)

# fileList <- list.files('Results/', full.names = TRUE)
# fContent <- lapply(fileList, function(X){
#   message(X)
#   try(read.csv(X), silent = TRUE)
# })
# fContent <- fContent[unlist(lapply(fContent, class)) != 'try-error']
# fContent <- do.call(rbind.data.frame, fContent)
# save(fContent, file = 'FileContent.RData')
load('FileContent.RData')

MT <- fContent
MT <- MT[MT$MTRATIO != 0,]
MT <- as.data.frame.array(MT)
MT <- MT[order(MT$TOTALCOUNT),]
MT$log10MTCOUNT <- log10(MT$MTCOUNT)
MT$log10TOTALCOUNT <- log10(MT$TOTALCOUNT)

expectedMT <- predict(lm(MT$log10MTCOUNT~MT$log10TOTALCOUNT), newdata = data.frame(TOTALCOUNT = MT$log10TOTALCOUNT), interval = 'prediction', level = 0.95)
expectedMT <- as.data.frame(expectedMT)

MT$expMTCOUNTlwr <- expectedMT$lwr
MT$expMTCOUNT <- expectedMT$fit
MT$expMTCOUNTupr <- expectedMT$upr
MT$CAT <- paste0(MT$PROTOCOL, ' - ', MT$SP)


p <- ggplot(MT, aes(x=log10TOTALCOUNT, y=log10MTCOUNT)) + 
  geom_point(cex = 0.1, alpha = 0.05) + geom_smooth(method=lm , color="red", se = FALSE) + 
  theme_bw() + geom_line(aes(y=expMTCOUNTlwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=expMTCOUNTupr), color = "red", linetype = "dashed")+
  xlab(parse(text = expression('log[10](Total~Counts)')))+
  ylab(parse(text = expression('log[10](Mitochondrial~Counts)')))

png('Figures/FigureS1.png', width = 86, height = 86, res = 200,  units = 'mm')
print(p)
dev.off()


MT <- MT[MT$log10MTCOUNT > MT$expMTCOUNTlwr & MT$log10MTCOUNT < MT$expMTCOUNTupr,]


p <- ggplot(MT, aes(x=SP, y=MTRATIO)) +  
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  coord_flip() + geom_jitter(shape=16, position=position_jitter(0.1), cex = 0.05, alpha = 0.01) +
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Specie') +
  ylab('Mitochondrial Ratio')

png('Figures/Figure1.png', width = 86, height = 40, res = 200,  units = 'mm')
print(p)
dev.off()




p <- ggplot(MT, aes(x=CAT, y=MTRATIO)) +  
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Technology and Specie') +
  ylab('Mitochondrial Ratio') 

png('Figures/Figure2.png', width = 86, height = 86, res = 200,  units = 'mm')
print(p)
dev.off()


C10X <- MT[MT$PROTOCOL == '10x chromium' & MT$SP == 'Homo sapiens' &  MT$CT != 'Unknown',]
C10X$CT[C10X$CT == 'Erythroid-like and erythroid precursor cells'] <- 'Erythroid precursor cells'
tCT <- table(C10X$CT)
C10X <- C10X[C10X$CT %in% names(tCT)[tCT > 1000],]
tCT <- table(C10X$CT)
C10X$CT <- as.factor(C10X$CT)
levels(C10X$CT) <- paste0(names(tCT), ' (', tCT, ')')

p <- ggplot(C10X, aes(x=CT, y=MTRATIO)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Cell Type') +
  ylab('Mitochondrial Ratio') 

png('Figures/Figure3.png', width = 86, height = 86*2, res = 200,  units = 'mm')
print(p)
dev.off()


C10X <- MT[MT$PROTOCOL == '10x chromium' & MT$SP == 'Mus musculus' &  MT$CT != 'Unknown',]
C10X$CT[C10X$CT == 'Erythroid-like and erythroid precursor cells'] <- 'Erythroid precursor cells'
tCT <- table(C10X$CT)
C10X <- C10X[C10X$CT %in% names(tCT)[tCT > 1000],]
tCT <- table(C10X$CT)
C10X$CT <- as.factor(C10X$CT)
levels(C10X$CT) <- paste0(names(tCT), ' (', tCT, ')')

p <- ggplot(C10X, aes(x=CT, y=MTRATIO)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Cell Type') +
  ylab('Mitochondrial Ratio') 

png('Figures/Figure4.png', width = 86, height = 86*2, res = 200,  units = 'mm')
print(p)
dev.off()

C10X <- MT[MT$PROTOCOL == '10x chromium' & MT$SP == 'Homo sapiens' &  MT$CT != 'Unknown',]
tTISSUE <- table(C10X$TISSUE)
C10X <- C10X[C10X$TISSUE %in% names(tTISSUE)[tTISSUE > 2000],]
C10X$TISSUE[C10X$TISSUE == 'Circulating tumor cells in hepatocellular carcinoma'] <- 'Hepatocellular carcinoma'

p <- ggplot(C10X, aes(x=TISSUE, y=MTRATIO)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Cell Type') +
  ylab('Mitochondrial Ratio')

png('Figures/Figure5.png', width = 86, height = 86*3, res = 200,  units = 'mm')
print(p)
dev.off()

C10X <- MT[MT$PROTOCOL == '10x chromium' & MT$SP == 'Mus musculus' &  MT$CT != 'Unknown',]
tTISSUE <- table(C10X$TISSUE)
C10X <- C10X[C10X$TISSUE %in% names(tTISSUE)[tTISSUE > 2000],]
C10X$TISSUE[C10X$TISSUE == 'Circulating tumor cells in hepatocellular carcinoma'] <- 'Hepatocellular carcinoma'

p <- ggplot(C10X, aes(x=TISSUE, y=MTRATIO)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Cell Type') +
  ylab('Mitochondrial Ratio') 

png('Figures/Figure6.png', width = 86, height = 86*3, res = 200,  units = 'mm')
print(p)
dev.off()
