# Size of the figures should be a single (86 mm) or double (178 mm) column width.
library(ggplot2)

# fileList <- list.files('Results/', full.names = TRUE)
# fContent <- lapply(fileList, function(X){
#   message(X)
#   try(read.csv(X), silent = TRUE)
# })
# fContent <- fContent[unlist(lapply(fContent, class)) != 'try-error']
# fContent <- do.call(rbind.data.frame, fContent)
# save(fContent, file = 'FileContent.RData')
load('FileContent.RData')
# Number of cells
nrow(fContent)
# Number of experiments
length(levels(fContent$SRS))
# Number of projects
length(levels(fContent$SRA))



MT <- fContent
table(MT$MTCOUNT != 0)
MT <- MT[MT$MTCOUNT != 0,]
MT <- as.data.frame.array(MT)
MT <- MT[order(MT$TOTALCOUNT),]
MT$log10MTCOUNT <- log10(MT$MTCOUNT)
MT$log10TOTALCOUNT <- log10(MT$TOTALCOUNT)
MT$log10NGENES <- log10(MT$NGENES)

expectedGENE <- predict(lm(MT$log10NGENES~poly(MT$log10TOTALCOUNT,2)), newdata = data.frame(TOTALCOUNT = MT$log10TOTALCOUNT), interval = 'prediction', level = 0.95)
expectedGENE <- as.data.frame(expectedGENE)
MT$expGENElwr <- expectedGENE$lwr
MT$expGENE <- expectedGENE$fit
MT$expGENEupr <- expectedGENE$upr

MT$CAT <- paste0(MT$PROTOCOL, ' - ', MT$SP)
CORVAL <- cor(MT$NGENES,MT$TOTALCOUNT, method = 'sp')
cor.test(MT$NGENES,MT$TOTALCOUNT, method = 'sp')

# p <- ggplot(MT, aes(x=log10TOTALCOUNT, y=log10NGENES)) + 
#   geom_point(cex = 0.1, alpha = 0.05) + xlab(parse(text = expression('log[10](Total~Counts)'))) +
#   geom_smooth(method=lm , color="red", se = FALSE, formula = y ~poly(x,2)) +
#   ylab(parse(text = expression('log[10](Total~Number~of~Genes)'))) +
#   geom_line(aes(y=expGENElwr), color = "red", linetype = "dashed")+
#   geom_line(aes(y=expGENEupr), color = "red", linetype = "dashed")+
#   theme_bw() + labs(subtitle = parse(text = paste0('rho==', round(CORVAL,2))))
# 
# png('Figures/FigureS1.png', width = 86, height = 86, res = 200,  units = 'mm')
# print(p)
# dev.off()

table(MT$log10NGENES > MT$expGENElwr & MT$log10NGENES < MT$expGENEupr)
table(MT$log10NGENES > MT$expGENElwr)
table(MT$log10NGENES < MT$expGENEupr)
MT <- MT[MT$log10NGENES > MT$expGENElwr & MT$log10NGENES < MT$expGENEupr,]
expectedMT <- predict(lm(MT$log10MTCOUNT~MT$log10TOTALCOUNT), newdata = data.frame(TOTALCOUNT = MT$log10TOTALCOUNT), interval = 'prediction', level = 0.95)
expectedMT <- as.data.frame(expectedMT)
MT$expMTCOUNTlwr <- expectedMT$lwr
MT$expMTCOUNT <- expectedMT$fit
MT$expMTCOUNTupr <- expectedMT$upr

CORVAL <- cor(MT$log10TOTALCOUNT, MT$log10MTCOUNT)
cor.test(MT$log10TOTALCOUNT, MT$log10MTCOUNT)
# p <- ggplot(MT, aes(x=log10TOTALCOUNT, y=log10MTCOUNT)) +
#   geom_point(cex = 0.1, alpha = 0.05) + geom_smooth(method=lm , color="red", se = FALSE) +
#   theme_bw() + geom_line(aes(y=expMTCOUNTlwr), color = "red", linetype = "dashed")+
#   geom_line(aes(y=expMTCOUNTupr), color = "red", linetype = "dashed")+
#   xlab(parse(text = expression('log[10](Total~Counts)')))+
#   ylab(parse(text = expression('log[10](Mitochondrial~Counts)')))+
#   labs(subtitle = parse(text = paste0('r ==', round(CORVAL,2))))
# 
# png('Figures/Figure1.png', width = 86, height = 86, res = 200,  units = 'mm')
# print(p)
# dev.off()

table(MT$log10MTCOUNT > MT$expMTCOUNTlwr & MT$log10MTCOUNT < MT$expMTCOUNTupr)
table(MT$log10MTCOUNT > MT$expMTCOUNTlwr)
table(MT$log10MTCOUNT < MT$expMTCOUNTupr)
MT <- MT[MT$log10MTCOUNT > MT$expMTCOUNTlwr & MT$log10MTCOUNT < MT$expMTCOUNTupr,]


p <- ggplot(MT, aes(x=SP, y=MTRATIO)) +  
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() + 
  coord_flip() + geom_jitter(shape=16, position=position_jitter(0.1), cex = 0.05, alpha = 0.01) +
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Specie') +
  ylab('Mitochondrial Ratio')

png('Figures/Figure2.png', width = 86, height = 40, res = 200,  units = 'mm')
print(p)
dev.off()

wilcox.test(MTRATIO~SP, MT, alternative = 'greater')
t.test(MTRATIO~SP, MT, alternative = 'greater')




p <- ggplot(MT, aes(x=CAT, y=MTRATIO)) +  
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Technology and Specie') +
  ylab('Mitochondrial Ratio') 

png('Figures/Figure3.png', width = 86, height = 86, res = 200,  units = 'mm')
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

png('Figures/FigureS2.png', width = 86*1.2, height = 86*2, res = 200,  units = 'mm')
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

png('Figures/FigureS3.png', width = 86*1.2, height = 86*3, res = 200,  units = 'mm')
print(p)
dev.off()

C10X <- MT[MT$PROTOCOL == '10x chromium' & MT$SP == 'Homo sapiens' &  MT$CT != 'Unknown',]
tTISSUE <- table(C10X$TISSUE)
C10X <- C10X[C10X$TISSUE %in% names(tTISSUE)[tTISSUE > 1000],]
C10X$TISSUE[C10X$TISSUE == 'Circulating tumor cells in hepatocellular carcinoma'] <- 'Hepatocellular carcinoma'
C10X$TISSUE[C10X$TISSUE == 'Peripheral blood mononuclear cell'] <- 'PBMC'
C10X$TISSUE[C10X$TISSUE == 'Peripheral blood mononuclear cells'] <- 'PBMC'
tT <- table(C10X$TISSUE)
C10X$TISSUE <- as.factor(C10X$TISSUE)
levels(C10X$TISSUE) <- paste0(names(tT), ' (', tT, ')')

p <- ggplot(C10X, aes(x=TISSUE, y=MTRATIO)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Tissue') +
  ylab('Mitochondrial Ratio')

png('Figures/Figure4.png', width = 86*1.3, height = 86*1.5, res = 200,  units = 'mm')
print(p)
dev.off()

pValues <- sapply(unique(C10X$TISSUE), function(X){
  t.test(C10X[C10X$TISSUE == X,]$MTRATIO, mu = 0.05, alternative = 'less')$p.value
})
sum(pValues > 0.05)
round(mean(pValues > 0.05)*100,1)
unique(C10X$TISSUE)[pValues > 0.05]

q75 <- sapply(unique(C10X$TISSUE), function(X){
  quantile(C10X[C10X$TISSUE == X,]$MTRATIO, 0.75)
})


pValues <- sapply(unique(C10X$TISSUE), function(X){
  t.test(C10X[C10X$TISSUE == X,]$MTRATIO, mu = 0.1, alternative = 'less')$p.value
})
sum(pValues > 0.05)
round(mean(pValues > 0.05)*100,1)
unique(C10X$TISSUE)[pValues > 0.05]

C10X <- MT[MT$PROTOCOL == '10x chromium' & MT$SP == 'Mus musculus' &  MT$CT != 'Unknown',]
tTISSUE <- table(C10X$TISSUE)
C10X <- C10X[C10X$TISSUE %in% names(tTISSUE)[tTISSUE > 1000],]
C10X$TISSUE[C10X$TISSUE == 'Circulating tumor cells in hepatocellular carcinoma'] <- 'Hepatocellular carcinoma'
C10X$TISSUE[C10X$TISSUE == 'Cortex 1'] <- 'Cortex'
C10X$TISSUE[C10X$TISSUE == 'Cortex 2'] <- 'Cortex'
C10X$TISSUE[C10X$TISSUE == 'Cortex 3'] <- 'Cortex'

tT <- table(C10X$TISSUE)
C10X$TISSUE <- as.factor(C10X$TISSUE)
levels(C10X$TISSUE) <- paste0(names(tT), ' (', tT, ')')

p <- ggplot(C10X, aes(x=TISSUE, y=MTRATIO)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw(base_size = 10) + coord_flip() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.05, alpha = 0.1) + 
  geom_hline(yintercept = 0.05, lty = 2, col = 'red') + 
  geom_hline(yintercept = 0.1, lty = 2, col = 'blue') +
  xlab('Tissue') +
  ylab('Mitochondrial Ratio') 

png('Figures/Figure5.png', width = 86*1.3, height = 86*3.2, res = 200,  units = 'mm')
print(p)
dev.off()


pValues <- sapply(unique(C10X$TISSUE), function(X){
  t.test(C10X[C10X$TISSUE == X,]$MTRATIO, mu = 0.05, alternative = 'less')$p.value
})
sum(pValues > 0.05)
round(mean(pValues > 0.05)*100,1)
unique(C10X$TISSUE)[pValues > 0.05]

round((table(MT$PROTOCOL)/nrow(MT))*100,1)


Hs <- density(MT[MT$SP == 'Homo sapiens',]$MTRATIO, from = 0)
Hs$y <- Hs$y/sum(Hs$y)
Mm <- density(MT[MT$SP == 'Mus musculus',]$MTRATIO, from = 0)
Mm$y <- Mm$y/sum(Mm$y)

png('Figures/FigureS0.png', width = 86, height = 86, res = 200,  units = 'mm')
par(mar=c(3,3,1,1), mgp = c(1.5,0.5,0))
plot(Hs, ylim = c(0,0.02), col = 'darkblue', xlab = 'Mitochondrial ratio', ylab = 'Probability', main = '')
lines(Mm, col = 'forestgreen')
abline(v = 0.05, col = 'red', lty = 2)
abline(v = 0.1, col = 'blue', lty = 2)
legend('topright', legend = c('Human', 'Mice'), col = c('darkblue', 'forestgreen'), lty = 1, bty = 'n', cex = 0.7)
dev.off()
