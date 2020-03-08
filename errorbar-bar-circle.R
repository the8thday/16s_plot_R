#久生

library(readxl)
library(tidyverse)
library(Hmisc)

dd <- read_excel('/Users/congliu/Work/GEOCHIP.xlsx', sheet = 2)
dim(dd)
dd[!complete.cases(dd),]
s1 <- group_by(dd, Subcategory1) %>% summarise(n())
t <- filter(dd, Subcategory1=='Aromatics')
cal_sd <- function(variables) {
  df <- variables
  df[is.na(df)] <- 0
  rf <- df[8:14]
  md <- df[15:21]
  for (i in seq(dim(rf)[1])) {
    rf1 <- rf[i,]
    md1 <- md[i,]
    hehe <- t.test(as.numeric(rf1), as.numeric(md1))
    print(hehe)
  }
}
cal_sd(t)
for (i in s1[,1]) {
  t <- filter(dd, Subcategory1==i)
  cal_sd(t)
}
#################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # 计算长度
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # 以 groupvars 为组,计算每组的长度,均值,以及标准差
  # ddply 就是 dplyr 中的 group_by + summarise
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # 重命名  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  # 计算标准偏差
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  # 计算置信区间
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
###画图，黑白柱状图，加bar
df <- read_tsv('df.tsv')
df <- select(df, -X1)
chao <- select(df, c('Sample','chao1','Type'))
shannon <- select(df, c('Sample','shannon','Type'))
simpson <- select(df, c('Sample','simpson','Type'))
aa <- summarySE(df, measurevar=c("chao1"), groupvars=c("Type"))
aa$Type <- factor(aa$Type)
simpson <- ggplot(data = aa, aes(x=Type, y=simpson)) + geom_bar(stat = "identity",position=position_dodge(), fill=c('white', 'gray', 'black'), color='black') + 
  geom_errorbar(aes(ymin=simpson-se, ymax=simpson+se), position=position_dodge(0.9), width=.2) +
  xlab('simpson')
shannon <- ggplot(data = aa, aes(x=Type, y=shannon)) + geom_bar(stat = "identity",position=position_dodge(), fill=c('white', 'gray', 'black'), color='black') + 
  geom_errorbar(aes(ymin=shannon-se, ymax=shannon+se), position=position_dodge(0.9), width=.2) +
  xlab('shannon')
chao_p <- ggplot(data = aa, aes(x=Type, y=chao1)) + geom_bar(stat = "identity",position=position_dodge(), fill=c('white', 'gray', 'black'), color='black') + 
  geom_errorbar(aes(ymin=chao1-se, ymax=chao1+se), position=position_dodge(0.9), width=.2) +
  xlab('chao1')

summarySE(chao_p, shannon, simpson, axis=2)
pp <- function(arg){
  aa <- summarySE(df, measurevar=c(arg), groupvars=c("Type"))
  aa$Type <- factor(aa$Type)
  ggplot(data = aa, aes(x=Type, y=arg)) + geom_bar(stat = "identity",position=position_dodge(), fill=c('white', 'gray', 'black'), color='black') + 
    geom_errorbar(aes(ymin=arg-se, ymax=arg+se), position=position_dodge(0.9), width=.2) +
    xlab(arg)
}
pp('')
##丰度图
zryh <- read_delim('/Users/congliu/prs/fenxi/zryh.tsv', delim = '\t')
zryh <- select(zryh, -c('X1','Description'))
p <- reshape2::melt(zryh, variable.name="Sample")
ggplot(data = p, aes(x=Sample, y=value, fill=Type)) + geom_bar(stat = "identity",position = "dodge")+
  ylab('Relative Abundance') + coord_flip() + theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + xlab('')

#PLDA图
library(tidyverse)
plda <- read_delim('/Users/congliu/prs/fenxi/meta/pls-da_result.txt', delim = '\t', col_names = F)
pld <- plda[c('X2','X3','X4')]
ggplot(data = pld) + geom_point(aes(x=pld$X3, y=pld$X4, color = pld$X2, shape=pld$X2)) +
  stat_ellipse(aes(x=pld$X3, y=pld$X4, color = pld$X2), type = "t") +
  theme(legend.position = 'right', legend.title = element_text(''))+
  xlab('COMPONENT 1') +
  ylab('COMPONENT 2') +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  scale_shape_manual(values=c(21, 22, 23, 24)) + labs(color='') + labs(shape='')
  #scale_color_manual()
  

library(ggord)
library(yyplot)


#OR值图
#Odds Ratio = (A / B)/(C / D) = A D/ B C
#95% CI of ln(OR) = ln(OR)±1.96(1/A + 1/B + 1/C + 1/D)0.5
#95% CI of OR = e95% CI of ln(OR)
or <- read_delim('/Users/congliu/prs/fenxi/meta/diff/H_CIA__vs__H_MTX_meta_diff.txt', delim = '\t', col_names = T)
odds <- or %>% select(OTU, oddsRatio, lower, upper, pvalues, `+samples in group 0`, `+samples in group 1`)
#odds <- or %>% select(OTU, oddsRatio)
odds_f <- odds %>% filter(oddsRatio != Inf, oddsRatio != 0, pvalues <= 0.05)
odds_f <- odds_f %>% arrange(desc(`+samples in group 0`), desc(`+samples in group 1`))
odds_f <- odds_f %>% filter(OTU %in% c('s__Bacteroides_finegoldi','f__Bacteroidaceae','s__Enterococcus_faecalis',
                                     's__Lactobacillus_kitasatonis','s__Lactococcus_garvieae','f__Streptococcaceae',
                                     'o__Lactobacillales','s__Clostridium_isatidis','s__Paeniclostridium_sordelli',
                                     'f__Helicobacteraceae', 'o__Campylobacterales'))
odds_l <- tidyr::gather(odds_f, num, value, -OTU)
odds_f <- slice(odds_f, 1:10)
ggplot(data = odds_f) +
  geom_bar(aes(y = oddsRatio, x = OTU), stat = "identity",position = "dodge") +
  xlab('Species') + ylab('or_value')
  geom_errorbar(aes(ymin=lower, ymax=upper, x = OTU), position=position_dodge(0.9), width=.2)







##丰度图
paths = list.files('/Users/congliu/prs/fenxi/huetu')
setwd('/Users/congliu/prs/fenxi/huetu')
foo <- function(path){
  name <- str_split(basename(path), '[.]')[[1]][1]
  zryh <- read_delim(path, delim = '\t')
  p <- select(zryh, -c('X1'))
  p <- reshape2::melt(p, variable.name="Sample")
  ggplot(data = p, aes(x=Sample, y=value, fill=group)) + geom_bar(stat = "identity",position = "dodge")+
    ylab('Relative Abundance') + coord_flip() + theme_bw()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + xlab('')
  #ggsave(pp, filename = paste('/Users/congliu/prs/fenxi/lefse/tupian/',name, '.pdf', sep = ''), width=8, height=8)
  #pdf(paste('/Users/congliu/prs/fenxi/huetu/',name, '.pdf', sep = ''))
}

for (i in paths) {
  foo(i)
}

foo("lefse_input_H_CIA_H_Control.txt")
foo("lefse_input_H_CIA_H_Control_H_MTX_H_HUMSCs.txt")
foo("lefse_input_H_CIA_H_MTX.txt")
foo("lefse_input_H_CIA_M_CIA.txt")
foo("lefse_input_H_HUMSCs_H_CIA.txt")
foo("lefse_input_H_HUMSCs_H_CIA_H_Control_H_MTX.txt")
foo("lefse_input_H_HUMSCs_H_Control.txt")
foo("lefse_input_H_HUMSCs_H_MTX.txt")
foo("lefse_input_H_HUMSCs_M_HUMSCs.txt")
foo("lefse_input_H_MTX_H_Control.txt")
foo("lefse_input_H_MTX_M_MTX.txt")
foo("lefse_input_M_Control_H_Control.txt")
foo("lefse_input_M_Control_M_CIA.txt")
foo("lefse_input_M_Control_M_HUMSCs.txt")
foo("lefse_input_M_Control_M_MTX.txt")
foo("lefse_input_M_HUMSCs_M_CIA.txt")
foo("lefse_input_M_HUMSCs_M_MTX.txt")
foo("lefse_input_M_MTX_M_CIA.txt")
foo("lefse_input_M_MTX_M_Control_M_HUMSCs_M_CIA.txt")

foo("lefse_input_H_CIA_H_Control_H_MTX_H_HUMSCs.txt")
foo("lefse_input_H_HUMSCs_H_CIA.txt")
foo("lefse_input_M_HUMSCs_M_CIA.txt")














