#配对图线段图
library(ggplot2)

paired_dat<-read.table('paired_ZNF280C_double.xls',sep = '\t',header = T)
paired_dat$ZNF280C_normal=log(paired_dat$ZNF280C_normal+1)
paired_dat$ZNF280C_tumor=log(paired_dat$ZNF280C_tumor+1)
paired_dat$normal<-2
paired_dat$tumor<-1
paired_dat$Normal=as.character(paired_dat$Normal)
paired_dat$Tumor=as.character(paired_dat$Tumor)

paired_dat2<- data.frame(group=c(paired_dat$Normal,paired_dat$Tumor),ZNF280C=c(paired_dat$ZNF280C_normal,paired_dat$ZNF280C_tumor))

paired_dat3<-data.frame(x=paired_dat$tumor,xend=paired_dat$normal,
                        y=paired_dat$ZNF280C_normal,yend=paired_dat$ZNF280C_tumor)

pdf("paired_ZNF280C.pdf",family="GB1")
ggplot(paired_dat2,aes(group,ZNF280C,color=group))+
  geom_point(size=3)+
  scale_color_manual(values=c("grey","red"))+
  geom_segment(data = paired_dat3,aes(x=x,xend=xend,y=y,yend=yend),colour= "black",size=1,linetype=1)+
  annotate(geom="text",x=1,y=3.5,label="p-value = 0.0004679\nn=41",size=5,alpha=1,hjust=0.5,vjust=1.5)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))

dev.off()


###配对线段图
paried <- data.table::fread("~/Desktop/mllt1_paired.txt",sep = "\t",header = T)
paried$adenoma <-1
paried$normal <-2
paired2<- data.frame(
  tissue=c(paried$tissuenormal,paried$tissueadenoma),
  mllt1 = c(paried$MLLT1_1normal,paried$MLLT1_1adenoma))

paired3 <- data.frame(
  x=paried$adenoma,
  xend=paried$normal,
  y= paried$MLLT1_1adenoma,
  yend=paried$MLLT1_1normal
)

library(ggplot2)
ggplot(paired2,aes(tissue,mllt1))+
  geom_point()+
  geom_segment(data=paired3,aes(x=x,xend=xend,y=y,yend=yend))



###library(ggalluvial)
library(patchwork)
library(d3Network)
library(Networkd3)
peidui <- read_delim('/Users/congliu/prs/fenxi/peidui/otu3.txt', delim = '\t')
#peidui <- read_delim('/Users/congliu/prs/fenxi/peidui/ff3.txt', delim = '\t')
input <- peidui %>% filter(Type %in% c('normal', 'cancer', '1000000'))
#input <- peidui %>% select(-X1)
map_file <- read_delim('/Users/congliu/prs/fenxi/peidui/map.txt', delim = '\t') %>% select(-c('#Sample','BarcodeSequence', 'LinkerPrimerSequence'))
#input <- reshape2::melt(input)
#input$Type <- factor(input$Type, levels = c('normal', 'cancer'))
#input <- within(input,
#                Type <- factor(Type,levels=names(sort(table(Type), decreasing=TRUE))))
input <- within(input,
                Type <- factor(Type, levels = c('normal', 'cancer', '1000000')))
#input2 <- filter(input, Description=='1000000')
#ggplot(data = input, aes(x = key, y = value, fill=Type)) + geom_bar(stat = 'identity')
cbPalette <- c("#8B0000","#E69F00", "#56B4E9", "#009E73","#800080", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#A0522D","#999999")

ggplot() + geom_bar(data = input, aes(x = Type, y = value, fill=variable),stat = 'identity',position="fill")+
  xlab('') + ylab('相对比例') + scale_x_discrete(labels = c('健康人群','结直肠肿瘤人群','您的样本')) +
  scale_fill_manual(values=cbPalette) +
  labs(fill="菌谱标记")





theme_grey(base_family = 'STKaiti')
theme(text = element_text(family='Kai'))
theme(legend.title=element_text('aa'))
  

tuiji <- function(x){
  peidui <- read_delim('/Users/congliu/prs/fenxi/peidui/otu3.txt', delim = '\t')
  map_file <- read_delim('/Users/congliu/prs/fenxi/peidui/map.txt', delim = '\t') %>% select(-c('#Sample','BarcodeSequence', 'LinkerPrimerSequence'))
  name <- as.character(filter(map_file, Description == x)[1,1])
  input <- peidui %>% filter(Type %in% c('normal', 'cancer', x))
  input <- within(input,
                  Type <- factor(Type, levels = c('normal', 'cancer', x)))
  pp <- ggplot() + geom_bar(data = input, aes(x = Type, y = value, fill=variable),stat = 'identity',position="fill")+
    xlab('') + ylab('相对比例') + scale_x_discrete(labels = c('normal'='健康人群','cancer'='结直肠肿瘤人群',x='您的样本')) +
    scale_fill_manual(values=cbPalette) +
    labs(fill="菌谱标记")
  ggsave(pp, filename = paste('/Users/congliu/prs/fenxi/peidui/',name, '_', x, '.png', sep = ''))
}

tuiji('1000000')


for(i in names(table(peidui$Type))){
  if(i %in% c('normal', 'cancer', 'adenoma')){
    next()
  } else {
    tuiji(i)
  }
}





##桑吉图
library(ggalluvial)
map_df()
sanky <- read_tsv('/Users/congliu/prs/fenxi/peidui/sanky.txt')
ggplot(sanky,
       aes(y = value,
           axis1 = normal, axis2 = adenima, axis3 = cancer,
           fill = variable)) +
  geom_alluvium() +
  scale_x_discrete(limits = c("normal", "adnoma", "cancer"))


























