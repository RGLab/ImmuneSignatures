library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra) # to combine plot and table

source("functions/inv_normal_transform.r")
source("functions/plot_size_bubbles.r")
source("functions/discretize.r")

fn = 'data/SDY404_combined_hai_titers.txt'
fn.y = 'data/SDY404_young_hai_titers.txt'
fn.o = 'data/SDY404_old_hai_titers.txt'

fn.out = sub(".tab$", "", fn)

# load data
df = read.table(fn, sep="\t", header=T,row.names=NULL, stringsAsFactors=F) %>%
  rename(IDs=subject)
df.y = read.table(fn.y, sep="\t", header=T,row.names=NULL, stringsAsFactors=F) %>%
  rename(IDs=subject)
df.o = read.table(fn.o, sep="\t", header=T,row.names=NULL, stringsAsFactors=F) %>%
  rename(IDs=subject)

df$ageGrp = NA
df$ageGrp[df$IDs %in% df.y$IDs] = "Young"
df$ageGrp[df$IDs %in% df.o$IDs] = "Older"

df.long = df %>%
  gather("virus_tp","hai",contains("_")) %>%
  separate(virus_tp, c("tp","virus"), sep="_", extra="merge")

# maximum across virus strains without normalization
df.max = df.long %>%
  group_by(IDs, tp) %>%
  summarise(hai_max = max(hai, na.rm=T)) %>%
  ungroup() %>%
  mutate(tp = paste0(tp,"_max")) %>%
  spread(tp,hai_max)

# response based on 4 fold change
df.4fc.resp = df.long %>%
  filter(tp=="fc") %>%
  group_by(IDs) %>%
  summarise(fc_4fc = cut(sum(hai>=4,na.rm=T),c(0:2,Inf), 0:2,right=F)) %>%
  ungroup()

# HAI normalization
df.norm = df.long %>%
  group_by(virus, tp) %>%
  mutate(hai_norm = (hai - median(hai, na.rm=T)) / sd(hai, na.rm=T)) %>%
  ungroup()

# maximum between normalized day0 and FC across virus strains 
# also calculate inver normal transformation and descretize to 20/80 and 30/70 quantiles
# fc_norm_max = MFC
df.norm.max = list()
df.norm.max$all = df.norm %>%
  group_by(IDs, tp) %>%
  summarise(hai_max = max(hai_norm, na.rm=T)) %>%
  ungroup() %>%
  mutate(tp = paste0(tp,"_norm_max")) %>%
  spread(tp,hai_max) %>%
  mutate(fc_norm_max_int = inv_normal_transform(fc_norm_max),
         fc_norm_max_d20 = discretize(fc_norm_max, 0.2),
         fc_norm_max_d30 = discretize(fc_norm_max, 0.3)) %>%
  inner_join(select(df, IDs,ageGrp), by="IDs")

# here we use all-normalized data to select young or older subjects  (but they can be processed individually)
df.norm.max$young = filter(df.norm.max$all, ageGrp=="Young")  
df.norm.max$older = filter(df.norm.max$all, ageGrp=="Older")


# decorrelation by binning along day0
# currently we bin manually (automatic binning is in development)
bins = list()
# SDY404
bins$all = c(1,5)
bins$young = c(1,5)
bins$older = c()

for (ag in c("all","young","older")) {
  
  # correlation before decorrelation
  cc = cor.test(df.norm.max[[ag]]$fc_norm_max_int, df.norm.max[[ag]]$d0_norm_max, method="spearman",exact=F)
  # axis limits
  xrng = range(df.norm.max[[ag]]$d0_norm_max, na.rm=T)
  yrng = range(df.norm.max[[ag]]$fc_norm_max_int, na.rm=T)
  xrng[1] = xrng[1]+diff(xrng)/2
  yrng[1] = yrng[1]+diff(yrng)/3*2
  
  # bin assignment and within-bin correlation
  df.norm.max[[ag]] = df.norm.max[[ag]] %>%
    mutate(bin = cut(d0_norm_max, breaks=c(-Inf, bins[[ag]], Inf), labels=1:(length(bins[[ag]])+1)) )
  df.bin = df.norm.max[[ag]] %>%
    group_by(bin) %>%
    summarise(n=n(), 
              fc_median = median(fc_norm_max_int, na.rm=T),
              cor = ifelse(n>2, cor(d0_norm_max, fc_norm_max_int, method="spearman", use="pairwise.complete.obs"), NA),
              cor.p = ifelse(n>2, cor.test(d0_norm_max, fc_norm_max_int, method="spearman", exact=F)$p.value, NA)) %>%
    ungroup()
  
  # table of bins on the plot
  tg = tableGrob(format(df.bin,digits=2), rows=NULL)
  
  # plot BEFORE decorrelation (one can go back and change bins)
  p = plot_size_bubbles(df.norm.max[[ag]]$d0_norm_max, df.norm.max[[ag]]$fc_norm_max_int) +
    scale_size_area(max_size = 8) + xlab("d0_norm_max") + ylab("fc_norm_max_ivt") +
    geom_step(data=data.frame(x=c(-Inf,bins[[ag]],Inf), y=c(df.bin$fc_median, df.bin$fc_median[length(df.bin$fc_median)])), col="red",lty=2) +
    ggtitle(sprintf("%s patients: n = %d, cor = %.2f, cor.p = %.2g",ag, nrow(df.norm.max[[ag]]), cc$estimate, cc$p.value)) +
    theme_bw() + theme(legend.key = element_blank())
  
  if (length(bins[[ag]])>0) {
    p = p + geom_vline(xintercept = bins[[ag]], col='red', lty=1) +
    annotation_custom(tableGrob(format(df.bin,digits=2), rows=NULL), xmin=xrng[1], xmax=xrng[2], 
                          ymin = yrng[1], ymax=yrng[2])
  }
  print(p)
  ggsave(sprintf("%s_%s_before_decorr.png", fn.out, ag), w=8,h=6)
}

for (ag in c("all","young","older")) {
  
  # decorrelated fold changes: fc_res_max = adjMFC
  # also descretized to 20/80 and 30/70 quantiles
  df.norm.max[[ag]] = df.norm.max[[ag]] %>%
    group_by(bin) %>%
    # limit number of subjects in the last bin (now only 1 will cause NA)
    mutate(fc_res_max = (fc_norm_max_int - median(fc_norm_max_int, na.rm=T)) / sd(fc_norm_max_int, na.rm=T)) %>%
    ungroup() %>%
    mutate(fc_res_max_d20 = discretize(fc_res_max, 0.2),
           fc_res_max_d30 = discretize(fc_res_max, 0.3))
  
  cc.dec = cor.test(df.norm.max[[ag]]$fc_res_max, df.norm.max[[ag]]$d0_norm_max, method="spearman",exact=F)
  df.dec.q = data.frame(y = quantile(df.norm.max[[ag]]$fc_res_max, c(0.2, 0.3, 0.7, 0.8), na.rm=T),
                        x = Inf, vjust = c(1.1,-0.1,1.1,-0.1), n = NA) %>%
    tibble::rownames_to_column("q")
  
  # compare responders
  df.dec.q$n[c(1,4)] = table(df.norm.max[[ag]]$fc_res_max_d20, exclude=1)
  df.dec.q$n[c(2,3)] = table(df.norm.max[[ag]]$fc_res_max_d30, exclude=1)
  
  # plot AFTER decorrelation 
  p = plot_size_bubbles(df.norm.max[[ag]]$d0_norm_max, df.norm.max[[ag]]$fc_res_max) +
    scale_size_area(max_size = 8) + xlab("d0_norm_max") + ylab("fc_norm_max_ivt") +
    geom_hline(yintercept = quantile(df.norm.max[[ag]]$fc_res_max, c(0.2, 0.8), na.rm=T), col="blue",lty=2) +
    geom_hline(yintercept = quantile(df.norm.max[[ag]]$fc_res_max, c(0.3, 0.7), na.rm=T), col="red",lty=2) +
    geom_text(data=df.dec.q, aes(x,y,label=sprintf("%s: n=%d",q,n), vjust=vjust), hjust=1) +
    ggtitle(sprintf("%s patients: n = %d, cor = %.2f, cor.p = %.2g",ag, nrow(df.norm.max[[ag]]), cc.dec$estimate, cc.dec$p.value)) +
    theme_bw() + theme(legend.key = element_blank())
  if (length(bins[[ag]])>0) {
    p = p + geom_vline(xintercept = bins[[ag]], col='red', lty=1)
  }
  print(p)
  ggsave(sprintf("%s_%s_after_decorr.png", fn.out, ag), w=8,h=6)
  
  # output all results
  df.out = inner_join(df, df.max, df.4fc.resp, by="IDs")
  df.out = inner_join(df.out, df.4fc.resp, by="IDs")
  df.out = inner_join(df.out, select(df.norm.max[[ag]], -ageGrp, -bin), by="IDs")
  
  write.table(df.out, file=sprintf("%s_%s_table.txt", fn.out, ag), sep="\t", row.names=F, col.names=T, quote=F)
}

# check if adjMFC correlates with baseline titer of indovodual virus strains
df.vir = df %>%
  select(-starts_with("fc")) %>%
  gather("virus","d0",-c(IDs,ageGrp)) %>%
  mutate(virus = sub("d0_","",virus)) %>%
  mutate(d0.log2 = log2(d0)) %>%
  inner_join(select(df.norm.max[[ag]], -ageGrp, -bin), by="IDs")

df.cc = df.vir %>%
  group_by(virus) %>%
  summarise(cc = cor.test(fc_res_max, d0.log2, method="spearman", exact = F)$estimate) %>%
  ungroup() %>%
  mutate(x=Inf, y=Inf, label=sprintf("rho = %.2g",cc)) %>%
  rename(f=virus)

p = plot_size_bubbles(df.vir$d0.log2, df.vir$fc_res_max, factor(df.vir$fc_res_max_d30), df.vir$virus) +
  scale_size_area(max_size = 8, name="count") +
  scale_color_manual(values=c("green","grey","red"), name="Response", labels=c("low","middle","high")) +
  geom_text(data=df.cc, aes(x,y,label=label), color="black", hjust=1.1, vjust=1.1) +
  xlab("day 0, log2") + ylab("adjMFC") +
  theme_bw() + theme(legend.key=element_blank())
p
ggsave(sprintf("%s_%s_viruses_check.png", fn.out, ag), w=12,h=5)
