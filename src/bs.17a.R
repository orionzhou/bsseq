source("functions.R")
cols6p = brewer.pal(6, "Paired")[2:6]
cols4 = viridis_pal(option='A', end=.7)(4)
cols2 = pal_startrek()(2)

yid = 'bs17a'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))
#
ref = t_cfg %>% filter(yid == !!yid) %>% pull(ref)
th = bsseq_sample_meta(yid) %>%
    mutate(Tissue = ifelse(Tissue=='V2','shoot',Tissue))

mid = 'm1'
tis = th$Tissue[th$MergeID == mid]
fi = sprintf("%s/raw/%s/%s.bed.gz", dird, yid, mid)
ti = read_tsv(fi, col_names=c('chrom','start','end','srd','cid','bin','ns','nm','nu','methy'))
ti2 = ti %>% separate('cid',c('type1','gid','type2'), sep="[\\|]")
tim = c('up'=1,'exon'=2,'down'=3,'intron'=4)
ti3 = ti2 %>% distinct(type2, bin) %>% mutate(gbin = bin + (tim[type2]-1) * 10)

t_methy = ti2 %>%
    filter(ns >= 2, nm+nu >= 20) %>%
    group_by(type1, gid, type2, bin, methy) %>%
    summarise(nm=sum(nm), nu=sum(nu), ns=sum(ns)) %>%
    ungroup() %>% mutate(ratio=nm/(nm+nu)) %>%
    inner_join(ti3, by=c('type2','bin')) %>%
    select(methy, type1,gid,bin=gbin,ratio)
tm = t_methy %>% select(-type1)

type = 'syn'
tg = read_gene_groups('gene')
tg %>% count(type)
fo = sprintf('%s/10.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(tm, tg, fo, 9,8,cols6p)

type = 'tf'
tg = read_gene_groups("TF")
tg %>% count(type)
fo = sprintf('%s/10.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(tm, tg, fo, 9,8, cols2, leg.nrow=1)

fd = '~/projects/grn/data/17_degree/degree.rds'
td = readRDS(fd)

type = 'grn'
tg = td$reg %>% filter(nid=='nc01', score>=3, deg > 0) %>%
    mutate(type = as.integer(cut_number(deg, 4))) %>%
    select(type, gid=reg.gid)
tg = td$tgt %>% filter(nid=='nc01', score==2, deg > 0) %>%
    mutate(type = cut(deg, c(1,5,10,Inf), include.lowest=T, right=F)) %>%
    select(type, gid=tgt.gid)
fo = sprintf('%s/10.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(tm, tg, fo, 9,8, cols4, leg.nrow=1)


