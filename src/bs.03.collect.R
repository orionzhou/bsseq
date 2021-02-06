source("functions.R")
#require(DMRcaller)
t_cfg
gcfg = read_genome_conf()
tsyn = read_syn(gcfg)
dirw = file.path(dird, 'test')

#{{{ create gene bins
gcfg$gene %>% select(ttype, len=size.rna) %>%
    group_by(ttype) %>% nest() %>%
    mutate(data2 = map(data, desc_stat)) %>% select(-data) %>%
    unnest() %>% spread(statK, statV)
#
tg = gcfg$gene.loc %>% filter(etype=='exon') %>%
    group_by(gid,tid,ttype) %>%
    summarise(chrom=chrom[1],start=min(start)-1,end=max(end),srd=srd[1]) %>%
    mutate(len= end - start) %>%
    mutate(start0 = start-len/2, end0 = end+len/2) %>%
    group_by(1) %>% mutate(idx=1:n()) %>% ungroup() %>% select(-`1`)
#
gr = with(tg, GRanges(seqnames=chrom, ranges=IRanges(start0,end0)))
nbin = 40
grl = tile(gr, n=nbin)

to = as_tibble(grl) %>% select(idx=group,chrom=seqnames,start,end) %>%
    inner_join(tg[c('gid','idx','srd')], by='idx') %>%
    select(-idx)
to1 = to %>% filter(srd == '+') %>%
    arrange(gid, chrom, start) %>%
    group_by(gid) %>%
    mutate(i = 1:n()) %>% ungroup()
to2 = to %>% filter(srd == '-') %>%
    arrange(gid, chrom, desc(start)) %>%
    group_by(gid) %>%
    mutate(i = 1:n()) %>% ungroup()
to = rbind(to1, to2)%>%
    mutate(start = start-1) %>%
    filter(start>0, end>0) %>%
    select(chrom, start, end, gid, i) %>%
    arrange(chrom, start, end)

fo = file.path(dirw, '10.gene.bed')
#write_tsv(to, fo, col_names=F)
#}}}

fi = '~/projects/bsseq/data/cache/bs13a/38_grange/s17.rds'
x = readRDS(fi)
x

go = as_tibble(x[1:10])


sum(x$readsN > 0)
fo = 'tmp.bed'

meta_plot_by_type <- function(t_methy, t_gene_group, fo, wd=8, ht=8) {
#{{{ mRNA / lncRNA
    tp0 = tp0 %>% mutate(type = factor(type))
    types = levels(tp0$type)
    tps = tp0 %>% mutate(type = as.character(type)) %>%
        count(type) %>% mutate(type2=sprintf("%s (N=%s)", type, number(n))) %>%
        mutate(type=factor(type, levels=types)) %>% arrange(type)
    level_key = tps$type2
    names(level_key) = tps$type
    tp0 = tp0 %>% mutate(type = recode_factor(type, !!!level_key))
    tp = t_methy %>% mutate(i=i) %>%
        inner_join(tp0, by = 'gid') %>%
        select(type, i, cg, chg, chh) %>%
        gather(methy, ratio, -type, -i) %>%
        mutate(methy = str_to_upper(methy)) %>%
        filter(!is.na(ratio)) %>%
        group_by(type,i,methy) %>%
        summarise(q5=quantile(ratio, .05), q25=quantile(ratio,.25),
            q50=quantile(ratio,.5), q75=quantile(ratio,.75), q95=quantile(ratio,.95)) %>%
        ungroup()
    #
    tss=11; tts=30
    p1 = ggplot(tp) +
        geom_ribbon(aes(x=i, ymin=q25, ymax=q75), fill='grey90') +
        geom_line(aes(x=i, y=q50)) +
        geom_vline(xintercept=tss, color=pal_npg()(2)[1], lty=2) +
        geom_vline(xintercept=tts, color=pal_npg()(2)[2], lty=2) +
        scale_x_continuous(name='bin', breaks=c(tss,tts), labels=c('TSS','TTS'), expand=c(0, 0)) +
        scale_y_continuous(name = 'methylation ratio', expand=expand_scale(mult=c(.05,.05))) +
        facet_grid(methy~type, scale='free_y') +
        otheme(legend.pos = 'top.left', legend.dir = 'v', legend.title = T,
            xtick=T, ytick=T, ygrid = T,
            ytitle=T, xtext=T, ytext=T)
    ggsave(p1, file=fo, width=wd, height=ht)
#}}}
}

type = 'ttype'
tp0 = tg %>% select(gid, type=ttype)
wd = 8; ht = 8
fo = sprintf('%s/27.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(t_methy, tp0, fo, wd, ht)
#
type = 'syn'
tp0 = tsyn %>% select(gid, type=ftype)
wd = 12; ht = 6
fo = sprintf('%s/27.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(t_methy, tp0, fo, wd, ht)
#
type = 'tf'
ft = '~/projects/genome/data/B73/61_functional/06.tf.tsv'
tt = read_tsv(ft) %>% mutate(tf = 'TF')
tp0 = tsyn %>% mutate(ftype = as.character(ftype)) %>%
    mutate(ftype = ifelse(ftype == 'non-syntenic', ftype, 'syntenic')) %>%
    left_join(tt, by = 'gid') %>% replace_na(list(tf = 'non-TF')) %>%
    mutate(type = sprintf("%s %s", ftype, tf))
wd = 10; ht = 6
fo = sprintf('%s/27.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(t_methy, tp0, fo, wd, ht)




