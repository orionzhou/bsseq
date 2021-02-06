source("functions.R")
t_cfg
gcfg = read_genome_conf()
tsyn = read_syn(gcfg)
dirw = file.path(dird, 'test')
tg = gcfg$gene.loc %>% filter(etype=='exon') %>%
    group_by(gid,tid,ttype) %>%
    summarise(chrom=chrom[1],start=min(start)-1,end=max(end),srd=srd[1]) %>%
    mutate(len= end - start) %>%
    mutate(start0 = start-len/2, end0 = end+len/2) %>%
    group_by(1) %>% mutate(idx=1:n()) %>% ungroup() %>% select(-`1`)

#{{{ create gene bin file
gcfg$gene %>% select(ttype, len=size.rna) %>%
    group_by(ttype) %>% nest() %>%
    mutate(data2 = map(data, desc_stat)) %>% select(-data) %>%
    unnest() %>% spread(statK, statV)

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
write_tsv(to, fo, col_names=F)
#}}}

# run intersectBed job

tis = 'shoot'
tis = 'root'
fi = sprintf('%s/21.%s.ovlp.bed', dirw, tis)
ti = read_tsv(fi, col_names=c('chrom','start','end','gid','i','chrom2','start2','end2','cg_m','cg_t','cg_n','chg_m','chg_t','chg_n','chh_m','chh_t','chh_n','bp'))
tm = ti %>% mutate(pct = bp/(end2-start2)) %>% filter(pct >= .5)
#
tm1 = tm %>% select(gid,i,nm=cg_m,nt=cg_t,ns=cg_n) %>% mutate(methy = 'CG')
tm2 = tm %>% select(gid,i,nm=chg_m,nt=chg_t,ns=chg_n) %>% mutate(methy = 'CHG')
tm3 = tm %>% select(gid,i,nm=chh_m,nt=chh_t,ns=chh_n) %>% mutate(methy = 'CHH')
t_methy = rbind(tm1, tm2, tm3) %>%
    filter(ns >= 2, nt >= 50) %>%
    group_by(gid, i, methy) %>%
    summarise(nm=sum(nm), nt=sum(nt), ns=sum(ns)) %>%
    ungroup() %>% mutate(ratio=nm/nt) %>%
    select(methy, gid,i,ratio)

meta_plot_by_type <- function(t_methy, t_gene_group, fo, wd=8, ht=8) {
#{{{
    tp0 = t_methy %>% #mutate(i=i) %>%
        inner_join(t_gene_group, by = 'gid') %>%
        mutate(methy = str_to_upper(methy)) %>%
        filter(!is.na(ratio)) %>%
        mutate(type = factor(type))
    types = levels(tp0$type)
    tps = tp0 %>% mutate(type = as.character(type)) %>%
        filter(i == 11, methy == 'CG') %>%
        count(type) %>% mutate(type2=sprintf("%s (N=%s)", type, number(n))) %>%
        mutate(type=factor(type, levels=types)) %>% arrange(type)
    level_key = tps$type2
    names(level_key) = tps$type
    tp = tp0 %>% mutate(type = recode_factor(type, !!!level_key)) %>%
        group_by(type,i,methy) %>%
        summarise(q5=quantile(ratio, .05), q25=quantile(ratio,.25),
                  q50=quantile(ratio,.5), avg = mean(ratio),
                  q75=quantile(ratio,.75), q95=quantile(ratio,.95)) %>%
        ungroup()
    #
    tss=11; tts=30
    p1 = ggplot(tp) +
        geom_ribbon(aes(x=i, ymin=q25, ymax=q75), fill='grey90') +
        geom_line(aes(x=i, y=q50, linetype='median')) +
        geom_line(aes(x=i, y=avg, linetype='mean')) +
        geom_vline(xintercept=tss, color=pal_npg()(2)[1], lty=2) +
        geom_vline(xintercept=tts, color=pal_npg()(2)[2], lty=2) +
        scale_x_continuous(name='bin', breaks=c(tss,tts), labels=c('TSS','TTS'), expand=c(0, 0)) +
        scale_y_continuous(name = 'methylation ratio', expand=expand_scale(mult=c(.05,.05))) +
        scale_linetype_manual(values=c(2,1)) +
        facet_grid(methy~type, scale='free_y') +
        otheme(legend.pos = 'bottom.right', legend.dir = 'v', legend.title = F,
            xtick=T, ytick=T, ygrid = T,
            ytitle=T, xtext=T, ytext=T)
    ggsave(p1, file=fo, width=wd, height=ht)
#}}}
}

type = 'ttype'
t_gene_group = tg %>% select(gid, type=ttype)
wd = 8; ht = 8
fo = sprintf('%s/27.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(t_methy, t_gene_group, fo, wd, ht)
#
type = 'syn'
t_gene_group = tsyn %>% select(gid, type=ftype)
wd = 12; ht = 6
fo = sprintf('%s/27.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(t_methy, t_gene_group, fo, wd, ht)
#
type = 'tf'
ft = '~/projects/genome/data/B73/61_functional/06.tf.tsv'
tt = read_tsv(ft) %>% mutate(tf = 'TF')
t_gene_group = tsyn %>% mutate(ftype = as.character(ftype)) %>%
    mutate(ftype = ifelse(ftype == 'non-syntenic', ftype, 'syntenic')) %>%
    left_join(tt, by = 'gid') %>% replace_na(list(tf = 'non-TF')) %>%
    mutate(type = sprintf("%s %s", ftype, tf))
wd = 10; ht = 6
fo = sprintf('%s/27.%s.%s.pdf', dirw, tis, type)
meta_plot_by_type(t_methy, t_gene_group, fo, wd, ht)




