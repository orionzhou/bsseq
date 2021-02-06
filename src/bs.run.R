source("functions.R")
require(DMRcaller)
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
#}}}

fi = '~/projects/bsseq/data/cache/bs13a/38_grange/s17.rds'
x = readRDS(fi)
x


