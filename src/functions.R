require(devtools)
load_all('~/git/rmaize')
require(GenomicRanges)
dird = '~/projects/bsseq/data'
dirc = '/scratch.global/zhoux379/bsseq'
f_cfg = '~/projects/master.xlsx'
t_cfg = read_xlsx(f_cfg, sheet='barn', col_names=T) %>%
    filter(libtype=='bsseq') %>%
    select(yid,author,study,genotype,tissue,n,ref)
gcfg = read_genome_conf()
tsyn = read_syn(gcfg)
#f_yml = file.path(dird, '11.cfg.yaml')
#Sys.setenv("R_CONFIG_FILE" = f_yml)

read_te <- function(genome='Zmays_B73', dirg='~/projects/genome/data') {
    #{{{
    f_te = file.path(dirg, genome, '50_annotation/30.TE.tsv')
    read_tsv(f_te)
    #}}}
}
read_tf <- function(ft='~/projects/genome/data/Zmays_B73/61_functional/06.tf.tsv')
    tt = read_tsv(ft)
te = read_te()
tf = read_tf()

read_gene_groups <- function(opt='gene+TE', gcfg.=gcfg, tsyn.=tsyn,
                             te.=te, tf.=tf) {
    #{{{
    if(opt == 'rna') {
        gcfg$gene %>% select(gid, type=ttype) %>%
            filter(type %in% c('mRNA','lnc_RNA'))
    } else if(opt == 'gene') {
        tsyn %>% select(gid, type=ftype)
    } else if(opt == 'gene+TE') {
        t_te = te %>% select(gid=id) %>% mutate(type='TE')
        tsyn %>% mutate(ftype=as.character(ftype)) %>%
            select(gid, type=ftype) %>%
            bind_rows(t_te) %>%
            mutate(type = fct_relevel(type, 'TE'))
    } else if(opt == 'TE') {
        te %>% select(gid=id) %>% mutate(type)
    } else if(opt == 'TF') {
        tt = tf %>% mutate(type2 = 'TF')
        tsyn %>% mutate(ftype = as.character(ftype)) %>%
            mutate(ftype = ifelse(ftype=='non-syntenic', ftype, 'syntenic')) %>%
            left_join(tt, by = 'gid') %>%
            replace_na(list(type2 = 'non-TF')) %>%
            filter(ftype == 'syntenic') %>%
            mutate(type = sprintf("%s %s", ftype, type2)) %>%
            select(gid, type)
    } else {
        stop(sprintf("unknown opt %s\n", opt))
    }
    #}}}
}

bsseq_sample_meta <- function(yid='bs17a',dird='~/projects/bsseq/data') {
    #{{{
    fh1 = sprintf("%s/05_meta/%s.tsv", dird, yid)
    read_tsv(fh1)
    #}}}
}

meta_plot_by_type <- function(t_methy, tg, fo, wd=8, ht=8, cols=cols6p, leg.nrow=2) {
#{{{
    tp0 = t_methy %>% #mutate(i=i) %>%
        inner_join(tg, by = 'gid') %>%
        mutate(methy = str_to_upper(methy)) %>%
        filter(!is.na(ratio)) %>%
        mutate(type = factor(type))
    types = levels(tp0$type)
    tps = tp0 %>% mutate(type = as.character(type)) %>%
        filter(bin == 11, methy == 'CG') %>%
        count(type) %>% mutate(type2=sprintf("%s (N=%s)", type, number(n))) %>%
        mutate(type=factor(type, levels=types)) %>% arrange(type)
    level_key = tps$type2
    names(level_key) = tps$type
    tp = tp0 %>% mutate(type = recode_factor(type, !!!level_key)) %>%
        group_by(type,bin,methy) %>%
        summarise(q5=quantile(ratio, .05), q25=quantile(ratio,.25),
                  q50=quantile(ratio,.5), avg = mean(ratio),
                  q75=quantile(ratio,.75), q95=quantile(ratio,.95)) %>%
        ungroup() %>% select(type,bin,methy,Mean=avg,Q25=q25,Q50=q50,Q75=q75) %>%
        gather(stat, v, -type, -bin, -methy)
    #
    tss=11; tts=30
    leg.vjust = ifelse(leg.nrow==2, -.1, -.4)
    p1 = ggplot(tp) +
        #geom_ribbon(aes(x=bin, ymin=q25, ymax=q75), fill='grey90') +
        #geom_line(aes(x=bin, y=v, linetype='median')) +
        #geom_line(aes(x=bin, y=avg, linetype='mean')) +
        #geom_vline(xintercept=tss, color=pal_npg()(2)[1], lty=2) +
        #geom_vline(xintercept=tts, color=pal_npg()(2)[2], lty=2) +
        geom_line(aes(x=bin, y=v, color=type)) +
        geom_vline(xintercept=tss, lty=2) +
        geom_vline(xintercept=tts, lty=2) +
        scale_x_continuous(name='bin', breaks=c(tss,tts), labels=c('TSS','TTS'), expand=c(0, 0)) +
        scale_y_continuous(name = 'methylation ratio', expand=expand_scale(mult=c(.05,.05))) +
        scale_linetype_manual(values=c(2,1)) +
        scale_color_manual(values=cols) +
        facet_grid(methy~stat, scale='free_y') +
        otheme(legend.pos='top.center.out', legend.dir='v', legend.title = F,
            xtick=T, ytick=T, ygrid = T, margin=c(1.5,.2,.2,.2),
            ytitle=T, xtext=T, ytext=T) +
        theme(legend.justification = c(.5,leg.vjust)) +
        guides(color = guide_legend(nrow = leg.nrow))
    ggsave(p1, file=fo, width=wd, height=ht)
#}}}
}


