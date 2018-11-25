source("me.fun.r")
source("sra.R")
t_cfg

sid = 'mec03'
sid = 'me13c'
#{{{ config
Sys.setenv(R_CONFIG_ACTIVE = sid)
dirw = file.path(dird, '11_qc', sid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))
#
cfg = t_cfg %>% filter(sid == !!sid)
stopifnot(nrow(cfg) == 1)
study = cfg %>% pull(study)
readtype = cfg %>% pull(readtype)
mapper = cfg %>% pull(mapper)
genome = cfg %>% pull(reference)
meta = cfg %>% pull(meta)
x = load(file.path(dirg, genome, '55.rda'))
#}}}

#{{{ [meta=F] mapping stats, raw read counts, normalize
th = get_read_list(dird, sid)
tiss = unique(th$Tissue); genos = unique(th$Genotype); treas = unique(th$Treatment)
reps = unique(th$Replicate)
paired = unique(th$paired)
if(length(paired) == 2) paired = 'both'

diri = file.path(dird, '08_raw_output', sid)
fi = file.path(diri, 'trimming.tsv')
tt1 = read_tsv(fi) %>% 
    separate(sid, c('SampleID','pe'), sep='[\\.]') %>%
    mutate(passed=passed_filter_reads,
    failed=low_quality_reads+too_many_N_reads+too_short_reads+too_long_reads) %>%
    mutate(passed = ifelse(pe=='pe', passed/2, passed)) %>%
    mutate(failed = ifelse(pe=='pe', failed/2, failed)) %>%
    mutate(total = passed+failed) %>%
    select(SampleID, total, passed, failed)
fi = file.path(diri, 'bamstats.tsv')
tt2 = read_tsv(fi) %>% select(SampleID=sid, everything())
fi = file.path(diri, 'multiqc_data/multiqc_featureCounts.txt')
tt3 = read_multiqc_featurecounts(fi)

tt = th %>% select(-paired) %>% 
    left_join(tt1, by = 'SampleID') %>%
    left_join(tt2, by = 'SampleID') %>%
    left_join(tt3, by = 'SampleID')
tt %>% mutate(nd = passed-pair-unpair) %>% pull(nd) %>% sum()
tt %>% mutate(nd = pair-pair_map-pair_orphan-pair_unmap) %>% pull(nd) %>% sum()
tt %>% mutate(nd = pair+unpair-Assigned-
              Unassigned_MultiMapping-
              Unassigned_NoFeatures-Unassigned_Ambiguity-Unassigned_Unmapped) %>% select(nd)#pull(nd) %>% sum()
#
tt %>% group_by(Tissue, Genotype, Treatment) %>%
    summarise(total = sum(total), Assigned = sum(Assigned)) %>%
    ungroup() %>% group_by(1) %>%
    summarise(total_median = median(total/1e6),
              total_mean = mean(total/1e6),
              assigned_median = median(Assigned/1e6),
              assigned_mean = mean(Assigned/1e6)) %>% print(n=1)

fo = file.path(dirw, '10.mapping.stat.tsv')
write_tsv(tt, fo)

fi = file.path(diri, 'featurecounts.tsv')
t_rc = read_tsv(fi)
tm = t_rc %>% gather(SampleID, ReadCount, -gid) %>%
    filter(SampleID %in% th$SampleID)
res = readcount_norm(tm, size.gene)
tl = res$tl; tm = res$tm
#}}}

#{{{ [meta=T] collect featurecounts data & normalize
sids = str_split(config::get("sids"), "[\\+]")[[1]] 
sids
#
th = tibble(); t_rc = tibble()
for (sid1 in sids) {
    diri = file.path(dird, '08_raw_output', sid1, 'multiqc_data')
    th1 = get_read_list(dird, sid1)
    if(sid1 == 'me12a') {
        th1 = th1 %>% filter(Treatment == 'WT')
    } else if(sid1 == 'me13b') {
        th1 = th1 %>% filter(!str_detect(Treatment, "ET"))
    } else if(sid1 == 'me18b') {
        th1 = th1 %>% filter(Genotype == 'B73')
    } else if(sid1 == 'me99b') {
        th1 = th1 %>% filter(Genotype == 'B73')
    }
    th1 = th1 %>% mutate(sid = sid1) %>% 
        select(sid, SampleID, Tissue, Genotype, Treatment, everything())
    th = rbind(th, th1)
    fi = file.path(diri, '../featurecounts.tsv')
    t_rc1 = read_tsv(fi) %>% select(one_of(c('gid', th1$SampleID)))
    stopifnot(ncol(t_rc1) == nrow(th1) + 1)
    if(nrow(t_rc) == 0) {
        t_rc = t_rc1
    } else {
        t_rc = t_rc %>% inner_join(t_rc1, by = 'gid')
    }
}
dim(th); dim(t_rc)
th %>% dplyr::count(sid)
th %>% dplyr::count(Tissue) %>% print(n=40)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

res = merge_reps(th, tm, sid)
th = res$th; tm = res$tm
th = th %>% mutate(Replicate = NA)
tiss = unique(th$Tissue); genos = unique(th$Genotype); treas = unique(th$Treatment)
reps = unique(th$Replicate)
#}}}

#{{{ [optional] fix/remove mis-labelled samples, save to mexx.c.tsv and 20.rc.norm.rda
#cls = cutree(ehc, h = .1)
#tcl = tibble(SampleID = names(cls), grp = as.integer(cls))
#th2 = th %>% inner_join(tcl, by = 'SampleID')
#th3 = th2 %>% group_by(Tissue, Genotype) %>% 
    #summarise(nrep = length(Replicate), ngrp = length(unique(grp))) %>%
    #ungroup() %>%
    #filter(nrep > 1, ngrp > 1)
#th2 %>% inner_join(th3, by = c("Tissue", "Genotype")) %>% print(n=40)

fh1 = sprintf("%s/05_read_list/%s.tsv", dird, sid)
th = read_tsv(fh1)
fh2 = sprintf("%s/05_read_list/%s.c.tsv", dird, sid)
if(sid == 'me13c') {
    th = th %>% 
        mutate(Genotype = ifelse(SampleID=='SRR767691','Oh43',Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='SRR651079','Oh7b',Genotype))
} else if(sid == 'me14c') {
    th = th %>% filter(! SampleID %in% c("SRR254169"))
} else if(sid == 'me14d') {
    th = th %>% filter(! SampleID %in% c("SRR1573518", 'SRR1573513'))
} else if(sid == 'me16b') {
    th = th %>% filter(!SampleID %in% c("SRR1620930","SRR1620929","SRR1620927",
                                        "SRR1620908","SRR1620913"))
} else if(sid == 'me17a') {
    th = th %>%
        mutate(Tissue = ifelse(SampleID == 'SRR445601', 'tassel', Tissue)) %>%
        mutate(Tissue = ifelse(SampleID == 'SRR445416', 'tassel', Tissue)) %>%
        mutate(Genotype = ifelse(SampleID == 'SRR426798', 'Mo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID == 'SRR426814', 'M37W', Genotype))
} else if(sid == 'me99b') {
    gts = c("B73", "Mo17", "B73xMo17")
    tissues = sort(unique(th$Tissue))
    th = th %>% 
        #filter(Genotype %in% gts) %>%
        filter(! SampleID %in% c('BR207', 'BR230', "BR235")) %>%
        mutate(Genotype = ifelse(SampleID=='BR003', 'Mo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR004', 'B73', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR006', 'B73xMo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR007', 'Mo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR029', 'B73', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR032', 'Mo17', Genotype))
    th = th %>% mutate(Replicate = '')
    th = sra_fill_replicate(th)
} else if(sid == 'me99c') {
    th = th %>% 
        mutate(Tissue = ifelse(SampleID == 'bm252', 'Leaf', Tissue))
}
write_tsv(th, fh2)#, na = '')

# re-normalize everything [if no sample is removed then no need to do this]
fi = file.path(diri, 'featurecounts.tsv')
t_rc = read_tsv(fi)
tm = t_rc %>% gather(SampleID, ReadCount, -gid) %>%
    filter(SampleID %in% th$SampleID)
res = readcount_norm(tm, size.gene)
tl = res$tl; tm = res$tm
#}}}

#{{{ save to 20.rc.norm.rda
fo = file.path(dirw, '20.rc.norm.rda')
save(th, tl, tm, file = fo)
#}}}

#{{{ read from 20.rc.norm.rda
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x
#}}}

#{{{ hclust
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)

cor_opt = "spearman"
cor_opt = "pearson"
hc_opt = "ward.D"
hc_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
edist <- as.dist(1-cor(e, method = cor_opt))
ehc <- hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
lnames = ehc$labels[ehc$order]
#
tp = th %>% mutate(taxa = SampleID, lab = SampleID) 
if(length(tiss)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Tissue), lab)
if(length(genos)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Genotype), lab)
if(length(treas)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Treatment), lab)
if(length(reps)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Replicate), lab)
t_hc = tp %>% select(taxa, everything())
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
plot_hclust_tree(tree, t_hc, fo, 
                 labsize = config::get("hc.labsize"), 
                 x.expand = config::get("hc.x.expand"),
                 x.off = config::get("hc.x.off"), 
                 wd = config::get("hc.wd"), ht = config::get("hc.ht"))
#}}}

#{{{ mec03-specific annotation
fh = file.path(dirw, '10.sample.tsv')
write_tsv(th, fh)
fh = file.path(dirw, '11.sample.curated.tsv')
th = read_tsv(fh)
th %>% dplyr::count(Tissue) %>% print(n = 23)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=10,
              pca = T, max_iter = 500)

#tissues_merge = th %>% count(Tissue) %>% filter(n < 5) %>% pull(Tissue)
tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID') %>%
    replace_na(list(Treatment='')) %>%
    #mutate(Treatment=ifelse(Tissue %in% tissues_merge, sprintf("%s_%s",Tissue,Treatment), Treatment)) %>%
    #mutate(Tissue=ifelse(Tissue %in% tissues_merge, 'misc', Tissue)) %>%
    mutate(txt = sprintf("%s_%s", Genotype, Replicate))
p_tsne = ggplot(tp) +
    geom_text_repel(data=tp, aes(x=V1,y=V2,label=txt), size=2, alpha=.8) +
    geom_point(aes(x=V1, y=V2, color=Tissue), size=2) +
    #stat_ellipse(aes(x=V1, y=V2, fill=Tissue), linetype=1, alpha=.4) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    #scale_shape_manual(values = c(15,4,16)) +
    scale_color_manual(values = pal_npg()(10)) +
    otheme(legend.pos = 'bottom.right', legend.dir = 'v',
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(color = guide_legend(ncol = 2, byrow = T))
fp = file.path(dirw, "25.tsne.pdf")
ggsave(p_tsne, filename = fp, width = 8, height = 8)
#}}}

me_output(sid, study, meta, dird)

#{{{ # ggtree + heatmap
is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips = tree$edge[is_tip,2]
lnames = tree$tip.label[ordered_tips]
df_mat = data.frame(as.matrix(1-edist))[,rev(lnames)]
#
labsize = config::get("hc.labsize")
x.expand = config::get("hc.x.expand")
x.off = config::get("hc.x.off")
cols1 = c('gray80','black','red','seagreen3', pal_d3()(5))
p1 = ggtree(tree) +
    #geom_tiplab(size = labsize, color = 'black') +
    scale_x_continuous(expand = expand_scale(mult=c(.02,x.expand))) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    geom_tiplab(aes(label = lab), size = labsize, offset = x.off, family = 'mono')
    #geom_text(aes(label = Rep), size = 2, nudge_x = .022, hjust = 0)
p2 = gheatmap(p1, df_mat, offset = 10, width = 10, colnames = F)
fo = file.path(dirw, '22.heatmap.pdf')
ggsave(p2, filename = fo, width = 30, height = 20)
#}}}
#{{{ # pheatmap
mat = 1 - as.matrix(edist)
mat = mat[lnames, lnames]
t_hm = t_hc %>% mutate(SampleID = factor(SampleID, levels = lnames)) %>%
    arrange(SampleID)
fo = sprintf("%s/22.heatmap.pdf", dirw)
pheatmap(
    mat               = mat,
    #color             = inferno(length(mat_breaks) - 1),
    #breaks            = mat_breaks,
    border_color      = NA,
    cluster_cols      = F,
    cluster_rows      = F,
    show_colnames     = F,
    show_rownames     = T,
    labels_row        = t_hm$lab,
    #annotation_row    = th[,c('SampleID','sid')],
    #annotation_colors = pal_d3()(10),
    drop_levels       = T,
    fontsize          = 6,
    main              = study,
    filename          = fo,
    width             = config::get("hm.wd"),
    height            = config::get("hm.ht")
)

  cellwidth = 30, cellheight = 30, scale = "none",
  treeheight_row = 200,
  kmeans_k = NA,
  show_rownames = T, show_colnames = F,
  clustering_method = "complete",
  cluster_rows = T, cluster_cols = T,
  #clustering_distance_rows = drows1, 
  #clustering_distance_cols = dcols1,
  #annotation_col = ta,
  #annotation_colors = ann_colors,
#}}}
#{{{ # PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(Treatment = factor(Treatment), Replicate = factor(Replicate))
fo = file.path(dirw, '25.pca.pdf')
plot_pca(tp, fo, opt = config::get("pca.opt"), labsize = config::get("pca.labsize"),
         wd = config::get("pca.wd"), ht = config::get("pca.ht"))
#}}}

#{{{ ##me99d - enders stress response 3' RNA-Seq
#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           lab = sprintf("%s %s %s", Tissue, Genotype, Treatment)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(-.02,3.3)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = lab), size = 2.5, nudge_x = .01, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 10)
#}}}
#}}}
#{{{ ##me99e - settles endosperm
#{{{ hclust tree
tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           lab = sprintf("%s %s %s", Tissue, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(-.02,5.5)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = lab), size = 2.5, nudge_x = .01, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 6, height = 8)
#}}}
#}}}


