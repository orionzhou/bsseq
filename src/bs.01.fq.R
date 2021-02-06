source("functions.R")
t_cfg %>% print(n=40)

diri = '~/projects/barn/data/15_read_list'
diro = '~/projects/bsseq/data/06_joinrep'

yid = 'bs15a'
fi = sprintf("%s/%s.tsv", diri, yid)
ti = read_tsv(fi)
ti %>% count(Tissue, Genotype, Treatment)
