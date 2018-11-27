source("functions.R")
source(file.path(dirr,"sra.R"))
t_cfg %>% select(sid, study, author) %>% print(n=40)

get_raw_read_list <- function(ti, sid) {
#{{{
if(sid == 'me10a') {
#{{{ Li2010
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = 'leaf',
                  Genotype = 'B73',
                  Treatment = SampleName,
                  Replicate = '',
                  paired = paired) %>%
        arrange(SampleID)
#}}}
}
#}}}
}

sid = 'sp056'
fi = file.path(dird, '03.raw.xlsx')
ti = read_xlsx(fi, sheet=sid)
#
diro1 = sprintf("%s/%s/09_fastq_raw", dirc, sid)
diro2 = sprintf("%s/%s/10_fastq", dirc, sid)
if(!dir.exists(diro1)) system(sprintf("mkdir -p %s", diro1))
if(!dir.exists(diro2)) system(sprintf("mkdir -p %s", diro2))

# create sym-links, write read list
ndig = floor(log10(nrow(ti))) + 1
fmt_sid = sprintf("s%%0%dd", ndig)
dir0 = '/home/springer/data_release/umgc/hiseq'
ti = ti %>% fill(Tissue, Genotype, directory) %>%
    mutate(SampleID = sprintf(fmt_sid, 1:nrow(ti))) %>%
    mutate(f1 = sprintf("%s/%s/%s_R1_001.fastq", dir0, directory, file),
           f2 = sprintf("%s/%s/%s_R2_001.fastq", dir0, directory, file)) %>%
    mutate(gz = ifelse(file.exists(f1), F, T)) %>%
    mutate(f1 = ifelse(gz, sprintf("%s.gz", f1), f1),
           f2 = ifelse(gz, sprintf("%s.gz", f2), f2)) %>%
    mutate(nf1 = sprintf("%s/%s_1.fq", diro1, SampleID),
           nf2 = sprintf("%s/%s_2.fq", diro1, SampleID)) %>%
    mutate(nf1 = ifelse(gz, sprintf("%s.gz", nf1), nf1),
           nf2 = ifelse(gz, sprintf("%s.gz", nf2), nf2)) %>%
    mutate(cmd1 = sprintf("ln -sf %s %s", f1, nf1),
           cmd2 = sprintf("ln -sf %s %s", f2, nf2)) %>%
    mutate(tag = file.exists(f1) & file.exists(f2))
sum(!ti$tag)

map_int(ti$cmd1, system)
map_int(ti$cmd2, system)

th
th = ti %>% select(SampleID, Tissue, Genotype) %>%
    mutate(Treatment='', Replicate = '', paired = T) %>%
    select(SampleID, Tissue, Genotype, Treatment, everything())
th = sra_fill_replicate(th)
th %>% count(Genotype, Tissue, Treatment) %>% print(n=50)
fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)



