#{{{ load required libraries, define common variables
require(grid)
require(tidyverse)
require(readxl)
require(gtable)
require(ggtree)
require(RColorBrewer)
require(viridis)
require(cluster)
require(Hmisc)
require(ggsignif)
require(cowplot)
require(GGally)
require(ggridges)
require(ggpubr)
require(ggsci)
require(ggrepel)
require(scales)
require(pheatmap)
require(yaml)
options(stringsAsFactors = FALSE)
dirr = '~/git/luffy/r'
source(file.path(dirr, 'plot.R'))
source('~/projects/genomes/src/ge.fun.r')
#
dird = '~/projects/maize.methylation/data'
dirc = '/scratch.global/zhoux379/maize.methylation'
f_cfg = file.path(dird, '01.cfg.xlsx')
t_cfg = read_xlsx(f_cfg, sheet=1, col_names=T)
#f_yml = file.path(dird, '11.cfg.yaml')
#Sys.setenv("R_CONFIG_FILE" = f_yml)
#}}}

get_read_list <- function(dird, sid) {
    #{{{
    fh1 = sprintf("%s/05_read_list/%s.tsv", dird, sid)
    fh2 = sprintf("%s/05_read_list/%s.c.tsv", dird, sid)
    fh = ifelse(file.exists(fh2), fh2, fh1)
    read_tsv(fh)
    #}}}
}

