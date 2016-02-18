require(ggplot2)
require(plyr)

dirw = '/home/youngn/zhoup/git/rgeneclust/test'
setwd(dirw)

### domain configuration
f21 = file.path(dirw, "21.tbl")
t21 = read.table(f21, header = T, sep = "\t", as.is = T)
fcoil = file.path(dirw, "ncoil/31.txt")
tcoil = read.table(fcoil, header = F, sep = "\t", as.is = T)

simp_nbs <- function(domstr) {
  doms = unlist(strsplit(domstr, split = ","))
  conf = c()
  if('TIR' %in% doms) {
    conf = c(conf, 'TIR')
  }
  if('NB-ARC' %in% doms) {
    conf = c(conf, 'NBS')
  }
  if(sum(c("LRR_2","LRR_3","LRR_4","LRR_6","LRR_8") %in% doms) > 0) {
    conf = c(conf, 'LRR')
  }
  paste(conf, sep = "-", collapse = "-")
}

doms = sapply(t21$doms, simp_nbs)

tcoil2 = cbind(tcoil, cc = 'CC')
t22 = merge(t21, tcoil2, by.x = 'id', by.y = 'V1', all.x = T)
t22$conf = doms
idxs = which(!is.na(t22$cc))
t22$conf[idxs] = paste(t22$cc[idxs], t22$conf[idxs], sep = '-')
table(t22$conf)

orgs = sapply(strsplit(t22$id, split = "[|]"), '[', 1)
gids = sapply(strsplit(t22$id, split = "[|]"), '[', 2)
t23 = cbind(t22, org = orgs, gid = gids)
table(t23[,c('conf','org')])

fo = file.path(dirw, "51.dom.stat.tbl")
write.table(t23, fo, sep = "\t", row.names = F, col.names = T, quote = F)

### sequence clustering
f32 = file.path(dirw, "32.tbl")
t32 = read.table(f32, header = T, sep = "\t", as.is = T)
orgs = sapply(strsplit(t32$id, split = "[|]"), '[', 1)
gids = sapply(strsplit(t32$id, split = "[|]"), '[', 2)
t33 = cbind(t32, org = orgs, gid = gids)

t41 = ddply(t33, .(grp), summarise, n = length(id), org = paste(sort(unique(org)), sep = ",", collapse = ","))

fo = file.path(dirw, "52.cluster.stat.tbl")
write.table(t33[order(t33$grp),], fo, sep = "\t", row.names = F, col.names = T, quote = F)
