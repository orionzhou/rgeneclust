require(ggplot2)
require(plyr)
require(VennDiagram)

dirw = '/home/youngn/zhoup/git/rgeneclust/rosa'
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

### venn diagram for RosaR3
tt = table(t41$org)
a1 = tt['Fv']; a2 = tt['Md']; a3 = tt['Pp']; a12 = tt['Fv,Md']; a23 = tt['Md,Pp']; a13 = tt['Fv,Pp']; a123 = tt['Fv,Md,Pp']
venn.plot <- draw.triple.venn(area1 = a1, area2 = a2, area3 = a3, 
  n12 = a12+a123, n23 = a23+a123, n13 = a13+a123, n123 = a123,
  category = c("Fv (Grape)", "Md (Apple)", "Pp (Pear)"),
  fill = c("dodgerblue", "firebrick", "green"), lty = "solid",
  cex = 1, cat.cex = 1, cat.dist = 0.04,
#  cat.just = list(c(-1, -1), c(1, 1)),
)
fo = sprintf("%s/54.share.pdf", dirw)
pdf(file = fo, width = 4, height = 4, bg = 'transparent')
grid.newpage()
grid.draw(venn.plot)
dev.off()

### AFS for RosaR12
t42 = ddply(t33, .(grp), summarise, no = length(sort(unique(org))))
tt = table(t42$no)
to = data.frame(norg = as.numeric(names(tt)), cnt = as.numeric(tt))
p = ggplot(to) +
  geom_bar(aes(x = norg, y = cnt), stat = 'identity', width = 0.6) +
  scale_y_continuous(name = "# NBS lineages") +
  scale_x_continuous(name = "# sharing species", breaks = 1:6)
ggsave(p, filename = file.path(dirw, "54.share.pdf"), width = 4, height = 4)
