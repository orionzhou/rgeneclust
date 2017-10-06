require(ggplot2)
require(plyr)
require(VennDiagram)
require(RColorBrewer)
require(ape)

dirw = '/home/springer/zhoux379/git/rgeneclust/rosa'
fcfg = file.path(dirw, "../rosa.csv")
dirw = '/home/springer/zhoux379/data/misc1/apple'
fcfg = file.path(dirw, "apple.csv")
setwd(dirw)
options(stringsAsFactors = FALSE)

## read gene location
tcfg = read.table(fcfg, header = F, sep = ",", as.is = T)

tl = data.frame()
for (i in 1:nrow(tcfg)) {
	org = tcfg$V1[i]; fl = tcfg$V3[i]
	if(fl == '') next
	tls = read.table(fl, header = T, sep = "\t", as.is = T)
	tl = rbind(tl, cbind(org = org, tls))
}
fl = file.path(dirw, "08.gene.loc.tsv")
write.table(tl, fl, sep = "\t", row.names = F, col.names = T, quote = F)

### domain configuration & cluster assignments
f08 = file.path(dirw, "08.gene.loc.tsv")
t08 = read.table(f08, header = T, sep = "\t", as.is = T)
f21 = file.path(dirw, "21.tsv")
t21 = read.table(f21, header = T, sep = "\t", as.is = T)
fcoil = file.path(dirw, "ncoil/31.txt")
tcoil = read.table(fcoil, header = F, sep = "\t", as.is = T)
f27 = file.path(dirw, "27.pfam.tsv")
t27 = read.table(f27, header = T, sep = "\t", as.is = T)
f32 = file.path(dirw, "32.tbl")
t32 = read.table(f32, header = T, sep = "\t", as.is = T)

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
t22 = merge(t21[,c(1:3,5)], tcoil2, by.x = 'id', by.y = 'V1', all.x = T)
t22$conf = doms
colnames(t22)[4] = "conf_NB_ARC"
idxs = which(!is.na(t22$cc))
t22$conf[idxs] = paste(t22$cc[idxs], t22$conf[idxs], sep = '-')
table(t22$conf)

ids = sub("[|]", "\01", t22$id)
orgs = sapply(strsplit(ids, split = "\01"), '[', 1)
gids = sapply(strsplit(ids, split = "\01"), '[', 2)
t23 = cbind(t22, org = orgs, gid = gids)
table(t23[,c('conf','org')])
stopifnot(nrow(t23) == nrow(21))
t23 = merge(t23, t32, by='id', all.x=T)

doms_simple = c("NB-ARC","TIR","LRR_2","LRR_3","LRR_4","LRR_6","LRR_8")
t28 = ddply(t27, .(qid), summarise, doms_extra = paste(unique(hid[!hid %in% doms_simple]), collapse=","), configuration = paste(sprintf("[%d-%d]%s:%g",qbeg,qend,hid,e), sep=" ", collapse=" "))

t08 = cbind(t08[,3:5], id = sprintf("%s|%s", t08$org, t08$id))
t24 = merge(t23, t08, by = 'id', all.x = T)
table(data.frame(org=t24$org, na=is.na(t24$chr)))

t51 = merge(t24, t28, by.x = 'id', by.y = 'qid', all.x = T)
t51 = t51[,c('id','org','gid','chr','beg','end','size', 'grp','conf','conf_NB_ARC','doms_extra','configuration')]
colnames(t51)[c(8:9,11:12)] = c("cluster","configuration_simple", "extra_domains", "configuration_detail")
fo = file.path(dirw, "51.dom.stat.tsv")
write.table(t51, fo, sep = "\t", row.names = F, col.names = T, quote = F, na = '')

### get NB seqs to Leon
f13 = file.path(dirw, "13.tsv")
t13 = read.table(f13, header = T, sep = "\t", as.is = T)

ids = sub("[|]", "\01", t13$id)
orgs = sapply(strsplit(ids, split = "\01"), '[', 1)
gids = sapply(strsplit(ids, split = "\01"), '[', 2)
idxs = which(orgs %in% c("Fa", "Fi", "Fni", "Fnu", "Fo", "Fv", "Ro"))
idxs = which(orgs %in% c("Md", "Mda", "Mdg", "Mdj"))

to = t13[idxs, c('id','nbs_qb','nbs_qe')]
to$nbs_qb = to$nbs_qb-1
fo = file.path(dirw, "13.bed")
write.table(to, fo, sep = "\t", row.names = F, col.names = F, quote = F)
#seqret.pl -d 05.pro.fas -b 13.bed -o 13.fas

### venn diagram for RosaR3
tt = table(t41$org)
ary = c(tt['Fv'], tt['Md'], tt['Pp'], tt['Fv,Md'], tt['Md,Pp'], tt['Fv,Pp'], tt['Fv,Md,Pp'])
ary[is.na(ary)] = 0
venn.plot <- draw.triple.venn(
  area1 = sum(ary[c(1,4,6,7)]), 
  area2 = sum(ary[c(2,4,5,7)]), 
  area3 = sum(ary[c(3,5,6,7)]), 
  n12 = ary[4]+ary[7], n23 = ary[5]+ary[7], n13 = ary[6]+ary[7], n123 = ary[7],
  category = c("Fv (Strawberry)", "Md (Apple)", "Pp (Pear)"),
  fill = c("dodgerblue", "firebrick", "green"), lty = "solid",
  cex = 1, cat.cex = 1, cat.dist = 0.00,
#  cat.just = list(c(-1, -1), c(1, 1)),
)
fo = sprintf("%s/54.share.pdf", dirw)
pdf(file = fo, width = 5, height = 5, bg = 'transparent')
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


### phylogeny
f_tree = file.path(dirw, "44.ft.nwk")
  tree = read.tree(f_tree)
	tl = t41
	rootlab = "Arabidopsis"
	tree = root(tree, outgroup = which(tree$tip.label==rootlab))
#	tree = rotate(tree, 1054)
  
  grps = tree$tip.label
  tip.cols = rep('black', length(grps))
  tip.cols[grps %in% tl$grp[tl$norg == 1 & tl$org == 'Fv']] = '#E41A1C'
  tip.cols[grps %in% tl$grp[tl$norg == 1 & tl$org == 'Md']] = '#377EB8'
  tip.cols[grps %in% tl$grp[tl$norg == 1 & tl$org == 'Pp']] = '#4DAF4A'
  tip.cols[grps %in% tl$grp[tl$norg == 2]] = 'grey'
  tip.cols[grps %in% tl$grp[tl$norg == 3]] = 'black'
  

  scores = as.numeric(tree$node.label)
  node.bg = rep('white', tree$Nnode)
  node.bg[scores >= 0.9] = 'black'
  node.bg[scores >= 0.8 & scores < 0.9] = 'grey'
  
  fo = file.path(dirw, "45.pdf")
  pdf(fo, width = 15, height = 15)
  plot(tree, type = 'fan', show.node.label = F, show.tip.label = T,
    tip.color = tip.cols, label.offset = 0.08, 
    no.margin = T, cex = 0.6)
#  tiplabels(pch = 22, frame = 'none', adj = 0.55, bg = tip.cols)
  nodelabels(pch = 22, bg = node.bg, cex=0.5)
#  add.scale.bar(x = 0.02, y = tree$Nnode*0.9 , lcol = 'black')
  dev.off()

write.tree(tree, file.path(dirw, "x.txt"))
tree = read.tree(file.path(dirw, "x.txt"))
labs = tree$tip.label
write.table(rev(labs), file.path(dirw, "44.ft.txt"), sep = "\t", row.names = F, col.names = F, quote = F, na = '')

### assign physical clusters based on chromosomal locations
source("/home/springer/zhoux379/git/luffy/r/Location.R")
cluster_gap_map = c(
	'Fv' = 270267,
	'Md' = 624739,
	'Pb' = 533904,
	'Pm' = 335848,
	'Pp' = 376109
)

require(xlsx)
fo = sprintf("%s/loccluster/11.xlsx", dirw)
for (org in names(cluster_gap_map)) {
	cluster_gap = cluster_gap_map[org]
	fi = sprintf("%s/loccluster/01.%s.tsv", dirw, org)
	tt = read.table(fi, header = T, sep = "\t")
	
	chrs = sort(unique(tt$Chr))
	chrmap = 1:length(chrs)
	names(chrmap) = chrs

	tt = cbind(tt, pos = chrmap[tt$Chr] * 100000000 + (tt$Start + tt$End) / 2)
	tt = tt[order(tt$pos), ]
	
	x = locCluster(tt$pos, cluster_gap)
	to = cbind(tt[,1:5], x[,3:4])
	#fo = sprintf("%s/loccluster/05.%s.tsv", dirw, org)
	#write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
	tag = ifelse(org == names(cluster_gap_map)[1], F, T)
	write.xlsx(to, file = fo, sheetName = org, row.names = F, col.names = T, append = tag)
}


dirw = '/home/springer/zhoux379/data/in/apple.gene.pos'
fi = file.path(dirw, "01.tsv")
ti = read.table(fi, sep = "\t", as.is = T, header = T)
to = data.frame()
for (org in names(cluster_gap_map)) {
	cluster_gap = cluster_gap_map[org]
	cluster_gap = 200000
	tt = ti[ti$org == org,]
	chrs = sort(unique(tt$chr))
	chrmap = 1:length(chrs)
	names(chrmap) = chrs

	tt = cbind(tt, pos = chrmap[tt$chr] * 100000000 + (tt$beg + tt$end) / 2)
	tt = tt[order(tt$pos), ]
	
	x = locCluster(tt$pos, cluster_gap)
	tos = cbind(tt[,-ncol(tt)], x[,3:4])
	to = rbind(to, tos)
}
fo = sprintf("%s/05.200kb.tsv", dirw)
write.table(to, fo, sep = "\t", row.names = F, col.names = T, quote = F)
