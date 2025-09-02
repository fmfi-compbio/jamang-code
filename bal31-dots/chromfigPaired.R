# work with 2 experiments (small, large) named based on effect size
# large is drawn in the background, bigger circles, darker colors

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop("Usage: genome.sizes densityS mergedS densityL mergedL output.n", call.=FALSE)
}

# color of points small, large
CPOINTS = "firebrick1"
CPOINTL = "firebrick4"
# color of density, small, large
CDENS = "gold2"
CDENL = "steelblue4"
# color of chrom
CCHR = "black"

N = 20
COLSIZE = 3e6
ROWS = 10
YCHROM1=0.0  # chromosome rectangle bottom
YCHROM2=0.01  # chromome rectangle top
YDENSE=YCHROM2+0.04
XDENSE=450
SDENSE=0.6   # scale - good for pcov without log
YMERGED=(YCHROM1+YCHROM2)/2
SMERGEDS=2  # circle size in mm for small
SMERGEDL=4  # circle size in mm for small
YCONTIG=-0.17
XCONTIG=3e4

chrom = read.table(args[1])
names(chrom)=c("contig","length")
chrom = subset(chrom, contig != "mtDNA")
chrom$x = ((0:(N-1)) %/% ROWS) * COLSIZE
chrom$y = ROWS - (0:(N-1)) %% ROWS

chrom2x = chrom$x
names(chrom2x) = chrom$contig
chrom2y = chrom$y
names(chrom2y) = chrom$contig

read_merged <- function(filename) {
  merged = try(read.table(filename))
  if (inherits(merged, "try-error"))
    merged = data.frame(matrix(ncol = 3, nrow = 0))
  names(merged) = c("contig", "start", "end")
  merged = subset(merged, contig != "mtDNA")
  merged$x = chrom2x[as.character(merged$contig)]
  merged$y = chrom2y[as.character(merged$contig)]
  return(merged)
}

mergedS = read_merged(args[3])
mergedL = read_merged(args[5])

read_density <- function(filename) {
  density = read.table(filename)
  names(density) = c("contig", "start", "end", "num", "cov", "len", "pcov")
  density = subset(density, contig != "mtDNA" & num > 0)
  density$log = log(density$pcov)
  density$val = density$log - min(density$log)
  density$x = chrom2x[as.character(density$contig)]
  density$y = chrom2y[as.character(density$contig)]
  return(density)
}

densityS = read_density(args[2])
densityL = read_density(args[4])

library(ggplot2)
ggplot() + geom_rect(data=chrom, mapping=aes(xmin=x, xmax=x+length, ymin=y+YCHROM1, ymax=y+YCHROM2), fill="blue", color=CCHR)
last_plot() + geom_rect(data=densityL, mapping=aes(xmin=x+XDENSE+start, xmax=x-XDENSE+end, ymin=y+YDENSE, ymax=y+YDENSE+pcov*SDENSE), fill=CDENL, color=CDENL)
last_plot() + geom_rect(data=densityS, mapping=aes(xmin=x+XDENSE+start, xmax=x-XDENSE+end, ymin=y+YDENSE, ymax=y+YDENSE+pcov*SDENSE), fill=CDENS, color=CDENS, linewidth=0)
last_plot() + geom_text(data=chrom, aes(x=x+XCONTIG,y=y+YCONTIG, label=contig), size=4, color="black", hjust=0)
if(nrow(mergedL)>0) 
  last_plot() + geom_point(data=mergedL, mapping=aes(x=x+(start+end)/2, y=y+YMERGED), size=SMERGEDL, stroke=0, fill=CPOINTL, color=CPOINTL, shape = "circle filled")
if(nrow(mergedS)>0) 
  last_plot() + geom_point(data=mergedS, mapping=aes(x=x+(start+end)/2, y=y+YMERGED), size=SMERGEDS, stroke=0, fill=CPOINTS, color=CPOINTS, shape = "circle filled")

LEGENDX = 0.4*COLSIZE
LEGENDY = 3
LEGENDYSKIP = 0.5
LEGENDW = XCONTIG/2
LEGENDH = 0.05
df = data.frame(x = LEGENDX, y=LEGENDY)
last_plot() + geom_point(data=df, mapping=aes(x=x, y=y), size=SMERGEDL, stroke=0, fill=CPOINTL, color=CPOINTL, shape = "circle filled")
last_plot() + geom_text(data=df, aes(x=x+XCONTIG,y=y, label=" detected, large dose/time"), size=4, color="black", hjust=0)
df = data.frame(x = LEGENDX, y=LEGENDY-LEGENDYSKIP)
last_plot() + geom_point(data=df, mapping=aes(x=x, y=y), size=SMERGEDS, stroke=0, fill=CPOINTS, color=CPOINTS, shape = "circle filled")
last_plot() + geom_text(data=df, aes(x=x+XCONTIG,y=y, label=" detected, small dose/time"), size=4, color="black", hjust=0)
df = data.frame(x = LEGENDX, y=LEGENDY-2*LEGENDYSKIP)
last_plot() + geom_rect(data=df, mapping=aes(xmin=x-LEGENDW, xmax=x+LEGENDW, ymin=y-LEGENDH, ymax=y+LEGENDH), fill=CDENL, color=CDENL) 
last_plot() + geom_text(data=df, aes(x=x+XCONTIG,y=y, label=" % significant, large dose/time"), size=4, color="black", hjust=0)
df = data.frame(x = LEGENDX, y=LEGENDY-3*LEGENDYSKIP)
last_plot() + geom_rect(data=df, mapping=aes(xmin=x-LEGENDW, xmax=x+LEGENDW, ymin=y-LEGENDH, ymax=y+LEGENDH), fill=CDENS, color=CDENS, linewidth=0) 
last_plot() + geom_text(data=df, aes(x=x+XCONTIG,y=y, label=" % significant, small dose/time"), size=4, color="black", hjust=0)

# smaller ratio means lower graph
last_plot() + theme_bw()+ coord_fixed(ratio = 0.3e6)

last_plot() + theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), 
  axis.ticks=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank(),
  panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.position="none",
)

ggsave(args[6])

