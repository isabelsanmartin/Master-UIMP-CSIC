source("scripts/plot_anc_range.util.R")

# file names
fp = "./" # edit to provide an absolute filepath
plot_fn = paste(fp, "output/epoch.range.pdf",sep="")
tree_fn = paste(fp, "output/epoch.ase.tre", sep="")
label_fn = paste(fp, "output/epoch.state_labels.txt", sep="")
color_fn = paste(fp, "range_colors.n4.txt", sep="")

#Attention: uncomment these commands if you haven't installed RevGadgets, devtools and dependencies. Alternatively, skip (comment out) these commands if you have already installed RevGadgets, devtools and dependencies

install.packages("devtools", dependencies=TRUE)
library(devtools)
install_github("GuangchuangYu/ggtree")
install_github("revbayes/RevGadgets")
library(RevGadgets)

# get state labels and state colors
labs <- c("1"  = "K",   "2"  = "O", 
          "3"  = "M",   "4"  = "H", 
          "5"  = "KO",  "6"  = "KM", 
          "7"  = "OM",  "8"  = "KH", 
          "9"  = "OH",  "10" = "MH", 
          "11" = "KOM", "12" = "KOH", 
          "13" = "KMH", "14" = "OMH", 
          "15" = "KOMH")
ancstates <- processAncStates(tree_fn, state_labels = labs)

# plot the ancestral states
pp=plotAncStatesPie(t = ancstates, 
                     # Include cladogenetic events
                     cladogenetic = T,
                     # adjust tip labels 
                     tip_labels_offset = 0.1)


# get plot dimensions
x_phy = max(pp$data$x)       # get height of tree
x_label = 3.5                # choose space for tip labels
x_start = 5.5                # choose starting age (greater than x_phy)
x0 = -(x_start - x_phy)      # determine starting pos for xlim
x1 = x_phy + x_label         # determine ending pos for xlim

# add axis
pp = pp + theme_tree2()
pp = pp + labs(x="Age (Ma)")

# change x coordinates
pp = pp + coord_cartesian(xlim=c(x0,x1), expand=TRUE)

# plot axis ticks
island_axis = sec_axis(~ ., breaks=x_phy-c(5.1, 2.95, 1.55, 0.5), labels=c("+K","+O","+M","+H") )
x_breaks = seq(0,x_start,1) + x0
x_labels = rev(seq(0,x_start,1))
pp = pp + scale_x_continuous(breaks=x_breaks, labels=x_labels, sec.axis=island_axis)

pp

# save 
ggsave(file=plot_fn, plot=pp, device="pdf", height=7, width=10, useDingbats=F)

