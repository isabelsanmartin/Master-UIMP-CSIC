source("scripts/plot_anc_range.util.R")

# file names
fp = "./" # edit to provide an absolute filepath
plot_fn = paste(fp, "output_5/simple.range.pdf",sep="")
tree_fn = paste(fp, "output_5/simple.ase.tre", sep="")
label_fn = paste(fp, "output_5/simple.state_labels.txt", sep="")
color_fn = paste(fp, "range_colors.n4.txt", sep="")

install.packages("devtools", dependencies=TRUE)
library(devtools)
install_github("GuangchuangYu/ggtree")
install_github("revbayes/RevGadgets")
library(RevGadgets)
library(ggtree)

# get state labels and state colors
labs <- c("1"  = "Wn",   "2"  = "En", 
          "3"  = "Wp",   "4"  = "Ep", 
          "5"  = "WnEn",  "6"  = "WnWp", 
          "7"  = "EnWp",  "8"  = "WnEp", 
          "9"  = "EnEp",  "10" = "WpEp", 
          "11" = "WnEnWp", "12" = "WnEnEp", 
          "13" = "WnWpEp", "14" = "EnWpEp", 
          "15" = "WnEnWpEp")
ancstates <- processAncStates(tree_fn, state_labels = labs)

# plot the ancestral states
pp=plotAncStatesPie(t = ancstates, 
                    # Include cladogenetic events
                    cladogenetic = T,
                    # Add text labels to the tip pie symbols
                    tip_labels_states = FALSE,
                    # tamaño del nombre de las sp.
                    tip_labels_size = 4,
                    # Offset those text labels slightly
                    tip_labels_states_offset = .05,
                    # Offset the tip labels to make room for tip pies
                    tip_labels_offset = 1, 
                    # Move tip pies right slightly 
                    tip_pie_nudge_x = .07,
                    # Change the size of node and tip pies  
                    tip_pie_size = 1,
                    node_pie_size = 1.5,
                    state_transparency = 1)


# get plot dimensions
x_phy = max(pp$data$x)       # get height of tree
x_label = 3.5                # choose space for tip labels
x_start = 45                  # choose starting age (greater than x_phy)
x0 = -(x_start - x_phy)      # determine starting pos for xlim
x1 = x_phy + x_label         # determine ending pos for xlim

# add axis
pp = pp + theme_tree2()
pp = pp + labs(x="Age (Ma)")

# change x coordinates
pp = pp + coord_cartesian(xlim=c(x0,x1), expand=TRUE)

# plot axis ticks
island_axis = sec_axis(~ ., breaks=x_phy-c(35, 25, 3.5), labels=c("-A","+T","+E") )
x_breaks = seq(0,x_start,5) + x0
x_labels = rev(seq(0,x_start,5))
pp = pp + scale_x_continuous(breaks=x_breaks, labels=x_labels, sec.axis=island_axis)

pp

# save 
ggsave(file=plot_fn, plot=pp, device="pdf", height=7, width=10, useDingbats=F)

