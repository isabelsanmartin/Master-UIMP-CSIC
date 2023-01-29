#PLOTTING STOCHASTIC CHARACTER MAPPING FROM REVBAYES OUTPUT
library(ggplot2)
library(ggsn)
library(ggtree)
library(RevGadgets)
library(ape)
library(phytools)


character_file = "output-stochastic/simple_marginal_character.tree"

sim = read.simmap(file=character_file, format="phylip")


# Define colours for each area or state
colors = vector()
for (i in 1:length( sim$maps ) ) { 
    colors = c(colors, names(sim$maps[[i]]) )
}
colors
colors = sort(as.numeric(unique(colors)))
colors

# We can use different two colour choices to plot the posterior tree as a "heatmap"
# cols = setNames( heat.colors(length(colors), rev=TRUE), colors)
# cols

# Or use a basic palette with red, yellow, green, blue, pink, etc. For marginal states, this works better.

cols = setNames( rainbow(length(colors), start=0.0, end=0.9), colors)
cols

# fsize is font size for tipe labels, lwd = line width for plotting, ftype = b (bold), i (italics)
# pts: whether to plot filled circles at each tree vertex, as well as transition points between mapped states: default is false

plotSimmap(sim, cols, fsize=1.0, lwd=2.0, split.vertical=TRUE, ftype="bi")

# Add legend


# To identify which colour corresponde to which state
leg = names(cols)
leg

add.simmap.legend(leg, colors=cols, cex=0.9, x=0.8, y=0.8, fsize=0.9)

# A message appears in console: "Click where you want to draw legend". Click and draw in RQuartz window to get the legend plotted.

# Save image using Save ----- RPlot

#PLOTTING POSTERIORS OF STOCHASTIC CHARACTER MAPPING BiSSE

posterior_file = "output-stochastic/simple_marginal_posterior.tree"

sim_p = read.simmap(file=posterior_file, format="phylip")


# Define colours for posterior probability 
colors = vector()
for (i in 1:length( sim_p$maps ) ) { 
    colors = c(colors, names(sim_p$maps[[i]]) )
}
colors = sort(as.numeric(unique(colors)))

# We can use different two colour choices to plot the posterior tree as a "heatmap". For posteriors, this works better.

cols = setNames( heat.colors(length(colors), rev=TRUE), colors)

# Or using a basic palette with red, yellow, blue, etc.

# cols = setNames( rainbow(length(colors), start=0.0, end=0.9, rev=TRUE), colors)



# fsize is font size for tipe labels, lwd = line width for plotting, ftype = b (bold), i (italics)
# pts: whether to plot filled circles at each tree vertex, as well as transition points between mapped states: default is false.

plotSimmap(sim_p, cols, fsize=1.0, lwd=2.0, split.vertical=TRUE, ftype="bi", pts=FALSE)

# Add legend
# To identify which colour corresponde to which value of the posterior probability
leg = names(cols)
leg

add.simmap.legend(leg, colors=cols, cex=0.2, x=0.2, y=0.2, fsize=0.3)

# A message appears in console: "Click where you want to draw legend". Click and draw in RQuartz window to get the legend plotted.

# Save image using Save ----- RPlot



