# library
library(rgl)

# This is to output a rgl plot in a rmarkdown document.
setupKnitr()

data = read.csv("C:/GitHub/OMP/data/cluster_output_8_threads_dataset.txt", 
  header = FALSE, sep = " ")

# Add a new column with color
mycolors <- c('royalblue1', 'darkcyan', 'oldlace', 'purple')
data$color <- mycolors[ as.numeric(data$V4) ]

# Plot
plot3d( 
  x=data$`V1`, y=data$`V2`, z=data$`V3`, 
  col = data$color, 
  type = 's', 
  radius = .1,
  xlab="X", ylab="Y", zlab="Z")

# To display in an R Markdown document:
rglwidget()

# To save to a file:
htmlwidgets::saveWidget(rglwidget(width = 520, height = 520), 
  file = "./3dscatter.html",
  libdir = "libs",
  selfcontained = FALSE
)

