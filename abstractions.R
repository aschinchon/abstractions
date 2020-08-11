######################################################################
# Load in libraries
######################################################################
library(Rcpp) # to iterate fast
library(tidyverse) # to plot and transform data
library(colourlovers) # to color drawings with nice colors
library(reshape2) # to convert matrix into data frames

# Import C++ code
sourceCpp('abstractions_funs.cpp')

# Default aesthetics of the ggplot
opt <-  theme(panel.border = element_rect(color="black", fill = NA),
              legend.position = "none",
              axis.ticks       = element_blank(),
              panel.grid       = element_blank(),
              axis.title       = element_blank(),
              axis.text        = element_blank())

######################################################################
# Parameters
######################################################################
# Framework environment
imageW <- 800 # image width (pixels)
imageH <- 800 # image heigth (pixels)
decayT <- 0.1 # Trail-map chemoattractant diffusion decay factor

# Agent
FL <-  22.5 * pi / 180 # FL sensor angle from forward position (degrees, 22.5 - 45)
FR <- -22.5 * pi / 180 # FR sensor angle from forward position (degrees, 22.5 - 45)
RA <-  45 * pi / 180 # Agent rotation angle (degrees)
SO <- 6 # Sensor offset distance (pixels)
SS <- 1 # Step size—how far agent moves per step (pixels) 
depT <- 15 # Chemoattractant deposition per step

iters <- 2000 # Number if iterations
agents <- 1000 # Number of agents

######################################################################
# Initialization of environment layer and particle positions
######################################################################
# Environment matrix, initialized with zeros
envM <- matrix(0 , imageH, imageW)

# Create a magnetic disc
for (i in 1:nrow(envM)){
  for (j in 1:ncol(envM)){
    if(sqrt((i-imageH/2)^2+(j-imageH/2)^2)>imageH/8 &
       sqrt((i-imageH/2)^2+(j-imageH/2)^2)<imageH/6) envM[i,j]=5
  }
}

# Place agents in a circle
tibble(h = seq(from = 0, to = 2*pi, length.out = agents)) %>% 
  mutate(x = (imageH/20)*cos(h)+imageH/2,
         y = (imageH/20)*sin(h)+imageH/2,
         h = jitter(h+pi, amount = 0) ) -> parF


# Make agents dance
envM <- physarum(envM, parF, decayT, FL, FR, 
                 RA, SO, SS, depT, iters)

# Transform resulting environment matrix into data frame
df <- melt(envM)
colnames(df) <- c("x","y","v") # to name columns

# Pick a top palette from colourlovers
palette <- sample(clpalettes('top'), 1)[[1]] 
colors <- palette %>% swatch %>% .[[1]]

# Do the plot
ggplot(data = df %>% filter(v>0), aes(x = x, y = y, fill = log(v))) + 
  geom_raster(interpolate = TRUE) +
  coord_equal() +
  scale_fill_gradientn(colours = colors) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0)) +
  opt -> plot

plot

# Do you like it? Save it!
ggsave("choose_a_name.png", plot, width  = 3, height = 3)