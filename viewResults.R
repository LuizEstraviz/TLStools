# install and load necessary packages
if(!require(rgl)){
  install.packages('rgl')
  require(rgl)
}

if(!require(lidR)){
  install.packages('lidR')
  require(lidR)
}


# set working directory
setwd('')


# set point cloud name
cloudName = 'lcer.las'

# set optional parameters for the las2rings application
optionalParameters = ''

# assemble the command and run las2rings over the point cloud
las2rings = paste('las2rings -i', cloudName, optionalParameters)
system(las2rings)


# read the output files generated in the process
cloud = readLAS('lcer.las')
stemsCloud = readLAS('lcer_stems.laz')
layerStack = readLAS('lcer_layerStack.laz')
treeStats = read.table('lcer_result.txt', head=T)


# plot the different data
# --open a new rgl device with black background
open3d()
bg3d('black')

# --load the points from the detected tree positions on multiple stacked layers
rgl.points(layerStack@data, col='blue')

# --load the stem points
rgl.points(stemsCloud@data, col='darkred')

# load the original point cloud
rgl.points(cloud@data, col='white', size=.5)

# plot spheres with the corresponding diameters for every detected tree segment, at the right coordinates
with(treeStats, spheres3d(x, y, (z_min + z_max) / 2, rad, col='yellow') )

# plot the axes
axes3d()

