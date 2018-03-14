# install and load necessary packages
if(!require(rgl)){
  install.packages('rgl')
  require(rgl)
}

if(!require(lidR)){
  install.packages('lidR')
  require(lidR)
}

# function for plotting 3d cylinders
plotCylinder = function(xCenter = 0, yCenter = 0, hBase = 0, hTop = 1, radius = 0.5, col = 'green'){
  
  axis = matrix(c(
    rep(xCenter, 2),
    rep(yCenter, 2),
    seq(hBase, hTop, length.out = 2)
  ), ncol = 3, byrow = F)
  
  cyl = cylinder3d(axis, radius = radius)
  
  # shade3d(addNormals(subdivision3d(cyl, depth = 0)), col = col)
  mesh = shade3d(cyl, col=col)
}


# set working directory
setwd('C://Work/TLStools/bin/Release/')

# set point cloud name
cloudName = 'square.las'

# set optional parameters for the las2rings application
optionalParameters = ''

# assemble the command and run las2rings over the point cloud
las2rings = paste('TLStrees -i', cloudName, optionalParameters)
system(las2rings)


# read the output files generated in the process
cloud = readLAS(cloudName)
stemsCloud = readLAS('square_trees.laz')
layerStack = readLAS('square_segmt.laz')
treeStats = read.table('square_reslt.txt', head=T)


# plot the different data
# --open a new rgl device with black background
clear3d()
bg3d('black')

# --load the points from the detected tree positions on multiple stacked layers
rgl.points(layerStack@data, col='blue')

# --load the stem points
rgl.points(stemsCloud@data, col='darkred')

# load the original point cloud
rgl.points(cloud@data, col='white', size=.5)

# fast - plot spheres with the corresponding diameters for every detected tree segment, at the right coordinates
# with(treeStats, spheres3d(x, y, (z_min + z_max) / 2, rad, col='yellow') )

# slow - plot cylinders instead of spheres
apply(treeStats, 1, function(row) plotCylinder(row['x'], row['y'], row['z_min'], row['z_max'], row['rad'], 'yellow'))

# plot the axes
axes3d()

