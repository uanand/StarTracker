import os
import sys

sys.path.append(os.path.abspath('./lib'))
import starTracker
import utils

catalogueFileName = './dataset/scHIP4to6p5.txt'
FOVDegree = 20

catalogue = starTracker.readCatalogue(catalogueFileName)
catalogue = starTracker.addCartesianCoordinates(catalogue)
starTracker.createDistanceTable(catalogue,FOVDegree)


