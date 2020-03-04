import os
import sys
import pandas

sys.path.append(os.path.abspath('./lib'))
import starTracker
import utils

catalogueFileName = './dataset/scHIP4to6p5.txt'
FOVDegree = 20

# catalogue = starTracker.readCatalogue(catalogueFileName)
# catalogue = starTracker.addCartesianCoordinates(catalogue)
# starTracker.createDistanceTable(catalogue,FOVDegree)
# starTracker.createKvector()
start,end = starTracker.kVectorSearch(0.450,0.455)

# kVector = pandas.read_csv('./dataset/kVector.dat')
# dTable = pandas.read_csv('./dataset/dTable.dat',delimiter='\t')
# dTable = dTable.head(1000)

# eps = 1e-10
# n = dTable.shape[0]
# yMin,yMax = dTable['distance'][0],dTable['distance'][n-1]
# m = (yMax-yMin+2*eps)/(n-1)
# q = yMin-m-eps

