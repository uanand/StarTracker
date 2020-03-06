import os
import sys
import pandas

sys.path.append(os.path.abspath('./lib'))
from starTracker import starTracker
import utils

catalogueFileName = './dataset/scHIP4to6p5.txt'
FOVDegree = 20

ST = starTracker(catalogueFileName,FOVDegree)
ST.simCam(1,1,1,5,0.00001)

numStars = ST.stars.shape[0]
if (numStars>=3):
    ST.createTTable()
    ST.kVectorSearch()
    print (ST.tTable)
# start,end = starTracker.kVectorSearch(0.450,0.455)
