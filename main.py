import os
import sys
import pandas

sys.path.append(os.path.abspath('./lib'))
from starTracker import starTracker
import utils

catalogueFileName = './dataset/scHIP4to6p5.txt'
FOVDegree = 5

ST = starTracker(catalogueFileName,FOVDegree)
ST.simCam(1,1,1,5,0.0)

numStars = ST.stars.shape[0]
if (numStars>=3):
    ST.createTTable()
    ST.kVectorSearch()
    ST.pyramid()
    ST.geometricPyramid()
    
