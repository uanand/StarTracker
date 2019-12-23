import os
import sys
import numpy
import pandas
from tqdm import tqdm

############################################################
# READ THE STAR CATALOGUE AS A NUMPY ARRAY
############################################################
def readCatalogue(fileName,skiprows=1):
    catalogue = numpy.loadtxt(fileName,dtype='float64',skiprows=skiprows)
    catalogue = pandas.DataFrame(data=catalogue,columns=['RA(deg)','DEC(deg)','HIPNum','Mag'])
    catalogue['HIPNum'] = catalogue['HIPNum'].astype('int32')
    return catalogue
############################################################

############################################################
# CONVERT THE CATALOGUE STAR ORIENTATION TO UNIT VECTOR
############################################################
def addCartesianCoordinates(catalogue):
    [row,col] = catalogue.shape
    alphaR,deltaR = numpy.deg2rad(catalogue['RA(deg)']),numpy.deg2rad(catalogue['DEC(deg)'])
    catalogue['x'] = numpy.cos(deltaR)*numpy.cos(alphaR)
    catalogue['y'] = numpy.cos(deltaR)*numpy.sin(alphaR)
    catalogue['z'] = numpy.sin(deltaR)
    return catalogue
############################################################

############################################################
# CALCULATE AND SAVE DISTANCE TABLE
############################################################
def createDistanceTable(catalogue,FOVDegree):
    outFile = open('./dataset/dTable.dat','w')
    outFile.write('HIPNum1\tHIPNum2\tdistance\n')
    [row,col] = catalogue.shape
    for i in tqdm(range(row-1)):
        v1 = [catalogue['x'][i],catalogue['y'][i],catalogue['z'][i]]
        for j in range(i+1,row):
            v2 = [catalogue['x'][j],catalogue['y'][j],catalogue['z'][j]]
            dotProduct = numpy.dot(v1,v2)
            if (dotProduct<-1.0):
                dotProduct = -1.0
            elif (dotProduct>1.0):
                dotProduct = 1.0
            theta = numpy.rad2deg(numpy.arccos(dotProduct))
            if (theta<=FOVDegree):
                outFile.write('%d\t%d\t%f\n' %(catalogue['HIPNum'][i],catalogue['HIPNum'][j],theta))
    outFile.close()
    
    dTable = pandas.read_csv('./dataset/dTable.dat',delimiter='\t')
    dTable.sort_values(by='distance',inplace=True)
    outFile = open('./dataset/dTable.dat','w')
    outFile.write('HIPNum1\tHIPNum2\tdistance\n')
    for i in tqdm(dTable.index):
        outFile.write('%d\t%d\t%f\n' %(dTable['HIPNum1'][i],dTable['HIPNum2'][i],dTable['distance'][i]))
    outFile.close()
############################################################

############################################################
# CALCULATE K-VECTOR
############################################################
def createKvector(self):
    pass
############################################################
