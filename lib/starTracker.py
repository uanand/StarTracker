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
def createKvector():
    dTable = pandas.read_csv('./dataset/dTable.dat',delimiter='\t')
    
    eps = 1e-10
    n = dTable.shape[0]
    yMin,yMax = dTable['distance'][0],dTable['distance'][n-1]
    m = (yMax-yMin+2*eps)/(n-1)
    c = yMin-eps
    jMin = 0
    
    outFile = open('./dataset/kVector.dat','w')
    outFile.write('Index\n0\n')
    for i in tqdm(range(1,n-1)):
        y = m*i+c
        for j in range(jMin,n-1):
            if (dTable['distance'][j]<=y and dTable['distance'][j+1]>y):
                outFile.write('%d\n' %(j+1))
                jMin = j
                break
    outFile.write('%d\n' %(n))
    outFile.close()
############################################################

############################################################
# K-VECTOR SEARCH
############################################################
def kVectorSearch(dMin,dMax):
    kVector = pandas.read_csv('./dataset/kVector.dat',delimiter='\t')
    dTable = pandas.read_csv('./dataset/dTable.dat',delimiter='\t')
    
    eps = 1e-10
    n = dTable.shape[0]
    yMin,yMax = dTable['distance'][0],dTable['distance'][n-1]
    m = (yMax-yMin+2*eps)/(n-1)
    c = yMin-eps
    
    startIndex = max(int(numpy.floor((dMin-c)/m)),0)
    endIndex = min(int(numpy.ceil((dMax-c)/m)),n-1)
    start = kVector['Index'][startIndex]
    end = kVector['Index'][endIndex]
    
    return start,end
############################################################
