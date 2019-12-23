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
    for i in tqdm(range(50)):#tqdm(range(row-1)):
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
    for i in dTable.index:
        outFile.write('%d\t%d\t%f\n' %(dTable['HIPNum1'][i],dTable['HIPNum2'][i],dTable['distance'][i]))
    outFile.close()
    return dTable
############################################################

############################################################
# CALCULATE K-VECTOR
############################################################
def createKvector(self):
    distanceTable = []
    for i in range(self.catalogueRow-1):
        for j in range(i+1,self.catalogueRow):
            distance = arccos(dot(self.catalogue[i,0:3], self.catalogue[j,0:3]))
            if (distance <= self.fovR):
                distanceTable.append([self.catalogue[i,3], self.catalogue[j,3], distance])
    distanceTable = numpy.asarray(distanceTable)
    distanceTable = numpy.sort(distanceTable.view('f8,f8,f8'), order=['f2'], axis=0).view(numpy.float64)
    mkdir('preProcess')
    numpy.save('./preProcess/catalogueDistanceTable.npy', distanceTable)
############################################################






# def readCatalogue(self, skiprows=1):
        # print "READING THE CATALOGUE AND CONVERTING TO VECTORS"
        # catalogue = numpy.loadtxt('./catalogue/'+self.catalogueName, dtype='float64', skiprows=skiprows)
        # alphaR, deltaR = deg2rad(catalogue[:,0]), deg2rad(catalogue[:,1])
        # [row,col] = catalogue.shape
        # self.catalogue = numpy.column_stack(((cos(deltaR)*cos(alphaR)), cos(deltaR)*sin(alphaR), sin(deltaR), catalogue[:,2], catalogue[:,3]))
        # [self.catalogueRow, self.catalogueCol] = self.catalogue.shape


############################################################
