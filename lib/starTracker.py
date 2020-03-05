import os
import sys
import numpy
import pandas
from tqdm import tqdm
from numpy import sin,cos,deg2rad,rad2deg,cross,dot

class starTracker:
    ############################################################
    # INITIALIZATION AFTER OBJECT IS DEFINED
    ############################################################
    def __init__(self,catalogueFileName,FOVDegree):
        self.catalogueFileName = catalogueFileName
        self.FOVDegree = FOVDegree
        self.readCatalogue(catalogueFileName)
        self.addCartesianCoordinates()
        [self.catalogueRow,self.catalogueCol] = self.catalogue.shape
        if (os.path.exists('./dataset/dTable.dat')==False):
            self.createDTable()
        if (os.path.exists('./dataset/kVector.dat')==False):
            self.createKvector()
        self.dTable = pandas.read_csv('./dataset/dTable.dat',delimiter='\t')
        self.kVector = pandas.read_csv('./dataset/kVector.dat',delimiter='\t')
    ############################################################
    
    ############################################################
    # READ THE CATALOGUE, SAVE RA AND DEC ANGLES IN RADIANS AND
    # SAVE AS A GLOBAL CLASS DATAFRAME 
    ############################################################
    def readCatalogue(self,fileName,skiprows=1):
        catalogue = numpy.loadtxt(fileName,dtype='float64',skiprows=skiprows)
        self.catalogue = pandas.DataFrame(data=catalogue,columns=['RA(deg)','DEC(deg)','HIPNum','Mag'])
        self.catalogue['HIPNum'] = self.catalogue['HIPNum'].astype('int32')
        self.catalogue['RA(rad)'] = deg2rad(self.catalogue['RA(deg)'])
        self.catalogue['DEC(rad)'] = deg2rad(self.catalogue['DEC(deg)'])
    ############################################################
    
    ############################################################
    # ADD CARTESIAN COORDINATES OF STARS TO CATALOGUE
    ############################################################
    def addCartesianCoordinates(self):
        [row,col] = self.catalogue.shape
        alphaR,deltaR = self.catalogue['RA(rad)'],self.catalogue['DEC(rad)']
        self.catalogue['x'] = cos(deltaR)*cos(alphaR)
        self.catalogue['y'] = cos(deltaR)*sin(alphaR)
        self.catalogue['z'] = sin(deltaR)
    ############################################################
    
    ############################################################
    # CREATE THE DISTANCE TABLE BASED ON FOV
    ############################################################
    def createDTable(self):
        print ('Inter-star distance table not available. Creating now.')
        outFile = open('./dataset/dTable.dat','w')
        outFile.write('HIPNum1\tHIPNum2\tdistance\n')
        [row,col] = self.catalogue.shape
        for i in tqdm(range(row-1)):
            v1 = [self.catalogue['x'][i],self.catalogue['y'][i],self.catalogue['z'][i]]
            for j in range(i+1,row):
                v2 = [self.catalogue['x'][j],self.catalogue['y'][j],self.catalogue['z'][j]]
                dotProduct = dot(v1,v2)
                if (dotProduct<-1.0):
                    dotProduct = -1.0
                elif (dotProduct>1.0):
                    dotProduct = 1.0
                theta = rad2deg(numpy.arccos(dotProduct))
                if (theta<=self.FOVDegree):
                    outFile.write('%d\t%d\t%f\n' %(self.catalogue['HIPNum'][i],self.catalogue['HIPNum'][j],theta))
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
    # CREATE K-VECTOR FOR SINGLE SHOT SEARCHING
    ############################################################
    def createKvector(self):
        eps = 1e-10
        n = self.dTable.shape[0]
        yMin,yMax = self.dTable['distance'][0],self.dTable['distance'][n-1]
        m = (yMax-yMin+2*eps)/(n-1)
        c = yMin-eps
        jMin = 0
        
        outFile = open('./dataset/kVector.dat','w')
        outFile.write('Index\n0\n')
        for i in tqdm(range(1,n-1)):
            y = m*i+c
            for j in range(jMin,n-1):
                if (self.dTable['distance'][j]<=y and self.dTable['distance'][j+1]>y):
                    outFile.write('%d\n' %(j+1))
                    jMin = j
                    break
        outFile.write('%d\n' %(n))
        outFile.close()
    ############################################################
    
    ############################################################
    # SIMULATE CAMERA USING THE THREE ANGLES FOR ORIENTATION
    # NUMBER OF BRIGHTEST STARS WE WANT TO SIMULATE
    # UNIFORMLY DISTRIBUTED ERROR IN CENTROID MEASUREMENT OF STARS
    ############################################################
    def simCam(self,psi,phi,theta,num,error):
        SF2 = sin(deg2rad(self.FOVDegree/2))
        
        X = [cos(theta)*cos(psi) - sin(theta)*sin(phi)*sin(psi),\
             cos(theta)*sin(psi) + sin(theta)*sin(phi)*cos(psi),\
             -sin(theta)*cos(phi)\
            ]
        Y = [-cos(phi)*sin(psi),\
             cos(phi)*cos(psi),\
             sin(phi)\
            ]
        Z = cross(X,Y)
        
        counter = 0
        for i in range(self.catalogueRow):
            alpha,delta = self.catalogue['RA(rad)'][i],self.catalogue['DEC(rad)'][i]
            V = [cos(delta)*cos(alpha),\
                 cos(delta)*sin(alpha),\
                 sin(delta)\
                ]
            VX = abs(dot(V,X))
            VY = abs(dot(V,Y))
            VZ = dot(V,Z)
            if (VZ>0 and VX<=SF2 and VY<=SF2):
                counter+=1
                if (counter==1):
                    stars = numpy.array([self.catalogue['RA(deg)'][i],self.catalogue['DEC(deg)'][i],self.catalogue['HIPNum'][i],self.catalogue['Mag'][i]])
                else:
                    stars = numpy.row_stack((stars,numpy.array([self.catalogue['RA(deg)'][i],self.catalogue['DEC(deg)'][i],self.catalogue['HIPNum'][i],self.catalogue['Mag'][i]])))
        stars = pandas.DataFrame(data=stars,columns=['RA(deg)','DEC(deg)','HIPNum','Mag'])
        stars['HIPNum'] = stars['HIPNum'].astype('int32')
        stars.sort_values('Mag',ascending=False,inplace=True)
        stars.reset_index(drop=True,inplace=True)
        [row,col] = stars.shape
        
        ############################################################
        # SELECTING NUM BRIGHTEST STARS
        if (row>num):
            stars = stars.head(num)
        else:
            print ('Found only %d stars.' %(row))
        [row,col] = stars.shape
        
        ############################################################
        # INTRODUCING ERROR IN STAR POSITION
        errRA = (numpy.random.rand(row)-0.5)*error
        errDEC = (numpy.random.rand(row)-0.5)*error
        stars['RA(deg)'] = stars['RA(deg)']+errRA
        stars['DEC(deg)'] = stars['DEC(deg)']+errDEC
        self.stars = stars
    ############################################################
    


# ############################################################
# # K-VECTOR SEARCH
# ############################################################
# def kVectorSearch(dMin,dMax):
    # kVector = pandas.read_csv('./dataset/kVector.dat',delimiter='\t')
    # dTable = pandas.read_csv('./dataset/dTable.dat',delimiter='\t')
    
    # eps = 1e-10
    # n = dTable.shape[0]
    # yMin,yMax = dTable['distance'][0],dTable['distance'][n-1]
    # m = (yMax-yMin+2*eps)/(n-1)
    # c = yMin-eps
    
    # startIndex = max(int(numpy.floor((dMin-c)/m)),0)
    # endIndex = min(int(numpy.ceil((dMax-c)/m)),n-1)
    # start = kVector['Index'][startIndex]
    # end = kVector['Index'][endIndex]
    
    # return start,end
# ############################################################
