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
        [self.catalogueRow,self.catalogueCol] = self.catalogue.shape
        if (os.path.exists('./dataset/dTable.dat')==False):
            self.createDTable()
        self.dTable = pandas.read_csv('./dataset/dTable.dat',delimiter='\t')    
        if (os.path.exists('./dataset/kVector.dat')==False):
            self.createKvector()
        self.kVector = pandas.read_csv('./dataset/kVector.dat',delimiter='\t')
    ############################################################
    
    ############################################################
    # READ THE CATALOGUE, SAVE RA AND DEC ANGLES IN RADIANS AND
    # SAVE AS A GLOBAL CLASS DATAFRAME 
    ############################################################
    def readCatalogue(self,fileName,skiprows=1):
        catalogue = numpy.loadtxt(fileName,dtype='float64',skiprows=skiprows)
        catalogue = pandas.DataFrame(data=catalogue,columns=['RA(deg)','DEC(deg)','HIPNum','Mag'])
        catalogue['HIPNum'] = catalogue['HIPNum'].astype('int32')
        catalogue['RA(rad)'] = deg2rad(catalogue['RA(deg)'])
        catalogue['DEC(rad)'] = deg2rad(catalogue['DEC(deg)'])
        alphaR,deltaR = catalogue['RA(rad)'],catalogue['DEC(rad)']
        catalogue['x'] = cos(deltaR)*cos(alphaR)
        catalogue['y'] = cos(deltaR)*sin(alphaR)
        catalogue['z'] = sin(deltaR)
        self.catalogue = catalogue
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
                    outFile.write('%d\t%d\t%.10f\n' %(self.catalogue['HIPNum'][i],self.catalogue['HIPNum'][j],theta))
        outFile.close()
        
        dTable = pandas.read_csv('./dataset/dTable.dat',delimiter='\t')
        dTable.sort_values(by='distance',inplace=True)
        outFile = open('./dataset/dTable.dat','w')
        outFile.write('HIPNum1\tHIPNum2\tdistance\n')
        for i in tqdm(dTable.index):
            outFile.write('%d\t%d\t%.10f\n' %(dTable['HIPNum1'][i],dTable['HIPNum2'][i],dTable['distance'][i]))
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
        row = stars.shape[0]
        stars['id'] = range(1,row+1)
        
        ############################################################
        # INTRODUCING ERROR IN STAR POSITION
        errRA = (numpy.random.rand(row)-0.5)*error
        errDEC = (numpy.random.rand(row)-0.5)*error
        stars['RA(deg)'] = stars['RA(deg)']+errRA
        stars['DEC(deg)'] = stars['DEC(deg)']+errDEC
        stars['RA(rad)'] = deg2rad(stars['RA(deg)'])
        stars['DEC(rad)'] = deg2rad(stars['DEC(deg)'])
        
        ############################################################
        # CONVERTING TO CARTESIAN COORDINATES
        alphaR,deltaR = stars['RA(rad)'],stars['DEC(rad)']
        stars['x'] = cos(deltaR)*cos(alphaR)
        stars['y'] = cos(deltaR)*sin(alphaR)
        stars['z'] = sin(deltaR)
        
        ############################################################
        # CALCULATING THE MAXIMUM POSSIBLE ERROR
        angle1,angle2 = deg2rad(error),-deg2rad(error)
        vec1 = [cos(angle1)*cos(angle1),cos(angle1)*sin(angle1),sin(angle1)]
        vec2 = [cos(angle2)*cos(angle2),cos(angle2)*sin(angle2),sin(angle2)]
        maxError = numpy.arccos(dot(vec1,vec2))
        stars['Error'] = maxError
        self.stars = stars
    ############################################################
    
    ############################################################
    # CREATE THE DISTANCE TABLE FOR SIMULATED CAMERA IMAGE
    ############################################################
    def createTTable(self):
        [row,col] = self.stars.shape
        counter = 0
        tTable = numpy.zeros([int(row*(row-1)/2),8])
        for i in range(row-1):
            for j in range(i+1,row):
                vector1 = [self.stars['x'][i],self.stars['y'][i],self.stars['z'][i]]
                vector2 = [self.stars['x'][j],self.stars['y'][j],self.stars['z'][j]]
                distance = rad2deg(numpy.arccos(dot(vector1,vector2)))
                error = rad2deg(self.stars['Error'][i])+rad2deg(self.stars['Error'][j])
                tTable[counter][0] = i+1
                tTable[counter][1] = j+1
                tTable[counter][2] = distance
                tTable[counter][3] = error
                tTable[counter][4] = distance-error
                tTable[counter][5] = distance+error
                tTable[counter][6] = self.stars['HIPNum'][i]
                tTable[counter][7] = self.stars['HIPNum'][j]
                counter += 1
        tTable = pandas.DataFrame(data=tTable,columns=['id1','id2','distance','error','distance-error','distance+error','HIPNum1','HIPNum2'])
        tTable['id1'] = tTable['id1'].astype('int32')
        tTable['id2'] = tTable['id2'].astype('int32')
        tTable['HIPNum1'] = tTable['HIPNum1'].astype('int32')
        tTable['HIPNum2'] = tTable['HIPNum2'].astype('int32')
        self.tTable = tTable
    ############################################################
    
    ############################################################
    # K-VECTOR SEARCH FOR ALL THE DISTANCE IN T-TABLE
    ############################################################
    def kVectorSearch(self):
        eps = 1e-10
        n = self.dTable.shape[0]
        yMin,yMax = self.dTable['distance'][0],self.dTable['distance'][n-1]
        m = (yMax-yMin+2*eps)/(n-1)
        c = yMin-eps
        
        startList,endList = [],[]
        row = self.tTable.shape[0]
        for i in range(row):
            startIndex = max(int(numpy.floor((self.tTable['distance-error'][i]-c)/m)),0)
            endIndex = min(int(numpy.ceil((self.tTable['distance-error'][i]-c)/m)),n-1)
            start = self.kVector['Index'][startIndex]
            end = self.kVector['Index'][endIndex]
            for j in range(start,-1,-1):
                if (self.dTable['distance'][j]<=self.tTable['distance-error'][i]):
                    start = j
                    break
            for j in range(start,n):
                if (self.dTable['distance'][j]>=self.tTable['distance-error'][i]):
                    start = j-1
                    break
            for j in range(end,n):
                if (self.dTable['distance'][j]>=self.tTable['distance+error'][i]):
                    end = j
                    break
            for j in range(end,-1,-1):
                if (self.dTable['distance'][j]<=self.tTable['distance+error'][i]):
                    end = j+1
                    break
            startList.append(start)
            endList.append(end)
        self.tTable['start'] = startList
        self.tTable['end'] = endList
    ############################################################
    
    ############################################################
    # CHECK IF THE ENTRIES OF  T TABLE ARE CORRECT
    ############################################################
    def checkTTable(self):
        row = self.tTable.shape[0]
        for i in range(row):
            start,end = self.tTable['start'][i],self.tTable['end'][i]
            print (self.dTable.loc[start:end])
    ############################################################
    
    ############################################################
    # RUN THE PYRAMID STAR IDENTIFICATION ALGORITHM
    ############################################################
    def pyramid(self):
        candidates = {}
        numStars = self.stars.shape[0]
        for i in range(1,numStars+1):
            candidates[i] = {}
            for j in range(1,numStars+1):
                candidates[i][j] = []
                
        row = self.tTable.shape[0]
        for i in range(row):
            start,end = self.tTable['start'][i],self.tTable['end'][i]
            star1,star2 = self.tTable['id1'][i],self.tTable['id2'][i]
            for j in range(start,end+1):
                candidates[star1][star2].append([self.dTable['HIPNum1'][j],self.dTable['HIPNum2'][j]])
                candidates[star1][star2].append([self.dTable['HIPNum2'][j],self.dTable['HIPNum1'][j]])
                candidates[star2][star1].append([self.dTable['HIPNum1'][j],self.dTable['HIPNum2'][j]])
                candidates[star2][star1].append([self.dTable['HIPNum2'][j],self.dTable['HIPNum1'][j]])
                
        for i in range(1,numStars-1):
            for j in range(i+1,numStars):
                for k in range(j+1,numStars+1):
                    star1,star2,star3,flag = findTriangle(candidates[i][j],candidates[j][k],candidates[k][i])
                    if (flag==1):
                        break
                if (flag==1):
                    break
            if (flag==1):
                break
                
        for i in range(1,numStars+1):
            if (i==star1 or i==star2 or i==star3):
                pass
            else:
                star4,flag = findPyramid(candidates[star1][i],candidates[star2][i],candidates[star3][i])
    ############################################################
    
    ############################################################
    # RUN THE GEOMETRIC-VOTING-PYRAMID STAR IDENTIFICATION ALGORITHM
    ############################################################
    def geometricPyramid(self):
        pass
    ############################################################
