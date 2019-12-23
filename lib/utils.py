import os
import sys
import numpy
import pandas
from tqdm import tqdm

############################################################
# MAKE DIRECTORY
############################################################
def mkdir(dirName):
    if (os.path.exists(dirName) == False):
        os.makedirs(dirName)
############################################################
