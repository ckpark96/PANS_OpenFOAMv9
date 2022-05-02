import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate

''' Dataset without duplpicates '''
# Obtain general dataset from Paraview:
# import msh file by creating .foam file -> Save Data
data = pd.read_csv('dataarrays.csv', header=0) # Not required if the cell centers values are obtained (as done below)
epsData = data['epsilon']
kData = data['k']
volData = data['V']

centerData = pd.read_csv('cellCenters.csv', header=0)
xcoord = centerData['Points_0']
ycoord = centerData['Points_1']

''' Dataset with duplicates '''
# Obtain cell length dataset from Paraview - these contain duplicates:
# import msh file by creating .foam file -> filter 'Extract Edges' -> filter 'Cell Size' -> Save Data
dupData = pd.read_csv('lengthData0.0.csv', header=0)
lenData = dupData['Length']
dupVolData = dupData['V']


##### Constants #####
Cmu = 0.09

''' Extract k and epsilon data from OpenFOAM files '''
##### Constants #####
dt = 0.01
endT = 0.3
rowskip = 21

def clean_data(filename, rowstoskip):
    '''
    For OpenFOAM files

    Cleans k or epsilon data, leaving only the relevant internalField values
    Thus, exclude boundaryField values and wall values
    '''
    data = pd.read_csv(filename, header=None, skiprows = rowstoskip)
    data = pd.to_numeric(data[0], errors='coerce') # convert non-numerics to NaN
    dataNaN = data.isnull()
    idx = np.where(dataNaN)[0][0]
    data = data[:idx]
    data = data.to_frame()

    return data.to_numpy()


def time_avg_data(timeStep, latestTime, rowstoskip): # Not required if temporal statistics filter is used in Paraview
    '''
    For OpenFOAM files
    
    Computes the time averages of k and epsilon
    '''
    variousTime = np.arange(timeStep, latestTime, timeStep)
    for count, t in enumerate(variousTime):
        if '9' in str(count): # Every 9th count does the t only have 1 decimal point
            t = np.round(t,1)
        else:
            t = np.round(t,2)
        kfile = clean_data(str(t)+'/k', rowstoskip)
        epsfile = clean_data(str(t)+'/epsilon', rowstoskip)
        
        if count == 0:
            kBigData = np.zeros((kfile.shape[0],variousTime.shape[0]))
            epsBigData = np.zeros((epsfile.shape[0],variousTime.shape[0]))
            kBigData[:,count] = kfile.reshape((kfile.shape[0]))
            epsBigData[:,count] = epsfile.reshape(epsfile.shape[0])
    
    kAvg = np.mean(kBigData, axis=1)
    epsAvg = np.mean(epsBigData, axis=1)

    return kAvg, epsAvg


''' Functions needed to extract minimum, maximum and average cell lengths '''

def drop_successive(df):
    '''
    Instead of dropping all the duplicates, this function only drops the consecutive duplicates

    Input:
    df should be a pandas.DataFrame of a a pandas.Series
    Output:
    df of ts with masked or dropped values
    '''
    #### Method 1 #####
    # df = df.mask(df.shift(1) == df) # Mask keeping the first occurrence
    # return df.dropna(axis=0, how='all') # Drop the values (e.g. rows are deleted)

    #### Method 2 ####
    dupBool = df.shift() != df # Mark those successive duplicates as False (non-duplicates as True)
    newdf = df.loc[dupBool] # Extract the values with True
    return newdf, dupBool

# cleanVolData, dupBoolean = drop_successive(dupVolData) # Only needed if taking from lengthData csv
# firstData = dupBoolean[dupBoolean].index # Extracts indices of True values
# invDupBoolean = ~dupBoolean



def consecutive_ranges(data, val):
    '''
    Get ranges for CONSECUTIVE duplicates
    '''
    isval = np.concatenate(([0], np.equal(data, val).view(np.int8), [0]))
    absdiff = np.abs(np.diff(isval))
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)

    return ranges

# test = np.array([np.nan, np.nan, 1, 0, 0, 2, 0, 0])


def get_length_data(duplicate_data, length_data, whichstat):
    '''
    Takes length_data and produces a statistical array with minimum, median and maximum length for each cell
    ### Not the most efficient method (but does the job) ###

    Input:
    duplicate_data: raw (volume) data that is inclusive of duplicates
    length_data: length data
    whichstat: 'max', 'min' or 'avg' for maximum, minimum and average values of each cell respectively
    '''
    non_duplicate_data = drop_successive(duplicate_data)[0]
    lenDataCopy = copy.deepcopy(length_data) # to not dirty the original length data
    lenStat = np.zeros((non_duplicate_data.shape[0],3))

    for count, val in enumerate(non_duplicate_data):
        _range = consecutive_ranges(duplicate_data, val)
        maxVal = np.nanmax(lenDataCopy[_range[0,0]:_range[0,1]])
        minVal = np.nanmin(lenDataCopy[_range[0,0]:_range[0,1]])
        midVal = np.nanmedian(lenDataCopy[_range[0,0]:_range[0,1]])

        if midVal==np.nan or minVal==np.nan or maxVal==np.nan:
            print(count, val)
        lenStat[count] = np.array([minVal, midVal, maxVal])

        for i in np.arange(_range[0,0],_range[0,1]):
            lenDataCopy[i] = np.random.rand() # Replace already processed values with random values to not include in range of next occurrences 

    if whichstat == 'min':
        return lenStat[:,0]
    elif whichstat == 'max':
        return lenStat[:,2]
    elif whichstat == 'avg':
        return np.mean(lenStat, axis=1)
    elif whichstat == 'all':
        return lenStat


kTimeAvg, epsTimeAvg = time_avg_data(dt, endT, rowskip)
lenStat = get_length_data(dupVolData, lenData, 'avg')
taylorLen = kTimeAvg**(3/2) / epsTimeAvg # Turbulence Taylor scale length
optfk = 1/np.sqrt(Cmu)*(lenStat / taylorLen)**(2/3) # (One of the many) formula for optimal f_k

def visualise_fk(x_coords, y_coords, fk):
    '''
    x_coords: x-coordinates of each cell (Obtained from Paraview via 'Cell center' filter)
    y_coords: y-coordinates of each cell (")
    # mask_cond: tuple of masking conditions - x first then y. e.g. (x>1, x<1) TO BE IMPLEMENTED
    '''
    dataLen = x_coords.shape[0]
    x = np.linspace(np.min(x_coords), np.max(x_coords), dataLen)
    y = np.linspace(np.min(y_coords), np.max(y_coords), dataLen)
    xx, yy = np.meshgrid(x,y)
    Z = interpolate.griddata((x_coords, y_coords), fk, (xx,yy),method='linear')
    maskWhere = np.add((xx<0)*1, (yy<0)*1) >1 ### TO BE MANUALLY SET
    Z = np.ma.array(Z, mask=maskWhere)
    plt.figure()
    plt.imshow(Z, origin='lower',
            extent=(xcoord.min(), xcoord.max(), ycoord.min(), ycoord.max()),
            aspect=(xcoord.max() - xcoord.min()) / (ycoord.max() - ycoord.min()))
    plt.colorbar()
    plt.show()

# visualise_fk(xcoord, ycoord, optfk)

# plt.figure()
# plt.contour(xx, yy, Z)
# plt.colorbar()
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(xcoord, ycoord, optfk, c='red')
# plt.show()