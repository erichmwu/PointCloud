#####################################################################################################################
#
# (C) 2017 - Eric Wu, Laboratory for Advanced Visualization & Applications, University of Hawaii at Manoa.
#
# Version 03/10/2017
#
# An octree data structure used to find an appropriate splat size for any given data set.
#
#####################################################################################################################

import sys
import numpy as np
import math
from scipy.spatial.distance import pdist

# Computes the distance between 2 points
def distance(p0, p1):
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 + (p0[2] - p1[2])**2)

# Finds the most common distance between all the points
def find_closest(bin):
    if len(bin) < 2:
        return 10000000000

    # np.set_printoptions(precision=6)
    distances = pdist(bin, 'euclidean')[np.nonzero(pdist(bin, 'euclidean'))]
    sortedDistance = np.asarray(np.unique(distances.round(decimals=20), return_counts=True)).T
    sortedDistance = sortedDistance[sortedDistance[:,1].argsort()[::-1]]
    # print sortedDistance
    
    ideal = np.where(sortedDistance[:,1] == sortedDistance[0][1])[0]
    maxAppearances = sortedDistance[ideal]
    
    mostCommonDistances = maxAppearances[maxAppearances[:,0].argsort()]
    
    return mostCommonDistances[0][0]

# Recursive octree traversal down the branch with the most points
def octree(bin):
    if len(bin) < 10000:
        return find_closest(bin)
        
    # print len(bin)    
        
    bins =[[], [], [], [], [], [], [], []]
        
    xyzmin = np.min(np.array(bin), axis=0)
    xyzmax = np.max(np.array(bin), axis=0)

    origin = [(xyzmax[0]-xyzmin[0])/2+xyzmin[0], (xyzmax[1]-xyzmin[1])/2+xyzmin[1], (xyzmax[2]-xyzmin[2])/2+xyzmin[2]]
    # print xyzmin
    # print xyzmax
    # print origin
    
    for point in bin:
        # print point
        if (point[0] < origin[0]) and (point[1] > origin[1]) and (point[2] < origin[2]):
            #print 1
            bins[0].append(point)
        if (point[0] > origin[0]) and (point[1] > origin[1]) and (point[2] < origin[2]):
            #print 2
            bins[1].append(point)
        if (point[0] < origin[0]) and (point[1] < origin[1]) and (point[2] < origin[2]):
            #print 3
            bins[2].append(point)
        if (point[0] > origin[0]) and (point[1] < origin[1]) and (point[2] < origin[2]):
            #print 4
            bins[3].append(point)
        if (point[0] < origin[0]) and (point[1] > origin[1]) and (point[2] > origin[2]):
            #print 5
            bins[4].append(point)
        if (point[0] > origin[0]) and (point[1] > origin[1]) and (point[2] > origin[2]):
            #print 6
            bins[5].append(point)
        if (point[0] < origin[0]) and (point[1] < origin[1]) and (point[2] > origin[2]):
            #print 7
            bins[6].append(point)
        if (point[0] > origin[0]) and (point[1] < origin[1]) and (point[2] > origin[2]):
            #print 8
            bins[7].append(point)
            
    densestBin = [len(bins[x]) for x in range(0, 8)]
    
    # print densestBin
    
    densestBin = [i[0] for i in sorted(enumerate(densestBin), key=lambda x:x[1], reverse=True)]
    
    return octree(bins[densestBin[0]])

# Main function which handles the reading of the data file and finding the splat size
def main(arg1):    
    xyz = []
    rgba = []
        
    with open(arg1, 'r') as infile:
        for line in infile:
            value = line.split()
            if len(value) > 1:
                xyz_val = [float(i) for i in value[0:3]]
                rgba_val = [float(i) for i in value[3:]]
                xyz.append(xyz_val)
                rgba.append(rgba_val)

    xyzmin = np.min(np.array(xyz), axis=0)
    xyzmax = np.max(np.array(xyz), axis=0)
    numPoints = len(xyz)
    
    print numPoints, xyzmin[0], xyzmin[1], xyzmin[2], xyzmax[0], xyzmax[1], xyzmax[2], octree(xyz),
    
if __name__=='__main__':
    sys.exit(main(sys.argv[1]))
