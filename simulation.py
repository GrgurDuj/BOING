import math as m
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import shapely.geometry as sg
from shapely.geometry.polygon import Polygon
from Point import Point

def orbit_spherical(inclination, right_ascension, time):  # time is in seconds
    clockangle = (time / period) * 2 * m.pi
    latitude = m.sin(clockangle) * inclination
    longitude = clockangle + right_ascension - (
            (time / day) * 2 * m.pi)  # rotation of the Earth beneath a prograde orbit
    if (latitude > 0.5 * m.pi):
        latitude = m.pi - latitude
        longitude += m.pi
    if (latitude < -0.5 * m.pi):
        latitude = -(m.pi + latitude)
        longitude += m.pi
    coordinate = [latitude, longitude]
    return coordinate

def STXY(latitude, longitude):  # lo = longitude is perpendicular to equator, la = latitude is parallel to equator
    # Gall-Peters projection
    x = longitude % (2 * m.pi)  # longitude in radians from central meridian
    y = 2 * m.sin(latitude)  # latitude in radians from equator
    coordinate = Point(x, y)
    return coordinate

def SunPosition(time):
    Earth_tilt = (23.5 / 180) * m.pi  # Earth's tilt in radians
    solar_inclination = Earth_tilt * m.sin(
        (time / (day * 365.25)) * 2 * m.pi)  # declination is 0 on March 21st (80th day of the year)
    solar_right_ascension = (time / day) * 2 * m.pi  # right ascension of sun rising
    position = Point(solar_inclination, solar_right_ascension)
    return position

def sunCircle(time):
    Earth_tilt = (23.5 / 180) * m.pi  # Earth's tilt in radians
    solar_inclination = Earth_tilt * m.sin(
        (time / (day * 365.25)) * 2 * m.pi)  # declination is 0 on March 21st (80th day of the year)
    solar_right_ascension = (time / day) * 2 * m.pi  # right ascension of sun rising
    r = 0.5 * m.pi
    n = 64
    pos = np.empty(n, dtype=Point)
    for i in range(n):
        theta = i * ((2 * m.pi) / n)
        pos[i] = Point(m.sin(theta) * r + solar_inclination, m.cos(theta) * r + solar_right_ascension)
        if (pos[i].latitude > 0.5 * m.pi):
            pos[i].latitude = m.pi - pos[i].latitude
            pos[i].longitude += m.pi
        if (pos[i].latitude < -0.5 * m.pi):
            pos[i].latitude = -(m.pi + pos[i].latitude)
            pos[i].longitude += m.pi
    return pos

def calcPeriod(a):
    mu = 3.986E14
    period = 2 * m.pi * m.sqrt((a ** 3) / mu)
    return period

def visibleAreaPitchAngle(altitude):
    minPitch = (0.7/180)*m.pi  # rad
    maxPitch = (6.08/180)*m.pi  # rad5
    minAltitude = 450000  # m
    maxAltitude = 700000  # m
    interpolateLength = interp1d([minAltitude, maxAltitude], [minPitch, maxPitch])
    return interpolateLength(altitude)

def visibleAreaRollAngle(altitude):
    # corresponds to 30deg of roll from sat perspective (15 each side)
    roll = ((0.0048310913 * (altitude / 1000) + 0.00623608)/180)*m.pi # rad
    #print(roll, "roll")
    # roll gives observable angle with respect to the center of earth
    #print(roll, altitudes)
    return roll

def satellitePrism(inclination, right_ascension, time):
    # satellite position
    clockangle = (time / period) * 2 * m.pi
    latitude = m.sin(clockangle) * inclination
    longitude = clockangle + right_ascension - ((time / day) * 2 * m.pi)  # rotation of the Earth beneath a prograde orbit
    
    coordinate = Point(latitude, longitude)
    # visible area prism
    #print(altitudes)
    height = visibleAreaPitchAngle(altitudes)  # function of altitude and observing time
    width = visibleAreaRollAngle(altitudes)  # function of altitude and observing time
    #print(height,width)
    # height = (6/180)*m.pi
    # width = (3.4/180)*m.pi
    # rotation
    direction = m.cos(clockangle) * inclination  # angle of the direction of the orbit, 0 = horizontal
    points = np.empty(4, dtype=Point)
    points.put(0, Point(width, height))
    points.put(1, Point(-width,-height))
    points.put(2, Point(-width,height))
    points.put(3, Point(width, -height))
    rotation_matrix = [[m.cos(direction), m.sin(direction)],
                       [-m.sin(direction), m.cos(direction)]]  # CW rotation matrix
    pointsRotated = np.empty(4, dtype=Point)
    for index in range(4):
        rotatedCoords = np.matmul(rotation_matrix, points[index].asArray())
        pointsRotated.put(index, Point(rotatedCoords[0], rotatedCoords[1]))
    coords = np.empty(4, dtype=Point)
    for index in range(4):
        coords.put(index, Point(pointsRotated[index].latitude + coordinate.latitude,
                                pointsRotated[index].longitude + coordinate.longitude))
        #print(coords[index].latitude)
        #overflowing
        #if(coordspindex.latitude > 0.5*m.pi):
        # if (coords[index].latitude > 0.5 * m.pi):
        #     coords[index].latitude = m.pi - latitude
        #     coords[index].longitude += m.pi
        # if (coords[index].latitude < -0.5 * m.pi):
        #     coords[index].latitude = -(m.pi + latitude)
        #     coords[index].longitude += m.pi
        
    return coords

def isInRectangle(vertices, originalPoint):
    point = sg.Point(originalPoint.latitude, originalPoint.longitude)
    polygon = Polygon([(vertices[0].latitude, vertices[0].longitude),
                      (vertices[2].latitude, vertices[2].longitude),
                      (vertices[1].latitude, vertices[1].longitude),
                       (vertices[3].latitude, vertices[3].longitude)])
    return polygon.contains(point)

def calculateObservableArea(visibleRegionInput):
    vertical = m.pi
    horizontal = 2 * m.pi
    n_vertical = int(vertical / res)
    n_horizontal = int(horizontal / res)
    
    mincol = n_horizontal
    maxcol = 0
    minrow = n_vertical
    maxrow = 0
    
    observableAreaMatrix = np.zeros([n_vertical, n_horizontal])
    matrixPosition = []

    for vertex in visibleRegionInput:
        latitude = (vertex.latitude * -1) + (0.5*m.pi)
        matrixPosition.append(Point(int((latitude)/res), int(((vertex.longitude % (2 * m.pi))/res))))
    for position in matrixPosition:
        if(position.latitude < minrow):
            minrow = position.latitude
        if(position.latitude > maxrow):
            maxrow = position.latitude
        if(position.longitude < mincol):
            mincol = position.longitude
        if(position.longitude > maxcol):
            maxcol = position.longitude
    rowheight = maxrow - minrow
    colwidth = maxcol - mincol
    
    if(colwidth > 0.5*observableAreaMatrix.shape[1]):
        mincol = 2*n_horizontal
        maxcol = 0
        for position in matrixPosition:
            if(position.longitude < 0.5*observableAreaMatrix.shape[1]):
                position.longitude += observableAreaMatrix.shape[1]
            if(position.longitude < mincol):
                mincol = position.longitude
            if(position.longitude > maxcol):
                maxcol = position.longitude
        colwidth = maxcol - mincol
    for i in range(colwidth):
        for j in range(rowheight):
            row = int(j + minrow)
            col = int(i + mincol)
            if isInRectangle(matrixPosition, Point(row, col)):
                if(row < 0):
                    row = int(-row - 1)
                    col = int(col + 0.5*n_horizontal)
                if(row > (observableAreaMatrix.shape[0]-1)):
                    row = int(observableAreaMatrix.shape[0] - row%observableAreaMatrix.shape[0] - 1)
                    col = int(col + 0.5*n_horizontal)
                if(col > (observableAreaMatrix.shape[1]-1)):
                    col = col%observableAreaMatrix.shape[1]
                observableAreaMatrix[row][col] = 1
                    
    counter = 0
    points = np.empty([20000,2])
    for i in range(observableAreaMatrix.shape[0]):
        for j in range(observableAreaMatrix.shape[1]):
            if observableAreaMatrix[i, j] != 0:
                pointFirst = Point(((0.5 * m.pi) - ((m.pi / observableAreaMatrix.shape[0]) * i)), (((2 * m.pi) / observableAreaMatrix.shape[1]) * j))
                pointProjected = STXY(pointFirst.latitude, pointFirst.longitude)
                points[counter,0] = pointProjected.latitude
                points[counter,1] = pointProjected.longitude
                counter += 1
    #print(counter)
    plt.scatter(points[:counter,0],points[:counter,1],color='green', s = 0.5)
                

    #np.savetxt("observableAreaMatrix.csv", observableAreaMatrix, delimiter=',')

    return observableAreaMatrix

def ellipse(inclination, right_ascension):
    n = 4
    pos = np.empty([n, 2])
    width = 1
    height = 1
    for i in range(n):
        theta = i * ((2 * m.pi) / n)
        pos[i] = [m.sin(theta) * height + inclination, m.cos(theta) * width + right_ascension]
    return pos

def dayNightMatrix(time):
    #Sunposition
    Earth_tilt = (23.5/180)*m.pi #Earth's tilt in radians
    solar_inclination = Earth_tilt*m.sin((time/(day*365.25))*2*m.pi) #declination is 0 on March 21st (80th day of the year)
    solar_right_ascension = (time/day)*2*m.pi #right ascension of sun rising
    r = 0.5*m.pi
    #Matrix
    DayNightMatrix = np.zeros([n_vertical, n_horizontal])
    n = int(2*r/res)
    start_ra = solar_right_ascension%(2*m.pi) - r
    start_column = int(start_ra/res)
    for i in range(n): #for each column in the shape
        col = i + start_column
        colsfromcenter = int(col - start_column - 0.5*n)
        distancefromcenter = colsfromcenter*res/(0.5*m.pi)
        theta = np.arccos(distancefromcenter)
        centerrowheight = int((-solar_inclination + 0.5*m.pi)/res)
        absoluterowsheight = int((m.sin(theta)*0.5*m.pi)/res)
        for j in range(2*absoluterowsheight):
            row = int(centerrowheight - absoluterowsheight + j)
            #Overflowing
            if(row >= DayNightMatrix.shape[0]):
                row = int(DayNightMatrix.shape[0] - 1 - (row%DayNightMatrix.shape[0]))
                col = int(col+n)
            if(row < 0):
                row = int(-row-1)
                col = int(col+n)
            if(col < 0):
                col = int(DayNightMatrix.shape[1] + col + 1)
            if(col >= DayNightMatrix.shape[1]-1):
                col = int(col%DayNightMatrix.shape[1])
            DayNightMatrix[row,col] = 1      
            col = i + start_column 
    return DayNightMatrix

def dataMatrix(daynightmatrix, cloudcoveragematrix):
    DataMatrix = daynightmatrix*cloudcoveragematrix
    return DataMatrix

def requirementMatrix(DayNightMatrix, ObservableAreaMatrix):
    RequirementMatrix = DayNightMatrix*ObservableAreaMatrix
    return RequirementMatrix

# Matrices importing
AOI_matrix = np.loadtxt("AOI_matrix.csv", delimiter=',', dtype=int)
CLOUD_matrix = ((np.delete((np.loadtxt("clouds_integer_yearly.csv", delimiter=',', skiprows = 1, dtype=str)), 0, 1)).astype(float)).astype(int)
#DATA_matrix = AOI_matrix*CLOUD_matrix
#print(DATA_matrix, DATA_matrix.shape)


counter = 0
for i in range(AOI_matrix.shape[0]):
    for j in range(AOI_matrix.shape[1]):
        if (AOI_matrix[i, j] != 0):
            point = Point((0.5 * m.pi - (m.pi / AOI_matrix.shape[0] * i)), (2 * m.pi / AOI_matrix.shape[1] * j))
            point_proj = STXY(point.latitude, point.longitude)
            # plt.scatter(point_proj[0],point_proj[1],color='blue')
            counter += 1
print(counter)

#Matrices importing
# AOI_matrix = np.loadtxt("AOI_matrix_corrected.csv", delimiter=',', dtype=int)
# cloudcoverage_matrix = np.loadtxt("cloudcoverage_map_weighted.csv", delimiter=',') #skiprows=1
# cloudcoverage_matrix = cloudcoverage_matrix[:,1:]

#inputs for multiple satellites
inclinations = [(50/180)*m.pi, (10/180)*m.pi] #radians
right_ascensions = [(20/180)*m.pi, (40/180)*m.pi] #radians
altitudes = 700*10**3 #meters
timelength = timestamp_length #seconds
res = (0.1/180)*m.pi #radians height and width per image

a = 6371000 + altitudes  # meters
period = int(calcPeriod(a))  # seconds
n = int(timelength / timestamp_length)
n_horizontal = int(horizontal/res)
n_vertical= int(vertical/res)
orbits = timelength / period
days = timelength / day
numSatellites = len(inclinations)

# figure
image = plt.imread("Gall–Peters_projection.jpg")
fig = plt.figure()
fig, ax = plt.subplots()
# for i in range(int(days)):
#     string = str(i)
#     string = image
#     string = ax.imshow(image, extent = [i*2*m.pi,(1+i)*2*m.pi,-2,2])
img = image
img = ax.imshow(image, extent=[0, 2 * m.pi, -2, 2])
plt.xlim(0, 2 * m.pi)
plt.ylim(-2, 2)

# Main loop

#positions_sun = np.empty(n, dtype=Point)

for i in range(n):
    visibleAreaPitchAngle(450000)
    visibleAreaRollAngle(450000)
    time = i * timestamp_length + 2601.06 * day
    
    #DayNightMatrix
    counter = 0
    DayNightMatrix = dayNightMatrix(time)
    points = np.empty([DayNightMatrix.shape[0]*DayNightMatrix.shape[1],2])
    for p in range(DayNightMatrix.shape[0]):
        for j in range(DayNightMatrix.shape[1]):
            if(DayNightMatrix[p,j] != 0):
                point = [(0.5*m.pi -(m.pi/DayNightMatrix.shape[0] * p)), (2*m.pi/DayNightMatrix.shape[1] * j)]
                point_proj = STXY(point[0],point[1])
                points[counter,0] = point_proj.latitude
                points[counter,1] = point_proj.longitude       
                counter += 1
    plt.scatter(points[:counter,0],points[:counter,1],color='yellow', s = 0.5)
    
    for j in range(numSatellites):
        print("main loop")
        # Satellite
        point_spherical = orbit_spherical(inclinations[j], right_ascensions[j], time)
        point_projection = STXY(point_spherical[0], point_spherical[1])
        plt.scatter(point_projection.latitude, point_projection.longitude, color='red', s=1)

        #Visible area
        visibleRegion = satellitePrism(inclinations[j], right_ascensions[j], time)
        ObservableAreaMatrix = calculateObservableArea(visibleRegion)
        projections = np.empty(4, dtype=Point)
        OAprojections = np.empty(4, dtype=Point)
        for index in range(4):
            #overflowing
            if (visibleRegion[index].latitude > 0.5 * m.pi):
                visibleRegion[index].latitude = m.pi - visibleRegion[index].latitude
                visibleRegion[index].longitude += m.pi
            if (visibleRegion[index].latitude < -0.5 * m.pi):
                visibleRegion[index].latitude = -(m.pi + visibleRegion[index].latitude)
                visibleRegion[index].longitude += m.pi
            projections.put(index, STXY(visibleRegion[index].latitude, visibleRegion[index].longitude))
        # for projection in projections:
        #     plt.scatter(projection.latitude, projection.longitude, color='red', s=1)
        
        #Requirement matrix
        counter = 0
        RequirementMatrix = requirementMatrix(DayNightMatrix, ObservableAreaMatrix)
        points = np.empty([DayNightMatrix.shape[0]*DayNightMatrix.shape[1],2])
        for q in range(RequirementMatrix.shape[0]):
            for t in range(RequirementMatrix.shape[1]):
                if(RequirementMatrix[q,t] != 0):
                    point = [(0.5*m.pi -(m.pi/RequirementMatrix.shape[0] * q)), (2*m.pi/RequirementMatrix.shape[1] * t)]
                    point_proj = STXY(point[0],point[1])
                    points[counter,0] = point_proj.latitude
                    points[counter,1] = point_proj.longitude       
                    counter += 1
        plt.scatter(points[:counter,0],points[:counter,1],color='blue', s = 0.5)
                    
    
#AOI_matrix
# counter = 0
# points = np.empty([AOI_matrix.shape[0]*AOI_matrix.shape[1],2])
# for i in range(AOI_matrix.shape[0]):
#     for j in range(AOI_matrix.shape[1]):
#         if(AOI_matrix[i,j] != 0):
#             point = [(0.5*m.pi -(m.pi/AOI_matrix.shape[0] * i)), ((2*m.pi/AOI_matrix.shape[1] * j) + m.pi)]
#             point_proj = STXY(point[0],point[1])
#             points[counter,0] = point_proj.latitude
#             points[counter,1] = point_proj.longitude
#             counter += 1
# plt.scatter(points[:,0], points[:,1], color='blue', s = 0.01)

#Cloud coverage matrix
# projection_red = np.empty([60000,2])
# projection_orange = np.empty([5400000,2])
# projection_yellow = np.empty([6200000,2])
# projection_green = np.empty([2000000,2])
# projection_blue = np.empty([3030000,2])
# counter_red = 0
# counter_orange = 0
# counter_yellow = 0
# counter_green = 0
# counter_blue = 0
# for i in range(cloudcoverage_matrix.shape[0]):
#     for j in range(cloudcoverage_matrix.shape[1]):
#         if(cloudcoverage_matrix[i,j] != 0):
#             if(cloudcoverage_matrix[i,j] == 1):
#                 position = [(0.5*m.pi -(m.pi/cloudcoverage_matrix.shape[0] * i)), (2*m.pi/cloudcoverage_matrix.shape[1] * j)]
#                 projection = STXY(position[0],position[1])
#                 projection_red[counter_red,0] = projection.latitude
#                 projection_red[counter_red,1] = projection.longitude
#                 counter_red += 1
#             elif(cloudcoverage_matrix[i,j] <= 2):
#                 position = [(0.5*m.pi -(m.pi/cloudcoverage_matrix.shape[0] * i)), (2*m.pi/cloudcoverage_matrix.shape[1] * j)]
#                 projection = STXY(position[0],position[1])
#                 projection_orange[counter_orange,0] = projection.latitude
#                 projection_orange[counter_orange,1] = projection.longitude
#                 counter_orange += 1
#             elif(cloudcoverage_matrix[i,j] <= 3):
#                 position = [(0.5*m.pi -(m.pi/cloudcoverage_matrix.shape[0] * i)), (2*m.pi/cloudcoverage_matrix.shape[1] * j)]
#                 projection = STXY(position[0],position[1])
#                 projection_yellow[counter_yellow,0] = projection.latitude
#                 projection_yellow[counter_yellow,1] = projection.longitude
#                 counter_yellow += 1
#             elif(cloudcoverage_matrix[i,j] <= 4):
#                 position = [(0.5*m.pi -(m.pi/cloudcoverage_matrix.shape[0] * i)), (2*m.pi/cloudcoverage_matrix.shape[1] * j)]
#                 projection = STXY(position[0],position[1])
#                 projection_green[counter_green,0] = projection.latitude
#                 projection_green[counter_green,1] = projection.longitude
#                 counter_green += 1
#             else:
#                 position = [(0.5*m.pi -(m.pi/cloudcoverage_matrix.shape[0] * i)), (2*m.pi/cloudcoverage_matrix.shape[1] * j)]
#                 projection = STXY(position[0],position[1])
#                 projection_blue[counter_blue,0] = projection.latitude
#                 projection_blue[counter_blue,1] = projection.longitude
#                 counter_blue += 1
# #print(counter_red,counter_orange,counter_yellow,counter_green,counter_blue)
# plt.scatter(projection_red[:counter_red,0], projection_red[:counter_red,1],color='red', s = 0.1, marker='.', linewidths=0)
# plt.scatter(projection_orange[:counter_orange,0], projection_orange[:counter_orange,1],color='orange', s = 0.1, marker='.', linewidths=0)
# plt.scatter(projection_yellow[:counter_yellow,0], projection_yellow[:counter_yellow,1],color='yellow', s = 0.1, marker='.', linewidths=0)
# plt.scatter(projection_green[:counter_green,0], projection_green[:counter_green,1],color='green', s = 0.1, marker='.', linewidths=0)
# plt.scatter(projection_blue[:counter_blue,0], projection_blue[:counter_blue,1],color='blue', s = 0.1, marker='.', linewidths=0)
  
#Datamatrix
# counter = 0
# Data_matrix = dataMatrix(AOI_matrix, cloudcoverage_matrix)
# points = []

# interpolateGrayscape = interp1d([30, 0], [0, 255])
# for i in range(Data_matrix.shape[0]):
#         for j in range(Data_matrix.shape[1]):
#             if(Data_matrix[i,j] == 0 or Data_matrix[i,j] > 30):
#                 continue
#             position = [(0.5*m.pi -(m.pi/Data_matrix.shape[0] * i)), (2*m.pi/Data_matrix.shape[1] * j + m.pi)]
#             projection = STXY(position[0],position[1])
#             points.append([projection.latitude, projection.longitude, interpolateGrayscape(Data_matrix[i,j])])
#             #print(projection.latitude,projection.longitude)
#             counter += 1
# points = np.array(points)
# plt.scatter(points[:,0], points[:,1], c = points[:,2], cmap='inferno', s=0.2, marker=',', linewidths=0)
# #print(counter)

print("Saving")
# fig.savefig("Datamatrix.jpg", dpi=2000)   


# dataframe = pd.DataFrame(Data_matrix)
# dataframe.to_csv("Data_matrix.csv", header = False)

plt.show()
print("Done")