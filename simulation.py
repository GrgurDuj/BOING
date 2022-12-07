import math as m
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import astropy.constants as c
#np.set_printoptions(threshold=np.inf)
def orbit_cartesian(orbit_inclination, orbit_right_ascension, time): #time is in radians
    x = m.sin(time+orbit_right_ascension)*m.cos(orbit_inclination)
    y = m.sin(time-orbit_right_ascension)*m.sin(orbit_inclination)
    z = m.cos(time+orbit_right_ascension)
    coordinate = [x,y,z]
    return coordinate

def orbit_spherical(inclination, right_ascension, time): #time is in seconds
    clockangle = (time/period)*2*m.pi
    latitude = m.sin(clockangle)*inclination
    longitude = clockangle+right_ascension - ((time/day)*2*m.pi) #rotation of the Earth beneath a prograde orbit
    if(latitude > 0.5*m.pi):
        latitude = m.pi - latitude
        longitude += m.pi
    if(latitude < -0.5*m.pi):
        latitude = -(m.pi + latitude)
        longitude += m.pi
    coordinate = [latitude, longitude]
    return coordinate

def STXY(latitude, longitude): #lo = longitude is perpendicular to equator, la = latitude is parallel to equator
    #Gall-Peters projection
    x = longitude%(2*m.pi) #longitude in radians from central meridian
    y = 2*m.sin(latitude) #latitude in radians from equator
    coordinate = [x,y]
    return coordinate

def SunPosition(time):
    Earth_tilt = (23.5/180)*m.pi #Earth's tilt in radians
    solar_inclination = Earth_tilt*m.sin((time/(day*365.25))*2*m.pi) #declination is 0 on March 21st (80th day of the year)
    solar_right_ascension = (time/day)*2*m.pi #right ascension of sun rising
    position = [solar_inclination, solar_right_ascension]
    return position

def sunCircle(time):
    Earth_tilt = (23.5/180)*m.pi #Earth's tilt in radians
    solar_inclination = Earth_tilt*m.sin((time/(day*365.25))*2*m.pi) #declination is 0 on March 21st (80th day of the year)
    solar_right_ascension = (time/day)*2*m.pi #right ascension of sun rising
    r = 0.5*m.pi
    n = 64
    pos = np.empty([n,2])
    for i in range(n):
        theta = i*((2*m.pi)/n)
        pos[i] = [m.sin(theta)*r + solar_inclination, m.cos(theta)*r + solar_right_ascension] 
        if(pos[i,0] > 0.5*m.pi):
            pos[i,0] = m.pi - pos[i,0]
            pos[i,1] += m.pi
        if(pos[i,0] < -0.5*m.pi):
            pos[i,0] = -(m.pi + pos[i,0])
            pos[i,1] += m.pi
    return pos

def circle(inclination, right_ascension):
    r = 0.5*m.pi
    n = 4
    pos = np.empty([n,2])
    for i in range(n):
        theta = i*((2*m.pi)/n)
        pos[i] = [m.sin(theta)*r + inclination, m.cos(theta)*r + right_ascension]
    return pos

def period(a):
    mu = 3.986E14
    period = 2*m.pi*m.sqrt((a**3)/mu)
    return period

def satellitePrism(inclination, right_ascension, time):
    #satellite position
    clockangle = (time/period)*2*m.pi
    latitude = m.sin(clockangle)*inclination
    longitude = clockangle+right_ascension - ((time/day)*2*m.pi) #rotation of the Earth beneath a prograde orbit
    if(latitude > 0.5*m.pi):    #if the inclination overflows 0.5*pi, the angle is reversed and rotated to the other side of the globe
        latitude = m.pi - latitude
        longitude += m.pi
    if(latitude < -0.5*m.pi):
        latitude = -(m.pi + latitude)
        longitude += m.pi
    coordinate = [latitude, longitude]
    
    #visible area prism    
    width = (6.08/180)*m.pi  #function of altitude and observing time
    height = (30/180)*m.pi     #function of altitude and observing time
    
    #rotation
    direction = m.cos(clockangle)*inclination   #angle of the direction of the orbit, 0 = horizontal
    point1 = [height,width]
    point2 = [(2*m.pi*(measure_time/period) - height),-width]
    point3 = [(2*m.pi*(measure_time/period) - height),width]
    point4 = [height,-width]
    rotation_matrix = [[m.cos(direction), m.sin(direction)],[-m.sin(direction), m.cos(direction)]] #CW rotation matrix
    point1_rotated = np.matmul(rotation_matrix,point1)
    point2_rotated = np.matmul(rotation_matrix,point2)
    point3_rotated = np.matmul(rotation_matrix,point3)
    point4_rotated = np.matmul(rotation_matrix,point4)
    coord1 = [point1_rotated[0] + coordinate[0], point1_rotated[1] + coordinate[1]]
    coord2 = [point2_rotated[0] + coordinate[0], point2_rotated[1] + coordinate[1]]
    coord3 = [point3_rotated[0] + coordinate[0], point3_rotated[1] + coordinate[1]]
    coord4 = [point4_rotated[0] + coordinate[0], point4_rotated[1] + coordinate[1]]
    coords = [coord1,coord2,coord3,coord4]
    
    #overflowing
    for i in range(len(coords)):
        if(coords[i][0] > 0.5*m.pi):
            coords[i][0] = m.pi - coords[i][0]
            coords[i][1] += m.pi
        elif(coords[i][0] < -0.5*m.pi):
            coords[i][0] = -(m.pi + coords[i][0])
            coords[i][1] += m.pi
    return coords

def ellipse(inclination, right_ascension):
    n = 4
    pos = np.empty([n,2])
    width = 1
    height = 1
    for i in range(n):
        theta = i*((2*m.pi)/n)
        pos[i] = [m.sin(theta)*height + inclination,m.cos(theta)*width + right_ascension]
    return pos

def DayNightMatrix(time):
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
            if(col >= DayNightMatrix.shape[1]):
                col = int(col%DayNightMatrix.shape[1])
            if(col < 0):
                col = int(DayNightMatrix.shape[1] + col + 1)
            DayNightMatrix[row,col] = 1          
            col = i + start_column 
    return DayNightMatrix

#parameters
timestamp_length = 180 #seconds
day = 86164 #seconds, 365.25/366.25 * 24 * 3600, seconds in a sidereal day
r_earth = 6371*1000 #meters
measure_time = 120 #seconds
vertical = m.pi #radians
horizontal = 2*m.pi #radians


#inputs for multiple satellites
inclinations = [(20/180)*m.pi] #radians
right_ascensions = [(20/180)*m.pi] #radians
altitudes = 700*10**3 #meters
timelength = 1*timestamp_length #seconds
res = (3/180)*m.pi #radians height and width per image

a = 6371000 + altitudes #meters
period = int(period(a)) #seconds
n = int(timelength/timestamp_length) 
orbits = timelength/period
days = timelength/day
n_vertical = int(vertical/res)
n_horizontal = int(horizontal/res)

#figure
image = plt.imread("Gall–Peters_projection.jpg")
fig = plt.figure()
fig, ax = plt.subplots()
# for i in range(int(days)):
#     string = str(i)
#     string = image
#     string = ax.imshow(image, extent = [i*2*m.pi,(1+i)*2*m.pi,-2,2])
img = image
img = ax.imshow(image, extent = [0,2*m.pi,-2,2])
plt.xlim(0,2*m.pi)
plt.ylim(-2,2) 

#Main loop
positions = np.empty([n,2])
positions_sun = np.empty([n,2])
for i in range(n):
    time = i*timestamp_length + 100.32*day
    for j in range(len(inclinations)):
        #Satellite
        point_spherical = orbit_spherical(inclinations[j],right_ascensions[j],time)
        point_projection = STXY(point_spherical[0],point_spherical[1])
        positions[j] = point_projection
        plt.scatter(positions[j,0], positions[j,1], color = 'red', s = 1)
        
        #Visible area
        visibleRegion = satellitePrism(inclinations[j], right_ascensions[j], time)
        pos1 = visibleRegion[0]
        pos2 = visibleRegion[1]
        pos3 = visibleRegion[2]
        pos4 = visibleRegion[3]
        pos1_proj = STXY(pos1[0],pos1[1])
        pos2_proj = STXY(pos2[0],pos2[1])
        pos3_proj = STXY(pos3[0],pos3[1])
        pos4_proj = STXY(pos4[0],pos4[1])
        plt.scatter(pos1_proj[0],pos1_proj[1], color='red', s=1)
        plt.scatter(pos2_proj[0],pos2_proj[1], color='red', s=1)
        plt.scatter(pos3_proj[0],pos3_proj[1], color='red', s=1)
        plt.scatter(pos4_proj[0],pos4_proj[1], color='red', s=1)
    #Sun
    point_sun = SunPosition(time)
    point_sun_projection = STXY(point_sun[0], point_sun[1])
    positions_sun[i] = point_sun_projection
    plt.scatter(positions_sun[i,0],positions_sun[i,1], color='yellow', s=5)
    
    #Illuminated area
    illuminatedRegion = sunCircle(time)
    for i in range(illuminatedRegion.shape[0]):
        point = STXY(illuminatedRegion[i,0], illuminatedRegion[i,1])
        plt.scatter(point[0], point[1], color='yellow', s=5)
    #DayNightMatrix
    DayNightMatrix = DayNightMatrix(time)
    counter = 0
    for i in range(DayNightMatrix.shape[0]):
        for j in range(DayNightMatrix.shape[1]):
            if(DayNightMatrix[i,j] != 0):
                point = [(0.5*m.pi -(m.pi/DayNightMatrix.shape[0] * i)), (2*m.pi/DayNightMatrix.shape[1] * j)]
                point_proj = STXY(point[0],point[1])
                plt.scatter(point_proj[0],point_proj[1],color='blue')
                counter += 1
    print(counter)
#AOI matrix
# vertical = m.pi
# horizontal = 2*m.pi
# res = (5/180)*m.pi
# n_vertical = int(vertical/res)
# n_horizontal = int(horizontal/res)
# AOI_input = np.ones([n_vertical, n_horizontal])
# counter = 0

# for i in range(AOI_input.shape[0]):
#     for j in range(AOI_input.shape[1]):
#         if(AOI_input[i,j] != 0):
#             counter +=1
# AOI_interest = np.empty([counter,2])
# print(counter)
# counter = 0
# for i in range(AOI_input.shape[0]):
#     for j in range(AOI_input.shape[1]):
#         if(AOI_input[i,j] != 0):
#             AOI_interest[counter,0] = (0.5*m.pi -(m.pi/AOI_input.shape[0] * i)) #latitude
#             AOI_interest[counter,1] = (2*m.pi/AOI_input.shape[1] * j) #longitude
#             counter += 1

# for i in range(counter):
#     point = STXY(AOI_interest[i,0], AOI_interest[i,1])
#     plt.scatter(point[0], point[1], color = 'blue', s = 1)
#fig.savefig("DayNightMatrix.jpg", dpi=1000)
plt.show()
#DayNightMatrix = DayNightMatrix(time)
#dataframe = pd.DataFrame(DayNightMatrix)
#dataframe.to_csv("DayNightMatrix.csv")
print("Done")