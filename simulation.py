import math as m
import matplotlib.pyplot as plt
import numpy as np
# import astropy.constants as c

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
def satelliteEllipse(inclination, right_ascension, time):
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
    
    #visible area ellipse
    #width = (10/180)*m.pi #longitude
    #height = (45/180)*m.pi #latitude
    
    width = (6.08/180)*m.pi  #function of altitude and observing time
    height = (30/180)*m.pi     #function of altitude and observing time
    
    direction = m.cos(clockangle)*inclination   #angle of the direction of the orbit, 0 = horizontal
    x1 = height
    x2 = (2*m.pi*(measure_time/period) - height)
    x3 = x2
    x4 = x1
    y1 = width
    y2 = -width
    y3 = y1
    y4 = y2
    x1rotated = x1*m.cos(direction) + y1*m.sin(direction)      #rotation matrix for 2d, CCW rotation
    y1rotated = -x1*m.sin(direction) + y1*m.cos(direction)
    x2rotated = x2*m.cos(direction) + y2*m.sin(direction)      #rotation matrix for 2d, CCW rotation
    y2rotated = -x2*m.sin(direction) + y2*m.cos(direction)
    x3rotated = x3*m.cos(direction) + y3*m.sin(direction)
    y3rotated = -x3*m.sin(direction) + y3*m.cos(direction)
    x4rotated = x4*m.cos(direction) + y4*m.sin(direction)
    y4rotated = -x4*m.sin(direction) + y4*m.cos(direction)
    coord1 = [x1rotated + coordinate[0], y1rotated + coordinate[1]]
    coord2 = [x2rotated + coordinate[0], y2rotated + coordinate[1]]
    coord3 = [x3rotated + coordinate[0], y3rotated + coordinate[1]]
    coord4 = [x4rotated + coordinate[0], y4rotated + coordinate[1]]
    if(coord1[0] > 0.5*m.pi):
        coord1[0] = m.pi - coord1[0]
        coord1[1] += m.pi
    if(coord1[0] < -0.5*m.pi):
        coord1[0] = -(m.pi + coord1[0])
        coord1[0] += m.pi
    if(coord2[0] > 0.5*m.pi):
        coord2[0] = m.pi - coord2[0]
        coord2[1] += m.pi
    if(coord2[0] < -0.5*m.pi):
        coord2[0] = -(m.pi + coord2[0])
        coord2[0] += m.pi
    if(coord3[0] > 0.5*m.pi):
        coord3[0] = m.pi - coord3[0]
        coord3[1] += m.pi
    if(coord3[0] < -0.5*m.pi):
        coord3[0] = -(m.pi + coord3[0])
        coord3[0] += m.pi
    if(coord4[0] > 0.5*m.pi):
        coord4[0] = m.pi - coord4[0]
        coord4[1] += m.pi
    if(coord4[0] < -0.5*m.pi):
        coord4[0] = -(m.pi + coord4[0])
        coord4[0] += m.pi
    # if(pos[i,0] > 0.5*m.pi):
    #     pos[i,0] = m.pi - pos[i,0]
    #     pos[i,1] += m.pi
    # if(pos[i,0] < -0.5*m.pi):
    #     pos[i,0] = -(m.pi + pos[i,0])
    #     pos[i,1] += m.pi
    return [coord1,coord2,coord3,coord4]
def ellipse(inclination, right_ascension):
    n = 4
    pos = np.empty([n,2])
    width = 1
    height = 1
    for i in range(n):
        theta = i*((2*m.pi)/n)
        pos[i] = [m.sin(theta)*height + inclination,m.cos(theta)*width + right_ascension]
    return pos


#parameters
timestamp_length = 180 #seconds
day = 86164 #seconds, 365.25/366.25 * 24 * 3600, seconds in a sidereal day

#inputs for multiple satellites
inclinations = [(50/180)*m.pi]#, (80/180)*m.pi]
right_ascensions = [(20/180)*m.pi]#, (50/180)*m.pi]
altitudes = 700*10**3 #meters
a = 6371000 + altitudes #meters
period = int(period(a)) #seconds
timelength = 20*timestamp_length
n = int(timelength/timestamp_length)
orbits = timelength/period
days = timelength/day
r_earth = 6371*1000 #meters
measure_time = 120 #seconds

#orbit plot 2d
# n = timelength/timestamp_length
# inclination = (45/180)*m.pi
# right_ascension = (10/180)*m.pi
# positions = np.empty([int(n),3])
# for i in range(int(n)):
#     time = (i*(timestamp_length/period)*2*m.pi)
#     point = orbit_cartesian(inclination, right_ascension, time)
#     positions[i] = point
# plt.figure()
# plt.xlim(-1,1)
# plt.ylim(-1,1)
# plt.plot(positions[:,0], positions[:,1])
# plt.show()

#2d projection
image = plt.imread("Gall–Peters_projection.jpg")
positions = np.empty([n,2])
positions_sun = np.empty([n,2])

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

#main loop
for i in range(n):
    time = i*timestamp_length + 2505.82*day
    for j in range(len(inclinations)):
        point_spherical = orbit_spherical(inclinations[j],right_ascensions[j],time)
        point_projection = STXY(point_spherical[0],point_spherical[1])
        positions[j] = point_projection
        print(positions)
        #plt.plot(positions[:,0], positions[:,1], color = 'red')
        
        visibleRegion = satelliteEllipse(inclinations[j], right_ascensions[j], time)
        pos1 = visibleRegion[0]
        pos2 = visibleRegion[1]
        pos3 = visibleRegion[2]
        pos4 = visibleRegion[3]
        print(pos1, pos2, pos3, pos4)
        pos1_proj = STXY(pos1[0],pos1[1])
        pos2_proj = STXY(pos2[0],pos2[1])
        pos3_proj = STXY(pos3[0],pos3[1])
        pos4_proj = STXY(pos4[0],pos4[1])
        plt.scatter(pos1_proj[0],pos1_proj[1], color='red', s=1)
        plt.scatter(pos2_proj[0],pos2_proj[1], color='red', s=1)
        plt.scatter(pos3_proj[0],pos3_proj[1], color='red', s=1)
        plt.scatter(pos4_proj[0],pos4_proj[1], color='red', s=1)
        # for k in range(visibleRegion.shape[0]):
        #     point = STXY(visibleRegion[k,0], visibleRegion[k,1])
        #     plt.scatter(point[0], point[1], color='red', s = 10)
        
    point_sun = SunPosition(time)
    point_sun_projection = STXY(point_sun[0], point_sun[1])
    positions_sun[i] = point_sun_projection
    #illuminated area
    illuminatedRegion = sunCircle(time)
    for i in range(illuminatedRegion.shape[0]):
        point = STXY(illuminatedRegion[i,0], illuminatedRegion[i,1])
        plt.scatter(point[0], point[1], color='yellow', s = 5)
    
for i in range(illuminatedRegion.shape[0]):
        point = STXY(illuminatedRegion[i,0], illuminatedRegion[i,1])
        plt.scatter(point[0], point[1], color='yellow', s = 5)
    

#plt.scatter(positions[:,0], positions[:,1], color = 'red', s = 1)
#plt.plot(positions[:,0], positions[:,1], color = 'red')
plt.scatter(positions_sun[n-1,0], positions_sun[n-1,1], color='yellow', s=5)

# position_custom = [(51.838035/180)*m.pi, (184.391580/180)*m.pi]
# position_custom = STXY(position_custom[0], position_custom[1])
# plt.scatter(position_custom[0],position_custom[1], color = 'red', s = 0.5)



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
print("Done")
fig.savefig("test_square_observable.jpg", dpi=1000)
plt.show()