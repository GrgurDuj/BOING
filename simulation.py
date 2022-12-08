import math as m
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from Point import Point


# import astropy.constants as c
def orbit_cartesian(orbit_inclination, orbit_right_ascension, time):  # time is in radians
    x = m.sin(time + orbit_right_ascension) * m.cos(orbit_inclination)
    y = m.sin(time - orbit_right_ascension) * m.sin(orbit_inclination)
    z = m.cos(time + orbit_right_ascension)
    coordinate = [x, y, z]
    return coordinate


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


def circle(inclination, right_ascension):
    r = 0.5 * m.pi
    n = 4
    pos = np.empty([n, 2])
    for i in range(n):
        theta = i * ((2 * m.pi) / n)
        pos[i] = [m.sin(theta) * r + inclination, m.cos(theta) * r + right_ascension]
    return pos


def period(a):
    mu = 3.986E14
    period = 2 * m.pi * m.sqrt((a ** 3) / mu)
    return period

def visibleAreaPitchAngle(altitude):
    minPitch = 0.7  #deg
    maxPitch = 6.08 #deg
    minAltitude = 450000 #m
    maxAltitude = 700000 #m
    interpolateLength = interp1d([minAltitude, maxAltitude], [minPitch, maxPitch])
    print(interpolateLength(altitude), "pitch")
    return interpolateLength(altitude)

def visibleAreaRollAngle(altitude):
    # corresponds to 30deg of roll from sat perspective (15 each side)
    roll = 0.0048310913 * (altitude/1000) + 0.00623608
    print(roll, "roll")
    # roll gives observable angle with respect to the center of earth
    return roll

def satellitePrism(inclination, right_ascension, time):
    # satellite position
    clockangle = (time / period) * 2 * m.pi
    latitude = m.sin(clockangle) * inclination
    longitude = clockangle + right_ascension - (
            (time / day) * 2 * m.pi)  # rotation of the Earth beneath a prograde orbit
    if latitude > 0.5 * m.pi:  # if the inclination overflows 0.5*pi, the angle is reversed and rotated to the other side of the globe
        latitude = m.pi - latitude
        longitude += m.pi
    if latitude < -0.5 * m.pi:
        latitude = -(m.pi + latitude)
        longitude += m.pi

    coordinate = Point(latitude, longitude)

    # visible area prism
    height = (6.08 / 180) * m.pi  # function of altitude and observing time
    width = (30 / 180) * m.pi  # function of altitude and observing time

    # rotation
    direction = m.cos(clockangle) * inclination  # angle of the direction of the orbit, 0 = horizontal
    '''
    point1 = [width,height]
    point2 = [(2*m.pi*(measure_time/period) - width),-height]
    point3 = [(2*m.pi*(measure_time/period) - width),height]
    point4 = [width,-height]
    '''
    points = np.empty(4, dtype=Point)
    points.put(0, Point(width, height))
    points.put(1, Point(2 * m.pi * (measure_time / period) - width, -height))
    points.put(2, Point(2 * m.pi * (measure_time / period) - width, height))
    points.put(3, Point(width, -height))
    rotation_matrix = [[m.cos(direction), m.sin(direction)],
                       [-m.sin(direction), m.cos(direction)]]  # CW rotation matrix
    pointsRotated = np.empty(4, dtype=Point)
    for index in range(4):
        rotatedCoords = np.matmul(rotation_matrix, points[index].asArray())
        pointsRotated.put(index, Point(rotatedCoords[0], rotatedCoords[1]))
    '''    
    pointsRotated.put(0, np.matmul(rotation_matrix, points[0]))
    pointsRotated.put(1, np.matmul(rotation_matrix, points[1]))
    pointsRotated.put(2, np.matmul(rotation_matrix, points[2]))
    pointsRotated.put(3, np.matmul(rotation_matrix, points[3]))
    '''
    coords = np.empty(4, dtype=Point)
    for index in range(4):
        coords.put(index, Point(pointsRotated[index].latitude + coordinate.latitude,
                                pointsRotated[index].longitude + coordinate.longitude))
    '''
    coord1 = [pointsRotated[0].latitude + coordinate.latitude, pointsRotated[0].longitude + coordinate.longitude]
    coord2 = [pointsRotated[1].latitude + coordinate.latitude, pointsRotated[1].longitude + coordinate.longitude]
    coord3 = [pointsRotated[2].latitude + coordinate.latitude, pointsRotated[2].longitude + coordinate.longitude]
    coord4 = [pointsRotated[3].latitude + coordinate.latitude, pointsRotated[3].longitude + coordinate.longitude]
    coords = [coord1,coord2,coord3,coord4]
    '''

    # overflowing
    for i in range(len(coords)):
        if (coords[i].latitude > 0.5 * m.pi):
            coords[i].latitude = m.pi - coords[i].latitude
            coords[i].longitude += m.pi
        elif (coords[i].latitude < -0.5 * m.pi):
            coords[i].latitude = -(m.pi + coords[i].latitude)
            coords[i].longitude += m.pi
    return coords


def calculateObservableArea(visibleRegionInput, res):
    xycoords = np.empty(4, dtype=Point)
    count = 0
    for point in visibleRegionInput:
        xycoords.put(count, STXY(point.latitude, point.longitude))
        count += 1

    for point2 in xycoords:
        print(str(point2.latitude) + ", " + str(point2.longitude))
        plt.scatter(point2.latitude, point2.longitude, color='blue', s=3)

    # dayNightMatrix = np.zeros([n_vertical, n_horizontal])
    #for column in range(n_vertical):




def ellipse(inclination, right_ascension):
    n = 4
    pos = np.empty([n, 2])
    width = 1
    height = 1
    for i in range(n):
        theta = i * ((2 * m.pi) / n)
        pos[i] = [m.sin(theta) * height + inclination, m.cos(theta) * width + right_ascension]
    return pos


# parameters
timestamp_length = 180  # seconds
day = 86164  # seconds, 365.25/366.25 * 24 * 3600, seconds in a sidereal day
r_earth = 6371 * 1000  # meters
measure_time = 120  # seconds

# inputs for multiple satellites
inclinations = [(20 / 180) * m.pi]  # radians
right_ascensions = [(20 / 180) * m.pi]  # radians
altitudes = 500 * 10 ** 3  # meters
timelength = 5 * timestamp_length  # seconds

a = 6371000 + altitudes  # meters
period = int(period(a))  # seconds
n = int(timelength / timestamp_length)
orbits = timelength / period
days = timelength / day

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
positions = np.empty(n, dtype=Point)
positions_sun = np.empty(n, dtype=Point)
for i in range(n):
    visibleAreaPitchAngle(450000)
    visibleAreaRollAngle(450000)
    time = i * timestamp_length + 2505.82 * day
    for j in range(len(inclinations)):
        # Satellite
        point_spherical = orbit_spherical(inclinations[j], right_ascensions[j], time)
        point_projection = STXY(point_spherical[0], point_spherical[1])
        positions[j] = point_projection
        plt.scatter(positions[j].latitude, positions[j].longitude, color='red', s=1)

        # Visible area
        visibleRegion = satellitePrism(inclinations[j], right_ascensions[j], time)
        observableArea = calculateObservableArea(visibleRegion, (5 / 180) * m.pi)

        '''
        pos1 = visibleRegion[0]
        pos2 = visibleRegion[1]
        pos3 = visibleRegion[2]
        pos4 = visibleRegion[3]
        '''
        projections = np.empty(4, dtype=Point)
        OAprojections = np.empty(4, dtype=Point)
        for index in range(4):
            projections.put(index, STXY(visibleRegion[index].latitude, visibleRegion[index].longitude))
        '''
        pos1_proj = STXY(pos1.latitude,pos1.longitude)
        pos2_proj = STXY(pos2.latitude,pos2.longitude)
        pos3_proj = STXY(pos3.latitude,pos3.longitude)
        pos4_proj = STXY(pos4.latitude,pos4.longitude)
        '''
        for projection in projections:
            plt.scatter(projection.latitude, projection.longitude, color='red', s=1)
        '''
        plt.scatter(pos1_proj[0],pos1_proj[1], color='red', s=1)
        plt.scatter(pos2_proj[0],pos2_proj[1], color='red', s=1)
        plt.scatter(pos3_proj[0],pos3_proj[1], color='red', s=1)
        plt.scatter(pos4_proj[0],pos4_proj[1], color='red', s=1)
        '''
    # Sun
    point_sun = SunPosition(time)
    point_sun_projection = STXY(point_sun.latitude, point_sun.longitude)
    positions_sun[i] = point_sun_projection
    plt.scatter(positions_sun[i].latitude, positions_sun[i].longitude, color='yellow', s=5)

    # Illuminated area
    illuminatedRegion = sunCircle(time)
    for i in range(illuminatedRegion.shape[0]):
        point = STXY(illuminatedRegion[i].latitude, illuminatedRegion[i].longitude)
        plt.scatter(point.latitude, point.longitude, color='yellow', s=5)

# AOI matrix

vertical = m.pi
horizontal = 2 * m.pi
res = (5 / 180) * m.pi
n_vertical = int(vertical / res)
n_horizontal = int(horizontal / res)
AOI_input = np.ones([n_vertical, n_horizontal])
counter = 0
'''
for i in range(AOI_input.shape[0]):
    for j in range(AOI_input.shape[1]):
        if(AOI_input[i,j] != 0):
            counter +=1
AOI_interest = np.empty(counter, dtype=Point).fill(Point(0, 0))
print(counter)
counter = 0
for i in range(AOI_input.shape[0]):
    for j in range(AOI_input.shape[1]):
        if(AOI_input[i,j] != 0):
            AOI_interest[counter].latitude = (0.5*m.pi -(m.pi/AOI_input.shape[0] * i)) #latitude
            AOI_interest[counter].longitude = (2*m.pi/AOI_input.shape[1] * j) #longitude
            counter += 1

for i in range(counter):
    point = STXY(AOI_interest[i].latitude, AOI_interest[i].longitude)
    plt.scatter(point.latitude, point.longitude, color = 'blue', s = 1)
fig.savefig("test_square_observable.jpg", dpi=10)
'''

plt.show()
print("Done")
