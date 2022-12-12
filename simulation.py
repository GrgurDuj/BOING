import math as m
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import shapely.geometry as sg
from shapely.geometry.polygon import Polygon

from Point import Point


# import astropy.constants as c
# np.set_printoptions(threshold=np.inf)
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


def calcPeriod(a):
    mu = 3.986E14
    period = 2 * m.pi * m.sqrt((a ** 3) / mu)
    return period


def visibleAreaPitchAngle(altitude):
    minPitch = 0.7  # deg
    maxPitch = 6.08  # deg
    minAltitude = 450000  # m
    maxAltitude = 700000  # m
    interpolateLength = interp1d([minAltitude, maxAltitude], [minPitch, maxPitch])
    print(interpolateLength(altitude), "pitch")
    return interpolateLength(altitude)


def visibleAreaRollAngle(altitude):
    # corresponds to 30deg of roll from sat perspective (15 each side)
    roll = 0.0048310913 * (altitude / 1000) + 0.00623608
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


def isInRectangle(vertices, originalPoint):
    point = sg.Point(originalPoint.latitude, originalPoint.longitude)
    polygon = Polygon([(vertices[0].latitude, vertices[0].longitude),
                       (vertices[1].latitude, vertices[1].longitude),
                       (vertices[2].latitude, vertices[2].longitude),
                       (vertices[3].latitude, vertices[3].longitude)])
    return polygon.contains(point)


def calculateObservableArea(visibleRegionInput, res):
    xycoords = np.empty(4, dtype=Point)
    numberOfVertices = 0
    vertical = m.pi
    horizontal = 2 * m.pi
    n_vertical = int(vertical / res)
    n_horizontal = int(horizontal / res)
    observableAreaMatrix = np.zeros([n_vertical, n_horizontal])

    for point in visibleRegionInput:
        xycoords.put(numberOfVertices, STXY(point.latitude, point.longitude))
        numberOfVertices += 1

    for point in xycoords:
        print(str(point.latitude) + ", " + str(point.longitude))
        plt.scatter(point.latitude, point.longitude, color='blue', s=3)

    matrixPosition = []
    for point in xycoords:
        matrixPosition.append(Point(m.degrees(point.latitude), m.degrees(point.longitude)))

    for position in matrixPosition:
        observableAreaMatrix[int(position.latitude)][int(position.longitude)] = 1

    for i in range(n_horizontal):
        for j in range(n_vertical):
            if isInRectangle(matrixPosition, Point(i, j)):
                observableAreaMatrix[i][j] = observableAreaMatrix[i][j] + 2

    np.savetxt("observableAreaMatrix.csv", observableAreaMatrix, delimiter=',')
    '''
    poly = []
    for index in range(numberOfVertices):
        poly.append([xycoords[index].latitude, xycoords[index].longitude])

    observableAreaImage = Image.fromarray(observableAreaMatrix)
    draw = ImageDraw.Draw(observableAreaImage)
    draw.polygon([tuple(p) for p in poly], fill=1)
    observableAreaMatrix = np.asarray(observableAreaImage)
    plt.imshow(observableAreaMatrix)
    [(0.5 * m.pi - (m.pi / DayNightMatrix.shape[0] * i)), (2 * m.pi / DayNightMatrix.shape[1] * j)]
    np.savetxt("observableAreaMatrix.csv", observableAreaMatrix, delimiter=',')
    '''

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


def DayNightMatrix(time):
    # Sunposition
    Earth_tilt = (23.5 / 180) * m.pi  # Earth's tilt in radians
    solar_inclination = Earth_tilt * m.sin(
        (time / (day * 365.25)) * 2 * m.pi)  # declination is 0 on March 21st (80th day of the year)
    solar_right_ascension = (time / day) * 2 * m.pi  # right ascension of sun rising
    r = 0.5 * m.pi
    # Matrix
    res = (0.5 / 180) * m.pi
    n_vertical = int(vertical / res)
    n_horizontal = int(horizontal / res)
    DayNightMatrix = np.zeros([n_vertical, n_horizontal])
    n = int(2 * r / res)
    start_ra = solar_right_ascension % (2 * m.pi) - r
    start_column = int(start_ra / res)
    for i in range(n):  # for each column in the shape
        col = i + start_column
        colsfromcenter = int(col - start_column - 0.5 * n)
        distancefromcenter = colsfromcenter * res / (0.5 * m.pi)
        theta = np.arccos(distancefromcenter)
        centerrowheight = int((-solar_inclination + 0.5 * m.pi) / res)
        absoluterowsheight = int((m.sin(theta) * 0.5 * m.pi) / res)
        for j in range(2 * absoluterowsheight):
            row = int(centerrowheight - absoluterowsheight + j)
            # Overflowing
            if (row >= DayNightMatrix.shape[0]):
                row = int(DayNightMatrix.shape[0] - 1 - (row % DayNightMatrix.shape[0]))
                col = int(col + n)
            if (row < 0):
                row = int(-row - 1)
                col = int(col + n)
            if (col >= DayNightMatrix.shape[1]):
                col = int(col % DayNightMatrix.shape[1])
            if (col < 0):
                col = int(DayNightMatrix.shape[1] + col + 1)
            DayNightMatrix[row, col] = 1
            col = i + start_column
    return DayNightMatrix


# parameters
timestamp_length = 180  # seconds
day = 86164  # seconds, 365.25/366.25 * 24 * 3600, seconds in a sidereal day
r_earth = 6371 * 1000  # meters
measure_time = 120  # seconds
vertical = m.pi  # radians
horizontal = 2 * m.pi  # radians

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

# inputs for multiple satellites
inclinations = [(20 / 180) * m.pi]  # radians
right_ascensions = [(20 / 180) * m.pi]  # radians
altitudes = 700 * 10 ** 3  # meters
timelength = 1 * timestamp_length  # seconds
res = (3 / 180) * m.pi  # radians height and width per image

a = 6371000 + altitudes  # meters
period = int(calcPeriod(a))  # seconds
n = int(timelength / timestamp_length)
orbits = timelength / period
days = timelength / day

a = 6371000 + altitudes  # meters
period = int(calcPeriod(a))  # seconds
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
    time = i * timestamp_length + 250.52 * day
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
        observableArea = calculateObservableArea(visibleRegion, (0.5 / 180) * m.pi)

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
    # DayNightMatrix
    DayNightMatrix = DayNightMatrix(time)
    counter = 0
    for i in range(DayNightMatrix.shape[0]):
        for j in range(DayNightMatrix.shape[1]):
            if (DayNightMatrix[i, j] != 0):
                point = [(0.5 * m.pi - (m.pi / DayNightMatrix.shape[0] * i)), (2 * m.pi / DayNightMatrix.shape[1] * j)]
                point_proj = STXY(point[0], point[1])
                # plt.scatter(point_proj[0],point_proj[1],color='blue')
                counter += 1
    print(counter)

    # AOI_matrix
    # counter = 0
    # for i in range(AOI_matrix.shape[0]):
    #     for j in range(AOI_matrix.shape[1]):
    #         if(AOI_matrix[i,j] != 0):
    #             point = [(0.5*m.pi -(m.pi/DayNightMatrix.shape[0] * i)), (2*m.pi/DayNightMatrix.shape[1] * j)]
    #             point_proj = STXY(point[0],point[1])
    #             plt.scatter(point_proj[0],point_proj[1],color='blue')
    #             counter += 1

# AOI matrix
# vertical = m.pi
# horizontal = 2*m.pi
# res = (5/180)*m.pi
# n_vertical = int(vertical/res)
# n_horizontal = int(horizontal/res)
# AOI_input = np.ones([n_vertical, n_horizontal])
# counter = 0

vertical = m.pi
horizontal = 2 * m.pi
res = (0.5 / 180) * m.pi
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

# for i in range(counter):
#     point = STXY(AOI_interest[i,0], AOI_interest[i,1])
#     plt.scatter(point[0], point[1], color = 'blue', s = 1)
#fig.savefig("DayNightMatrix.jpg", dpi=1000)
plt.show()
#DayNightMatrix = DayNightMatrix(time)
#dataframe = pd.DataFrame(DayNightMatrix)
#dataframe.to_csv("DayNightMatrix.csv")
print("Done")
'''
