import math as m
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import shapely.geometry as sg
from shapely.geometry.polygon import Polygon
from Point import Point, fromArray
import time as tm

# parameters
timestamp_length = 180  # seconds
day = 86164  # seconds, 365.25/366.25 * 24 * 3600, seconds in a sidereal day
r_earth = 6371 * 1000  # meters
measure_time = 120  # seconds
vertical = m.pi  # radians
horizontal = 2 * m.pi  # radians
altitude = 700 * 10 ** 3  # meters
a = 6371000 + altitude  # meters
res = (0.1 / 180) * m.pi
n_horizontal = int(horizontal / res)
n_vertical = int(vertical / res)


def calcPeriod(alt):
    mu = 3.986E14
    return 2 * m.pi * m.sqrt((alt ** 3) / mu)


period = int(calcPeriod(a))  # seconds


def orbit_spherical(inclination, right_ascension, time, altitude):  # time is in seconds
    period = int(calcPeriod(altitude))
    clockangle = (time / period) * 2 * m.pi
    latitude = m.sin(clockangle) * inclination
    longitude = clockangle + right_ascension - (
            (time / day) * 2 * m.pi)  # rotation of the Earth beneath a prograde orbit
    if latitude > 0.5 * m.pi:
        latitude = m.pi - latitude
        longitude += m.pi
    if latitude < -0.5 * m.pi:
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


def visibleAreaPitchAngle(altitude):
    minPitch = (0.7 / 180) * m.pi  # rad
    maxPitch = (6.08 / 180) * m.pi  # rad5
    minAltitude = 450000  # m
    maxAltitude = 700000  # m
    interpolateLength = interp1d([minAltitude, maxAltitude], [minPitch, maxPitch])
    return interpolateLength(altitude)


def visibleAreaRollAngle(altitude):
    # corresponds to 30deg of roll from sat perspective (15 each side)
    roll = ((0.0048310913 * (altitude / 1000) + 0.00623608) / 180) * m.pi  # rad
    # roll gives observable angle with respect to the center of earth
    return roll


def satellitePrism(inclination, right_ascension, time):
    # satellite position
    clockangle = (time / period) * 2 * m.pi
    latitude = m.sin(clockangle) * inclination
    longitude = clockangle + right_ascension - (
            (time / day) * 2 * m.pi)  # rotation of the Earth beneath a prograde orbit

    coordinate = Point(latitude, longitude)
    # visible area prism
    height = visibleAreaPitchAngle(altitude)  # function of altitude and observing time
    width = visibleAreaRollAngle(altitude)  # function of altitude and observing time
    # rotation
    direction = m.cos(clockangle) * inclination  # angle of the direction of the orbit, 0 = horizontal
    points = np.empty(4, dtype=Point)
    points.put(0, Point(width, height))
    points.put(1, Point(-width, -height))
    points.put(2, Point(-width, height))
    points.put(3, Point(width, -height))
    rotation_matrix = [[m.cos(direction), m.sin(direction)],
                       [-m.sin(direction), m.cos(direction)]]  # CW rotation matrix
    pointsRotated = np.empty(4, dtype=Point)
    coords = np.empty(4, dtype=Point)
    for index in range(4):
        rotatedCoords = fromArray(np.matmul(rotation_matrix, points[index].asArray()))
        # pointsRotated.put(index, Point(rotatedCoords[0], rotatedCoords[1]))
        # pointsRotated.put(index, rotatedCoords)
        coords.put(index, Point(rotatedCoords.latitude + coordinate.latitude,
                                rotatedCoords.longitude + coordinate.longitude))
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
        latitude = (vertex.latitude * -1) + (0.5 * m.pi)
        matrixPosition.append(Point(int((latitude) / res), int(((vertex.longitude % (2 * m.pi)) / res))))
    for position in matrixPosition:
        if (position.latitude < minrow):
            minrow = position.latitude
        if (position.latitude > maxrow):
            maxrow = position.latitude
        if (position.longitude < mincol):
            mincol = position.longitude
        if (position.longitude > maxcol):
            maxcol = position.longitude
    rowheight = maxrow - minrow
    colwidth = maxcol - mincol

    if (colwidth > 0.5 * observableAreaMatrix.shape[1]):
        mincol = 2 * n_horizontal
        maxcol = 0
        for position in matrixPosition:
            if (position.longitude < 0.5 * observableAreaMatrix.shape[1]):
                position.longitude += observableAreaMatrix.shape[1]
            if (position.longitude < mincol):
                mincol = position.longitude
            if (position.longitude > maxcol):
                maxcol = position.longitude
        colwidth = maxcol - mincol
    for i in range(colwidth):
        for j in range(rowheight):
            row = int(j + minrow)
            col = int(i + mincol)
            if isInRectangle(matrixPosition, Point(row, col)):
                if (row < 0):
                    row = int(-row - 1)
                    col = int(col + 0.5 * n_horizontal)
                if (row > (observableAreaMatrix.shape[0] - 1)):
                    row = int(observableAreaMatrix.shape[0] - row % observableAreaMatrix.shape[0] - 1)
                    col = int(col + 0.5 * n_horizontal)
                if (col > (observableAreaMatrix.shape[1] - 1)):
                    col = col % observableAreaMatrix.shape[1]
                observableAreaMatrix[row][col] = 1

    counter = 0
    points = np.empty([20000, 2])
    for i in range(observableAreaMatrix.shape[0]):
        for j in range(observableAreaMatrix.shape[1]):
            if observableAreaMatrix[i, j] != 0:
                pointFirst = Point(((0.5 * m.pi) - ((m.pi / observableAreaMatrix.shape[0]) * i)),
                                   (((2 * m.pi) / observableAreaMatrix.shape[1]) * j))
                pointProjected = STXY(pointFirst.latitude, pointFirst.longitude)
                points[counter, 0] = pointProjected.latitude
                points[counter, 1] = pointProjected.longitude
                counter += 1
    plt.scatter(points[:counter, 0], points[:counter, 1], color='orange', s=0.5, alpha=0.03)

    return observableAreaMatrix


def calculateDayNightMatrix(time, graph=False):
    # Sunposition
    Earth_tilt = (23.5 / 180) * m.pi  # Earth's tilt in radians
    solar_inclination = Earth_tilt * m.sin(
        (time / (day * 365.25)) * 2 * m.pi)  # declination is 0 on March 21st (80th day of the year)
    solar_right_ascension = (time / day) * 2 * m.pi  # right ascension of sun rising
    r = 0.5 * m.pi
    # Matrix
    dayNightMatrix = np.zeros([n_vertical, n_horizontal])
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
            if (row >= dayNightMatrix.shape[0]):
                row = int(dayNightMatrix.shape[0] - 1 - (row % dayNightMatrix.shape[0]))
                col = int(col + n)
            if (row < 0):
                row = int(-row - 1)
                col = int(col + n)
            if (col < 0):
                col = int(dayNightMatrix.shape[1] + col + 1)
            if (col >= dayNightMatrix.shape[1] - 1):
                col = int(col % dayNightMatrix.shape[1])
            dayNightMatrix[row, col] = 1
            col = i + start_column
    if graph:
        points = pointsForGraph(dayNightMatrix)
        plt.scatter([point.latitude for point in points], [point.longitude for point in points], color='yellow', s=0.5,
                    alpha=0.02)
    return dayNightMatrix


def calculateDataMatrix(aoiMatrix, cloudCoverageMatrix, graph=False):
    dataMatrix = aoiMatrix * cloudCoverageMatrix
    points = []
    if graph:
        interpolateGrayscale = interp1d([30, 0], [0, 255])
        for i in range(dataMatrix.shape[0]):
            for j in range(dataMatrix.shape[1]):
                if dataMatrix[i, j] == 0 or dataMatrix[i, j] > 30:
                    continue
                position = Point(0.5 * m.pi - (m.pi / dataMatrix.shape[0] * i),
                                 2 * m.pi / dataMatrix.shape[1] * j + m.pi)
                projection = STXY(position.latitude, position.longitude)
                points.append([projection, interpolateGrayscale(dataMatrix[i, j])])
        fig, ax = plt.subplots()
        plt.scatter([entry[0].latitude for entry in points], [entry[0].longitude for entry in points],
                    c=[entry[1] for entry in points], cmap='inferno', s=0.2, marker=',', linewidths=0)
    return dataMatrix


def pointsForGraph(matrix):
    points = []
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if (matrix[i, j] != 0):
                point = Point(0.5 * m.pi - (m.pi / matrix.shape[0] * i), 2 * m.pi / matrix.shape[1] * j)
                point_proj = STXY(point.latitude, point.longitude)
                points.append(point_proj)
    return points


def calculateRequirementMatrix(dayNightMatrix, observableAreaMatrix, graph=False):
    requirementMatrix = dayNightMatrix * observableAreaMatrix
    if graph:
        points = pointsForGraph(requirementMatrix)
        plt.scatter([point.latitude for point in points], [point.longitude for point in points], color='blue', s=0.5,
                    alpha=0.03)
    return requirementMatrix


def calculateCloudCoverageMatrix(month=0, graph=False):
    if month == 0:
        cloudCoverageMatrix = np.loadtxt("cloudcoverage_map_weighted.csv", delimiter=',')
    else:
        cloudCoverageMatrix = np.loadtxt("cc_per_month/" + str(month) + ".csv", delimiter=',')
    if graph:
        points = []
        interpolateGrayscale = interp1d([30, 0], [0, 255])
        for i in range(cloudCoverageMatrix.shape[0]):
            for j in range(cloudCoverageMatrix.shape[1]):
                if cloudCoverageMatrix[i, j] == 0 or cloudCoverageMatrix[i, j] > 30:
                    continue
                position = Point(0.5 * m.pi - (m.pi / cloudCoverageMatrix.shape[0] * i),
                                 2 * m.pi / cloudCoverageMatrix.shape[1] * j)
                projection = STXY(position.latitude, position.longitude)
                points.append([projection, interpolateGrayscale(cloudCoverageMatrix[i, j])])
        plt.scatter([entry[0].latitude for entry in points], [entry[0].longitude for entry in points],
                    c=[entry[1] for entry in points], cmap='inferno', s=0.2, marker=',', linewidths=0)
    return cloudCoverageMatrix
