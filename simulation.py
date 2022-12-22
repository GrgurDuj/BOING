# import math as m
from math import pi, cos, sin, sqrt, ceil, trunc
import matplotlib.pyplot as plt
from numpy import arccos, empty, matmul, zeros, loadtxt, nonzero, minimum, argmin, where
# import numpy as np
from scipy.interpolate import interp1d
import shapely.geometry as sg
from shapely.geometry.polygon import Polygon
from Point import Point, fromArray
import time as tm

# parameters
timestamp_length = 180  # seconds
day = 86164  # seconds, 365.25/366.25 * 24 * 3600, seconds in a sidereal day
r_earth = 6371 * 1000  # meters
circumference = 2 * pi * r_earth  # meters
measure_time = 120  # seconds
vertical = pi  # radians
horizontal = 2 * pi  # radians
altitude = 700 * 10 ** 3  # meters
a = 6371000 + altitude  # meters
res = (0.1 / 180) * pi
n_horizontal = int(horizontal / res)
n_vertical = int(vertical / res)


def timeFunction(func):
    # This function shows the execution time of
    # the function object passed
    def wrap_func(*args, **kwargs):
        t1 = tm.time()
        result = func(*args, **kwargs)
        t2 = tm.time()
        print(f'Function {func.__name__!r} executed in {(t2 - t1):.4f}s')
        return result

    return wrap_func


def calcPeriod(alt):
    mu = 3.986E14
    return 2 * pi * sqrt((alt ** 3) / mu)

period = int(calcPeriod(a))  # seconds


def orbit_spherical(inclination, right_ascension, time, altitude):  # time is in seconds
    clockangle = (time / period) * 2 * pi
    latitude = sin(clockangle) * inclination
    longitude = clockangle + right_ascension - (
            (time / day) * 2 * pi)  # rotation of the Earth beneath a prograde orbit
    if latitude > 0.5 * pi:
        latitude = pi - latitude
        longitude += pi
    if latitude < -0.5 * pi:
        latitude = -(pi + latitude)
        longitude += pi
    coordinate = [latitude, longitude]
    return coordinate


def STXY(latitude, longitude):  # lo = longitude is perpendicular to equator, la = latitude is parallel to equator
    # Gall-Peters projection
    x = longitude % (2 * pi)  # longitude in radians from central meridian
    y = 2 * sin(latitude)  # latitude in radians from equator
    coordinate = Point(x, y)
    return coordinate


def visibleAreaPitchAngle(altitude):
    minPitch = (0.7 / 180) * pi  # rad
    maxPitch = (6.08 / 180) * pi  # rad5
    minAltitude = 450000  # m
    maxAltitude = 700000  # m
    interpolateLength = interp1d([minAltitude, maxAltitude], [minPitch, maxPitch])
    return interpolateLength(altitude)


def visibleAreaRollAngle(altitude):
    # corresponds to 30deg of roll from sat perspective (15 each side)
    roll = ((0.0048310913 * (altitude / 1000) + 0.00623608) / 180) * pi  # rad
    # roll gives observable angle with respect to the center of earth
    return roll


def satellitePrism(inclination, right_ascension, time):
    # satellite position
    clockangle = (time / period) * 2 * pi
    latitude = sin(clockangle) * inclination
    longitude = clockangle + right_ascension - (time / day) * 2 * pi  # rotation of the Earth beneath a prograde orbit

    coordinate = Point(latitude, longitude)
    # visible area prism
    height = visibleAreaPitchAngle(altitude)  # function of altitude and observing time
    width = visibleAreaRollAngle(altitude)  # function of altitude and observing time
    # rotation
    direction = cos(clockangle) * inclination  # angle of the direction of the orbit, 0 = horizontal
    points = empty(4, dtype=Point)
    points.put(0, Point(width, height))
    points.put(1, Point(-width, -height))
    points.put(2, Point(-width, height))
    points.put(3, Point(width, -height))
    rotation_matrix = [[cos(direction), sin(direction)],
                       [-sin(direction), cos(direction)]]  # CW rotation matrix
    coords = empty(4, dtype=Point)
    for index in range(4):
        rotatedCoords = fromArray(matmul(rotation_matrix, points[index].asArray()))
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

@timeFunction
def calculateObservableArea(visibleRegionInput, graph=False):
    mincol = n_horizontal
    maxcol = 0
    minrow = n_vertical
    maxrow = 0

    observableAreaMatrix = empty([n_vertical, n_horizontal])
    matrixPosition = []

    for vertex in visibleRegionInput:
        latitude = (vertex.latitude * -1) + (0.5 * pi)
        matrixPosition.append(Point(int((latitude) / res), int(((vertex.longitude % (2 * pi)) / res))))
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

    if (colwidth > 0.5 * n_horizontal):
        mincol = 2 * n_horizontal
        maxcol = 0
        for position in matrixPosition:
            if (position.longitude < 0.5 * n_horizontal):
                position.longitude += n_horizontal
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
                if (row > (n_vertical - 1)):
                    row = int(n_vertical - row % n_vertical - 1)
                    col = int(col + 0.5 * n_horizontal)
                if (col > (n_horizontal - 1)):
                    col = col % n_horizontal
                observableAreaMatrix[row][col] = 1
    if (graph):
        counter = 0
        points = empty([20000, 2])
        for i in range(n_vertical):
            for j in range(n_horizontal):
                if observableAreaMatrix[i, j] != 0:
                    pointFirst = Point(((0.5 * pi) - ((pi / n_vertical) * i)),
                                       (((2 * pi) / n_horizontal) * j))
                    pointProjected = STXY(pointFirst.latitude, pointFirst.longitude)
                    points[counter, 0] = pointProjected.latitude
                    points[counter, 1] = pointProjected.longitude
                    counter += 1
        plt.scatter(points[:counter, 0], points[:counter, 1], color='orange', s=0.5, alpha=0.03)

    return observableAreaMatrix


@timeFunction
def calculateDayNightMatrix(time, graph=False):
    # Sunposition
    Earth_tilt = (23.5 / 180) * pi  # Earth's tilt in radians
    solar_inclination = Earth_tilt * sin(
        (time / (day * 365.25)) * 2 * pi)  # declination is 0 on March 21st (80th day of the year)
    solar_right_ascension = (time / day) * 2 * pi  # right ascension of sun rising
    r = 0.5 * pi
    # Matrix
    dayNightMatrix = empty([n_vertical, n_horizontal])
    n = int(2 * r / res)
    start_ra = solar_right_ascension % (2 * pi) - r
    start_column = int(start_ra / res)
    for i in range(n):  # for each column in the shape
        col = i + start_column
        colsfromcenter = int(col - start_column - 0.5 * n)
        distancefromcenter = colsfromcenter * res / (0.5 * pi)
        theta = arccos(distancefromcenter)
        centerrowheight = int((-solar_inclination + 0.5 * pi) / res)
        absoluterowsheight = int((sin(theta) * 0.5 * pi) / res)
        for j in range(2 * absoluterowsheight):
            row = int(centerrowheight - absoluterowsheight + j)
            # Overflowing
            if (row >= n_vertical):
                row = int(n_vertical - 1 - (row % n_vertical))
                col = col + n
            elif (row < 0):
                row = -row - 1
                col = col + n
            if (col < 0):
                col = n_horizontal + col + 1
            elif (col >= n_horizontal - 1):
                col = int(col % n_horizontal)
            dayNightMatrix[row, col] = 1
            col = i + start_column
    if graph:
        points = pointsForGraph(dayNightMatrix)
        plt.scatter([point.latitude for point in points], [point.longitude for point in points], color='yellow', s=0.5,
                    alpha=0.02)
    return dayNightMatrix

@timeFunction
def calculateDataMatrix(aoiMatrix, cloudCoverageMatrix, graph=False):
    dataMatrix = aoiMatrix * cloudCoverageMatrix
    points = []
    if graph:
        interpolateGrayscale = interp1d([30, 0], [0, 255])
        for i in range(n_vertical):
            for j in range(n_horizontal):
                if dataMatrix[i, j] == 0 or dataMatrix[i, j] > 30:
                    continue
                position = Point(0.5 * pi - (pi / n_vertical * i),
                                 2 * pi / n_horizontal * j + pi)
                projection = STXY(position.latitude, position.longitude)
                points.append([projection, interpolateGrayscale(dataMatrix[i, j])])
        fig, ax = plt.subplots()
        plt.scatter([entry[0].latitude for entry in points], [entry[0].longitude for entry in points],
                    c=[entry[1] for entry in points], cmap='inferno', s=0.2, marker=',', linewidths=0)
    return dataMatrix


@timeFunction
def pointsForGraph(matrix):
    points = []
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if (matrix[i, j] != 0):
                point = Point(0.5 * pi - (pi / matrix.shape[0] * i), 2 * pi / matrix.shape[1] * j)
                point_proj = STXY(point.latitude, point.longitude)
                points.append(point_proj)
    return points


@timeFunction
def calculateRequirementMatrix(dayNightMatrix, observableAreaMatrix, graph=False):
    requirementMatrix = dayNightMatrix * observableAreaMatrix
    if graph:
        points = pointsForGraph(requirementMatrix)
        plt.scatter([point.latitude for point in points], [point.longitude for point in points], color='blue', s=0.5,
                    alpha=0.03)
    return requirementMatrix


@timeFunction
def calculateCloudCoverageMatrix(month=0, graph=False):
    if month == 0:
        cloudCoverageMatrix = loadtxt("cloudcoverage_map_weighted.csv", delimiter=',')
    else:
        cloudCoverageMatrix = loadtxt("cc_per_month/" + str(month) + ".csv", delimiter=',')
    if graph:
        points = []
        interpolateGrayscale = interp1d([30, 0], [0, 255])
        for i in range(n_vertical):
            for j in range(n_horizontal):
                if cloudCoverageMatrix[i, j] == 0 or cloudCoverageMatrix[i, j] > 30:
                    continue
                position = Point(0.5 * pi - (pi / cloudCoverageMatrix.shape[0] * i),
                                 2 * pi / cloudCoverageMatrix.shape[1] * j)
                projection = STXY(position.latitude, position.longitude)
                points.append([projection, interpolateGrayscale(cloudCoverageMatrix[i, j])])
        plt.scatter([entry[0].latitude for entry in points], [entry[0].longitude for entry in points],
                    c=[entry[1] for entry in points], cmap='inferno', s=0.2, marker=',', linewidths=0)
    return cloudCoverageMatrix


@timeFunction
def calculateFinalMatrix(DataMatrixDynamic, RequirementMatrix, graph=False):
    FinalMatrix = DataMatrixDynamic * RequirementMatrix
    size = 1  # indices, square
    # entries = nonzero(FinalMatrix)
    if (size == 1):
        if (len(nonzero(FinalMatrix)[0]) != 0):
            minval = min(FinalMatrix[nonzero(FinalMatrix)])
            index = where(FinalMatrix == minval)
            indices = [index[0][0], index[1][0]]  # This is the position that we are imaging
            DataMatrixDynamic[indices[0], indices[1]] += -1
            print("Observation done")
        else:
            print("No observation could be done")

    # Maximum number of observations, then least observations for resolving
    if graph:
        points = []
        for i in range(n_vertical):
            for j in range(n_horizontal):
                if (FinalMatrix[i, j] != 0):
                    point = Point(0.5 * pi - (pi / FinalMatrix.shape[0] * i), 2 * pi / FinalMatrix.shape[1] * j)
                    point_proj = STXY(point.latitude, point.longitude)
                    points.append(point_proj)
        plt.scatter([point.latitude for point in points], [point.longitude for point in points], color='red', s=0.2,
                    marker=',', linewidths=0)
    return FinalMatrix


# arclength = 2*pi*r_earth*(res*(180/pi)/360)
# fixedAltitude = (arclength/1000 + 8)/0.04
# print(fixedAltitude)
# print((arclength/20000)**-1)
print("You are running the wrong file idiot")
