import simulation as sim
import math as m
import matplotlib.pyplot as plt
import numpy as np
from Point import Point, fromArray
import time as tm


def main():
    # constant parameters
    timestamp_length = 180  # seconds
    day = 86164  # seconds, 365.25/366.25 * 24 * 3600, seconds in a sidereal day
    r_earth = 6371 * 1000  # meters
    measure_time = 120  # seconds
    vertical = m.pi  # radians
    horizontal = 2 * m.pi  # radians

    # inputs for multiple satellites

    inclinations = [(60 / 180) * m.pi, (50 / 180) * m.pi, (40 / 180) * m.pi, (10 / 180) * m.pi, (20 / 180) * m.pi,
                    (30 / 180) * m.pi, (55 / 180) * m.pi, (70 / 180) * m.pi, ]  # radians
    right_ascensions = [(20 / 180) * m.pi, (40 / 180) * m.pi, (60 / 180) * m.pi, (80 / 180) * m.pi, (100 / 180) * m.pi,
                        (120 / 180) * m.pi, (140 / 180) * m.pi, (160 / 180) * m.pi, ]  # radians
    altitude = 700 * 10 ** 3  # meters
    timelength = timestamp_length  # seconds
    # radians height and width per image

    a = 6371000 + altitude  # meters
    period = int(sim.calcPeriod(a))  # seconds
    n = int(timelength / timestamp_length)

    orbits = timelength / period
    days = timelength / day
    numSatellites = 1  # len(inclinations)

    # figure
    image = plt.imread("Gall–Peters_projection.jpg")
    fig, ax = plt.subplots()
    # for i in range(int(days)):
    #     string = str(i)
    #     string = image
    #     string = ax.imshow(image, extent = [i*2*m.pi,(1+i)*2*m.pi,-2,2])
    image = ax.imshow(image, extent=[0, 2 * m.pi, -2, 2])
    plt.xlim(0, 2 * m.pi)
    plt.ylim(-2, 2)

    '''
    # Main loop
    for i in range(12):
        cloudCoverageMatrix = np.loadtxt("cc_per_month/" + str(i) + ".csv", delimiter=',')
        dataMatrix = sim.calculateDataMatrix(aoiMatrix, cloudCoverageMatrix)
        np.savetxt(dataMatrix, "data_per_month/data" + str(i) + "1.csv", delimeter=',')

    # positions_sun = np.empty(n, dtype=Point)
    startTime = tm.time()
    for i in range(n):
        time = i * timestamp_length + 160.301 * day
        dayNightMatrix = sim.calculateDayNightMatrix(time)

        for j in range(numSatellites):
            time += 0.1 * day * j ** 2
            # Satellite
            point_spherical = sim.orbit_spherical(inclinations[j], right_ascensions[j], time)
            point_projection = sim.STXY(point_spherical[0], point_spherical[1])
            plt.scatter(point_projection.latitude, point_projection.longitude, color='red', s=1)

            # Visible area
            # visibleRegion = satellitePrism(inclinations[j], right_ascensions[j], time)
            # observableAreaMatrix = calculateObservableArea(visibleRegion)
            # projections = np.empty(4, dtype=Point)
            # OAprojections = np.empty(4, dtype=Point)
            # for index in range(4):
            #     #overflowing
            #     if (visibleRegion[index].latitude > 0.5 * m.pi):
            #         visibleRegion[index].latitude = m.pi - visibleRegion[index].latitude
            #         visibleRegion[index].longitude += m.pi
            #     if (visibleRegion[index].latitude < -0.5 * m.pi):
            #         visibleRegion[index].latitude = -(m.pi + visibleRegion[index].latitude)
            #         visibleRegion[index].longitude += m.pi
            #     projections.put(index, STXY(visibleRegion[index].latitude, visibleRegion[index].longitude))
            # # for projection in projections:
            # #     plt.scatter(projection.latitude, projection.longitude, color='red', s=1)

            # Requirement matrix

        # Cloud coverage matrix

    # Datamatrix
    # counter = 0
    # month = 0
    # Data_matrix = dataMatrix(AOI_matrix, cloudcoverage_matrix)
    # points = []

    print("Saving")
    # fig.savefig("Requirement_matrix_visualization.jpg", dpi=2000)

    # dataframe = pd.DataFrame(Data_matrix)
    # dataframe.to_csv("Data_matrix.csv", header = False)
    print("%.2f seconds for dayNightMatrix" % (tm.time() - startTime))
    '''
    plt.show()
    print("Done")


if __name__ == "__main__":
    main()
