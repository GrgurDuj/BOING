import simulation as sim
import math as m
import matplotlib.pyplot as plt
import numpy as np
from Point import Point, fromArray
import time as tm
from datetime import datetime


def main():
    print("Jk")
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
    altitude = 478 * 10 ** 3  # meters, AT 478KM, THE IMAGE SIZE IS 0.1 DEG 
    timelength = 100*timestamp_length  # seconds
    # radians height and width per image

    a = 6371000 + altitude  # meters
    period = int(sim.calcPeriod(a))  # seconds
    numTimestamps = int(timelength / timestamp_length)

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
    
    #Loading data matrix
    monthInit = datetime.now().month
    DataMatrixStaticMonthly = np.loadtxt("data_per_month/data" + str(monthInit) + "_new.csv", delimiter=',')
    DataMatrixDynamic = DataMatrixStaticMonthly
    print("file loaded after: " + str(tm.time()-startTime) + " seconds")
    #timeNow = datetime.timestamp(datetime.now())
    timeNow = 1671626098.291077 #for testing so time is constant
    # print(timeNow)
    #Main Loop
    for timestamp in range(numTimestamps):
        time = timeNow + timestamp*timestamp_length + 23*timestamp_length
        month = datetime.fromtimestamp(time).month
        if(month != monthInit):
            DataMatrixStaticMonthly = np.loadtxt("data_per_month/data" + str(month) + "_new.csv", delimiter=',')
            DataMatrixDynamic = DataMatrixStaticMonthly
            print("new file loaded")
            monthInit = month
        DayNightMatrix = sim.calculateDayNightMatrix(time)
        print("Main loop", (tm.time()-startTime))
        for satellite in range(numSatellites):
            VisibleRegion = sim.satellitePrism(inclinations[satellite], right_ascensions[satellite], time)
            ObservableAreaMatrix = sim.calculateObservableArea(VisibleRegion)
            RequirementMatrix = sim.calculateRequirementMatrix(DayNightMatrix, ObservableAreaMatrix)
            FinalMatrix = sim.calculateFinalMatrix(DataMatrixDynamic,RequirementMatrix, graph=True)
            print("Satellite", (tm.time()-startTime))
            continue
    
    print("Saving")
    fig.savefig("FinalMatrix", dpi=500)
    plt.show()
    print("Done")


if __name__ == "__main__":
    startTime = tm.time()
    main()
    print("Finished running after: " + str(tm.time()-startTime), " seconds")
