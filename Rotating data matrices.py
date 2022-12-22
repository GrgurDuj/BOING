import simulation as sim
import math as m
import matplotlib.pyplot as plt
import numpy as np
from Point import Point, fromArray
import time as tm
from datetime import datetime, timedelta

# DataMatrix = np.loadtxt("data_per_month/data12_new.csv", delimiter=',')
# DataMatrixNew = np.empty(DataMatrix.shape)
# for i in range(DataMatrix.shape[0]):
#     for j in range(DataMatrix.shape[1]):
#         col = int(j + 0.5*DataMatrix.shape[1])
#         if(col >= DataMatrix.shape[1]):
#             col = col%DataMatrix.shape[1]
#         DataMatrixNew[i,col] = DataMatrix[i,j]
# points = []

# image = plt.imread("Gall–Peters_projection.jpg")
# fig, ax = plt.subplots()
# # for i in range(int(days)):
# #     string = str(i)
# #     string = image
# #     string = ax.imshow(image, extent = [i*2*m.pi,(1+i)*2*m.pi,-2,2])
# image = ax.imshow(image, extent=[0, 2 * m.pi, -2, 2])
# plt.xlim(0, 2 * m.pi)
# plt.ylim(-2, 2)
    
# for i in range(DataMatrix.shape[0]):
#             for j in range(DataMatrix.shape[1]):
#                 if(DataMatrixNew[i,j] != 0):
#                     point = Point(0.5 * m.pi - (m.pi / DataMatrix.shape[0] * i), 2 * m.pi / DataMatrix.shape[1] * j)
#                     point_proj = sim.STXY(point.latitude, point.longitude)
#                     points.append(point_proj)
# plt.scatter([point.latitude for point in points], [point.longitude for point in points], color='red', s = 0.5)
# plt.show()
# #np.savetxt("data_per_month/data12_new.csv", DataMatrixNew, delimiter=',')