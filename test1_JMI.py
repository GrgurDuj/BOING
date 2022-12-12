import math as m
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#atest=np.loadtxt("cloudcoverage_map_yearly_average.csv", delimiter=',', dtype=float)
#print(atest.shape)
#cloudcoverage_matrix = (np.ceil((np.loadtxt("cloudcoverage_map_yearly_average.csv", delimiter=',', dtype=float))**(-1))).astype(int)
#print(cloudcoverage_matrix)
#print(np.average(cloudcoverage_matrix))
#dataframe = pd.DataFrame(cloudcoverage_matrix)
#dataframe.to_csv("clouds_integer_yearly.csv")

cloudcoverage_matrix = np.ceil((np.loadtxt("cloudcoverage_map_yearly_average.csv", delimiter=',', dtype=float))**(-1))
#print(cloudcoverage_matrix)
#print(np.average(cloudcoverage_matrix))
dataframe = pd.DataFrame(cloudcoverage_matrix)
dataframe.to_csv("clouds_integer_yearly.csv")

# cloudcoverage_matrix = np.transpose(cloudcoverage_matrix)
# cloudcoverage_matrix = np.flip(cloudcoverage_matrix)
# cloudcoverage_matrix = np.fliplr(cloudcoverage_matrix)
# np.savetxt("cloudcoverage_matrix_corrected.csv",cloudcoverage_matrix, delimiter=',')
#dataframe.to_csv("DayNightMatrix.csv")
#print(AOI_matrix.shape)



#from main simulation file:, around lines ~290
# AOI_matrix = np.loadtxt("AOI_matrix.csv", delimiter=',', dtype=int)
# CLOUD_matrix = ((np.delete((np.loadtxt("clouds_integer_yearly.csv", delimiter=',', skiprows = 1, dtype=str)), 0, 1)).astype(float)).astype(int)
# #DATA_matrix = AOI_matrix*CLOUD_matrix
# #print(DATA_matrix, DATA_matrix.shape)