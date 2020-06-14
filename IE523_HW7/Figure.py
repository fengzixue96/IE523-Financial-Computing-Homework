import numpy as np
import matplotlib.pyplot as plt
file=open('C:/Users/HP/Desktop/Lesson/523/Assignment/Assignment 7/time')
data=[]
for line in file.readlines():
    row = line.strip().split(',')
    data.append(row)
array=np.array(data)
print(array[:,0])
plt.plot(array[:,0],array[:,1])
plt.plot(array[:,0],array[:,2])
plt.legend(['Repeated Squares', 'Brute Force Multiplication'])
plt.xlabel('Exponent')
plt.ylabel("Computation Time")