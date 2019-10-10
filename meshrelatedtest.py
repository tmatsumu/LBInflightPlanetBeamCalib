import numpy as np

def meshgrid2array_XY(X,Y):
    num_X = len(X[0,:])
    num_Y = len(Y[:,0])
    x_arr = []; y_arr=[]
    for i in range(num_X):
		x_arr.append(X[0,i])
    for j in range(num_Y):
		y_arr.append(Y[j,0])
    return x_arr, y_arr

x = np.array([1,2,3])
y = np.array([5,6,7,8])

X, Y = np.meshgrid(x,y)

x_arr, y_arr = meshgrid2array_XY(X,Y)

print x_arr, y_arr

Z = X**2 + Y**2

def meshgrid2array_Z(Z):
    num_x = len(Z[:,0])
    num_y = len(Z[0,:])
    z_arr = []
    for i in range(num_x):
    	for j in range(num_y):
			z_arr.append(Z[i,j])
    return z_arr

print Z
print meshgrid2array_Z(Z)