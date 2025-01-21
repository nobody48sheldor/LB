# petit exemple python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import os

files = os.listdir("res")
FILES = {}
for i in files:
    FILES[int(i.split('res_')[-1].split('.')[0])] = i

n = 800
Nt = len(files)
L = list(FILES.keys())
L.sort()
files = [FILES[i] for i in L]
print(files)


L = 0.5

x = np.linspace(-L,L,n)
y = np.linspace(-L,L,n)

def read(file):
    with open("res/"+file, 'r') as f:
        content = f.readlines()

    i = 0
    j = 0
    Nx, Ny = content[0].split("|")[0:2]
    Nx, Ny = int(Nx), int(Ny)
    print("(nx, ny) = ", Nx, Ny)
    PSI = [[0 for j in range(Ny)] for i in range(Nx)]
    values_str = content[1].split("/")[:-1]
    values = [float(val) for val in values_str]
    for i in range(Nx):
        for j in range(Ny):
                PSI[i][j] = values[i + Ny*j]
    return(PSI,Nx,Ny)


# general func

fig = plt.figure()
plt.ion()

PSI = [[] for time in range(Nt)]
AVG = [[] for time in range(Nt)]
maxi = 0

for time in range(Nt):

    # read
    PSI[time], Nx, Ny = read(files[time])

    avg = np.mean(PSI[time])
    AVG[time] = avg
   # PSI[time] = [[PSI[time][i][j] - avg for i in range(Nx)] for j in range(Ny)]

    #maxi = 1.1*max(maxi, max([ max([ max([abs(PSI[time][i][j][k]) for i in range(Nx)]) for j in range(Ny)]) for k in range(Nz) ]))
   # maxi = 1.1*max([ max([abs(PSI[time][i][j]) for i in range(Nx)]) for j in range(Ny)])
   # PSI[time] = [[PSI[time][i][j]/maxi for i in range(Nx)] for j in range(Ny)]
    #print(PSI)
    print(avg)

#for time in range(Nt):
    #PSI[time] = [[[PSI[time][i][j][k]/maxi for i in range(Nx)] for j in range(Ny)] for k in range(Nz)]

for time in range(Nt):
    plt.cla()
    #plt.set_xlim([-L, L])
    #plt.set_ylim([-L, L])


    #arr = [[[PSI_C[i][j][k][l] for k in range(n_)] for j in range(n_)] for i in range(n_)]
    #print(len(arr))
    #print(len(arr[0]))
    #print(len(arr[0][0]))
    plt.title("frame : " + str(time))
    #print("done")
    plt.imshow(PSI[time])
    print("time : ", time)
    plt.pause(0.3)
    plt.savefig("renders/2D-{}.png".format(time), dpi=300)

print(AVG)
