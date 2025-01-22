import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import copy
import sys

def dist(x1,y1,x2,y2):
    return(np.sqrt((x2-x1)**2 + (y2-y1)**2))


def concentration(particules_x_coord, particules_y_coord, particule_n):
    psi = np.zeros(Ny,Nx)
    maxi = 0
    for j in range(Ny):
        for i in range(Nx):
            for k in range(particule_n):
                if particule_y_coord[k] == j:
                    if particule_x_coord[k] == i:
                        psi[j][i] += 1
                if psi[j][i] > maxi:
                    maxi = psi[j][i]
    return(psi/maxi)

def main():
    Nx = 799
    Ny = 499
    Nt = 40000
    tau = 1

    write_data = False
    calced_data = False

    if len(sys.argv) == 2:
        print(sys.argv[1])
        if sys.argv[1]=="-w":
            write_data = True
        if sys.argv[1]=="-c":
            calced_data = True
        if sys.argv[1]=="-cw" or sys.argv[1]=="-wc":
            write_data = True
            calced_data = True

    print("writing data at time 1000= ", write_data)
    print("using previously calculated data = ",calced_data)


    particule_n = 50
    particules_x = 26 + 10*np.random.randn(particule_n)
    particules_x_coord = [int(particules_x[i]) for i in range(particule_n)]
    particules_y = 240*np.ones(particule_n) + 60*np.random.randn(particule_n)
    particules_y_coord = [int(particules_y[i]) for i in range(particule_n)]
    mass = 0.2

    x = np.linspace(0, Nx-1, Nx)
    y = np.linspace(0, Ny-1, Ny)
    X, Y = np.meshgrid(x, y)

    my_cmap = copy.copy(plt.cm.get_cmap('gray')) # get a copy of the gray color map
    my_cmap.set_bad(alpha=0) # set how the colormap handles 'bad' values

    # vitesses
    Nl = 9 #nombre de vitesses discretes

    #vecteurs vitesse discrets
    cxs = np.array([0, 0, 1, 1,  1,  0, -1, -1, -1])
    cys = np.array([0, 1, 1, 0, -1, -1, -1,  0,  1])

    # poids
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])

    # conditions initiales
    V = np.ones((Ny,Nx,Nl)) + 0.03*np.random.randn(Ny,Nx,Nl)
    V[:,:,3] = 3
    if calced_data:
        with open("data/steady.dat",'r') as f:
            V_read = f.readlines()
            for j in range(Ny):
                for i in range(Nx):
                    V_l = V_read[j*Nx + i].split("/")
                    V[j][i] = [float(V_l[l]) for l in range(Nl)]

    V = V + 0.0001*np.random.randn(Ny,Nx,Nl)

    #on assigne à celle vers la droite (3eme noeud)

    obstacle = np.full((Ny,Nx), False) # False : libre, True : obstacle
    #center = (Nx//4, Ny//2)
    #radius = 40
    #for y in range(Ny):
    #    for x in range(Nx):
    #        if dist(center[0],center[1],x,y) < radius :
    #            Cylindre[y][x] = True

    #img = mpimg.imread('airless.png')
    #img = mpimg.imread('wing.png')
    img = mpimg.imread('obstacle.png')
    for y in range(Ny):
        for x in range(Nx):
            if img[y][x][3] != 0:
                obstacle[y][x] = True
    obstacle_shape = [[np.nan if obstacle[y][x] == False else 1 for x in range(Nx)] for y in range(Ny)]




    # solve
    plt.imshow(obstacle_shape, interpolation='nearest', cmap=my_cmap)
    plt.show()

    curl_mean=0
    rho = np.sum(V, 2) # on additionne les vitesses des noeuds
    momentum_x = np.sum(V*cxs, 2) / rho
    momentum_y = np.sum(V*cys, 2) / rho

    velocity_field = np.sqrt(momentum_x**2 + momentum_y**2)
    img = plt.imshow(velocity_field, cmap='rainbow')
    colorbar = plt.colorbar(img, label='velocity')

    for time in range(Nt):

        #boundaries conditions
        V[:,-1, [6,7,8]] = V[:,-2, [6,7,8]]
        V[:,0, [2,3,4]] = V[:,1, [2,3,4]]

        V[0,:, [8,1,2]] = V[1,:, [8,1,2]]
        V[-1,:, [6,5,4]] = V[-2,:, [6,5,4]]



        # roll velocities
        for i,cx,cy in zip(range(Nl),cxs,cys):
            V[:,:, i] = np.roll(V[:,:, i], cx, axis = 1) # pour x
            V[:,:, i] = np.roll(V[:,:, i], cy, axis = 0) # pour y


        # obstacle
        boundary = V[obstacle,:]
        # on inverse les vitesses
        boundary = boundary[:, [0,5,6,7,8,1,2,3,4]]


        # compute rho and momentum
        rho = np.sum(V, 2) # on additionne les vitesses des noeuds
        momentum_x = np.sum(V*cxs, 2) / rho
        momentum_y = np.sum(V*cys, 2) / rho


        # inverse velocities in obstacle
        V[obstacle, :] = boundary
        momentum_x[obstacle] = 0
        momentum_y[obstacle] = 0


        #collisions

        V_eq = np.zeros(V.shape)
        for i,cx,cy,w in zip(range(Nl), cxs, cys, weights):
            V_eq[:,:,i] = rho*w * ( 1 + (3 * (cx*momentum_x + cy*momentum_y)) + (4.5 * (cx*momentum_x + cy*momentum_y)**2) - (1.5 * (momentum_x**2 + momentum_y**2)) )

        V = V -(V-V_eq)/tau

        #plot
        if time==1000:
            if write_data:
                print("-w")
                with open("data/steady.dat", 'w+') as f:
                    for j in range(Ny):
                        for i in range(Nx):
                            for l in range(Nl):
                                f.write(str(V[j][i][l]) + "/")
                            f.write("\n")

        if time>200:
            plt.title("streaming")
            for k in range(len(particules_x)):
                particules_x[k] += momentum_x[ particules_y_coord[k] ][ particules_x_coord[k] ] * (tau/mass)
                particules_y[k] += momentum_y[ particules_y_coord[k] ][ particules_x_coord[k] ] * (tau/mass)
                particules_x_coord[k] = int(particules_x[k])
                particules_y_coord[k] = int(particules_y[k])





        # Plot the vector field

       
        step=20
        if time%step== 0:
            dfydx = momentum_x[2:, 1:-1] - momentum_x[:-2, 1:-1]
            dfxdy = momentum_y[1:-1, 2:] - momentum_y[1:-1, :-2]
            curl = dfydx - dfxdy
            curl_mean_out = []
            for y in range(2,Ny-2):
                for x in range(2,Nx-2):
                    if (x<(Nx//10)) or (abs(y-(Ny//2)) > Ny//(2*5)):
                       curl_mean_out.append(curl[y][x])
            curl_mean = np.mean(curl_mean_out)
                
            curl_normalized = curl - curl_mean
            curl_normalized = curl - np.mean(curl_normalized)
            print(np.mean(curl_normalized))
            print(np.mean(np.sqrt(momentum_x**2 + momentum_y**2)))
            #plt.imshow(curl_normalized, cmap="bwr", interpolation='nearest')
            #plt.imshow(curl_normalized , cmap="bwr", interpolation='nearest')
            velocity_field = np.sqrt(momentum_x**2 + momentum_y**2)
            img = plt.imshow(velocity_field, cmap='rainbow')
            img = plt.imshow(concentration(particules_x_coord,particules_y_coord,partiucle_n), cmap='Reds')
            colorbar.update_normal(img)
            plt.streamplot(X, Y, momentum_x, momentum_y, density=1, linewidth=2, arrowsize=2, arrowstyle='->', color='white')
            plt.imshow(obstacle_shape, interpolation='nearest', cmap=my_cmap)
            plt.scatter(particules_x, particules_y, color='black', marker='x', s=200)


            plt.savefig("res/res_stream"+str(time//step)+".png")
            plt.pause(0.01)
            plt.cla()

if __name__ == "__main__":
    main()
