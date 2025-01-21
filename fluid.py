import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import copy

def dist(x1,y1,x2,y2):
    return(np.sqrt((x2-x1)**2 + (y2-y1)**2))

def main():
    Nx = 799
    Ny = 499
    Nt = 40000
    tau = 1

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
    V = np.ones((Ny,Nx,Nl)) + 0.01*np.random.randn(Ny,Nx,Nl)

    # print(len(V), len(V[0]), len(V[0][0]))
    #on assigne à celle vers la droite (3eme noeud)
    #V[:,:,3] = 2.3
    V[:,:,3] = 3

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
    #plt.imshow(Cylindre)
    plt.imshow(obstacle_shape, interpolation='nearest', cmap=my_cmap)
    plt.show()
    curl_mean=0
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
           # V_eq[:,:,i] = rho*w * ( 1 + 3 * (cx*momentum_x + cy*momentum_y) + 9 * (cx*momentum_x + cy*momentum_y)**2 / 2 - 3 * (momentum_x**2 + momentum_y**2)/2 )
            V_eq[:,:,i] = rho*w * ( 1 + (3 * (cx*momentum_x + cy*momentum_y)) + (4.5 * (cx*momentum_x + cy*momentum_y)**2) - (1.5 * (momentum_x**2 + momentum_y**2)) )

        V = V -(V-V_eq)/tau

        #plot


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
            plt.imshow(np.sqrt(momentum_x**2 + momentum_y**2))
            print(np.mean(np.sqrt(momentum_x**2 + momentum_y**2)))
            #plt.imshow(curl_normalized, cmap="bwr", interpolation='nearest')
            #plt.imshow(curl_normalized , cmap="bwr", interpolation='nearest')
            plt.streamplot(X, Y, momentum_x, momentum_y, density=1, linewidth=2, arrowsize=2, arrowstyle='->', color='white')
            plt.imshow(obstacle_shape, interpolation='nearest', cmap=my_cmap)
            #plt.savefig("res/res_"+str(time//step)+".png")
            #plt.savefig("res/0res_airless"+str(time//step)+".png")


            plt.savefig("res/res_stream"+str(time//step)+".png")
            plt.pause(0.01)
            plt.cla()

if __name__ == "__main__":
    main()
