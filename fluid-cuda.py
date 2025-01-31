import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import Normalize
from multiprocessing import Process
import copy
import sys

def concentration(particules_x_coord, particules_y_coord, particule_n, Ny, Nx):
    psi = cp.zeros((Ny, Nx))
    for k in range(len(particules_x_coord)):
        psi[particules_y_coord[k], particules_x_coord[k]] += 1
    maxi = psi.max()
    return psi, psi / maxi, maxi

def main():
    Nx = 799
    Ny = 499
    Nt = 6000
    tau = 1
    time_save = 3000
    time_stream = 3000
    #particule_stream_active = True
    particule_stream_active = False

    
    write_data = False
    calced_data = False
    brownian = True
    
    diffusion = 0.1

    if len(sys.argv) == 2:
        print(sys.argv[1])
        if sys.argv[1]=="-w":
            write_data = True
        if sys.argv[1]=="-c":
            #Nt = 1800
            time_stream = 200
            calced_data = True
        if sys.argv[1]=="-cw" or sys.argv[1]=="-wc":
            #Nt = 2500
            time_stream = 3000
            write_data = True
            calced_data = True

    print("writing data at time 1000= ", write_data)
    print("using previously calculated data = ",calced_data)

    norm = Normalize(vmin=0, vmax=1)
    
    particule_n = 200000
    particules_x = 40 + 10 * cp.random.randn(particule_n)
    particules_x_coord = cp.asarray([int(x) for x in particules_x.get()])
    particules_y = 240 * cp.ones(particule_n) + 50 * cp.random.randn(particule_n)
    particules_y_coord = cp.asarray([int(y) for y in particules_y.get()])
    
    mass = 0.2
    
    x = np.linspace(0, Nx-1, Nx)
    y = np.linspace(0, Ny-1, Ny)
    X, Y = np.meshgrid(x, y)
    
    my_cmap = copy.copy(plt.colormaps['gray'])
    my_cmap_red = copy.copy(plt.colormaps['OrRd'])
    my_cmap.set_bad(alpha=0)
    my_cmap_red.set_bad(alpha=0)
    


    Nl = 9
    cxs = cp.array([0, 0, 1, 1, 1, 0, -1, -1, -1])
    cys = cp.array([0, 1, 1, 0, -1, -1, -1, 0, 1])
    weights = cp.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
    

    V = cp.ones((Ny, Nx, Nl)) + 0.03 * cp.random.randn(Ny, Nx, Nl)
    V[:, :, 3] = 3

    if calced_data:
        with open("data/steady.dat", 'r') as f:
            V_read = f.readlines()
            for j in range(Ny):
                for i in range(Nx):
                    V_l = V_read[j * Nx + i].split("/")
                    V[j, i] = cp.array([float(V_l[l]) for l in range(Nl)])
    
    V += 0.0001 * cp.random.randn(Ny, Nx, Nl)
    
    obstacle = cp.full((Ny, Nx), False)
    img = mpimg.imread('obstacle.png')
    for y in range(Ny):
        for x in range(Nx):
            if img[y, x, 3] != 0:
                obstacle[y, x] = True
    
    obstacle_shape = cp.where(obstacle, 1, cp.nan)
    
    plt.imshow(obstacle_shape.get(), interpolation='nearest', cmap=my_cmap)
    plt.show()
    
    rho = cp.sum(V, 2)
    momentum_x = cp.sum(V * cxs, 2) / rho
    momentum_y = cp.sum(V * cys, 2) / rho
    
    velocity_field = cp.sqrt(momentum_x**2 + momentum_y**2)

    if particule_stream_active:
        psi_concentration, alpha, maxi = concentration(particules_x_coord, particules_y_coord, particule_n, Ny, Nx)
        psi_concentration, alpha = psi_concentration.get(), alpha.get()
    
        img_ = plt.imshow(psi_concentration, cmap=my_cmap_red)
        colorbar_ = plt.colorbar(img_)
        colorbar_.set_label('Concentration', fontsize=20)
        colorbar_.ax.tick_params(labelsize=15)

    img = plt.imshow(velocity_field.get(), cmap='rainbow', alpha=1)
    colorbar = plt.colorbar(img)
    colorbar.set_label('Velocity', fontsize=20)
    colorbar.ax.tick_params(labelsize=15)

    for time in range(Nt):

        #boundaries conditions

        V[:,-1, [6,7,8]] = V[:,-2, [6,7,8]]
        V[:,0, [2,3,4]] = V[:,1, [2,3,4]]

        V[0,:, [8,1,2]] = V[1,:, [8,1,2]]
        V[-1,:, [6,5,4]] = V[-2,:, [6,5,4]]



        # roll velocities

        for i,cx,cy in zip(range(Nl),cxs,cys):
            V[:,:, i] = cp.roll(V[:,:, i], cx, axis = 1) # pour x
            V[:,:, i] = cp.roll(V[:,:, i], cy, axis = 0) # pour y


        # obstacle

        boundary = V[obstacle,:]


        # inverse velocities of boundary 

        boundary = boundary[:, [0,5,6,7,8,1,2,3,4]]


        # compute rho and momentum

        rho = np.sum(V, 2) # on additionne les vitesses des noeuds
        momentum_x = cp.sum(V*cxs, 2) / rho
        momentum_y = cp.sum(V*cys, 2) / rho


        # inverse velocities in obstacle

        V[obstacle, :] = boundary
        momentum_x[obstacle] = 0
        momentum_y[obstacle] = 0


        #collisions

        V_eq = cp.zeros(V.shape)
        for i,cx,cy,w in zip(range(Nl), cxs, cys, weights):
            V_eq[:,:,i] = rho*w * ( 1 + (3 * (cx*momentum_x + cy*momentum_y)) + (4.5 * (cx*momentum_x + cy*momentum_y)**2) - (1.5 * (momentum_x**2 + momentum_y**2)) )

        V = V -(V-V_eq)/tau



        # save computed velocities

        if time==time_save:
            if write_data:
                print("-w")
                V_cpu = V.get()
                with open("data/steady.dat", 'w+') as f:
                    for j in range(Ny):
                        for i in range(Nx):
                            for l in range(Nl):
                                f.write(str(V_cpu[j][i][l]) + "/")
                            f.write("\n")


        # particules streaming

        if (time > time_stream) and particule_stream_active:
            out_of_bounds_mask = (particules_x < 2) | (particules_x >= (Nx - 2)) | (particules_y >= (Ny - 2)) | (particules_y < 2)
            if not brownian:
                particules_x += momentum_x[particules_y_coord, particules_x_coord] * (tau/mass)
                particules_y += momentum_y[particules_y_coord, particules_x_coord] * (tau/mass)
            else:
                size = particules_x.size
                rand_x = cp.random.uniform(-1,1,size=size)
                rand_y = cp.random.uniform(-1,1,size=size)
                particules_x = momentum_x[particules_y_coord, particules_x_coord]*(tau/mass) +( diffusion * rand_x / mass )
                particules_y = momentum_y[particules_y_coord, particules_x_coord]*(tau/mass) +( diffusion * rand_y / mass )

            in_bounds_mask = ~out_of_bounds_mask
            particules_x_coord[in_bounds_mask] = cp.floor(particules_x[in_bounds_mask]).astype(int)
            particules_y_coord[in_bounds_mask] = cp.floor(particules_y[in_bounds_mask]).astype(int)

            # Remove particles that are out of bounds (use CuPy delete or boolean indexing)
            k_list = cp.nonzero(out_of_bounds_mask)[0]

            if k_list.size > 0:
                particules_x = cp.delete(particules_x, k_list)
                particules_y = cp.delete(particules_y, k_list)
                particules_x_coord = cp.delete(particules_x_coord, k_list)
                particules_y_coord = cp.delete(particules_y_coord, k_list)
    

           # Plot the vector field

       
        step=20
        if time%step== 0:
            #plt.title("time : "+str(time) + " | particules : " + str(len(particules_x)), size=50, pad=50)
            plt.title("time : "+str(time) + " | particules : ", size=50, pad=50)


            # curl

            #dfydx = momentum_x[2:, 1:-1] - momentum_x[:-2, 1:-1]
            #dfxdy = momentum_y[1:-1, 2:] - momentum_y[1:-1, :-2]
            #curl = dfydx - dfxdy
            #curl_mean_out = []
            #for y in range(2,Ny-2):
            #    for x in range(2,Nx-2):
            #        if (x<(Nx//10)) or (abs(y-(Ny//2)) > Ny//(2*5)):
            #           curl_mean_out.append(curl[y][x])
            #curl_mean = np.mean(curl_mean_out)
            #curl_normalized = curl - curl_mean
            #curl_normalized = curl - np.mean(curl_normalized)
            #print(np.mean(curl_normalized))
            #print(np.mean(np.sqrt(momentum_x**2 + momentum_y**2)))
            #plt.imshow(curl_normalized, cmap="bwr", interpolation='nearest')


            velocity_field = cp.sqrt(momentum_x**2 + momentum_y**2).get()
            img = plt.imshow(velocity_field, cmap='rainbow', alpha=0.2)
            plt.streamplot(X, Y, momentum_x.get(), momentum_y.get(), density=1.2, linewidth=1.3, arrowsize=2, arrowstyle='->', color='white')
            if particule_stream_active and time>time_stream :
                psi_concentration, alpha, maxi = concentration(particules_x_coord, particules_y_coord, particule_n, Ny, Nx)
                psi_concentration, alpha = psi_concentration.get(), alpha.get()

                img_ = plt.imshow(psi_concentration, cmap=my_cmap_red)

                colorbar_.update_normal(img_)
                img_.set_alpha(alpha)

            colorbar.update_normal(img)
            plt.imshow(obstacle_shape.get(), interpolation='nearest', cmap=my_cmap)
            #plt.scatter(particules_x, particules_y, color='black', marker='x', s=200)


            if len(str(time//step)) == 1:
                number = "00"+str(time//step)
            if len(str(time//step)) == 2:
                number = "0"+str(time//step)
            if len(str(time//step)) == 3:
                number = str(time//step)
            plt.savefig("res/res_stream"+number+".png")
            plt.pause(0.01)
            plt.cla()



if __name__ == "__main__":
    main()

