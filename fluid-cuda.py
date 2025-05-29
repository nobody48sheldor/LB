import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.ndimage import gaussian_filter
#from multiprocessing import Process
import copy
import sys

def concentration(particules_x_coord, particules_y_coord, particule_n, Ny, Nx):
    psi = np.zeros((Ny,Nx))
    maxi = 0
    particules_x_coord_CPU, particules_y_coord_CPU = particules_x_coord.get(), particules_y_coord.get()
    for k in range(particules_x_coord.size):
        psi[particules_y_coord_CPU[k]][ particules_x_coord_CPU[k]] += 1
    maxi = max( [max(psi[j,:]) for j in range(Ny)] )
    mean_val= np.mean( [np.mean(psi[j,:]) for j in range(Ny)] )
    return(psi, psi/maxi, maxi)

def norm_grad_psi(psi,dx,dy):
    partial_psi_x = np.zeros_like(psi)
    partial_psi_y = np.zeros_like(psi)

    partial_psi_x[:, 1:-1] = (psi[:, 2:] - psi[:, :-2]) / (2 * dx)
    partial_psi_y[1:-1, :] = (psi[2:, :] - psi[:-2, :]) / (2 * dy)

    # aux bords

    partial_psi_x[:, 0] = (psi[:, 1] - psi[:, 0]) / dx
    partial_psi_x[:, 0] = (psi[:, 1] - psi[:, 0]) / dx

    partial_psi_y[0, :] = (psi[1, :] - psi[0, :]) / dy
    partial_psi_y[0, :] = (psi[1, :] - psi[0, :]) / dy

    return(np.sqrt(partial_psi_x**2 + partial_psi_y**2 ))

def main():
    mode = "slow"


    img = mpimg.imread('obstacle-courrant-shape.png')
    my_cmap = copy.copy(plt.colormaps['gray'])
    my_cmap_red = copy.copy(plt.colormaps['OrRd'])
    my_cmap.set_bad(alpha=0)
    my_cmap_red.set_bad(alpha=0)


    Nx = len(img[0])
    Ny = len(img)
    time_stop = False
    #Nt = 10000
    Nt = 40000
    tau = 0.81
    time_save = 2000
    time_stream = 1000
    particule_stream_active = True
    flux = True
    #v_init = 3/9

    if mode == "slow":
        time_stop = True
        v_init = 0.125
        Nt = 40000
    if mode == "normal":
        time_stop = True
        v_init = 3/9
        Nt = 25000
    if mode == "fast":
        time_stop = True
        v_init = 0.6
        Nt = 20000
    
    write_data = False
    calced_data = False
    brownian = True

    L = 5
    l = L*Ny/Nx
    dx = L/Nx
    dy = 0.8
    dl = L/(100*Nx)
    steps_to_reach = time_stream*20
    D = 2.23*(10**(-6))
    #D = 2.3*(10**(-9))
    dt = dl*dl/(2*D)


    print("\n")
    print("\n")
    print(f"dx = {dx}, dy = {dy}, dt = {dt}, D = {D}, v_init = {v_init} ~ mode : {mode}")
    print("\n")
    print("\n")
    #mass = 0.005/(v_init*Nx/L)

    if len(sys.argv) == 2:
        if sys.argv[1]=="-w":
            write_data = True
        if sys.argv[1]=="-c":
            time_stream = 200
            calced_data = True
        if sys.argv[1]=="-cw" or sys.argv[1]=="-wc":
            time_stream = time_stream
            write_data = True
            calced_data = True

    x_offset = 20
    #y_offset = 245
    y_offset = 195
    x_range = 20
    #y_range = 25
    y_range = 15
    particule_n = 3000
    #particule_n = 1500
    particules_x = x_offset + x_range * cp.random.randn(particule_n)
    particules_x_coord = cp.asarray([int(x) for x in particules_x.get()])
    particules_y = y_offset + y_range * cp.random.randn(particule_n)
    particules_y_coord = cp.asarray([int(y) for y in particules_y.get()])
    

    x = np.linspace(0, Nx-1, Nx)
    y = np.linspace(0, Ny-1, Ny)
    X, Y = np.meshgrid(x, y)
    

    Nl = 9
    cxs = cp.array([0, 0, 1, 1, 1, 0, -1, -1, -1])
    cys = cp.array([0, 1, 1, 0, -1, -1, -1, 0, 1])
    weights = cp.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
    

    V = cp.ones((Ny, Nx, Nl)) + 0.03 * cp.random.randn(Ny, Nx, Nl)
    V[:, :, 3] = v_init*Nl


    if calced_data:
        with open("data/steady.dat", 'r') as f:
            V_read = f.readlines()
            for j in range(Ny):
                for i in range(Nx):
                    V_l = V_read[j * Nx + i].split("/")
                    V[j, i] = cp.array([float(V_l[l]) for l in range(Nl)])
    
    V += 0.0001 * cp.random.randn(Ny, Nx, Nl)
 

    obstacle = cp.full((Ny, Nx), False)
    for y in range(Ny):
        for x in range(Nx):
            if img[y, x, 3] != 0:
                obstacle[y, x] = True
    
    obstacle_shape = cp.where(obstacle, 1, cp.nan)
    #plt.imshow(obstacle_shape.get(), interpolation='nearest', cmap=my_cmap)
    #plt.show()
   

    fig = plt.figure(figsize=(38.40,21.60))


    rho = cp.sum(V, 2)
    momentum_x = cp.sum(V * cxs, 2) / rho
    momentum_y = cp.sum(V * cys, 2) / rho
    momentum_field = cp.sqrt(momentum_x**2 + momentum_y**2)


    if particule_stream_active:
        psi_concentration, alpha, maxi = concentration(particules_x_coord, particules_y_coord, particule_n, Ny, Nx)
        psi_concentration_gaussian = gaussian_filter(psi_concentration, sigma=3)
        #print("maxi_init = ", maxi)
    
        img_ = plt.imshow(psi_concentration_gaussian, cmap=my_cmap_red)
        img_.set_alpha(alpha)
        colorbar_ = plt.colorbar(img_)
        colorbar_.set_label('Concentration', fontsize=20)
        colorbar_.ax.tick_params(labelsize=15)


    img = plt.imshow(momentum_field.get(), cmap='rainbow', alpha=1)
    colorbar = plt.colorbar(img)
    colorbar.set_label('Momentum', fontsize=20)
    colorbar.ax.tick_params(labelsize=15)



    # MAIN LB COMPUTATIONS

    for time in range(Nt):
        prev_len = 0

        #boundaries conditions

        if time_stop:
            if (time <= time_save):

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

                rho = cp.sum(V, 2) # on additionne les vitesses des noeuds
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
        else: # fait la meme chose, mais meme si on a t>time_save

            V[:,-1, [6,7,8]] = V[:,-2, [6,7,8]]
            V[:,0, [2,3,4]] = V[:,1, [2,3,4]]

            V[0,:, [8,1,2]] = V[1,:, [8,1,2]]
            V[-1,:, [6,5,4]] = V[-2,:, [6,5,4]]


            for i,cx,cy in zip(range(Nl),cxs,cys):
                V[:,:, i] = cp.roll(V[:,:, i], cx, axis = 1) # pour x
                V[:,:, i] = cp.roll(V[:,:, i], cy, axis = 0) # pour y


            boundary = V[obstacle,:]


            boundary = boundary[:, [0,5,6,7,8,1,2,3,4]]


            rho = cp.sum(V, 2)
            momentum_x = cp.sum(V*cxs, 2) / rho
            momentum_y = cp.sum(V*cys, 2) / rho


            V[obstacle, :] = boundary
            momentum_x[obstacle] = 0
            momentum_y[obstacle] = 0


            V_eq = cp.zeros(V.shape)
            for i,cx,cy,w in zip(range(Nl), cxs, cys, weights):
                V_eq[:,:,i] = rho*w * ( 1 + (3 * (cx*momentum_x + cy*momentum_y)) + (4.5 * (cx*momentum_x + cy*momentum_y)**2) - (1.5 * (momentum_x**2 + momentum_y**2)) )

            V = V -(V-V_eq)/tau



    # OTHER COMPUTATIONS : pour les particules et autres


        # save computed velocities at initial time (ca sert a déterminer : mass)
        if time==0:
            print(x_offset,y_offset, "\n")
            print("momentum = ",np.sqrt( (momentum_y.get()[y_offset,x_offset] )**2 +  (momentum_x.get()[y_offset,x_offset] )**2 ), "\n")
            momentum_val = np.sqrt( (momentum_y.get()[y_offset,x_offset] )**2 +  (momentum_x.get()[y_offset,x_offset] )**2 )
            mass = np.sqrt( (momentum_y.get()[y_offset,x_offset] )**2 +  (momentum_x.get()[y_offset,x_offset] )**2 )/(2*v_init*Nx/L)


        if time==time_save: # sauvegarde le champ calculé pour la prochaine éxécution
            if write_data:
                V_cpu = V.get()
                with open("data/steady.dat", 'w+') as f:
                    for j in range(Ny):
                        for i in range(Nx):
                            for l in range(Nl):
                                f.write(str(V_cpu[j][i][l]) + "/")
                            f.write("\n")


        # particules streaming

        if (time > time_stream) and particule_stream_active:
            if not brownian:
                particules_x += momentum_x[particules_y_coord, particules_x_coord] * (dt/mass)
                particules_y += momentum_y[particules_y_coord, particules_x_coord] * (dt/mass)
            else:
                size = particules_x.size
                rand_x = cp.random.uniform(-1,1,size=size)
                rand_y = cp.random.uniform(-1,1,size=size)
                particules_x += momentum_x[particules_y_coord, particules_x_coord]*(dt/mass) + ( dl * rand_x )
                particules_y += momentum_y[particules_y_coord, particules_x_coord]*(dt/mass) + ( dl * rand_y )


            # Remove particles that are out of bounds (use CuPy delete or boolean indexing)

            if flux:
                momentum_x_avg = cp.mean( momentum_x[ (y_offset - y_range) : (y_offset + y_range), (x_offset - x_range) : (x_offset + x_range) ] ).get()

                # distance que la particule moyenne aura parcourrue
                distance = (momentum_x_avg / mass ) * dt

                # nombre moyen de particules à rajouter
                number_add = int( abs( particule_n * ( distance / x_range) ))

                particules_x_add = ( x_offset - (distance/2) ) + distance * cp.random.randn(number_add)
                particules_y_add = y_offset + y_range * cp.random.randn(number_add)
                particules_x = cp.concatenate( (particules_x, particules_x_add) )
                particules_y = cp.concatenate( (particules_y, particules_y_add) )
                particules_x_coord = cp.concatenate( (particules_x_coord, particules_x_add) )
                particules_y_coord = cp.concatenate( (particules_y_coord, particules_y_add) )


            out_of_bounds_mask = (particules_x < 1) | (particules_x >= (Nx - 2)) | (particules_y >= (Ny - 2)) | (particules_y < 1)
            in_bounds_mask = ~out_of_bounds_mask
            particules_x_coord = cp.floor(particules_x).astype(int)
            particules_y_coord = cp.floor(particules_y).astype(int)

            in_obstacle_mask = obstacle[particules_y_coord, particules_x_coord].astype(bool)

            full_mask = out_of_bounds_mask | in_obstacle_mask 

            k_list = cp.where(full_mask)[0]
            if k_list.size > 0:
                particules_x = cp.delete(particules_x, k_list)
                particules_y = cp.delete(particules_y, k_list)
                particules_x_coord = cp.delete(particules_x_coord, k_list)
                particules_y_coord = cp.delete(particules_y_coord, k_list)
   

    # PLOTTING

       
        step=20
        if time%step== 0:
            plt.title("time : "+str(dt*time)[:5] + "s" + " | particules : " + str(particules_x.size), size=50, pad=50)


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


            momentum_field = cp.sqrt(momentum_x**2 + momentum_y**2).get()

            img = plt.imshow(momentum_field, cmap='rainbow', alpha=0.2)

            plt.streamplot(X, Y, momentum_x.get(), momentum_y.get(), density=1.2, linewidth=1.3, arrowsize=2, arrowstyle='->', color='white')
            plt.xlabel("$x$", fontsize=30)
            plt.ylabel("$y$", fontsize=30)

            if particule_stream_active and time>=time_stream :
                psi_concentration, alpha, maxi = concentration(particules_x_coord, particules_y_coord, particule_n, Ny, Nx)
                psi_concentration_gaussian, alpha_gaussian = gaussian_filter(psi_concentration, sigma=5), gaussian_filter(alpha, sigma=5)**0.5
                #print("maxi_after = ",maxi)

                img_ = plt.imshow(psi_concentration_gaussian, cmap=my_cmap_red)

                colorbar_.update_normal(img_)
                img_.set_alpha(alpha_gaussian)

                #print("particle stream active; maxi = ", maxi,"\n")
   
            colorbar.update_normal(img)
            plt.imshow(obstacle_shape.get(), interpolation='nearest', cmap=my_cmap)
            #plt.scatter(particules_x, particules_y, color='black', marker='x', s=200)


            # saving images

            if len(str(time//step)) == 1:
                number = "000"+str(time//step)
            if len(str(time//step)) == 2:
                number = "00"+str(time//step)
            if len(str(time//step)) == 3:
                number = "0"+str(time//step)
            if len(str(time//step)) == 4:
                number = str(time//step)

            plt.xticks(ticks=[(i*Nx/10) for i in range(11)], labels=[str(i*L/10)[:5] for i in range(11)], fontsize=18)
            plt.yticks(ticks=[(i*Ny/10) for i in range(11)], labels=[str(i*l/10)[:5] for i in range(11)], fontsize=18)
            plt.savefig("res/res_stream"+number+".png", dpi=100)


            # message de progression

            if (time<time_save) and (time<time_stream):
                percent = f"  {str(int(100*time/Nt))}% ~ [ {number} / {Nt//step} ] ~ [ {time} / {Nt} ] ~ maxi = {maxi}"
            if (time<time_save) and (time>time_stream):
                percent = f"  {str(int(100*time/Nt))}% ~ [ {number} / {Nt//step} ] ~ [ {time} / {Nt} ] ~ maxi = {maxi} | number add : {number_add}"
            if (time>=time_save) and (time<time_stream):
                percent = f"  {str(int(100*time/Nt))}% ~ [ {number} / {Nt//step} ] ~ [ {time} / {Nt} ] ~ maxi = {maxi} | time-saved"
            if time>=time_save and (time>time_stream):
                percent = f"  {str(int(100*time/Nt))}% ~ [ {number} / {Nt//step} ] ~ [ {time} / {Nt} ] ~ maxi = {maxi} | time-saved | number add : {number_add}"

            print(" " * prev_len, end="\r")  # Clear previous output
            print(percent , end="\r", flush=True)  # Print new output
            prev_len = len(percent)
            plt.cla()


    # calcule et affiche D(x,y)

    if particule_stream_active:
        plt.clf()
        plt.title("D(x,y) at equilibrium, with a factor", size=20, pad=50)
        D_field = (momentum_field*psi_concentration_gaussian)/(norm_grad_psi(psi_concentration_gaussian,dx,dy))
        img_D = plt.imshow(momentum_field, interpolation='nearest', cmap=my_cmap)
        img_border = plt.imshow(obstacle_shape.get(), interpolation='nearest', cmap=my_cmap_red)
        colorbar_D = plt.colorbar(img_D)
        plt.savefig("res/D.png", dpi=100)
        plt.show()



if __name__ == "__main__":
    main()

