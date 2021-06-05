import numpy as np 
import os
from datetime import datetime

nVar = 5
dim = 2 

def plot_tecplot(in_file, out_file):
    
    infile = open(in_file, "r")
    
    # Read header lines 

    line = infile.readline()
    line = infile.readline()
    line = infile.readline()
    line = infile.readline()
    line = infile.readline()
    split_line = line.split()

    N_x = int(split_line[1])
    N_y = int(split_line[2])

    line = infile.readline()

    # Get number of cells in the mesh 

    split_line = line.split()
    n_cells = int(split_line[1])

    W = np.zeros((n_cells,nVar))

    coords = np.zeros((n_cells,dim))

    # Read the coordinates 

    for i in range(0, n_cells):
        line = infile.readline()
        split_line = line.split()
        coords[i, 0] = np.float64(split_line[0])
        coords[i, 1] = np.float64(split_line[1])

    # Read some more headers 
        
    line = infile.readline()
    line = infile.readline()
    line = infile.readline()

    for i in range(0, n_cells):
        line = infile.readline()
        split_line = line.split()
        for c in range (0, nVar):
            W[i,c] = np.float64(split_line[c])
        
    infile.close()

    # Reshape arrays to fit rectangular coordinates 

    coords = coords.reshape((N_y, N_x, dim))
    W = W.reshape((N_y, N_x, nVar))
    
    # Calculate density gradient 

    RHO = W[:,:,0]
    grad_RHO = np.zeros((N_y, N_x))
    
    for j in range(0, N_y):
        for i in range(0, N_x):
            if (j == 0):
                rho_y = RHO[j+1,i] - RHO[j,i]
            elif (j == N_y-1):
                rho_y = RHO[j,i] - RHO[j-1,i]
            else:
                rho_y = 0.5*(RHO[j+1,i] - RHO[j-1,i])

            if (i == 0):
                rho_x = RHO[j,i+1] - RHO[j,i]
            elif (i == N_x-1):
                rho_x = RHO[j,i] - RHO[j,i-1]
            else:
                rho_x = 0.5*(RHO[j,i+1] - RHO[j,i-1])

            grad_RHO[j,i] = np.hypot(rho_x, rho_y)

    grad_RHO = np.log(1.0 + np.abs(grad_RHO))

    # Start plotting 
            
    outfile = open(out_file, "w+")
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    
    print("TITLE = \"Multiphase Euler Equations. File created on: ", dt_string, "\"",  file = outfile)
    outfile.write("VARIABLES = \"x\", \"y\", \"Density\", \"x-Velocity\", \"y-Velocity\", \"Pressure\", \"Phi\",\"grad_RHO\" ")
    print('\nZone I = ', N_x, "J = ", N_y, file = outfile) 
    
    for j in range(0, N_y):
        for i in range (0, N_x):
            print(coords[j,i,0], coords[j,i,1], W[j,i,0], W[j,i,1], W[j,i,2], W[j,i,3], W[j,i,4], grad_RHO[j,i], file = outfile)

            
    outfile.close()
###############################################################################################################################

for in_file in os.listdir('.'):
    if in_file.endswith('.vtk'):
        print(in_file)
        out_file = in_file.replace("vtk", "dat")
        plot_tecplot(in_file, out_file)

