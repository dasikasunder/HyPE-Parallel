import numpy as np 
import matplotlib.pyplot as plt 
import vtk
from vtk.util import numpy_support
from vtk.util.numpy_support import vtk_to_numpy
import os
from datetime import datetime

dim = 2 # Number of dimensions in the simulation

##########################################################################################################################################
# MatPlotLib Plotter 
##########################################################################################################################################

def plotMatPlotLib(Name,fileName,fieldNames, fieldTypes, X, Y, sol):

    JMAX = sol.shape[0]
    IMAX = sol.shape[1]
    nVar = sol.shape[2]

    im_ratio = IMAX/JMAX # Reverse to get better size

    slice_interval = 10 # Increase for lesser arrow plots
    skip = (slice(None,None,slice_interval),slice(None,None,slice_interval))

    for c in range(0,nVar):
        # Set the title of the plot         
        if (fieldTypes[c] == 0):
            plt.figure(c)
            plt.contourf(X,Y,sol[:,:,c],cmap='turbo',levels=500)
            plt.axis('equal')
            plt.xlabel("x")
            plt.ylabel("y")
            plt.colorbar(fraction=0.047*im_ratio)
            outFile = fileName.replace("vts","png")
            outFile = outFile.replace("sol",fieldNames[c])
            plt.savefig(outFile,dpi=500)
            plt.clf()
        if (fieldTypes[c] == 1):
            U = sol[:,:,c]; V = sol[:,:,c+1]
            Mag = np.hypot(U,V)
            plt.figure(c)
            plt.contourf(X,Y,Mag,cmap='turbo',levels=500)
            plt.quiver(X[skip],Y[skip],U[skip],V[skip],units='xy',cmap='turbo')
            plt.axis('equal')
            plt.xlabel("x")
            plt.ylabel("y")
            plt.colorbar()
            plt.tight_layout()
            outFile = fileName.replace("vts","png")
            outFile = outFile.replace("sol",fieldNames[c])
            plt.savefig(outFile,dpi=500)
            plt.clf()
        


##########################################################################################################################################
# Tecplot Plotter 
##########################################################################################################################################

def plotTecPlot(Name,fileName,fieldNames, fieldTypes, X, Y, sol):
    
    JMAX = sol.shape[0]
    IMAX = sol.shape[1]
    nVar = sol.shape[2]
    
    # Print the header part 

    outFileName = fileName.replace("vts","dat")
    outFile = open(outFileName, "w+")
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")   
    title ='TITLE = "' + Name + '. File created on: ' + dt_string + '"'
    print(title,file=outFile)
    print('VARIABLES = "x", "y",', end = " ", file = outFile)

    for c in range(0,nVar):

        if (fieldTypes[c] == 0):
            text = fieldNames[c]
        else:
            text = fieldNames[c]+str(fieldTypes[c])

        if (c == nVar-1):
            print('"' + text + '"', file = outFile)
        else:
            print('"' + text + '",' , end = " ", file = outFile)

    print('Zone I = ', IMAX, 'J = ', JMAX, file = outFile)
    for j in range(0,JMAX):
        for i in range (0,IMAX):
            print(X[j,i], Y[j,i], end = " ", file = outFile)
            for c in range(0,nVar):
                if (c == nVar-1):
                    print(sol[j,i,c], file = outFile)
                else:
                    print(sol[j,i,c], end = " ", file = outFile)

##########################################################################################################################################
# VTK plotter
##########################################################################################################################################

def plotVTK(Name,fileName,fieldNames, fieldTypes, X, Y, sol):
    JMAX = sol.shape[0]
    IMAX = sol.shape[1]
    nVar = sol.shape[2]

    # Print the header part

    outFileName = fileName.replace("vts","vtk")
    outFile = open(outFileName, "w+")
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    outFile.write("# vtk DataFile Version 2.0\n")
    title = Name + ' Created On: ' + dt_string
    outFile.write(title+"\n")
    outFile.write("ASCII\n")
    outFile.write("DATASET STRUCTURED_GRID\n")
    print('DIMENSIONS', IMAX, JMAX, 1 , file = outFile)
    print('POINTS', IMAX*JMAX, "double" , file = outFile)

    # Plot Coordinates

    for j in range(0, JMAX):
        for i in range(0,IMAX):
            print(X[j,i], Y[j,i], 0.0, file = outFile)

    # Plot Various Fields

    print('POINT_DATA', IMAX*JMAX , file = outFile)

    for c in range(0,nVar):
        if (fieldTypes[c] == 0):
            print('SCALARS ', fieldNames[c],  ' double 1' , file = outFile)
            print('LOOKUP_TABLE default' , file = outFile)
            for j in range(0,JMAX):
                for i in range(0,IMAX):
                    print(sol[j,i,c], file = outFile)
        if (fieldTypes[c] == 1):
            print('VECTORS ', fieldNames[c],  ' double' , file = outFile)
            for j in range(0,JMAX):
                for i in range(0,IMAX):
                    print(sol[j,i,c], sol[j,i,c+1], 0.0, file = outFile)

##########################################################################################################################################
# Calculate magnitude of gradient
##########################################################################################################################################

def calcGrad(U):
    JMAX = U.shape[0]
    IMAX = U.shape[1]

    grad = np.zeros((JMAX, IMAX))

    for j in range(0, JMAX):
        for i in range(0, IMAX):
            if (j == 0):
                grad_y = U[j+1,i] - U[j,i]
            elif (j == JMAX-1):
                grad_y = U[j,i] - U[j-1,i]
            else:
                grad_y = 0.5*(U[j+1,i] - U[j-1,i])

            if (i == 0):
                grad_x = U[j,i+1] - U[j,i]
            elif (i == IMAX-1):
                grad_x = U[j,i] - U[j,i-1]
            else:
                grad_x = 0.5*(U[j,i+1] - U[j,i-1])

            grad[j,i] = np.hypot(grad_x, grad_y)

    return grad

##########################################################################################################################################
# VTS (XML) Reader 
##########################################################################################################################################

def vtsReader(fileName):
    # Read the source file.
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()  
    output = reader.GetOutput()

    # Get the Number of points in x, y, and z directions 
    IMAX = output.GetDimensions()[0]
    JMAX = output.GetDimensions()[1]
    KMAX = output.GetDimensions()[2]

    # Get the coordiantes of the points  
    grid_vtk = output.GetPoints().GetData()

    grid = vtk_to_numpy(grid_vtk)
    X = grid[:,0]; Y = grid[:,1]

    X = X.reshape((JMAX,IMAX))
    Y = Y.reshape((JMAX,IMAX))

    # Get solution arrays 
    sol_vtk = output.GetPointData()
    nVar = sol_vtk.GetNumberOfArrays()
    sol = np.zeros((JMAX,IMAX,nVar))

    for c in range(0,nVar):
        array = vtk_to_numpy(sol_vtk.GetArray(c))
        sol[:,:,c] = array.reshape((JMAX,IMAX))

    return X,Y,sol

##################################################################################################################################
# Variables to be plotted in the five equation model
##################################################################################################################################

def fiveEqnModel(fileName):
    Name = "5EqnModel"
    fieldNames = ["Density", "Velocity", "Velocity", "Pressure", "Vol-Fraction", "grad-RHO"]
    fieldTypes = [0, 1, 2, 0, 0, 0]
    X,Y,sol_raw = vtsReader(fileName)
    sol = np.zeros((X.shape[0],X.shape[1],6))
    sol[:,:,0] = sol_raw[:,:,5]*sol_raw[:,:,0] + (1.0 - sol_raw[:,:,5])*sol_raw[:,:,1]
    sol[:,:,1] = sol_raw[:,:,2]
    sol[:,:,2] = sol_raw[:,:,3]
    sol[:,:,3] = sol_raw[:,:,4]
    sol[:,:,4] = sol_raw[:,:,5]
    sol[:,:,5] = calcGrad(sol[:,:,0])
    return Name, fieldNames, fieldTypes, X, Y, sol

##################################################################################################################################
# Variables to be plotted in the four equation model
##################################################################################################################################

def fourEqnModel(fileName):
    Name = "5EqnModel"
    fieldNames = ["Density", "Velocity", "Velocity", "Pressure", "Vol-Fraction", "grad-RHO"]
    fieldTypes = [0, 1, 2, 0, 0, 0]
    X,Y,sol_raw = vtsReader(fileName)
    sol = np.zeros((X.shape[0],X.shape[1],6))
    sol[:,:,0] = sol_raw[:,:,0]
    sol[:,:,1] = sol_raw[:,:,1]
    sol[:,:,2] = sol_raw[:,:,2]
    sol[:,:,3] = sol_raw[:,:,3]
    sol[:,:,4] = sol_raw[:,:,4]
    sol[:,:,5] = calcGrad(sol[:,:,0])
    return Name, fieldNames, fieldTypes, X, Y, sol

##################################################################################################################################

for fileName in os.listdir('.'):
    if fileName.endswith('.vts'):
        print(fileName)
        Name, fieldNames, fieldTypes, X, Y, sol = fiveEqnModel(fileName)
        plotMatPlotLib(Name, fileName,fieldNames, fieldTypes, X, Y, sol)

##################################################################################################################################




