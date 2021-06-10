/*
 * hype.c
 *      Author: sunder
 */ 

static char help[] = "Parallel Hyperbolic PDE solver upto fourth order accuracy.\n\n";

#include "hype.h" 

//----------------------------------------------------------------------------
// Main function of the code 
//----------------------------------------------------------------------------

int main(int argc,char **argv) {

    // --------------------------------------------
    // Initialize MPI 
    //---------------------------------------------

    PetscErrorCode ierr;                    /* For catching PETSc errors */ 
    PetscLogDouble start_time, end_time;    /* For logging the time values */

    ierr = PetscInitialize(&argc, &argv, (char*)0, help);CHKERRQ(ierr);
    ierr =  PetscTime(&start_time);CHKERRQ(ierr); 

    // --------------------------------------------
    // Set important user defined parameters  
    //---------------------------------------------

    AppCtx Ctx; 

    Ctx.N               = 1;
    Ctx.x_min           = 0.0;
    Ctx.x_max           = 10.0;
    Ctx.y_min           = -2.5;
    Ctx.y_max           = 2.5;
    Ctx.N_x             = 400;
    Ctx.N_y             = 200;
    Ctx.CFL             = 0.9;                            
    Ctx.InitialStep     = 0;
    Ctx.InitialTime     = 0.0;
    Ctx.FinalTime       = 0.02;
    Ctx.WriteInterval   = 1000;
    Ctx.RestartInterval = 500;
    Ctx.ReconsPrimitive = PETSC_FALSE;
    Ctx.OutFormat       = vtk;
    Ctx.LeftBoundary    = transmissive;
    Ctx.RightBoundary   = transmissive;
    Ctx.BottomBoundary  = transmissive;
    Ctx.TopBoundary     = transmissive;
    Ctx.Restart         = PETSC_FALSE; 

    //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // No need to change anything beyond this point
    //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Ctx.nDOF            = (Ctx.N+1)*(Ctx.N+2)/2;
    Ctx.h = (Ctx.x_max - Ctx.x_min)/(PetscReal)(Ctx.N_x);  

    if (Ctx.N < 2) {
        Ctx.Ngp_Face = 1;
        Ctx.Ngp_Vol = 1;
    }
    else {
        Ctx.Ngp_Face = 2;
        Ctx.Ngp_Vol = 4;
    }

    Ctx.N_node = 4*Ctx.Ngp_Face;

    // --------------------------------------------
    // Data members  
    //---------------------------------------------

    Vec U;                           // Solution Vector (Conserved variables)
    Vec RHS;                         // RHS vector to update the solution 
    DM da;                           // Grid object 
    PetscInt time_steps;             // No. of time steps 
    TS ts;                           // Time stepping object 
    PetscMPIInt MyPID;               // Rank of the current processor 
    PetscMPIInt numProcs;            // Size of the communicator

    // --------------------------------------------
    // Obtain the rank of the process and size of 
    // the communicator 
    //---------------------------------------------

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcs);CHKERRQ(ierr); 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Code running with %d processes\n", numProcs);CHKERRQ(ierr); 

    // --------------------------------------------
    // Initialize the grid and set field names
    //---------------------------------------------

    DMBoundaryType x_boundary;
    DMBoundaryType y_boundary;

    if (Ctx.LeftBoundary == periodic || Ctx.RightBoundary == periodic)
        x_boundary = DM_BOUNDARY_PERIODIC;
    else
        x_boundary = DM_BOUNDARY_GHOSTED; 

    if (Ctx.BottomBoundary == periodic || Ctx.TopBoundary == periodic)
        y_boundary = DM_BOUNDARY_PERIODIC;
    else
        y_boundary = DM_BOUNDARY_GHOSTED; 

    ierr = DMDACreate2d(PETSC_COMM_WORLD, // Global communicator      
                        x_boundary,       // Boundary conditions in x-direction 
                        y_boundary,       // Boundary conditions in y-direction
                        DMDA_STENCIL_BOX, // Stencil type (other is star type)
                        Ctx.N_x,          // No. of cells in x-direction 
                        Ctx.N_y,          // No. of cells in y-direction
                        PETSC_DECIDE,     // Domain decomposition in x-direction 
                        PETSC_DECIDE,     // Domain decomposition in y-direction
                        nVar,             // No. of dofs per cell 
                        3,                // Width of the stencil
                        NULL,
                        NULL,
                        &da);CHKERRQ(ierr); // da object 

    ierr = DMSetUp(da);CHKERRQ(ierr);

    // Now create various global vectors 

    ierr = DMCreateGlobalVector(da, &U);CHKERRQ(ierr);
    ierr = VecDuplicate(U,&Ctx.W);CHKERRQ(ierr);
    ierr = VecDuplicate(U,&RHS);CHKERRQ(ierr);

    // Set coordinates of cell centers 

    ierr = DMDASetUniformCoordinates(da,
                                        Ctx.x_min + 0.5*Ctx.h, Ctx.x_max + 0.5*Ctx.h,
                                        Ctx.y_min + 0.5*Ctx.h, Ctx.y_max + 0.5*Ctx.h,
                                        0.0,0.0);CHKERRQ(ierr);
    // Set names of the fields
                                        
    ierr = PetscObjectSetName((PetscObject)U,"cons");CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)Ctx.W,"sol");CHKERRQ(ierr);

#ifdef EULER
    ierr = DMDASetFieldName(da,0,"density");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,1,"x-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,2,"y-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,3,"pressure");CHKERRQ(ierr);
#endif

#ifdef EFFECTIVE_GAMMA
    ierr = DMDASetFieldName(da,0,"density");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,1,"x-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,2,"y-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,3,"pressure");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,4,"phi");CHKERRQ(ierr);
#endif

#ifdef BAER_NUNZIATO
    ierr = DMDASetFieldName(da,0,"density_s");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,1,"x-velocity_s");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,2,"y-velocity_s");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,3,"pressure_s");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,4,"density_g");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,5,"x-velocity_g");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,6,"y-velocity_g");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,7,"pressure_g");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,8,"phi_s");CHKERRQ(ierr);
#endif

#ifdef KAPILA_5EQN
    ierr = DMDASetFieldName(da,0,"density1");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,1,"density2");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,2,"x-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,3,"y-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,4,"pressure");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,5,"phi");CHKERRQ(ierr);
#endif

#ifdef DIFFUSE_INTERFACE
    ierr = DMDASetFieldName(da,0,"Density");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,1,"x-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,2,"y-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,3,"Pressure");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,4,"alpha");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,5,"x-velocity_s");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,6,"y-velocity_s");CHKERRQ(ierr);
#endif

    // --------------------------------------------
    // Allocate memory for boundary values and 
    // upwind fluxes
    //---------------------------------------------

    // Set quadrature points

    ierr = PetscMalloc1(Ctx.Ngp_Face,&Ctx.xGP_Face);CHKERRQ(ierr);
    ierr = PetscMalloc1(Ctx.Ngp_Face,&Ctx.wGP_Face);CHKERRQ(ierr);
    ierr = PetscMalloc1(Ctx.Ngp_Vol,&Ctx.xGP_Vol);CHKERRQ(ierr);
    ierr = PetscMalloc1(Ctx.Ngp_Vol,&Ctx.yGP_Vol);CHKERRQ(ierr);
    ierr = PetscMalloc1(Ctx.Ngp_Vol,&Ctx.wGP_Vol);CHKERRQ(ierr);
    ierr = PetscMalloc1(Ctx.N_node,&Ctx.x_node);CHKERRQ(ierr);
    ierr = PetscMalloc1(Ctx.N_node,&Ctx.y_node);CHKERRQ(ierr);


    if (Ctx.N < 2) {
        Ctx.xGP_Face[0] = 0.0; Ctx.wGP_Face[0] = 1.0;
        Ctx.xGP_Vol[0] = 0.0; Ctx.yGP_Vol[0] = 0.0; Ctx.wGP_Vol[0] = 1.0;

        Ctx.x_node[0] =  0.5; Ctx.y_node[0] = 0.0;
        Ctx.x_node[1] = -0.5; Ctx.y_node[1] = 0.0;
        Ctx.x_node[2] =  0.0; Ctx.y_node[2] = 0.5;
        Ctx.x_node[3] =  0.0; Ctx.y_node[3] =-0.5;
    }

    else {
        Ctx.xGP_Face[0] = -0.28867513459481287; Ctx.wGP_Face[0] = 0.5;
        Ctx.xGP_Face[1] =  0.28867513459481287; Ctx.wGP_Face[1] = 0.5;

        Ctx.xGP_Vol[0] = -0.28867513459481287; Ctx.yGP_Vol[0] = -0.28867513459481287; Ctx.wGP_Vol[0] = 0.25;
        Ctx.xGP_Vol[1] = -0.28867513459481287; Ctx.yGP_Vol[1] =  0.28867513459481287; Ctx.wGP_Vol[1] = 0.25;
        Ctx.xGP_Vol[2] =  0.28867513459481287; Ctx.yGP_Vol[2] = -0.28867513459481287; Ctx.wGP_Vol[2] = 0.25;
        Ctx.xGP_Vol[3] =  0.28867513459481287; Ctx.yGP_Vol[3] =  0.28867513459481287; Ctx.wGP_Vol[3] = 0.25;

        Ctx.x_node[0] =  0.5;                 Ctx.y_node[0] = -0.28867513459481287;
        Ctx.x_node[1] =  0.5;                 Ctx.y_node[1] =  0.28867513459481287;
        Ctx.x_node[2] = -0.5;                 Ctx.y_node[2] = -0.28867513459481287;
        Ctx.x_node[3] = -0.5;                 Ctx.y_node[3] =  0.28867513459481287;
        Ctx.x_node[4] = -0.28867513459481287; Ctx.y_node[4] =  0.5;
        Ctx.x_node[5] =  0.28867513459481287; Ctx.y_node[5] =  0.5;
        Ctx.x_node[6] = -0.28867513459481287; Ctx.y_node[6] = -0.5;
        Ctx.x_node[7] =  0.28867513459481287; Ctx.y_node[7] = -0.5;
    }

    PetscInt xs,ys,xm,ym,q,k;
    PetscReal value, grad_x, grad_y;

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    Ctx.u_bnd        = allocate5d(ym+2, xm+2, nVar, 4, Ctx.Ngp_Face);  // 4->number of faces in a cell
#ifdef VISCOUS
    Ctx.u_bnd_grad   = allocate6d(ym+2, xm+2, nVar, 4, Ctx.Ngp_Face, DIM); // 4->number of faces, 2->number of quadraure points
#endif
    Ctx.F            = allocate3d(ym, xm+1, nVar);
    Ctx.G            = allocate3d(ym+1, xm, nVar);
    Ctx.D            = allocate3d(ym, xm+1, nVar);
    Ctx.E            = allocate3d(ym+1, xm, nVar);
    Ctx.phiFace      = allocate3d(4, Ctx.Ngp_Face, Ctx.nDOF);     // 4->number of faces in a cell
    Ctx.phiNode      = allocate2d(Ctx.N_node, Ctx.nDOF);
    Ctx.phiVol       = allocate2d(Ctx.Ngp_Vol, Ctx.nDOF);
    Ctx.gradphiVol_x = allocate2d(Ctx.Ngp_Vol, Ctx.nDOF);
    Ctx.gradphiVol_y = allocate2d(Ctx.Ngp_Vol, Ctx.nDOF);
    Ctx.gradphiFace_x = allocate3d(4, Ctx.Ngp_Face,Ctx.nDOF); // 4 -> Number of faces in a cell
    Ctx.gradphiFace_y = allocate3d(4, Ctx.Ngp_Face,Ctx.nDOF); // 4 -> Number of faces in a cell
    

    // Find the value of basis functions and gradients at face quadrature points

    for (q = 0; q < Ctx.Ngp_Face; ++q) {
        for (k = 0; k < Ctx.nDOF; ++k) {

            // Left face
            value = basis(-0.5, Ctx.xGP_Face[q], k);
            basis_grad(-0.5, Ctx.xGP_Face[q], k, &grad_x, &grad_y);
            set_element_3d(Ctx.phiFace, 0, q, k, value);
            set_element_3d(Ctx.gradphiFace_x, 0, q, k, grad_x);
            set_element_3d(Ctx.gradphiFace_y, 0, q, k, grad_y);

            // Right face
            value = basis( 0.5, Ctx.xGP_Face[q], k);
            basis_grad( 0.5, Ctx.xGP_Face[q], k, &grad_x, &grad_y);
            set_element_3d(Ctx.phiFace, 1, q, k, value);
            set_element_3d(Ctx.gradphiFace_x, 1, q, k, grad_x);
            set_element_3d(Ctx.gradphiFace_y, 1, q, k, grad_y);

            // Bottom face
            value = basis(Ctx.xGP_Face[q], -0.5, k);
            basis_grad(Ctx.xGP_Face[q], -0.5, k, &grad_x, &grad_y);
            set_element_3d(Ctx.phiFace, 2, q, k, value);
            set_element_3d(Ctx.gradphiFace_x, 2, q, k, grad_x);
            set_element_3d(Ctx.gradphiFace_y, 2, q, k, grad_y);

            // Top face
            value = basis(Ctx.xGP_Face[q], 0.5, k);
            basis_grad(Ctx.xGP_Face[q], 0.5, k, &grad_x, &grad_y);
            set_element_3d(Ctx.phiFace, 3, q, k, value);
            set_element_3d(Ctx.gradphiFace_x, 3, q, k, grad_x);
            set_element_3d(Ctx.gradphiFace_y, 3, q, k, grad_y);
        }
    }
    
    // Find the value of basis functions on interior nodes 
    
    for (q = 0; q < Ctx.N_node; ++q) {
        for (k = 0; k < Ctx.nDOF; ++k) {
            set_element_2d(Ctx.phiNode, q, k, basis(Ctx.x_node[q], Ctx.y_node[q], k));
        }
    }
    
    // Find the value of basis functions and gradients on interior quadrature points 
    
    for (q = 0; q < Ctx.Ngp_Vol; ++q) {
        for (k = 0; k < Ctx.nDOF; ++k) {
            set_element_2d(Ctx.phiVol, q, k, basis(Ctx.xGP_Vol[q], Ctx.yGP_Vol[q], k));
            basis_grad(Ctx.xGP_Vol[q], Ctx.yGP_Vol[q], k, &grad_x, &grad_y);
            set_element_2d(Ctx.gradphiVol_x, q, k, grad_x);
            set_element_2d(Ctx.gradphiVol_y, q, k, grad_y);
        }
    }

    

    ierr = DMCreateLocalVector(da,&Ctx.localU);CHKERRQ(ierr);

    // --------------------------------------------
    // Initialize the solution (either with initial
    // condition or restart file)
    //---------------------------------------------

    if (Ctx.Restart) {
        
        // Initialize by reading the restart file 
        
        PetscViewer    viewer_binary;
        ierr = PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from restart1.bin ...\n");CHKERRQ(ierr);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"restart1.bin",FILE_MODE_READ,&viewer_binary);CHKERRQ(ierr);
        ierr = VecLoad(U,viewer_binary);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer_binary);CHKERRQ(ierr);
    }

    else {
        
        // Initialize by initial condition 
        ierr = InitializeSolution(U, da, Ctx);CHKERRQ(ierr);
    }

    // --------------------------------------------
    // Advance solution in time   
    //---------------------------------------------
    
    ierr = TSCreate(PETSC_COMM_SELF, &ts);CHKERRQ(ierr);              
    ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);           
    ierr = TSSetDM(ts,da);CHKERRQ(ierr);                              
    
    if (Ctx.ReconsPrimitive) {
        
        ierr = RHSFunctionPrim(ts, Ctx.InitialTime, U, RHS, &Ctx);CHKERRQ(ierr);
        ierr = TSSetRHSFunction(ts,NULL,RHSFunctionPrim, &Ctx);CHKERRQ(ierr);
    
    }
    
    else {
        
        ierr = RHSFunction(ts, Ctx.InitialTime, U, RHS, &Ctx);CHKERRQ(ierr);
        ierr = TSSetRHSFunction(ts,NULL,RHSFunction, &Ctx);CHKERRQ(ierr);
    }
    
    ierr = TSSetStepNumber(ts,Ctx.InitialStep);
    ierr = TSSetTime(ts, Ctx.InitialTime);
    ierr = TSSetTimeStep(ts, Ctx.dt);CHKERRQ(ierr);                     
    ierr = TSMonitorSet(ts,MonitorFunction,&Ctx,NULL);CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts, Ctx.FinalTime);CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts);CHKERRQ(ierr); 

    if (Ctx.N == 0) {
        ierr = TSSetType(ts,TSEULER);CHKERRQ(ierr);
    }

    else {
        ierr = TSSetType(ts, TSSSP);CHKERRQ(ierr);

        if (Ctx.N == 1) {
            ierr = TSSSPSetType(ts, TSSSPRKS2);CHKERRQ(ierr);
            ierr = TSSSPSetNumStages(ts,2);CHKERRQ(ierr);
        }

        else if (Ctx.N == 2) {
            ierr = TSSSPSetType(ts, TSSSPRKS3);CHKERRQ(ierr);
            ierr = TSSSPSetNumStages(ts,4);CHKERRQ(ierr);
        }

        else {
            ierr = TSSSPSetType(ts, TSSSPRK104);CHKERRQ(ierr);
        }
    }

    ierr = TSSolve(ts, U);CHKERRQ(ierr);
    ierr = TSGetStepNumber(ts,&time_steps);CHKERRQ(ierr); 
    
    // --------------------------------------------
    // Output solution in vtk format   
    //--------------------------------------------
    
    ierr = ComputePrimitiveVariables(U, Ctx.W, da);CHKERRQ(ierr);
    char filename[20]; 
    PetscViewer viewer;

    if (Ctx.OutFormat == vts) {
        sprintf(filename, "sol-%08d.vts", time_steps);
        ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
    }

    if (Ctx.OutFormat == vtk) {
        sprintf(filename, "sol-%08d.vtk", time_steps);
        ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);CHKERRQ(ierr);
        ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK);CHKERRQ(ierr);
    }

    ierr = DMView(da, viewer);CHKERRQ(ierr);
    ierr = VecView(Ctx.W, viewer);CHKERRQ(ierr);
    
    // --------------------------------------------
    // Get the norms of errors (only for periodic
    // test cases)
    //---------------------------------------------

    PetscReal nrm_2, nrm_inf;
    ierr = ErrorNorms(U, da, Ctx, &nrm_2, &nrm_inf);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Norm2 = %.7e, NormMax = %.7e\n", nrm_2, nrm_inf);CHKERRQ(ierr);

    // --------------------------------------------
    // Print the time taken for simulation       
    //---------------------------------------------

    ierr =  PetscTime(&end_time);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time taken =  %g\n",(double)(end_time - start_time));CHKERRQ(ierr);

    // --------------------------------------------
    // Free all the memory, finalize MPI and exit   
    //---------------------------------------------

    ierr = VecDestroy(&U);CHKERRQ(ierr);
    ierr = VecDestroy(&RHS);CHKERRQ(ierr);
    ierr = DMDestroy(&da);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = TSDestroy(&ts);CHKERRQ(ierr);
    ierr = VecDestroy(&Ctx.localU);CHKERRQ(ierr);
    ierr = VecDestroy(&Ctx.W);CHKERRQ(ierr);

    free5d(Ctx.u_bnd);
    free3d(Ctx.F);
    free3d(Ctx.G); 
    free3d(Ctx.D);
    free3d(Ctx.E);
    free3d(Ctx.phiFace);
    free2d(Ctx.phiNode);
    free2d(Ctx.phiVol);
    free2d(Ctx.gradphiVol_x);
    free2d(Ctx.gradphiVol_y);
    free3d(Ctx.gradphiFace_x);
    free3d(Ctx.gradphiFace_y);
#ifdef VISCOUS
    free6d(Ctx.u_bnd_grad);
#endif

    ierr = PetscFree(Ctx.xGP_Face);CHKERRQ(ierr);
    ierr = PetscFree(Ctx.wGP_Face);CHKERRQ(ierr);
    ierr = PetscFree(Ctx.xGP_Vol);CHKERRQ(ierr);
    ierr = PetscFree(Ctx.yGP_Vol);CHKERRQ(ierr);
    ierr = PetscFree(Ctx.wGP_Vol);CHKERRQ(ierr);
    ierr = PetscFree(Ctx.x_node);CHKERRQ(ierr);
    ierr = PetscFree(Ctx.y_node);CHKERRQ(ierr);

    ierr = PetscFinalize();CHKERRQ(ierr);

    return ierr; 
}
