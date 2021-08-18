/*
 * post_process.c
 *      Author: sunder
 */ 

#include "hype.h"  

//----------------------------------------------------------------------------
// Monitor function for additional processing in the intermediate time steps 
//----------------------------------------------------------------------------

PetscErrorCode MonitorFunction (TS ts,PetscInt step, PetscReal time, Vec U, void *ctx) {

    PetscErrorCode ierr; 
    AppCtx *Ctx = (AppCtx*)ctx;

    // Set the time step based on CFL condition 

    //ierr = TSSetTimeStep(ts, Ctx->dt);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"%d t = %.5e\n", step, time);CHKERRQ(ierr);

    // Plot the solution at the required time interval 

    if (Ctx->WriteInterval != 0) {

        if(step%Ctx->WriteInterval == 0) {
            
            DM da;
            ierr = TSGetDM(ts,&da);CHKERRQ(ierr);
            
            ierr = ComputePrimitiveVariables(U, Ctx->W, da);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data at t = %.5e, step = %d\n", time, step);CHKERRQ(ierr);
            char filename[20];
            PetscViewer viewer;

            sprintf(filename, "sol-%08d.vts", step); // 8 is the padding level, increase it for longer simulations
            ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);

            ierr = DMView(da, viewer);
            VecView(Ctx->W, viewer);
            ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
            
        }
    }

    if (Ctx->RestartInterval != 0) {
        
        PetscViewer viewer_binary;

        if(step%Ctx->RestartInterval == 0) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data in binary to restart1.bin at t = %.7e\n", time);CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"restart1.bin",FILE_MODE_WRITE, &viewer_binary);CHKERRQ(ierr);
            ierr = VecView(U,viewer_binary);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer_binary); CHKERRQ(ierr);
        }
        
        if((step+10)%(Ctx->RestartInterval) == 0) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data in binary to restart2.bin at t = %.7e\n", time);CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"restart2.bin",FILE_MODE_WRITE, &viewer_binary);CHKERRQ(ierr);
            ierr = VecView(U,viewer_binary);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer_binary);CHKERRQ(ierr);
        }
    }
    
    return ierr; 
}

//----------------------------------------------------------------------------
// Find L2 and L_inf errors for periodic test cases with square domain and 
// where final solution conicides with the initial condition
//----------------------------------------------------------------------------

PetscErrorCode ErrorNorms(Vec U, DM da, AppCtx Ctx, PetscReal* l2, PetscReal* linf) {

    PetscErrorCode ierr;

    DM          coordDA;
    Vec         coordinates;
    DMDACoor2d  **coords;
    Field   **u;
    Field   **u_exact;
    PetscInt    xs, ys, xm, ym, i, j, c, l , m;
    PetscReal integral[nVar]; 
    PetscReal xc, yc, xGP, yGP;
    PetscReal h = Ctx.h;
    PetscReal Q0[nVar];
    Vec U_exact;

    const PetscInt N_gp5 = 5;
    const PetscReal x_gp5[] = { 0.0000000000000000, -0.2692346550528416,  0.2692346550528416, -0.4530899229693320,  0.4530899229693320};
    const PetscReal w_gp5[] = {0.28444444444444444, 0.23931433524968326, 0.23931433524968326, 0.11846344252809456, 0.11846344252809456};
    
    ierr = VecDuplicate(U,&U_exact);CHKERRQ(ierr);
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(da, &coordDA);CHKERRQ(ierr);
    ierr = DMGetCoordinates(da, &coordinates);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U_exact, &u_exact);CHKERRQ(ierr);

    // Use five point gauss quadrature

    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {
            
            // Get coordinates of center of the cell 
            
            xc = coords[j][i].x; 
            yc = coords[j][i].y;
            
            for (c = 0; c < nVar; ++c)
                integral[c] = 0.0;
            
            for(l = 0; l < N_gp5; ++l) {
                for (m = 0; m < N_gp5; ++m) {

                    xGP = xc + h*x_gp5[l];
                    yGP = yc + h*x_gp5[m];
                    
                    InitialCondition(xGP,yGP,Q0);
                    
                    for (c = 0; c < nVar; ++c) 
                        integral[c] += w_gp5[l]*w_gp5[m]*Q0[c];
                }
            }
            
            for (c = 0; c < nVar; ++c) {
                if (c == 5) {
                    u_exact[j][i].comp[c] = integral[c]; 
                }
                
                else {
                    u[j][i].comp[c] = 0.0;
                    u_exact[j][i].comp[c] = 0.0; 
                }
            }
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, U_exact, &u_exact);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);

    PetscReal nrm_inf, nrm2; 

    ierr = VecAXPY(U_exact, -1.0, U);CHKERRQ(ierr);
    ierr = VecNorm(U_exact, NORM_INFINITY, &nrm_inf);CHKERRQ(ierr);
    ierr = VecNorm(U_exact, NORM_2, &nrm2);CHKERRQ(ierr);

    nrm2 = nrm2/((PetscReal)Ctx.N_x);

    *l2 = nrm2; 
    *linf = nrm_inf; 

    ierr = VecDestroy(&U_exact);CHKERRQ(ierr);

    return ierr; 
}

//----------------------------------------------------------------------------
// Find the cell averages of primitive variables in each cell 
//----------------------------------------------------------------------------

PetscErrorCode ComputePrimitiveVariables(Vec U, Vec W, DM da) {
    
    PetscErrorCode ierr;           // For catching PETSc errors 
    PetscInt c,i,j,xs,ys,xm,ym;    // Corners of the grid on the given solution 
    Field  **u;                    // Local array of the conserved variables 
    Field  **w;                    // Local array of the primitive variables 
    PetscReal Q[nVar], V[nVar];  

    // Read the local solution to the array u  

    ierr = DMDAVecGetArrayRead(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, W, &w);CHKERRQ(ierr);

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);


    for (j=ys; j<ys+ym; ++j) {
        for (i=xs; i<xs+xm; ++i) {
        
            // Select stencils for each component of the solution 
    
            for (c = 0 ; c < nVar; ++c)
                Q[c] = u[j][i].comp[c];
    
            PDECons2Prim(Q, V);
            
            for (c = 0 ; c < nVar; ++c)
                w[j][i].comp[c] = V[c];
        }
    }  

    ierr = DMDAVecRestoreArrayRead(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,W,&w);CHKERRQ(ierr);

    return ierr; 
}

