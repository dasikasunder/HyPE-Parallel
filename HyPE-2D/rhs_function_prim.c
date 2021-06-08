/*
 * rhs_function_prim.c
 *      Author: sunder
 */ 
#include "hype.h" 

//----------------------------------------------------------------------------
// Compute the value of RHS for each cell in the domain using 
// primitive variables
//----------------------------------------------------------------------------

PetscErrorCode RHSFunctionPrim(TS ts, PetscReal t, Vec U, Vec RHS, void* ctx) {

    PetscErrorCode ierr;           
    AppCtx *Ctx = (AppCtx*)ctx; 
    DM da;                         
    PetscInt c,i,j,k,f,q,xs,ys,xm,ym,xs_g,ys_g,xm_g,ym_g,oned_begin,oned_end,irhs;                
    Field   **w;                   
    Field   **rhs;                 
    PetscReal s_c, s_max_c = 0.0;
    PetscReal w_x_loc[s_width], w_y_loc[s_width], w_xy_loc[s_width];  
    PetscReal dt; 
    PetscReal *coeffs, *u_coeffs;
    PetscReal **sol;
    PetscReal value, grad_x, grad_y; 
    PetscReal r1_h = 1./(Ctx->h); 
    PetscReal x_loc, y_loc, nx, ny; 
    PetscInt local_i, local_j;
    PetscBool PAD; 
    
    PetscReal **Vgp, **grad_Vgp_x, **grad_Vgp_y, **VNode;
    PetscReal V[nVar], gradV_x[nVar], gradV_y[nVar], BgradQ[nVar];
    PetscReal VL[nVar], VR[nVar];
    PetscReal Flux_conv[nVar], Flux[nVar], Fluc[nVar], NCP[nVar]; 

    ierr = TSGetDM(ts,&da);CHKERRQ(ierr);
    
    // Calculate primitive variables in each cell 
    
    ierr = ComputePrimitiveVariables(U, Ctx->W, da);CHKERRQ(ierr);
    
    // Scatter global->local to have access to the required ghost values 

    ierr = DMGlobalToLocalBegin(da, Ctx->W, INSERT_VALUES, Ctx->localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, Ctx->W, INSERT_VALUES,   Ctx->localU);CHKERRQ(ierr);

    // Read the local solution to the array u  

    ierr = DMDAVecGetArrayRead(da, Ctx->localU, &w);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, RHS, &rhs);CHKERRQ(ierr);

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners(da, &xs_g, &ys_g, NULL, &xm_g, &ym_g, NULL);CHKERRQ(ierr);

#ifdef VISCOUS
    PetscInt iDim;
    PetscReal s_v, grad_VL[nVar][DIM]; PetscReal grad_VR[nVar][DIM], s_max_v = 0.0;
    PetscReal Flux_visc[nVar];
#endif

    // Allocate memory

    Malloc1D(&coeffs,Ctx->nDOF);
    Malloc1D(&u_coeffs,Ctx->nDOF);
    Malloc2D(&sol, Ctx->nDOF, nVar);
    Malloc2D(&Vgp,Ctx->Ngp_Vol,nVar);
    Malloc2D(&grad_Vgp_x,Ctx->Ngp_Vol,nVar);
    Malloc2D(&grad_Vgp_y,Ctx->Ngp_Vol,nVar);
    Malloc2D(&VNode,Ctx->N_node,nVar);

    //--------------------------------------------------------------
    // Apply Boundary Conditions 
    //--------------------------------------------------------------

    oned_begin = 0; oned_end = Ctx->N_x-1; 

    for (j = ys_g; j < ys_g+ym_g; ++j) {
        for (i = xs_g; i < xs_g+xm_g; ++i) {
            
            if (i < 0)  { // Left boundary 
                
                // Transmissive/Outflow Boundary 
                
                if (Ctx->LeftBoundary == transmissive) {
                    
                    irhs = oned_begin; 
                
                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = w[j][irhs].comp[c];
                }
                
                // Reflective Boundary 
                
                if (Ctx->LeftBoundary == reflective) {
                    
                    irhs = oned_begin - 1 - i; 
                
                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = w[j][irhs].comp[c];

#ifdef EULER
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component
#ifdef VISCOUS
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component
#endif
#endif

#ifdef EFFECTIVE_GAMMA
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component
#ifdef VISCOUS
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component
#endif
#endif

#ifdef BAER_NUNZIATO
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component (Phase 1)
                    w[j][i].comp[5] = -w[j][i].comp[5]; // Reflect x-velocity component (Phase 2)
#ifdef VISCOUS
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component (Phase 1)
                    w[j][i].comp[6] = -w[j][i].comp[6]; // Reflect y-velocity component (Phase 2)
#endif
#endif

#ifdef KAPILA_5EQN
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect x-velocity component
#ifdef VISCOUS
                    w[j][i].comp[3] = -w[j][i].comp[3]; // Reflect y-velocity component
#endif
#endif

                }

                // Inlet boundary

                if (Ctx->LeftBoundary == inlet) {

                    x_loc = Ctx->x_min + ((PetscReal)(i) + 0.5)*Ctx->h;
                    y_loc = Ctx->y_min + ((PetscReal)(j) + 0.5)*Ctx->h;

                    InletBCPrim(x_loc,y_loc,t,VR);

                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = VR[c];
                }
                
            }
            
            if (i >= Ctx->N_x) { // Right Boundary 
                
                // Transmissive/Outflow Boundary 
                
                if (Ctx->RightBoundary == transmissive) {
                    
                    irhs = oned_end; 
                    
                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = w[j][irhs].comp[c];
                }
                
                // Reflective Boundary
                
                if (Ctx->RightBoundary == reflective) {
                    
                    irhs = 2*oned_end - i + 1; 
                    
                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = w[j][irhs].comp[c];

#ifdef EULER
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component
#ifdef VISCOUS
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component
#endif
#endif

#ifdef EFFECTIVE_GAMMA
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component
#ifdef VISCOUS
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component
#endif
#endif

#ifdef BAER_NUNZIATO
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component (Phase 1)
                    w[j][i].comp[5] = -w[j][i].comp[5]; // Reflect x-velocity component (Phase 2)
#ifdef VISCOUS
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component (Phase 1)
                    w[j][i].comp[6] = -w[j][i].comp[6]; // Reflect y-velocity component (Phase 2)
#endif
#endif

#ifdef KAPILA_5EQN
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect x-velocity component
#ifdef VISCOUS
                    w[j][i].comp[3] = -w[j][i].comp[3]; // Reflect y-velocity component
#endif
#endif

                }
            }
        }
    }

    oned_begin = 0; oned_end = Ctx->N_y-1;

    for (j = ys_g; j < ys_g+ym_g; ++j) {
        for (i = xs_g; i < xs_g+xm_g; ++i) {
            
            if (j < 0) { // Bottom boundary 
                
                // Transmissive/Outflow Boundary
                
                if (Ctx->BottomBoundary == transmissive) {
                
                    irhs = oned_begin;
                    
                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = w[irhs][i].comp[c]; 
                }
                
                // Reflective Boundary
                
                if (Ctx->BottomBoundary == reflective) {
                
                    irhs = oned_begin - 1 - j;
                    
                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = w[irhs][i].comp[c]; 

#ifdef EULER
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component
#ifdef VISCOUS
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component
#endif
#endif

#ifdef EFFECTIVE_GAMMA
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component
#ifdef VISCOUS
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component
#endif
#endif

#ifdef BAER_NUNZIATO
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component (Phase 1)
                    w[j][i].comp[6] = -w[j][i].comp[6]; // Reflect y-velocity component (Phase 2)
#ifdef VISCOUS
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component (Phase 1)
                    w[j][i].comp[5] = -w[j][i].comp[5]; // Reflect x-velocity component (Phase 2)
#endif
#endif

#ifdef KAPILA_5EQN
                    w[j][i].comp[3] = -w[j][i].comp[3]; // Reflect y-velocity component
#ifdef VISCOUS
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect x-velocity component
#endif
#endif
                }
            }
            
            if (j >= Ctx->N_y) { // Top boundary 
                
                // Transmissive/Outflow Boundary
                
                if (Ctx->TopBoundary == transmissive) {
                
                    irhs = oned_end;
                    
                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = w[irhs][i].comp[c];
                }
                
                // Reflective Boundary
                
                if (Ctx->TopBoundary == reflective) {
                
                    irhs = 2*oned_end - j + 1;;
                    
                    for (c = 0; c < nVar; ++c)
                        w[j][i].comp[c] = w[irhs][i].comp[c];

#ifdef EULER
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect x-velocity component
#ifdef VISCOUS
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component
#endif
#endif

#ifdef EFFECTIVE_GAMMA
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component
#ifdef VISCOUS
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component
#endif
#endif

#ifdef BAER_NUNZIATO
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect y-velocity component (Phase 1)
                    w[j][i].comp[6] = -w[j][i].comp[6]; // Reflect y-velocity component (Phase 2)
#ifdef VISCOUS
                    w[j][i].comp[1] = -w[j][i].comp[1]; // Reflect x-velocity component (Phase 1)
                    w[j][i].comp[5] = -w[j][i].comp[5]; // Reflect x-velocity component (Phase 2)
#endif
#endif

#ifdef KAPILA_5EQN
                    w[j][i].comp[3] = -w[j][i].comp[3]; // Reflect y-velocity component
#ifdef VISCOUS
                    w[j][i].comp[2] = -w[j][i].comp[2]; // Reflect x-velocity component
#endif
#endif
                }
            }
        }
    }

    //--------------------------------------------------------------
    // Do WENO reconstruction for each cell
    //--------------------------------------------------------------
    
    for (j=ys-1; j<ys+ym+1; j++) {
        for (i=xs-1; i<xs+xm+1; i++) {
            
            local_j = j - (ys-1); 
            local_i = i - (xs-1);
            PAD = PETSC_TRUE; 
            
            // Select stencils for each component of the solution 
        
            for (c = 0 ; c < nVar; ++c) {
                
                w_x_loc[0] = w[j][i-2].comp[c]; 
                w_x_loc[1] = w[j][i-1].comp[c]; 
                w_x_loc[2] = w[j][i].comp[c]; 
                w_x_loc[3] = w[j][i+1].comp[c]; 
                w_x_loc[4] = w[j][i+2].comp[c];
                
                w_y_loc[0] = w[j-2][i].comp[c]; 
                w_y_loc[1] = w[j-1][i].comp[c]; 
                w_y_loc[2] = w[j][i].comp[c]; 
                w_y_loc[3] = w[j+1][i].comp[c]; 
                w_y_loc[4] = w[j+2][i].comp[c];

                w_xy_loc[0] = w[j][i].comp[c];
                w_xy_loc[1] = w[j+1][i+1].comp[c];
                w_xy_loc[2] = w[j-1][i+1].comp[c];
                w_xy_loc[3] = w[j+1][i-1].comp[c];
                w_xy_loc[4] = w[j-1][i-1].comp[c];
                
                Reconstruct(w_x_loc, w_y_loc, w_xy_loc, Ctx->N, coeffs, u_coeffs);

#ifdef VISCOUS
                // Calculate boundary extrpolated gradients

                for (f = 0; f < 4; ++f) {

                    for (q = 0; q < Ctx->Ngp_Face; ++q) {

                        grad_x = 0.0; grad_y = 0.0;

                        for (k = 0; k < Ctx->nDOF; ++k) {
                            grad_x += u_coeffs[k]*get_element_3d(Ctx->gradphiFace_x,f,q,k);
                            grad_y += u_coeffs[k]*get_element_3d(Ctx->gradphiFace_y,f,q,k);
                        }

                        grad_x = r1_h*grad_x; grad_y = r1_h*grad_y;

                        set_element_6d(Ctx->u_bnd_grad, local_j, local_i, c, f, q, 0, grad_x);
                        set_element_6d(Ctx->u_bnd_grad, local_j, local_i, c, f, q, 1, grad_y);
                    }
                }
#endif
                
                // Store the coefficients 
                
                for (k = 0; k < Ctx->nDOF; ++k) {
                    sol[k][c] = coeffs[k]; 
                }
                
                // Evaluate solution at various nodes to check physical admissibility of the solution  
                
                for (q = 0; q < Ctx->N_node; ++q) {
                    
                    VNode[q][c] = 0.0;
                    
                    for (k = 0; k < Ctx->nDOF; ++k)
                        VNode[q][c] += coeffs[k]*get_element_2d(Ctx->phiNode,q,k);
                }
            }
            
            // Check physical admissibility of the solution  and reduce to TVD if necessary 
            
            for (q = 0; q < Ctx->N_node; ++q) {
                
                for (c = 0; c < nVar; ++c) {
                    V[c] = VNode[q][c];
                }
                
                PAD = PDECheckPADPrim(V);
                
                    
                if (PAD == PETSC_FALSE)
                    break; 
            }
            
            
            if (PAD == PETSC_FALSE) {
                
                
                for (c = 0 ; c < nVar; ++c) {
                    
                    sol[0][c] = w[j][i].comp[c];
                    sol[1][c] = minmod(w[j][i+1].comp[c]-w[j][i].comp[c],w[j][i].comp[c]-w[j][i-1].comp[c]);
                    sol[2][c] = minmod(w[j+1][i].comp[c]-w[j][i].comp[c],w[j][i].comp[c]-w[j-1][i].comp[c]);
                    
                    for (k = 3; k < Ctx->nDOF; ++k) {
                        sol[k][c] = 0.0; 
                    }
                }
            
            }
            
            // Find the values of conserved variables at face quadrature points 
            
            for (c = 0 ; c < nVar; ++c) {
            
                // Get coefficients 
                
                for (f = 0; f < 4; ++f) {
                    for (q = 0; q < Ctx->Ngp_Face; ++q) {
                        
                        value = 0.0; 
                        
                        for (k = 0; k < Ctx->nDOF; ++k)
                            value += sol[k][c]*get_element_3d(Ctx->phiFace,f,q,k);
                        
                        set_element_5d(Ctx->u_bnd, local_j, local_i, c, f, q, value);
                    
                    }
                }
                
                for (q = 0; q < Ctx->Ngp_Vol; ++q) {
                    
                    value = 0.0; grad_x = 0.0; grad_y = 0.0; 
                    
                    for (k = 0; k < Ctx->nDOF; ++k) {
                        value  += sol[k][c]*get_element_2d(Ctx->phiVol,q,k);
                        grad_x += sol[k][c]*get_element_2d(Ctx->gradphiVol_x,q,k);
                        grad_y += sol[k][c]*get_element_2d(Ctx->gradphiVol_y,q,k);
                    }
                    
                    Vgp[q][c] = value;
                    grad_Vgp_x[q][c] = r1_h*grad_x; 
                    grad_Vgp_y[q][c] = r1_h*grad_y;
                }
            }
        

            // Add smooth part of the non-conservative product 
            
            if (i >= xs && i < xs+xm && j >= ys && j < ys+ym) {
                
                for (c = 0 ; c < nVar; ++c) 
                    rhs[j][i].comp[c] = 0.0; 
                
                for (q = 0; q < Ctx->Ngp_Vol; ++q) {
                
                    for (c = 0 ; c < nVar; ++c) {
                        V[c] = Vgp[q][c];
                        gradV_x[c] = grad_Vgp_x[q][c];
                        gradV_y[c] = grad_Vgp_y[q][c];
                    }
                    
                    PDENCPPrim(V, gradV_x, gradV_y, BgradQ);
                    
                    for (c = 0 ; c < nVar; ++c) 
                        rhs[j][i].comp[c] += -Ctx->wGP_Vol[q]*BgradQ[c];
                }
            }
        }
    } // End of cell loop
    
    // Find the upwind flux on each face in x-dirction 
    
    nx = 1.0; ny = 0.0; 

    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm+1; ++i) {
            
            local_j = j - (ys-1); 
            local_i = i - (xs-1);
            
            x_loc = Ctx->x_min + (PetscReal)(i)*Ctx->h;
            y_loc = Ctx->y_min + (PetscReal)(j)*Ctx->h;
            
            for (c = 0; c < nVar; ++c) {
                Flux[c] = 0.0; 
                NCP[c] = 0.0; 
            }

            for (q = 0; q < Ctx->Ngp_Face; ++q) {
                
                for (c = 0; c < nVar; ++c) {
                    VL[c] = get_element_5d(Ctx->u_bnd, local_j, local_i-1, c, 1, q);
                    VR[c] = get_element_5d(Ctx->u_bnd, local_j, local_i,   c, 0, q);
#ifdef VISCOUS
                    for (iDim = 0; iDim < DIM; ++iDim) {
                        grad_VL[c][iDim] = get_element_6d(Ctx->u_bnd_grad, local_j, local_i-1, c, 1, q, iDim);
                        grad_VR[c][iDim] = get_element_6d(Ctx->u_bnd_grad, local_j, local_i,   c, 0, q, iDim);
                    }
#endif
                }
                
                s_c = RiemannSolverPrim(VL, VR, nx, ny, x_loc, y_loc, Flux_conv, Fluc); if (s_c>s_max_c) s_max_c = s_c;
#ifdef VISCOUS
                s_v = ViscRiemannSolverPrim(VL, grad_VL, VR, grad_VR, nx, ny, Flux_visc);  if (s_v>s_max_v) s_max_v = s_v;
#endif
                
                for (c = 0; c < nVar; ++c) {
                    Flux[c] += Ctx->wGP_Face[q]*Flux_conv[c];
                    NCP[c]  += Ctx->wGP_Face[q]*Fluc[c];
#ifdef VISCOUS
                    Flux[c] += Ctx->wGP_Face[q]*Flux_visc[c];
#endif
                }
            }
        
            for (c = 0; c < nVar; ++c) {
                set_element_3d(Ctx->F, j-ys, i-xs, c, Flux[c]);
                set_element_3d(Ctx->D, j-ys, i-xs, c, NCP[c]);
            }
        }
    }
            
    // Find the upwind flux on each face in y-dirction 
    
    nx = 0.0; ny = 1.0;

    for (j = ys; j < ys+ym+1; ++j) {
        for (i = xs; i < xs+xm; ++i) {
            
            local_j = j - (ys-1); 
            local_i = i - (xs-1); 
            
            x_loc = Ctx->x_min + (PetscReal)(i)*Ctx->h;
            y_loc = Ctx->y_min + (PetscReal)(j)*Ctx->h;
            
            for (c = 0; c < nVar; ++c) {
                Flux[c] = 0.0; 
                NCP[c] = 0.0; 
            }

            for (q = 0; q < Ctx->Ngp_Face; ++q) {
                
                for (c = 0; c < nVar; ++c) {
                    VL[c] = get_element_5d(Ctx->u_bnd, local_j-1, local_i, c, 3, q);
                    VR[c] = get_element_5d(Ctx->u_bnd, local_j,   local_i, c, 2, q);
#ifdef VISCOUS
                    for (iDim = 0; iDim < DIM; ++iDim) {
                        grad_VL[c][iDim] = get_element_6d(Ctx->u_bnd_grad, local_j-1, local_i, c, 3, q, iDim);
                        grad_VR[c][iDim] = get_element_6d(Ctx->u_bnd_grad, local_j,   local_i, c, 2, q, iDim);
                    }
#endif
                }
                
                s_c = RiemannSolverPrim(VL, VR, nx, ny, x_loc, y_loc, Flux_conv, Fluc); if (s_c>s_max_c) s_max_c = s_c;
#ifdef VISCOUS
                s_v = ViscRiemannSolverPrim(VL, grad_VL, VR, grad_VR, nx, ny, Flux_visc);  if (s_v>s_max_v) s_max_v = s_v;
#endif
                
                for (c = 0; c < nVar; ++c) {
                    Flux[c] += Ctx->wGP_Face[q]*Flux_conv[c];
                    NCP[c]  += Ctx->wGP_Face[q]*Fluc[c];
#ifdef VISCOUS
                    Flux[c] += Ctx->wGP_Face[q]*Flux_visc[c];
#endif
                }
            }
            
            for (c = 0; c < nVar; ++c) {
                set_element_3d(Ctx->G, j-ys, i-xs, c, Flux[c]);
                set_element_3d(Ctx->E, j-ys, i-xs, c, NCP[c]);
            }
        }
    }

    // Now find the rhs in each cell 

    for (j=ys; j<ys+ym; ++j) {
        for (i=xs; i<xs+xm; ++i) {

            for (c = 0 ; c < nVar; ++c) {
        
                rhs[j][i].comp[c]  += -r1_h*(get_element_3d(Ctx->F, j-ys, i+1-xs, c) - get_element_3d(Ctx->F, j-ys, i-xs, c) + 
                                       0.5*(get_element_3d(Ctx->D, j-ys, i+1-xs, c) + get_element_3d(Ctx->D, j-ys, i-xs, c)))
                                      -r1_h*(get_element_3d(Ctx->G, j+1-ys, i-xs, c) - get_element_3d(Ctx->G, j-ys, i-xs, c) + 
                                       0.5*(get_element_3d(Ctx->E, j+1-ys, i-xs, c) + get_element_3d(Ctx->E, j-ys, i-xs, c)));
                                        
            }
            
        }
    }

    // Free all the memory

    Free1D(&coeffs);
    Free1D(&u_coeffs);
    Free2D(&sol,Ctx->nDOF);
    Free2D(&Vgp,Ctx->Ngp_Vol);
    Free2D(&grad_Vgp_x,Ctx->Ngp_Vol);
    Free2D(&grad_Vgp_y,Ctx->Ngp_Vol);
    Free2D(&VNode,Ctx->N_node);

    ierr = DMDAVecRestoreArray(da,Ctx->localU,&w);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,RHS,&rhs);CHKERRQ(ierr);

#ifdef VISCOUS
     dt = (Ctx->CFL*Ctx->h)/( 2.0*(s_max_c + (r1_h*s_max_v)*2.0) ); // Outer 2.0 corresponds to two dimensions
#else
    dt = (Ctx->CFL*Ctx->h)/(2.0*s_max_c); // 2.0 corresponds to 2 dimensions
#endif

    ierr = MPI_Allreduce(&dt, &Ctx->dt, 1, MPIU_REAL,MPIU_MIN, PetscObjectComm((PetscObject)da));CHKERRQ(ierr);

    return ierr;
}

