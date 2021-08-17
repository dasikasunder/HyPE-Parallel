/*
 * kapila_5eqn.c
 *      Author: sunder
 */

#include "hype.h"

#ifdef KAPILA_5EQN

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable
//----------------------------------------------------------------------------

void PDECons2Prim(const PetscReal *Q, PetscReal *V) {

    PetscReal alpha_1 = Q[5];
    PetscReal alpha_2 = 1.0 - alpha_1;
    PetscReal GAMMA   = 1.0 + ((GAMMA_1-1.0)*(GAMMA_2-1.0))/((GAMMA_2-1.0)*alpha_1 + (GAMMA_1-1.0)*alpha_2);
    PetscReal PI      = ((GAMMA -1.0)/GAMMA)*( GAMMA_1*PI_1*alpha_1/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*alpha_2/(GAMMA_2 - 1.0) );


    PetscReal rho     = Q[0] + Q[1];
    PetscReal u       = Q[2]/rho;
    PetscReal v       = Q[3]/rho;
    PetscReal E       = Q[4];
    PetscReal p       = (GAMMA-1.0)*(E - 0.5*rho*(u*u + v*v)) - GAMMA*PI;

    V[0] = Q[0]/alpha_1;
    V[1] = Q[1]/alpha_2;
    V[2] = u;
    V[3] = v;
    V[4] = p;
    V[5] = alpha_1;

}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
//----------------------------------------------------------------------------

void PDEPrim2Cons(const PetscReal *V, PetscReal *Q) {

    PetscReal alpha_1 = V[5];
    PetscReal alpha_2 = 1.0 - alpha_1;
    PetscReal GAMMA   = 1.0 + ((GAMMA_1-1.0)*(GAMMA_2-1.0))/((GAMMA_2-1.0)*alpha_1 + (GAMMA_1-1.0)*alpha_2);
    PetscReal PI      = ((GAMMA -1.0)/GAMMA)*( GAMMA_1*PI_1*alpha_1/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*alpha_2/(GAMMA_2 - 1.0) );

    PetscReal rho     = alpha_1*V[0] + alpha_2*V[1];
    PetscReal u       = V[2];
    PetscReal v       = V[3];
    PetscReal p       = V[4];
    PetscReal E       = (p + GAMMA*PI)/(GAMMA-1.0) + 0.5*rho*(u*u + v*v);

    Q[0] = alpha_1*V[0];
    Q[1] = alpha_2*V[1];
    Q[2] = rho*V[2];
    Q[3] = rho*V[3];
    Q[4] = E;
    Q[5] = alpha_1;
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFlux(const PetscReal *Q,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    // First extract the primitive variables
    PetscReal alpha_1 = Q[5];
    PetscReal alpha_2 = 1.0 - alpha_1;
    PetscReal GAMMA   = 1.0 + ((GAMMA_1-1.0)*(GAMMA_2-1.0))/((GAMMA_2-1.0)*alpha_1 + (GAMMA_1-1.0)*alpha_2);
    PetscReal PI      = ((GAMMA -1.0)/GAMMA)*( GAMMA_1*PI_1*alpha_1/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*alpha_2/(GAMMA_2 - 1.0) );

    PetscReal rho_1   = Q[0]/alpha_1;
    PetscReal rho_2   = Q[1]/alpha_2;
    PetscReal rho     = Q[0] + Q[1];
    PetscReal irho    = 1.0/rho;
    PetscReal u       = irho*Q[2];
    PetscReal v       = irho*Q[3];
    PetscReal E       = Q[4];
    PetscReal p       = (GAMMA-1.0)*(E - 0.5*rho*(u*u + v*v)) - GAMMA*PI;
    PetscReal un = u*nx + v*ny;

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        printf("Negative density = %f\n", rho);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ( p + PI < prs_floor) {
        printf("Negative pressure = %f\n", p + PI);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    F[0] = alpha_1*rho_1*un;
    F[1] = alpha_2*rho_2*un;
    F[2] = rho*u*un + p*nx;
    F[3] = rho*v*un + p*ny;
    F[4] = (E+p)*un;
    F[5] = 0.0;

    // Obtain the maximum eigenvalue

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(GAMMA*(p+PI)*irho);

    return s_max;
}

//----------------------------------------------------------------------------
// Viscous part of the flux in the normal direction
//----------------------------------------------------------------------------

PetscReal PDEViscFlux(const PetscReal* Q, const PetscReal grad_Q[nVar][DIM], PetscReal nx, PetscReal ny, PetscReal* F) {

    PetscReal r2_3 = 2./3.; PetscReal r4_3 = 4./3.;

    // Find the phase fractions

    PetscReal phi = Q[5];

    // Effective viscosity

    PetscReal mu = phi*MU_1 + (1.0-phi)*MU_2;

    PetscReal irho = 1.0/(Q[0] + Q[1]);


    PetscReal u = irho*Q[2];
    PetscReal v = irho*Q[3];
    PetscReal rho_x = grad_Q[0][0] + grad_Q[1][0];
    PetscReal rho_y = grad_Q[0][1] + grad_Q[1][1];
    PetscReal u_x = irho*(grad_Q[2][0] - u*rho_x);
    PetscReal u_y = irho*(grad_Q[2][1] - u*rho_y);
    PetscReal v_x = irho*(grad_Q[3][0] - v*rho_x);
    PetscReal v_y = irho*(grad_Q[3][1] - v*rho_y);
    PetscReal div_v = u_x + v_y;

        // Stress tensor

    PetscReal tau_xx = mu*(2.0*u_x - r2_3*div_v);
    PetscReal tau_xy = mu*(u_y + v_x) ;
    PetscReal tau_yy = mu*(2.0*v_y - r2_3*div_v);

    F[0] = 0.0;
    F[1] = 0.0;
    F[2] = -nx*tau_xx - ny*tau_xy;
    F[3] = -nx*tau_xy - ny*tau_yy;
    F[4] = -nx*(u*tau_xx + v*tau_xy) - ny*(u*tau_xy + v*tau_yy);
    F[5] = 0.0;

    return r4_3*mu*irho;
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCP(const PetscReal *Q, const PetscReal *grad_Q_x, const PetscReal *grad_Q_y, PetscReal *BgradQ) {

    PetscReal alpha_1 = Q[5];
    PetscReal alpha_2 = 1.0 - alpha_1;
    PetscReal GAMMA   = 1.0 + ((GAMMA_1-1.0)*(GAMMA_2-1.0))/((GAMMA_2-1.0)*alpha_1 + (GAMMA_1-1.0)*alpha_2);
    PetscReal PI      = ((GAMMA -1.0)/GAMMA)*( GAMMA_1*PI_1*alpha_1/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*alpha_2/(GAMMA_2 - 1.0) );

    PetscReal rho_1   = Q[0]/alpha_1;
    PetscReal rho_2   = Q[1]/alpha_2;
    PetscReal rho     = Q[0] + Q[1];
    PetscReal irho    = 1.0/rho;
    PetscReal u       = irho*Q[2];
    PetscReal v       = irho*Q[3];
    PetscReal E       = Q[4];
    PetscReal p       = (GAMMA-1.0)*(E - 0.5*rho*(u*u + v*v)) - GAMMA*PI;

    PetscReal a_1_sq  = GAMMA_1*(p + PI_1)/rho_1;
    PetscReal a_2_sq  = GAMMA_2*(p + PI_2)/rho_2;

    PetscReal K = alpha_1*alpha_2*(rho_2*a_2_sq - rho_1*a_1_sq)/(alpha_1*rho_2*a_2_sq + alpha_2*rho_1*a_1_sq);

    PetscReal alpha1_x = grad_Q_x[5]; PetscReal alpha1_y = grad_Q_y[5];
    PetscReal u_x = irho*irho*(grad_Q_x[2]*rho - Q[2]*(grad_Q_x[0] + grad_Q_x[1]));
    PetscReal v_y = irho*irho*(grad_Q_y[3]*rho - Q[3]*(grad_Q_y[0] + grad_Q_y[1]));

    K = 0.0;

    // Now find the fluxes

    BgradQ[0] = 0.0;
    BgradQ[1] = 0.0;
    BgradQ[2] = 0.0;
    BgradQ[3] = 0.0;
    BgradQ[4] = 0.0;
    BgradQ[5] = u*alpha1_x + v*alpha1_y - K*(u_x + v_y);
}

//----------------------------------------------------------------------------
// Non-conservative matrix B
//----------------------------------------------------------------------------

void PDEmatrixB(const PetscReal *Q, PetscReal nx, PetscReal ny, PetscReal B[nVar][nVar]) {

    PetscInt i, j;
    PetscReal alpha_1 = Q[5];
    PetscReal alpha_2 = 1.0 - alpha_1;
    PetscReal GAMMA   = 1.0 + ((GAMMA_1-1.0)*(GAMMA_2-1.0))/((GAMMA_2-1.0)*alpha_1 + (GAMMA_1-1.0)*alpha_2);
    PetscReal PI      = ((GAMMA -1.0)/GAMMA)*( GAMMA_1*PI_1*alpha_1/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*alpha_2/(GAMMA_2 - 1.0) );

    PetscReal rho_1   = Q[0]/alpha_1;
    PetscReal rho_2   = Q[1]/alpha_2;
    PetscReal rho     = Q[0] + Q[1];
    PetscReal irho    = 1.0/rho;
    PetscReal u       = irho*Q[2];
    PetscReal v       = irho*Q[3];
    PetscReal E       = Q[4];
    PetscReal p       = (GAMMA-1.0)*(E - 0.5*rho*(u*u + v*v)) - GAMMA*PI;

    PetscReal a_1_sq  = GAMMA_1*(p + PI_1)/rho_1;
    PetscReal a_2_sq  = GAMMA_2*(p + PI_2)/rho_2;

    PetscReal K = alpha_1*alpha_2*(rho_2*a_2_sq - rho_1*a_1_sq)/(alpha_1*rho_2*a_2_sq + alpha_2*rho_1*a_1_sq);

    PetscReal un = u*nx + v*ny;

    K = 0.0;

    for (i = 0; i < nVar; ++i)
        for (j = 0; j < nVar; ++j)
            B[i][j] = 0.0;

    B[5][0] =  irho*K*un;
    B[5][1] =  irho*K*un;
    B[5][2] =  -irho*K*nx;
    B[5][3] =  -irho*K*ny;
    B[5][5] =  un;
}



//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPAD(const PetscReal *Q) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal phi_1 = Q[5];
    PetscReal phi_2 = 1.0 - phi_1;

    PetscReal temp = phi_1/(GAMMA_1-1.0) + phi_2/(GAMMA_2 - 1.0);
    PetscReal gamma = 1.0 + 1.0/temp;
    PetscReal pi = (phi_1*GAMMA_1*PI_1/(GAMMA_1-1.0) + phi_2*GAMMA_2*PI_2/(GAMMA_2 - 1.0))/(1.0 + temp);

    PetscReal rho = Q[0] + Q[1];
    PetscReal irho = 1.0/rho;
    PetscReal u = irho*Q[2];
    PetscReal v = irho*Q[3];
    PetscReal E = Q[4];
    PetscReal p = (gamma - 1.0)*(E - 0.5*rho*(u*u + v*v)) - gamma*pi;

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p + pi)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if (phi_1 < 1.0e-4 || phi_1 > 1.0-1.0e-4) {
        PAD = PETSC_FALSE;
    }

    return PAD;
}

//----------------------------------------------------------------------------
// Inlet boundary conditions (in conservative form)
//----------------------------------------------------------------------------

void InletBC(PetscReal x, PetscReal y, PetscReal t, PetscReal* Q) {
    const PetscReal M = 2.0;
    const PetscReal rho_air = 1.0;
    const PetscReal rho_water = 1000.0;
    const PetscReal P = 1.013;
    const PetscReal U_avg = M*PetscSqrtReal(GAMMA_2*P/rho_air);
    const PetscReal phi_R = 1.0 - 1.0e-6;
    const PetscReal phi_L = 1.0e-6;
    const PetscReal PHI = 0.5*( phi_R + phi_L );
    const PetscReal dif_phi = (phi_L - phi_R)/2.0;
    const PetscReal h = 30.0/500.0;
    const PetscReal epsilon = 3.0*h;
    const PetscReal r_0 = 1.0;
    PetscReal V[nVar];
    PetscReal r;

    r = PetscAbsReal(y);

    V[0] = rho_water;
    V[1] = rho_air;

    if (r < 1.0)
        V[2] = 2.0*U_avg*(1.0 - y*y);
    else
        V[2] = 0.0;

    V[3] = 0.0;
    V[4] = P;
    V[5] = PHI + dif_phi*PetscTanhReal(-(r-r_0)/epsilon);

    PDEPrim2Cons(V,Q);
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFluxPrim(const PetscReal *V,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    PetscReal phi_1 = V[5];
    PetscReal phi_2 = 1.0 - phi_1;

    PetscReal temp = phi_1/(GAMMA_1-1.0) + phi_2/(GAMMA_2 - 1.0);
    PetscReal gamma = 1.0 + 1.0/temp;
    PetscReal pi = (phi_1*GAMMA_1*PI_1/(GAMMA_1-1.0) + phi_2*GAMMA_2*PI_2/(GAMMA_2 - 1.0))/(1.0 + temp);

    PetscReal rho_1 = V[0];
    PetscReal rho_2 = V[1];

    PetscReal rho = phi_1*rho_1 + phi_2*rho_2;
    PetscReal irho = 1.0/rho;
    PetscReal u = V[2];
    PetscReal v = V[3];
    PetscReal p = V[4];
    PetscReal E = (p + gamma*pi)/(gamma - 1.0) + 0.5*rho*(u*u + v*v);

    PetscReal un = u*nx+v*ny;

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        printf("Negative density = %f\n", rho);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ( p + pi < prs_floor) {
        printf("Negative pressure = %f\n", p + pi);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    F[0] = phi_1*rho_1*un;
    F[1] = phi_2*rho_2*un;
    F[2] = rho*u*un + nx*p;
    F[3] = rho*v*un + ny*p;
    F[4] = un*(E + p);
    F[5] = 0.0;

    // Obtain the maximum eigenvalue

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(gamma*(p+pi)*irho);

    return s_max;
}

//----------------------------------------------------------------------------
// Viscous part of the flux in the normal direction
//----------------------------------------------------------------------------

PetscReal PDEViscFluxPrim(const PetscReal* V, const PetscReal grad_V[nVar][DIM], PetscReal nx, PetscReal ny, PetscReal* F) {

    PetscReal r2_3 = 2./3.; PetscReal r4_3 = 4./3.;

    // Find the phase fractions

    PetscReal phi = V[5];

    // Effective viscosity

    PetscReal mu = phi*MU_1 + (1.0-phi)*MU_2;

    PetscReal u = V[2];
    PetscReal v = V[3];
    PetscReal u_x = grad_V[2][0];
    PetscReal u_y = grad_V[2][1];
    PetscReal v_x = grad_V[3][0];
    PetscReal v_y = grad_V[3][1];
    PetscReal div_v = u_x + v_y;

    // Stress tensor

    PetscReal tau_xx = mu*(2.0*u_x - r2_3*div_v);
    PetscReal tau_xy = mu*(u_y + v_x) ;
    PetscReal tau_yy = mu*(2.0*v_y - r2_3*div_v);

    F[0] = 0.0;
    F[1] = 0.0;
    F[2] = -nx*tau_xx - ny*tau_xy;
    F[3] = -nx*tau_xy - ny*tau_yy;
    F[4] = -nx*(u*tau_xx + v*tau_xy) - ny*(u*tau_xy + v*tau_yy);
    F[5] = 0.0;

    return r4_3*mu/V[0];
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCPPrim(const PetscReal *V, const PetscReal *grad_V_x, const PetscReal *grad_V_y, PetscReal *BgradQ) {

    PetscReal phi_1 = V[5];
    PetscReal phi_2 = 1.0 - phi_1;

    PetscReal rho_1 = V[0];
    PetscReal rho_2 = V[1];

    PetscReal u = V[2];
    PetscReal v = V[3];
    PetscReal p = V[4];

    PetscReal phi_x = grad_V_x[5]; PetscReal phi_y = grad_V_y[5];
    PetscReal u_x = grad_V_x[2];
    PetscReal v_y = grad_V_y[3];

    PetscReal a_1_sq = GAMMA_1*(p + PI_1)/rho_1;
    PetscReal a_2_sq = GAMMA_2*(p + PI_2)/rho_2;
    PetscReal K = (phi_1*phi_2*(rho_2*a_2_sq - rho_1*a_1_sq))/(phi_1*rho_2*a_2_sq + phi_2*rho_1*a_1_sq);

    K = 0.0;

    // Find the Non-Conservative product

    BgradQ[0] = 0.0;
    BgradQ[1] = 0.0;
    BgradQ[2] = 0.0;
    BgradQ[3] = 0.0;
    BgradQ[4] = 0.0;
    BgradQ[5] = u*phi_x + v*phi_y - K*(u_x + v_y);
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPADPrim(const PetscReal *V) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal phi_1 = V[5];
    PetscReal phi_2 = 1.0 - phi_1;

    PetscReal temp = phi_1/(GAMMA_1-1.0) + phi_2/(GAMMA_2 - 1.0);
    PetscReal pi = (phi_1*GAMMA_1*PI_1/(GAMMA_1-1.0) + phi_2*GAMMA_2*PI_2/(GAMMA_2 - 1.0))/(1.0 + temp);

    PetscReal rho_1 = V[0];
    PetscReal rho_2 = V[1];

    PetscReal rho = phi_1*rho_1 + phi_2*rho_2;
    PetscReal p = V[4];

    // Check if the input state is physically admissible


    if (rho < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p + pi)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if (phi_1 < 0.0 || phi_1 > 1.0) {
        PAD = PETSC_FALSE;
    }

    return PAD;
}

//----------------------------------------------------------------------------
// Inlet boundary conditions (in primitive form)
//----------------------------------------------------------------------------

void InletBCPrim(PetscReal x, PetscReal y, PetscReal t, PetscReal* V) {
    const PetscReal M = 2.0;
    const PetscReal rho_air = 1.0;
    const PetscReal rho_water = 1000.0;
    const PetscReal P = 1.013;
    const PetscReal U_avg = M*PetscSqrtReal(GAMMA_2*P/rho_air);
    const PetscReal phi_R = 1.0 - 1.0e-6;
    const PetscReal phi_L = 1.0e-6;
    const PetscReal PHI = 0.5*( phi_R + phi_L );
    const PetscReal dif_phi = (phi_L - phi_R)/2.0;
    const PetscReal h = 30.0/500.0;
    const PetscReal epsilon = 3.0*h;
    const PetscReal r_0 = 1.0;
    PetscReal r;

    r = PetscAbsReal(y);

    V[0] = rho_water;
    V[1] = rho_air;

    if (r < 1.0)
        V[2] = 2.0*U_avg*(1.0 - y*y);
    else
        V[2] = 0.0;

    V[3] = 0.0;
    V[4] = P;
    V[5] = PHI + dif_phi*PetscTanhReal(-(r-r_0)/epsilon);
}

//---------------------------------------------------------------------------------------------------------------------------------------
// Test cases
//---------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Smooth Vortex Advection
// [x,y] \in [-3,3] x [-3,3]
// Final Time: 6.0
// BC: L-P, R-P, B-P, T-P
// GAMMA_1 = 4.0; GAMMA_2 = 1.4; PI_1 = 20.0; PI_2 = 0.0
//----------------------------------------------------------------------------

void smoothVortex_KP5(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    PetscReal alpha_min = 0.25;
    PetscReal alpha_max = 0.75;
    PetscReal omega = PETSC_PI/3.0;
    PetscReal p_atm = 1.0;
    PetscReal k_eps = 0.3;
    PetscReal R = 1.0;
    PetscReal delta = 0.1;
    PetscReal rho1 = 1000.0;
    PetscReal rho2 = 1.0;
    PetscReal u0 = 1.0;
    PetscReal v0 = 1.0;
    PetscReal r = PetscSqrtReal(x*x + y*y);

    V0[0] = rho1 + delta*rho1*PetscSinReal(omega*(2.0*x+y))*PetscCosReal(omega*(x-2.0*y));
    V0[1] = rho2 + delta*rho2*PetscSinReal(omega*(x-2.0*y))*PetscCosReal(omega*(2.0*x+y));
    V0[2] = u0;
    V0[3] = v0;
    V0[4] = p_atm;
    V0[5] = alpha_min + 0.5*(alpha_max-alpha_min)*erfc((r/R - 1.0)/k_eps);

    /*
    V0[0] = 1.0;
    V0[1] = 0.001;
    V0[2] = 1.0;
    V0[3] = 1.0;
    V0[4] = 1.01325;
    V0[5] = 0.5 + 0.25*PetscSinReal(PETSC_PI*(x+y));
    */

    PDEPrim2Cons(V0, Q0);

}

//----------------------------------------------------------------------------
// Advection of an interface
// [x,y] \in [0,1] x [-0.05,0.05]
// Final Time: 0.26/150.0
// BC: L-T, R-T, B-T, T-T
// GAMMA_1 = 4.4; GAMMA_2 = 1.4; PI_1 = 6.0e8; PI_2 = 0.0
//----------------------------------------------------------------------------

void interface_advection_kp5(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    if (x < 0.5) {
        V0[0] = 1000.0;
        V0[1] = 1.2;
        V0[2] = 150.0;
        V0[3] = 0.0;
        V0[4] = 101325.0;
        V0[5] = 0.01;
    }

    else {
        V0[0] = 1000.0;
        V0[1] = 1.2;
        V0[2] = 150.0;
        V0[3] = 0.0;
        V0[4] = 101325.0;
        V0[5] = 0.99;
    }

    PDEPrim2Cons(V0, Q0);

}

//----------------------------------------------------------------------------
// Air Helium-Bubble Shock Interaction
// [x,y] \in [0.0,0.356] x [0.0,0.089]
// Final Time: 674.0e-6
// BC: L-T, R-T, B-R, T-R
// GAMMA_1 = 1.4; GAMMA_2 = 1.648; PI_1 = 0.0; PI_2 = 0.0
//----------------------------------------------------------------------------

void AirHelium_KP5(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    PetscReal x0 = 0.245, y0 = 0.0455, R = 0.025;

    // (Air) Shock

    if (x >= 0.275) {
        V0[0] = 1.92691;
        V0[1] = 0.25463;
        V0[2] = -114.42;
        V0[3] = 0.0;
        V0[4] = 1.56980e5;
        V0[5] = 1.0-1.0e-6;
    }

    else {
        V0[0] = 1.4;
        V0[1] = 0.25463;
        V0[2] = 0.0;
        V0[3] = 0.0;
        V0[4] = 1.0e5;
        V0[5] = 1.0-1.0e-6;
    }

    // (Helium) Bubble

    if ((x-x0)*(x-x0) + (y-y0)*(y-y0) <= R*R) {
        V0[0] = 1.4;
        V0[1] = 0.25463;
        V0[2] = 0.0;
        V0[3] = 0.0;
        V0[4] = 1.0e5;
        V0[5] = 1.0e-6;
    }

    PDEPrim2Cons(V0, Q0);
}


//----------------------------------------------------------------------------
// Shock in water hitting air coloumn
// [x,y] \in [0,10] x [-2.5,2.5]
// Final Time: 0.02
// BC: L-T, R-T, B-T, T-T
// GAMMA_1 = 4.4; GAMMA_2 = 1.4; PI_1 = 6000.0; PI_2 = 0.0
//----------------------------------------------------------------------------

void WaterAir_KP5(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    PetscReal x0 = 4.375; PetscReal y0 = 0.0;
    PetscReal R = 1.0;

    // (Water) Shock

    if (x < 1.0) {
        V0[0] = 1.325;
        V0[1] = 0.001;
        V0[2] = 68.52;
        V0[3] = 0.0;
        V0[4] = 1.915e4;
        V0[5] = 1.0-1.0e-6;
    }

    else {
        V0[0] = 1.0;
        V0[1] = 0.001;
        V0[2] = 0.0;
        V0[3] = 0.0;
        V0[4] = 1.0;
        V0[5] = 1.0-1.0e-6;
    }

    // (Air) Bubble

    if ((x-x0)*(x-x0) + (y-y0)*(y-y0) < R*R) {

        V0[0] = 1.0;
        V0[1] = 0.001;
        V0[2] = 0.0;
        V0[3] = 0.0;
        V0[4] = 1.0;
        V0[5] = 1.0e-6;
    }

    PDEPrim2Cons(V0, Q0);
}

//----------------------------------------------------------------------------
// Mach 6 Shock in air hitting water cylinder
// [x,y] \in [0,8] x [-1,1]
// Final Time: 0.896
// BC: L-T, R-T, B-R, T-R
// g1 = 4.4; g2 = 1.4; p1 = 6000.0; g1 = 0.0
//----------------------------------------------------------------------------

void WaterCylinder_KP5(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    const PetscReal x0 = 2.0; const PetscReal y0 = 0.0;
    const PetscReal R = 0.562; const PetscReal eps = 1.0e-2;

    // Shock

    if (x < 1.0) {

        V0[0] = 1000.0;
        V0[1] = 5.26829;
        V0[2] = 5.79;
        V0[3] = 0.0;
        V0[4] = 42.39;
        V0[5] = eps;
    }

    else {
        V0[0] = 1000.0;
        V0[1] = 1.0;
        V0[2] = 0.0;;
        V0[3] = 0.0;
        V0[4] = 1.013;
        V0[5] = eps;
    }

    // Water Cylinder

    if ((x-x0)*(x-x0) + (y-y0)*(y-y0) < R*R) {

        V0[0] = 1000.0;
        V0[1] = 1.0;
        V0[2] = 0.0;;
        V0[3] = 0.0;
        V0[4] = 1.013;
        V0[5] = 1.0-eps;
    }

    PDEPrim2Cons(V0, Q0);
}

//----------------------------------------------------------------------------
// Air-Jet Entering water reservoir
// [x,y] \in [0,30] x [-15,15]
// Final Time: 1000.0
// BC: L-I, R-T, B-R, T-R
// GAMMA_1 = 4.4; GAMMA_2 = 1.4; PI_1 = 6000.0; PI_2 = 0.0
//----------------------------------------------------------------------------

void AirJet_KP5(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    V0[0] = 1000.0;
    V0[1] = 1.0;
    V0[2] = 0.0;
    V0[3] = 0.0;
    V0[4] = 1.013;
    V0[5] = 1.0 - 1.0e-6;

    PDEPrim2Cons(V0, Q0);
}

# endif


