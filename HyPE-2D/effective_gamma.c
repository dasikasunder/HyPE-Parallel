/*
 * effective_gamma.c
 *      Author: sunder
 */

#include "hype.h"

#ifdef EFFECTIVE_GAMMA

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable
//----------------------------------------------------------------------------

void PDECons2Prim(const PetscReal *Q, PetscReal *V) {

    PetscReal phi = Q[4];
    PetscReal g = 1.0 + (GAMMA_1-1.0)*(GAMMA_2-1.0)/((1.0-phi)*(GAMMA_1 -1.0) + phi*(GAMMA_2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( GAMMA_1*PI_1*phi/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*(1.0 - phi)/(GAMMA_2 - 1.0) );
    PetscReal irho = 1.0/Q[0];

    V[0] = Q[0];
    V[1] = irho*Q[1];
    V[2] = irho*Q[2];
    V[3] = (g -1.0)*( Q[3] - 0.5*irho*(Q[1]*Q[1] + Q[2]*Q[2]) )  - g*P_inf;
    V[4] = phi;

}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
//----------------------------------------------------------------------------

void PDEPrim2Cons(const PetscReal *V, PetscReal *Q) {

    PetscReal g = 1.0 + (GAMMA_1-1.0)*(GAMMA_2-1.0)/((1.0-V[4])*(GAMMA_1 -1.0) + V[4]*(GAMMA_2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( GAMMA_1*PI_1*V[4]/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*(1.0 - V[4])/(GAMMA_2 - 1.0) );

    PetscReal e = (V[3] + g*P_inf)/(g - 1.0);
    PetscReal k = 0.5*V[0]*(V[1]*V[1] + V[2]*V[2]);

    Q[0] = V[0];
    Q[1] = V[0]*V[1];
    Q[2] = V[0]*V[2];
    Q[3] = k + e;
    Q[4] = V[4];
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFlux(const PetscReal *Q,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    PetscReal phi = Q[4];
    PetscReal g = 1.0 + (GAMMA_1-1.0)*(GAMMA_2-1.0)/((1.0-phi)*(GAMMA_1 -1.0) + phi*(GAMMA_2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( GAMMA_1*PI_1*phi/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*(1.0 - phi)/(GAMMA_2 - 1.0) );
    PetscReal rho = Q[0];
    PetscReal irho = 1.0/rho;
    PetscReal u = irho*Q[1];
    PetscReal v = irho*Q[2];
    PetscReal E = Q[3];
    PetscReal p = (g -1.0)*( E - 0.5*rho*(u*u + v*v) )  - g*P_inf;
    PetscReal un = u*nx + v*ny;
    PetscReal rhoun = rho*un;

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        printf("Negative density = %f\n", rho);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ((p + P_inf)  < prs_floor) {
        printf("Negative pressure, p + p_inf = %f\n", p + P_inf);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Now find the fluxes

    F[0] = rhoun;
    F[1] = rhoun*u + p*nx;
    F[2] = rhoun*v + p*ny;
    F[3] = un*(E + p);
    F[4] = 0.0;

    // Also obtain the maximum eigen value

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(g*irho*(p + P_inf));

    return s_max;
}

//----------------------------------------------------------------------------
// Viscous part of the flux in the normal direction
//----------------------------------------------------------------------------

PetscReal PDEViscFlux(const PetscReal* Q, const PetscReal grad_Q[nVar][DIM], PetscReal nx, PetscReal ny, PetscReal* F) {

    PetscReal r2_3 = 2./3.; PetscReal r4_3 = 4./3.;

    // Find the phase fractions

    PetscReal phi = Q[4];

    // Effective viscosity

    PetscReal mu = phi*MU_1 + (1.0-phi)*MU_2;


    PetscReal irho = 1.0/Q[0];

    PetscReal u = irho*Q[1];
    PetscReal v = irho*Q[2];
    PetscReal u_x = irho*(grad_Q[1][0] - u*grad_Q[0][0]);
    PetscReal u_y = irho*(grad_Q[1][1] - u*grad_Q[0][1]);
    PetscReal v_x = irho*(grad_Q[2][0] - v*grad_Q[0][0]);
    PetscReal v_y = irho*(grad_Q[2][1] - v*grad_Q[0][1]);
    PetscReal div_v = u_x + v_y;

        // Stress tensor

    PetscReal tau_xx = mu*(2.0*u_x - r2_3*div_v);
    PetscReal tau_xy = mu*(u_y + v_x) ;
    PetscReal tau_yy = mu*(2.0*v_y - r2_3*div_v);

    F[0] = 0.0;
    F[1] = -nx*tau_xx - ny*tau_xy;
    F[2] = -nx*tau_xy - ny*tau_yy;
    F[3] = -nx*(u*tau_xx + v*tau_xy) - ny*(u*tau_xy + v*tau_yy);
    F[4] = 0.0;

    return r4_3*mu*irho;
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCP(const PetscReal *Q, const PetscReal *grad_Q_x, const PetscReal *grad_Q_y, PetscReal *BgradQ) {

    PetscReal phi_x = grad_Q_x[4]; PetscReal phi_y = grad_Q_y[4];

    // Now find the fluxes

    BgradQ[0] = 0.0;
    BgradQ[1] = 0.0;
    BgradQ[2] = 0.0;
    BgradQ[3] = 0.0;
    BgradQ[4] = (Q[1]*phi_x + Q[2]*phi_y)/Q[0];
}

//----------------------------------------------------------------------------
// Non-conservative matrix B
//----------------------------------------------------------------------------

void PDEmatrixB(const PetscReal *Q, PetscReal nx, PetscReal ny, PetscReal B[nVar][nVar]) {

    PetscInt i, j;

    for (i = 0; i < nVar; ++i)
        for (j = 0; j < nVar; ++j)
            B[i][j] = 0.0;

    B[4][4] =  (Q[1]*nx + Q[2]*ny)/Q[0];
}



//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPAD(const PetscReal *Q) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal phi = Q[4];
    PetscReal g = 1.0 + (GAMMA_1-1.0)*(GAMMA_2-1.0)/((1.0-phi)*(GAMMA_1 -1.0) + phi*(GAMMA_2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( GAMMA_1*PI_1*phi/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*(1.0 - phi)/(GAMMA_2 - 1.0) );
    PetscReal rho = Q[0];
    PetscReal irho = 1.0/rho;
    PetscReal u = irho*Q[1];
    PetscReal v = irho*Q[2];
    PetscReal E = Q[3];
    PetscReal p = (g -1.0)*( E - 0.5*rho*(u*u + v*v) )  - g*P_inf;

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p + P_inf)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if (phi < 0.0 || phi > 1.0) {
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

    if (r < 1.0) {
        V[0] = rho_air;
        V[1] = 2.0*U_avg*(1.0 - y*y);
    }
    else {
        V[0] = rho_water;
        V[1] = 0.0;
    }

    V[2] = 0.0;
    V[3] = P;
    V[4] = PHI + dif_phi*PetscTanhReal(-(r-r_0)/epsilon);

    PDEPrim2Cons(V,Q);
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFluxPrim(const PetscReal *V,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    PetscReal phi = V[4];
    PetscReal g = 1.0 + (GAMMA_1-1.0)*(GAMMA_2-1.0)/((1.0-phi)*(GAMMA_1 -1.0) + phi*(GAMMA_2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( GAMMA_1*PI_1*phi/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*(1.0 - phi)/(GAMMA_2 - 1.0) );
    PetscReal rho = V[0];
    PetscReal u = V[1];
    PetscReal v = V[2];
    PetscReal p = V[3];
    PetscReal E = (p + g*P_inf)/(g - 1.0) + 0.5*rho*(u*u + v*v);
    PetscReal un = u*nx + v*ny;
    PetscReal rhoun = rho*un;

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        printf("Negative density = %f\n", rho);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ((p + P_inf)  < prs_floor) {
        printf("Negative pressure, p + p_inf = %f\n", p + P_inf);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Now find the fluxes

    F[0] = rhoun;
    F[1] = rhoun*u + p*nx;
    F[2] = rhoun*v + p*ny;
    F[3] = un*(E + p);
    F[4] = 0.0;

    // Also obtain the maximum eigen value

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(g*(p + P_inf)/rho);

    return s_max;
}

//----------------------------------------------------------------------------
// Viscous part of the flux in the normal direction
//----------------------------------------------------------------------------

PetscReal PDEViscFluxPrim(const PetscReal* V, const PetscReal grad_V[nVar][DIM], PetscReal nx, PetscReal ny, PetscReal* F) {

    PetscReal r2_3 = 2./3.; PetscReal r4_3 = 4./3.;

    // Find the phase fractions

    PetscReal phi = V[4];

    // Effective viscosity

    PetscReal mu = phi*MU_1 + (1.0-phi)*MU_2;

    PetscReal u = V[1];
    PetscReal v = V[2];
    PetscReal u_x = grad_V[1][0];
    PetscReal u_y = grad_V[1][1];
    PetscReal v_x = grad_V[2][0];
    PetscReal v_y = grad_V[2][1];
    PetscReal div_v = u_x + v_y;

    // Stress tensor

    PetscReal tau_xx = mu*(2.0*u_x - r2_3*div_v);
    PetscReal tau_xy = mu*(u_y + v_x) ;
    PetscReal tau_yy = mu*(2.0*v_y - r2_3*div_v);

    F[0] = 0.0;
    F[1] = -nx*tau_xx - ny*tau_xy;
    F[2] = -nx*tau_xy - ny*tau_yy;
    F[3] = -nx*(u*tau_xx + v*tau_xy) - ny*(u*tau_xy + v*tau_yy);
    F[4] = 0.0;

    return r4_3*mu/V[0];
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCPPrim(const PetscReal *V, const PetscReal *grad_V_x, const PetscReal *grad_V_y, PetscReal *BgradQ) {

    PetscReal phi_x = grad_V_x[4]; PetscReal phi_y = grad_V_y[4];

    // Now find the fluxes

    BgradQ[0] = 0.0;
    BgradQ[1] = 0.0;
    BgradQ[2] = 0.0;
    BgradQ[3] = 0.0;
    BgradQ[4] = V[1]*phi_x + V[2]*phi_y;
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPADPrim(const PetscReal *V) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal phi = V[4];
    PetscReal g = 1.0 + (GAMMA_1-1.0)*(GAMMA_2-1.0)/((1.0-phi)*(GAMMA_1 -1.0) + phi*(GAMMA_2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( GAMMA_1*PI_1*phi/(GAMMA_1 - 1.0) + GAMMA_2*PI_2*(1.0 - phi)/(GAMMA_2 - 1.0) );
    PetscReal rho = V[0];
    PetscReal p = V[3];

    // Check if the input state is physically admissible


    if (rho < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p + P_inf)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if (phi < 0.0 || phi > 1.0) {
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

    if (r < 1.0) {
        V[0] = rho_air;
        V[1] = 2.0*U_avg*(1.0 - y*y);
    }
    else {
        V[0] = rho_water;
        V[1] = 0.0;
    }

    V[2] = 0.0;
    V[3] = P;
    V[4] = PHI + dif_phi*PetscTanhReal(-(r-r_0)/epsilon);
}

//---------------------------------------------------------------------------------------------------------------------------------------
// Test cases
//---------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Air Helium-Bubble Shock Interaction
// [x,y] \in [0.0,0.356] x [0.0,0.089]
// Final Time: 674.0e-6
// BC: L-T, R-T, B-R, T-R
// GAMMA_1 = 1.4; GAMMA_2 = 1.648; PI_1 = 0.0; PI_2 = 0.0
//----------------------------------------------------------------------------

void AirHelium_EG(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    PetscReal x0 = 0.245, y0 = 0.0455, R = 0.025;

    // (Air) Shock

    if (x >= 0.275) {
        V0[0] = 1.92691;
        V0[1] = -114.42;
        V0[2] = 0.0;
        V0[3] = 1.56980e5;
        V0[4] = 1.0-1.0e-6;
    }

    else {
        V0[0] = 1.4;
        V0[1] = 0.0;
        V0[2] = 0.0;
        V0[3] = 1.0e5;
        V0[4] = 1.0-1.0e-6;
    }

    // (Helium) Bubble

    if ((x-x0)*(x-x0) + (y-y0)*(y-y0) <= R*R) {
        V0[0] = 0.25463;
        V0[1] = 0.0;
        V0[2] = 0.0;
        V0[3] = 1.0e5;
        V0[4] = 1.0e-6;
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

void WaterAir_EG(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    PetscReal x0 = 4.375; PetscReal y0 = 0.0;
    PetscReal R = 1.0;

    // (Water) Shock

    if (x < 1.0) {
        V0[0] = 1.325;
        V0[1] = 68.52;
        V0[2] = 0.0;
        V0[3] = 1.915e4;
        V0[4] = 1.0-1.0e-6;
    }

    else {
        V0[0] = 1.0;
        V0[1] = 0.0;
        V0[2] = 0.0;
        V0[3] = 1.0;
        V0[4] = 1.0-1.0e-6;
    }

    // (Air) Bubble

    if ((x-x0)*(x-x0) + (y-y0)*(y-y0) < R*R) {

        V0[0] = 0.001;
        V0[1] = 0.0;
        V0[2] = 0.0;
        V0[3] = 1.0;
        V0[4] = 1.0e-6;
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

void AirJet_EG(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    V0[0] = 1000.0;
    V0[1] = 0.0;
    V0[2] = 0.0;
    V0[3] = 1.013;
    V0[4] = 1.0 - 1.0e-6;

    PDEPrim2Cons(V0, Q0);
}

# endif
