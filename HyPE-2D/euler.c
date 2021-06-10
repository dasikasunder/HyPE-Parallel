/*
 * euler.c
 *      Author: sunder
 */

#include "hype.h"

#ifdef EULER

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable
//----------------------------------------------------------------------------

void PDECons2Prim(const PetscReal *Q, PetscReal *V) {

    PetscReal irho = 1.0/Q[0];

    V[0] = Q[0];
    V[1] = irho*Q[1];
    V[2] = irho*Q[2];
    V[3] = (GAMMA_1 -1.0)*( Q[3] - 0.5*irho*(Q[1]*Q[1] + Q[2]*Q[2]) );

}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
//----------------------------------------------------------------------------

void PDEPrim2Cons(const PetscReal *V, PetscReal *Q) {

    PetscReal e = (V[3])/(GAMMA_1 - 1.0);
    PetscReal k = 0.5*V[0]*(V[1]*V[1] + V[2]*V[2]);

    Q[0] = V[0];
    Q[1] = V[0]*V[1];
    Q[2] = V[0]*V[2];
    Q[3] = k + e;
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFlux(const PetscReal *Q,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    PetscReal rho = Q[0];
    PetscReal irho = 1.0/rho;
    PetscReal u = irho*Q[1];
    PetscReal v = irho*Q[2];
    PetscReal un = nx*u + ny*v;
    PetscReal rhoun = rho*un;
    PetscReal E = Q[3];
    PetscReal p = (GAMMA_1 -1.0)*( E - 0.5*rho*(u*u + v*v) );

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        printf("Negative density = %f\n", rho);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if (p  < prs_floor) {
        printf("Negative pressure, p = %f\n", p);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Now find the fluxes

    F[0] = rhoun;
    F[1] = rhoun*u + nx*p;
    F[2] = rhoun*v + ny*p;
    F[3] = un*(E + p);

    // Also obtain the maximum eigen value

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(GAMMA_1*p*irho);

    return s_max;
}

//----------------------------------------------------------------------------
// Find the Viscous flux in the normal direction F_v.n
//----------------------------------------------------------------------------

PetscReal PDEViscFlux(const PetscReal* Q, const PetscReal grad_Q[nVar][DIM], PetscReal nx, PetscReal ny, PetscReal* F) {

    const PetscReal r2_3 = 2./3;
    const PetscReal r4_3 = 4./3;

    const PetscReal C_P   = 1.4;     // Specific heat at constant pressure
    const PetscReal C_V   = 1.0;     // Specific heat at constant volume
    const PetscReal PR    = 0.7;     // Prandtl number

    PetscReal rho = Q[0];
    PetscReal irho = 1.0/rho;

    PetscReal u = irho*Q[1];
    PetscReal v = irho*Q[2];
    PetscReal u_x = irho*(grad_Q[1][0] - u*grad_Q[0][0]);
    PetscReal u_y = irho*(grad_Q[1][1] - u*grad_Q[0][1]);
    PetscReal v_x = irho*(grad_Q[2][0] - v*grad_Q[0][0]);
    PetscReal v_y = irho*(grad_Q[2][1] - v*grad_Q[0][1]);
    PetscReal div_v = u_x + v_y;

    PetscReal mu = MU_1;     // Viscosity (you can use Sutherlands law here)
    PetscReal k = C_P*mu/PR; // Heat conductivity
    PetscReal R = C_P - C_V; // Gas constant

    // Viscous flux (for derivations see the sympy notebook)

    PetscReal tau_xx = mu*(2.0*u_x - r2_3*div_v);
    PetscReal tau_xy = mu*(u_y + v_x) ;
    PetscReal tau_yy = mu*(2.0*v_y - r2_3*div_v);

    PetscReal q_x = -k*(GAMMA_1 - 1.0)*(-2.0*(Q[1]*grad_Q[1][0] + Q[2]*grad_Q[2][0])*Q[0] +
            (Q[1]*Q[1] + Q[2]*Q[2])*grad_Q[0][0] +
            (-2.0*Q[3]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2])*grad_Q[0][0] +
            2*Q[0]*Q[0]*grad_Q[3][0])/(2.0*R*Q[0]*Q[0]*Q[0]);

    PetscReal q_y = -k*(GAMMA_1 - 1.0)*(-2*(Q[1]*grad_Q[1][1] + Q[2]*grad_Q[2][1])*Q[0] +
            (Q[1]*Q[1] + Q[2]*Q[2])*grad_Q[0][1] +
            (-2.0*Q[3]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2])*grad_Q[0][1] +
            2*Q[0]*Q[0]*grad_Q[3][1])/(2.0*R*Q[0]*Q[0]*Q[0]);

    F[0] =  0.0;
    F[1] = -nx*tau_xx - ny*tau_xy;
    F[2] = -nx*tau_xy - ny*tau_yy;
    F[3] = -nx*(u*tau_xx + v*tau_xy - 0.0*q_x) - ny*(u*tau_xy + v*tau_yy - 0.0*q_y);

    return PetscMax(r4_3*mu*irho, GAMMA_1*mu*irho/PR);
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCP(const PetscReal *Q, const PetscReal *grad_Q_x, const PetscReal *grad_Q_y, PetscReal *BgradQ) {

    BgradQ[0] = 0.0;
    BgradQ[1] = 0.0;
    BgradQ[2] = 0.0;
    BgradQ[3] = 0.0;
}


//----------------------------------------------------------------------------
// Non-conservative matrix B
//----------------------------------------------------------------------------

void PDEmatrixB(const PetscReal *Q, PetscReal nx, PetscReal ny, PetscReal B[nVar][nVar]) {

    PetscInt i, j;

    for (i = 0; i < nVar; ++i)
        for (j = 0; j < nVar; ++j)
            B[i][j] = 0.0;
}


//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPAD(const PetscReal *Q) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal rho = Q[0];
    PetscReal irho = 1.0/rho;
    PetscReal u = irho*Q[1];
    PetscReal v = irho*Q[2];
    PetscReal E = Q[3];
    PetscReal p = (GAMMA_1 -1.0)*( E - 0.5*rho*(u*u + v*v) );

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if (p  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    return PAD;
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFluxPrim(const PetscReal *V,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    PetscReal rho = V[0];
    PetscReal u = V[1];
    PetscReal v = V[2];
    PetscReal un = nx*u + ny*v;
    PetscReal rhoun = rho*un;
    PetscReal p = V[3];
    PetscReal E = p/(GAMMA_1 - 1.0) + 0.5*rho*(u*u + v*v);

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        printf("Negative density = %f\n", rho);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if (p  < prs_floor) {
        printf("Negative pressure, p = %f\n", p);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Now find the fluxes

    F[0] = rhoun;
    F[1] = rhoun*u + nx*p;
    F[2] = rhoun*v + ny*p;
    F[3] = un*(E + p);

    // Also obtain the maximum eigen value

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(GAMMA_1*p/rho);

    return s_max;
}

//----------------------------------------------------------------------------
// Viscous part of the flux in the normal direction
//----------------------------------------------------------------------------

PetscReal PDEViscFluxPrim(const PetscReal* V, const PetscReal grad_V[nVar][DIM], PetscReal nx, PetscReal ny, PetscReal* F) {

    PetscReal r2_3 = 2./3.; PetscReal r4_3 = 4./3.;


    // viscosity

    PetscReal mu = MU_1;

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

    // Find the Non-Conservative product

    BgradQ[0] = 0.0;
    BgradQ[1] = 0.0;
    BgradQ[2] = 0.0;
    BgradQ[3] = 0.0;
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPADPrim(const PetscReal *V) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal rho = V[0];
    PetscReal p = V[3];

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if (p  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    return PAD;
}

//----------------------------------------------------------------------------
// Test Cases
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Lid driven cavity
// [x,y] \in [0,1] x [0,1]
// Final Time: 10.0
// BC: L-wall, R-wall, B-wall, T-moving wall
// GAMMA = 1.4; MU = 1.0e-2
//----------------------------------------------------------------------------

void LidDrivenCavity_NS(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    V0[0] = 1.0;
    V0[1] = 0.0;
    V0[2] = 0.0;
    V0[3] = 100.0/GAMMA_1;

    PDEPrim2Cons(V0,Q0);
}

//----------------------------------------------------------------------------
// Viscous Shock Tube
// [x,y] \in [0,1] x [0,1]
// Final Time: 1.0
// BC: L-wall, R-wall, B-wall, T-wall
// GAMMA = 1.4; MU =1/Re (200,1000)
//----------------------------------------------------------------------------

void ViscousShockTube_NS(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    if (x <= 0.5) {
        V0[0] = 120.0;
        V0[1] = 0.0;
        V0[2] = 0.0;
        V0[3] = V0[0]/GAMMA_1;
    }

    else {
        V0[0] = 1.2;
        V0[1] = 0.0;
        V0[2] = 0.0;
        V0[3] = V0[0]/GAMMA_1;
    }

    PDEPrim2Cons(V0,Q0);
}

#endif


