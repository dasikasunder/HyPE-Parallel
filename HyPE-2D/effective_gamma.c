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
        printf("Negative pressure, p = %f\n", p + P_inf);
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

# endif