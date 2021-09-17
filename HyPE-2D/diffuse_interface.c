/*
 * diffuse_interface.c
 *      Author: sunder
 */

#include "hype.h"

#ifdef DIFFUSE_INTERFACE

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable
//----------------------------------------------------------------------------

void PDECons2Prim(const PetscReal *Q, PetscReal *V) {

    PetscReal alpha = Q[4];
    PetscReal u_s   = Q[5];
    PetscReal v_s   = Q[6];
    PetscReal rho   = Q[0]/alpha;
    PetscReal u     = Q[1]/Q[0];
    PetscReal v     = Q[2]/Q[0];
    PetscReal E     = Q[3]/alpha;
    PetscReal p     = (GAMMA_1 - 1.0)*( E - 0.5*rho*(u*u + v*v) ) - GAMMA_1*PI_1;

    // Now assign them to vector V

    V[0] = rho;    // Density
    V[1] = u;      // x-velocity
    V[2] = v;      // y-velocity
    V[3] = p;      // pressure
    V[4] = alpha;  // Volume fraction
    V[5] = u_s;    // x-velocity of solid
    V[6] = v_s;    // y-velocity of solid
}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
//----------------------------------------------------------------------------

void PDEPrim2Cons(const PetscReal *V, PetscReal *Q) {

    Q[0] = V[4]*V[0]; // alpha*rho
    Q[1] = Q[0]*V[1]; // alpha*rho*u
    Q[2] = Q[0]*V[2]; // alpha*rho*v
    Q[3] = Q[0]*( 0.5*(V[1]*V[1] + V[2]*V[2]) + (V[3] + GAMMA_1*PI_1)/(V[0]*(GAMMA_1-1.0)) ); // alpha*E
    Q[4] = V[4]; // alpha
    Q[5] = V[5]; // u_s
    Q[6] = V[6]; // v_s
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFlux(const PetscReal *Q,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {


    PetscReal alpha, u_s, v_s, rho, irho, u, v, E, p, un, aun, rhoun;

    alpha = Q[4];
    u_s   = Q[5];
    v_s   = Q[6];
    rho   = Q[0]/alpha;
    irho  = 1.0/rho;
    u     = irho*Q[1];
    v     = irho*Q[2];
    E     = Q[3]/alpha;
    p     = (GAMMA_1 - 1.0)*( E - 0.5*rho*(u*u + v*v) ) - GAMMA_1*PI_1;
    un    = nx*u + ny*v;
    aun   = alpha*un;
    rhoun = rho*aun;


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

    // Find the fluxes

    F[0] = rhoun;
    F[1] = rhoun*u + nx*alpha*p;
    F[2] = rhoun*v + ny*alpha*p;
    F[3] = (E+p)*aun;
    F[4] = 0.0;
    F[5] = 0.0;
    F[6] = 0.0;

    PetscReal s_g = PetscAbsReal(un) + PetscSqrtReal(GAMMA_1*irho*(p + PI_1));
    PetscReal s_s = PetscAbsReal(u_s*nx + v_s*ny);

    return PetscMax(s_s, s_g);
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCP(const PetscReal *Q, const PetscReal *grad_Q_x, const PetscReal *grad_Q_y, PetscReal *BgradQ) {

    PetscReal alpha, u_s, v_s, rho, irho, u, v, E, p, alpha_x, alpha_y;

    alpha = Q[4];
    u_s   = Q[5];
    v_s   = Q[6];
    rho   = Q[0]/alpha;
    irho  = 1.0/rho;
    u     = irho*Q[1];
    v     = irho*Q[2];
    E     = Q[3]/alpha;
    p     = (GAMMA_1 - 1.0)*( E - 0.5*rho*(u*u + v*v) ) - GAMMA_1*PI_1;

    // Now find the fluxes

    alpha_x = grad_Q_x[4];
    alpha_y = grad_Q_y[4];

    BgradQ[0] = 0.0;
    BgradQ[1] = -p*alpha_x;
    BgradQ[2] = -p*alpha_y;
    BgradQ[3] = -p*(u_s*alpha_x + v_s*alpha_y);
    BgradQ[4] = u_s*alpha_x + v_s*alpha_y;
    BgradQ[5] = 0.0;
    BgradQ[6] = 0.0;
}

//----------------------------------------------------------------------------
// Non-conservative matrix B
//----------------------------------------------------------------------------

void PDEmatrixB(const PetscReal *Q, PetscReal nx, PetscReal ny, PetscReal B[nVar][nVar]) {

    PetscReal alpha, u_s, v_s, rho, irho, u, v, E, p;

    alpha = Q[4];
    u_s   = Q[5];
    v_s   = Q[6];
    rho   = Q[0]/alpha;
    irho  = 1.0/rho;
    u     = irho*Q[1];
    v     = irho*Q[2];
    E     = Q[3]/alpha;
    p     = (GAMMA_1 - 1.0)*( E - 0.5*rho*(u*u + v*v) ) - GAMMA_1*PI_1;

    PetscInt i, j;

    for (i = 0; i < nVar; ++i)
        for (j = 0; j < nVar; ++j)
            B[i][j] = 0.0;

    B[1][4] = -nx*p;
    B[2][4] = -ny*p;
    B[3][4] = -nx*p*u_s - ny*p*v_s;

    B[4][4] =  nx*u_s + ny*v_s;
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPAD(const PetscReal *Q) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal alpha, rho, irho, u, v, E, p;

    alpha = Q[4];
    rho   = Q[0]/alpha;
    irho  = 1.0/rho;
    u     = irho*Q[1];
    v     = irho*Q[2];
    E     = Q[3]/alpha;
    p     = (GAMMA_1 - 1.0)*( E - 0.5*rho*(u*u + v*v) ) - GAMMA_1*PI_1;

    // Check if the input state is physically admissible

    if (rho < rho_floor)
        PAD = PETSC_FALSE;

    if ((p + PI_1)  < prs_floor)
        PAD = PETSC_FALSE;

    return PAD;
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFluxPrim(const PetscReal *V,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    PetscReal alpha, u_s, v_s, rho, u, v, E, p, un, aun, rhoun;

    rho   = V[0];
    u     = V[1];
    v     = V[2];
    p     = V[3];
    alpha = V[4];
    u_s   = V[5];
    v_s   = V[6];
    E     = (p + GAMMA_1*PI_1)/(GAMMA_1 - 1.0) + 0.5*rho*(u*u + v*v);
    un    = nx*u + ny*v;
    aun   = alpha*un;
    rhoun = rho*aun;


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

    // Find the fluxes

    F[0] = rhoun;
    F[1] = rhoun*u + nx*alpha*p;
    F[2] = rhoun*v + ny*alpha*p;
    F[3] = (E+p)*aun;
    F[4] = 0.0;
    F[5] = 0.0;
    F[6] = 0.0;

    PetscReal s_g = PetscAbsReal(un) + PetscSqrtReal(GAMMA_1*(p + PI_1)/rho);
    PetscReal s_s = PetscAbsReal(u_s*nx + v_s*ny);

    return PetscMax(s_s, s_g);
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCPPrim(const PetscReal *V, const PetscReal *grad_V_x, const PetscReal *grad_V_y, PetscReal *BgradQ) {

    PetscReal u_s, v_s, p, alpha_x, alpha_y;

    p     = V[3];
    u_s   = V[5];
    v_s   = V[6];

    alpha_x = grad_V_x[4];
    alpha_y = grad_V_y[4];

    BgradQ[0] = 0.0;
    BgradQ[1] = -p*alpha_x;
    BgradQ[2] = -p*alpha_y;
    BgradQ[3] = -p*(u_s*alpha_x + v_s*alpha_y);
    BgradQ[4] = u_s*alpha_x + v_s*alpha_y;
    BgradQ[5] = 0.0;
    BgradQ[6] = 0.0;
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPADPrim(const PetscReal *V) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal alpha, rho, p;

    rho   = V[0];
    p     = V[3];
    alpha = V[4];

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p + PI_1)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if (alpha < 0.0 || alpha > 1.0) {
        PAD = PETSC_FALSE;
    }

    return PAD;
}


void shockDiffractionCylinder_DIM(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    PetscReal M_s = 4.0;
    PetscReal rho_1 = 1.4;
    PetscReal p_1 = 1.0;
    PetscReal a_1 = PetscSqrtReal(GAMMA_1*p_1/rho_1);
    PetscReal prs_ratio = 1.0 + (2.0*GAMMA_1)/(GAMMA_1 + 1.0)*(M_s*M_s - 1.0);
    PetscReal tempa = (GAMMA_1 + 1.0)/(GAMMA_1 - 1.0);
    PetscReal rho_ratio = (1.0 + tempa*prs_ratio)/(tempa + prs_ratio);
    tempa = (GAMMA_1-1.0)/(GAMMA_1+1.0);
    PetscReal tempb = ((2.0*GAMMA_1)/(GAMMA_1+1.0))/(prs_ratio + tempa);
    PetscReal u_p = (a_1/GAMMA_1)*(prs_ratio -1.0)*PetscSqrtReal(tempb);

    // Gas

    if (x < -0.5) {
        V0[0] = rho_1*rho_ratio;
        V0[1] = u_p;
        V0[2] = 0.0;
        V0[3] = p_1*prs_ratio;
    }

    else {

        V0[0] = rho_1;
        V0[1] = 0.0;
        V0[2] = 0.0;
        V0[3] = p_1;
    }

    // Solid

    PetscReal x0 = 0.0; PetscReal y0 = 0.0; PetscReal R0 = 1.0;

    if ((x-x0)*(x-x0) + (y-y0)*(y-y0) <= R0*R0)
        V0[4] = 1.0e-2;
    else
        V0[4] = 1.0 - 1.0e-2;

    V0[5] = 0.0;
    V0[6] = 0.0;

    PDEPrim2Cons(V0, Q0);
}

#endif
