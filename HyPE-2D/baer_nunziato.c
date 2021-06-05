/*
 * baer_nunziato.c
 *      Author: sunder
 */

#include "hype.h"

#ifdef BAER_NUNZIATO

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable
//----------------------------------------------------------------------------

void PDECons2Prim(const PetscReal *Q, PetscReal *V) {

    PetscReal phi_1 = Q[8];
    PetscReal phi_2 = 1.0 - Q[8];

    PetscReal rho_1 = Q[0]/phi_1;
    PetscReal rho_2 = Q[4]/phi_2;

    PetscReal u_1 = Q[1]/Q[0]; PetscReal v_1 = Q[2]/Q[0];
    PetscReal u_2 = Q[5]/Q[4]; PetscReal v_2 = Q[6]/Q[4];

    PetscReal E_1 = Q[3]/phi_1;
    PetscReal E_2 = Q[7]/phi_2;

    PetscReal p_1 = (GAMMA_1 - 1.0)*( E_1 - 0.5*rho_1*(u_1*u_1 + v_1*v_1) ) - GAMMA_1*PI_1;
    PetscReal p_2 = (GAMMA_2 - 1.0)*( E_2 - 0.5*rho_2*(u_2*u_2 + v_2*v_2) ) - GAMMA_2*PI_2;

    V[0] = rho_1;
    V[1] = u_1;
    V[2] = v_1;
    V[3] = p_1;
    V[4] = rho_2;
    V[5] = u_2;
    V[6] = v_2;
    V[7] = p_2;
    V[8] = phi_1;

}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
//----------------------------------------------------------------------------

void PDEPrim2Cons(const PetscReal *V, PetscReal *Q) {

    Q[0] = V[8]*V[0]; // phi_1*rho_1
    Q[1] = Q[0]*V[1]; // phi_1*rho_1*vx_1
    Q[2] = Q[0]*V[2]; // phi_1*rho_1*vy_1
    Q[3] = Q[0]*( 0.5*(V[1]*V[1] + V[2]*V[2]) + (V[3] + GAMMA_1*PI_1)/(V[0]*(GAMMA_1-1.0)) ); // phi_1*E_1

    Q[4] = (1.0-V[8])*V[4]; // phi_2*rho_2
    Q[5] = Q[4]*V[5];       // phi_2*rho_2*vx_2
    Q[6] = Q[4]*V[6];       // phi_2*rho_2*vy_2
    Q[7] = Q[4]*( 0.5*(V[5]*V[5] + V[6]*V[6]) + (V[7] + GAMMA_2*PI_2)/(V[4]*(GAMMA_2-1.0)) ); // phi_2*E_2

    Q[8] = V[8]; // phi_1
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

PetscReal PDEFlux(const PetscReal *Q,
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    PetscReal phi_1 = Q[8];
    PetscReal phi_2 = 1.0 - Q[8];

    PetscReal rho_1 = Q[0]/phi_1;
    PetscReal rho_2 = Q[4]/phi_2;

    PetscReal u_1 = Q[1]/Q[0]; PetscReal v_1 = Q[2]/Q[0];
    PetscReal u_2 = Q[5]/Q[4]; PetscReal v_2 = Q[6]/Q[4];

    PetscReal E_1 = Q[3]/phi_1;
    PetscReal E_2 = Q[7]/phi_2;

    PetscReal p_1 = (GAMMA_1 - 1.0)*( E_1 - 0.5*rho_1*(u_1*u_1 + v_1*v_1) ) - GAMMA_1*PI_1;
    PetscReal p_2 = (GAMMA_2 - 1.0)*( E_2 - 0.5*rho_2*(u_2*u_2 + v_2*v_2) ) - GAMMA_2*PI_2;

    PetscReal un_1 = nx*u_1 + ny*v_1;
    PetscReal un_2 = nx*u_2 + ny*v_2;

    // Check if the input state is physically admissible

    if (rho_1 < rho_floor) {
        printf("Negative density (solid) = %f\n", rho_1);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if (rho_2 < rho_floor) {
        printf("Negative density (gas) = %f\n", rho_2);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ((p_1 + PI_1)  < prs_floor) {
        printf("Negative pressure (solid), p = %f\n", p_1 + PI_1);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ((p_2 + PI_2)  < prs_floor) {
        printf("Negative pressure (gas), p = %f\n", p_2 + PI_2);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    F[0] = phi_1*rho_1*un_1;
    F[1] = phi_1*(rho_1*u_1*un_1 + nx*p_1);
    F[2] = phi_1*(rho_1*v_1*un_1 + ny*p_1);
    F[3] = phi_1*un_1*(E_1 + p_1);
    F[4] = phi_2*rho_2*un_2;
    F[5] = phi_2*(rho_2*u_2*un_2 + nx*p_2);
    F[6] = phi_2*(rho_2*v_2*un_2 + ny*p_2);
    F[7] = phi_2*un_2*(E_2 + p_2);
    F[8] = 0.0;

    // Obtain the maximum eigenvalue

    PetscReal s_1 = PetscAbsReal(un_1) + PetscSqrtReal(GAMMA_1*(p_1 + PI_1)/rho_1);
    PetscReal s_2 = PetscAbsReal(un_2) + PetscSqrtReal(GAMMA_2*(p_2 + PI_2)/rho_2);

    return PetscMax(s_1, s_2);
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCP(const PetscReal *Q, const PetscReal *grad_Q_x, const PetscReal *grad_Q_y, PetscReal *BgradQ) {

    // u_I = u_1; p_I = p_2 (standard approximation)

    PetscReal phi_2 = 1.0-Q[8];

    PetscReal rho_2 = Q[4]/phi_2;

    PetscReal u_1 = Q[1]/Q[0]; PetscReal v_1 = Q[2]/Q[0];
    PetscReal u_2 = Q[5]/Q[4]; PetscReal v_2 = Q[6]/Q[4];

    PetscReal E_2 = Q[7]/phi_2;

    PetscReal p_I = (GAMMA_2 - 1.0)*( E_2 - 0.5*rho_2*(u_2*u_2 + v_2*v_2) ) - GAMMA_2*PI_2;

    PetscReal phi_x = grad_Q_x[8]; PetscReal phi_y = grad_Q_y[8];

    // Find the Non-Conservative product

    BgradQ[0] = 0.0;
    BgradQ[1] = -p_I*phi_x;
    BgradQ[2] = -p_I*phi_y;
    BgradQ[3] = -p_I*(u_1*phi_x + v_1*phi_y);
    BgradQ[4] = 0.0;
    BgradQ[5] = p_I*phi_x;
    BgradQ[6] = p_I*phi_y;
    BgradQ[7] = p_I*(u_1*phi_x + v_1*phi_y);
    BgradQ[8] = u_1*phi_x + v_1*phi_y;
}

//----------------------------------------------------------------------------
// Non-conservative matrix B
//----------------------------------------------------------------------------

void PDEmatrixB(const PetscReal *Q, PetscReal nx, PetscReal ny, PetscReal B[nVar][nVar]) {

    // u_I = u_1; p_I = p_2 (standard approximation)

    PetscInt i, j;

    PetscReal phi_2 = 1.0 - Q[8];

    PetscReal rho_2 = Q[4]/phi_2;

    PetscReal u_1 = Q[1]/Q[0]; PetscReal v_1 = Q[2]/Q[0];
    PetscReal u_2 = Q[5]/Q[4]; PetscReal v_2 = Q[6]/Q[4];

    PetscReal E_2 = Q[7]/phi_2;

    PetscReal p_I = (GAMMA_2 - 1.0)*( E_2 - 0.5*rho_2*(u_2*u_2 + v_2*v_2) ) - GAMMA_2*PI_2;

    for (i = 0; i < nVar; ++i)
        for (j = 0; j < nVar; ++j)
            B[i][j] = 0.0;

    B[1][8] = -nx*p_I;
    B[2][8] = -ny*p_I;
    B[3][8] = -nx*p_I*u_1 - ny*p_I*v_1;

    B[5][8] =  nx*p_I;
    B[6][8] =  ny*p_I;
    B[7][8] =  nx*p_I*u_1 +  ny*p_I*v_1;

    B[8][8] =  nx*u_1 + ny*v_1;
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPAD(const PetscReal *Q) {

    PetscBool PAD = PETSC_TRUE;


    PetscReal phi_1 = Q[8];
    PetscReal phi_2 = 1.0 - Q[8];

    PetscReal rho_1 = Q[0]/phi_1;
    PetscReal rho_2 = Q[4]/phi_2;

    PetscReal u_1 = Q[1]/Q[0]; PetscReal v_1 = Q[2]/Q[0];
    PetscReal u_2 = Q[5]/Q[4]; PetscReal v_2 = Q[6]/Q[4];

    PetscReal E_1 = Q[3]/phi_1;
    PetscReal E_2 = Q[7]/phi_2;

    PetscReal p_1 = (GAMMA_1 - 1.0)*( E_1 - 0.5*rho_1*(u_1*u_1 + v_1*v_1) ) - GAMMA_1*PI_1;
    PetscReal p_2 = (GAMMA_2 - 1.0)*( E_2 - 0.5*rho_2*(u_2*u_2 + v_2*v_2) ) - GAMMA_2*PI_2;

    // Check if the input state is physically admissible

    if (rho_1 < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if (rho_2 < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p_1 + PI_1)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p_2 + PI_2)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if (phi_1 < 0.0 || phi_1 > 1.0) {
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

    PetscReal phi_1 = V[8];
    PetscReal phi_2 = 1.0 - V[8];

    PetscReal rho_1 = V[0];
    PetscReal rho_2 = V[4];

    PetscReal u_1 = V[1]; PetscReal v_1 = V[2];
    PetscReal u_2 = V[5]; PetscReal v_2 = V[6];

    PetscReal p_1 = V[3];
    PetscReal p_2 = V[8];

    PetscReal E_1 = (p_1 + GAMMA_1*PI_1)/(GAMMA_1 - 1.0) + 0.5*rho_1*(u_1*u_1 + v_1*v_1);
    PetscReal E_2 = (p_2 + GAMMA_2*PI_2)/(GAMMA_2 - 1.0) + 0.5*rho_2*(u_2*u_2 + v_2*v_2);

    PetscReal un_1 = nx*u_1 + ny*v_1;
    PetscReal un_2 = nx*u_2 + ny*v_2;

    // Check if the input state is physically admissible

    if (rho_1 < rho_floor) {
        printf("Negative density (solid) = %f\n", rho_1);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if (rho_2 < rho_floor) {
        printf("Negative density (gas) = %f\n", rho_2);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ((p_1 + PI_1)  < prs_floor) {
        printf("Negative pressure (solid), p = %f\n", p_1 + PI_1);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ((p_2 + PI_2)  < prs_floor) {
        printf("Negative pressure (gas), p = %f\n", p_2 + PI_2);
        printf("At x = %f, y = %f\n", x, y);
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    F[0] = phi_1*rho_1*un_1;
    F[1] = phi_1*(rho_1*u_1*un_1 + nx*p_1);
    F[2] = phi_1*(rho_1*v_1*un_1 + ny*p_1);
    F[3] = phi_1*un_1*(E_1 + p_1);
    F[4] = phi_2*rho_2*un_2;
    F[5] = phi_2*(rho_2*u_2*un_2 + nx*p_2);
    F[6] = phi_2*(rho_2*v_2*un_2 + ny*p_2);
    F[7] = phi_2*un_2*(E_2 + p_2);
    F[8] = 0.0;

    // Obtain the maximum eigenvalue

    PetscReal s_1 = PetscAbsReal(un_1) + PetscSqrtReal(GAMMA_1*(p_1 + PI_1)/rho_1);
    PetscReal s_2 = PetscAbsReal(un_2) + PetscSqrtReal(GAMMA_2*(p_2 + PI_2)/rho_2);

    return PetscMax(s_1, s_2);
}

//----------------------------------------------------------------------------
// Non-conservative product BgradQ
//----------------------------------------------------------------------------

void PDENCPPrim(const PetscReal *V, const PetscReal *grad_V_x, const PetscReal *grad_V_y, PetscReal *BgradQ) {

    PetscReal u_I = V[1]; PetscReal v_I = V[2]; // u_I = u_1
    PetscReal p_I = V[7]; // p_I = p_2
    PetscReal phi_x = grad_V_x[8]; PetscReal phi_y = grad_V_y[8];

    // Find the Non-Conservative product

    BgradQ[0] = 0.0;
    BgradQ[1] = -p_I*phi_x;
    BgradQ[2] = -p_I*phi_y;
    BgradQ[3] = -p_I*(u_I*phi_x + v_I*phi_y);
    BgradQ[4] = 0.0;
    BgradQ[5] = p_I*phi_x;
    BgradQ[6] = p_I*phi_y;
    BgradQ[7] = p_I*(u_I*phi_x + v_I*phi_y);
    BgradQ[8] = u_I*phi_x + v_I*phi_y;
}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state
//----------------------------------------------------------------------------

PetscBool PDECheckPADPrim(const PetscReal *V) {

    PetscBool PAD = PETSC_TRUE;

    PetscReal phi_1 = V[8];

    PetscReal rho_1 = V[0];
    PetscReal rho_2 = V[4];

    PetscReal p_1 = V[3];
    PetscReal p_2 = V[7];

    // Check if the input state is physically admissible

    if (rho_1 < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if (rho_2 < rho_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p_1 + PI_1)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if ((p_2 + PI_2)  < prs_floor) {
        PAD = PETSC_FALSE;
    }

    if (phi_1 < 0.0 || phi_1 > 1.0) {
        PAD = PETSC_FALSE;
    }

    return PAD;
}

#endif
