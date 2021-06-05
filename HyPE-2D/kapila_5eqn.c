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

# endif

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

void smooth_vortex_kp5(PetscReal x, PetscReal y, PetscReal* Q0) {

    PetscReal V0[nVar];

    PetscReal alpha_min = 0.01;
    PetscReal alpha_max = 0.99;
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

void air_helium_kp5(PetscReal x, PetscReal y, PetscReal* Q0) {

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

void water_air_kp5(PetscReal x, PetscReal y, PetscReal* Q0) {

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
        V0[5] = 0.99;
    }

    else {
        V0[0] = 1.0;
        V0[1] = 0.001;
        V0[2] = 0.0;
        V0[3] = 0.0;
        V0[4] = 1.0;
        V0[5] = 0.99;
    }

    // (Air) Bubble

    if ((x-x0)*(x-x0) + (y-y0)*(y-y0) < R*R) {

        V0[0] = 1.0;
        V0[1] = 0.001;
        V0[2] = 0.0;
        V0[3] = 0.0;
        V0[4] = 1.0;
        V0[5] = 0.01;
    }

    PDEPrim2Cons(V0, Q0);
}