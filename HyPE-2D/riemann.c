/*
 * riemann.c
 *      Author: sunder
 */

#include "hype.h"

//----------------------------------------------------------------------------
// The Roe matrix B̃(Ua,Ub) between two generic states Qa and Qb is numerically
// via Gaussian quadrature, by the function RoeMatrix
//----------------------------------------------------------------------------

void RoeMatrix(const PetscReal *Qa, const PetscReal *Qb, PetscReal nx, PetscReal ny, PetscReal BRoe[nVar][nVar]) {

    // 3 - point quadrature points in [0,1] for path integrals

    const PetscInt N_gps = 3;
    const PetscReal s_gp[] = {0.1127016653792583, 0.5, 0.8872983346207417};
    const PetscReal w_gp[] = {0.2777777777777778, 0.4444444444444444, 0.2777777777777778};

    PetscInt i, j, c, q;

    PetscReal Q[nVar];
    PetscReal B[nVar][nVar];

    // First make the BRoe matrix zero

    for (i = 0; i < nVar; ++i)
        for (j = 0; j < nVar; ++j)
            BRoe[i][j] = 0.0;

    for (q = 0; q < N_gps; ++q) {

        for (c = 0; c < nVar; ++c)
            Q[c] = Qa[c] + s_gp[q]*(Qb[c] - Qa[c]);

        PDEmatrixB(Q, nx, ny, B);

        for (i = 0; i < nVar; ++i) {
            for (j = 0; j < nVar; ++j) {
                BRoe[i][j] += w_gp[q]*B[i][j];
            }
        }
    }
}

//----------------------------------------------------------------------------
// Given the left and right state values, find the upwind flux in the
// direction (nx, ny) using LLF/Rusanov riemann solver
//----------------------------------------------------------------------------

PetscReal RiemannSolver(const PetscReal *QL, const PetscReal *QR,
                           const PetscReal nx, const PetscReal ny,
                           const PetscReal x,  const PetscReal y,
                           PetscReal *F, PetscReal *D) {

    int i, j, c;
    PetscReal FL[nVar], FR[nVar], Q_jump[nVar];
    PetscReal B[nVar][nVar];

    PetscReal s_max_l = PDEFlux(QL, nx, ny, x, y, FL);
    PetscReal s_max_r = PDEFlux(QR, nx, ny, x, y, FR);

    PetscReal s_max = PetscMax(s_max_l, s_max_r);

    for (c = 0; c < nVar; ++c) {
        Q_jump[c] = QR[c] - QL[c];
        F[c] = 0.5*( FR[c] + FL[c] - s_max*Q_jump[c] );
    }

    RoeMatrix(QL, QR, nx, ny, B);

    // Multiply B with the jump term

    for (i = 0; i < nVar; ++i) {

        D[i] = 0.0;

        for (j = 0; j < nVar; ++j) {
            D[i] += B[i][j]*Q_jump[j];
        }
    }

    return s_max;
}

//----------------------------------------------------------------------------
// Given the left and right state values, find the upwind flux in the
// direction (nx, ny) using LLF/Rusanov riemann solver
//----------------------------------------------------------------------------

PetscReal RiemannSolverPrim(const PetscReal *VL, const PetscReal *VR,
                           const PetscReal nx, const PetscReal ny,
                           const PetscReal x,  const PetscReal y,
                           PetscReal *F, PetscReal *D) {


    PetscInt c;
    PetscReal FL[nVar], FR[nVar], grad_V_x[nVar], grad_V_y[nVar];
    PetscReal QL[nVar], QR[nVar];
    PetscReal V_av[nVar];

    PDEPrim2Cons(VL,QL); PDEPrim2Cons(VR,QR);

    PetscReal s_max_l = PDEFluxPrim(VL, nx, ny, x, y, FL);
    PetscReal s_max_r = PDEFluxPrim(VR, nx, ny, x, y, FR);

    PetscReal s_max = PetscMax(s_max_l, s_max_r);

    for (c = 0; c < nVar; ++c) {
        V_av[c] = 0.5*(VR[c]+VL[c]);
        grad_V_x[c] = nx*(VR[c]-VL[c]);
        grad_V_y[c] = ny*(VR[c]-VL[c]);
        F[c] = 0.5*( FR[c] + FL[c] - s_max*(QR[c]-QL[c]) );
    }

    PDENCPPrim(V_av, grad_V_x, grad_V_y, D);

    return s_max;
}

//----------------------------------------------------------------------------
// Viscous Riemann Solver (Does average of the two fluxes)
//----------------------------------------------------------------------------

PetscReal ViscRiemannSolver(const PetscReal* QL, PetscReal grad_QL[nVar][DIM],
                               const PetscReal* QR, PetscReal grad_QR[nVar][DIM],
                               const PetscReal nx, const PetscReal ny,
                               PetscReal* Flux) {

    PetscReal FL[nVar], FR[nVar];
    PetscInt c;

    PetscReal s_max_l = PDEViscFlux(QL, grad_QL, nx, ny, FL);
    PetscReal s_max_r = PDEViscFlux(QR, grad_QR, nx, ny, FR);

    PetscReal s_max = PetscMax(s_max_l, s_max_r);

    for (c = 0; c < nVar; ++c)
        Flux[c] = 0.5*(FR[c] + FL[c]);

    return s_max;
}

//----------------------------------------------------------------------------
// Viscous Riemann Solver (Does average of the two fluxes)
//----------------------------------------------------------------------------

PetscReal ViscRiemannSolverPrim(const PetscReal* VL, PetscReal grad_VL[nVar][DIM],
                               const PetscReal* VR, PetscReal grad_VR[nVar][DIM],
                               const PetscReal nx, const PetscReal ny,
                               PetscReal* Flux) {

    PetscReal FL[nVar], FR[nVar];
    PetscInt c;

    PetscReal s_max_l = PDEViscFluxPrim(VL, grad_VL, nx, ny, FL);
    PetscReal s_max_r = PDEViscFluxPrim(VR, grad_VR, nx, ny, FR);

    PetscReal s_max = PetscMax(s_max_l, s_max_r);

    for (c = 0; c < nVar; ++c)
        Flux[c] = 0.5*(FR[c] + FL[c]);

    return s_max;
}
