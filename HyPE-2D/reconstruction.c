/*
 * reconstruction.c
 *      Author: sunder
 */ 
#include "hype.h"

//----------------------------------------------------------------------------
// Value of Nth Order basis functions 
//----------------------------------------------------------------------------

PetscReal basis(PetscReal x, PetscReal y, PetscInt n) {
    
    switch (n) {
        case 0:
            return 1.0;
            break; 
        case 1:
            return x; 
            break;
        case 2:
            return y;
            break;
        case 3:
            return x*x - 1./12.;
            break;
        case 4:
            return y*y - 1./12.;
            break;
        case 5:
            return x*y;
            break; 
        case 6:
            return x*(x*x - 3./20.);
            break;
        case 7:
            return y*(y*y - 3./20.);
            break; 
        case 8:
            return y*(x*x - 1./12.);
            break; 
        case 9:
            return x*(y*y - 1./12.);
            break;
        default:
            return 0.0;
    }
}

//----------------------------------------------------------------------------
// Gradients of Nth Order basis functions 
//----------------------------------------------------------------------------

void basis_grad(PetscReal x, PetscReal y, PetscInt n, PetscReal* grad_x, PetscReal* grad_y) {
    
    switch (n) {
        case 0:
            *grad_x = 0.0;
            *grad_y = 0.0; 
            break; 
        case 1:
            *grad_x = 1.0;
            *grad_y = 0.0; 
            break;
        case 2:
            *grad_x = 0.0;
            *grad_y = 1.0; 
            break;
        case 3:
            *grad_x = 2.0*x;
            *grad_y = 0.0; 
            break;
        case 4:
            *grad_x = 0.0;
            *grad_y = 2.0*y; 
            break;
        case 5:
            *grad_x = y;
            *grad_y = x; 
            break; 
        case 6:
            *grad_x = 3.0*x*x - 3./20.;
            *grad_y = 0.0; ;
            break;
        case 7:
            *grad_x = 0.0;
            *grad_y = 3.0*y*y - 3./20.;
            break; 
        case 8:
            *grad_x = 2.0*x*y;
            *grad_y = (x*x - 1./12.); 
            break; 
        case 9:
            *grad_x = (y*y - 1./12.);
            *grad_y = 2.0*x*y;
            break;
        default:
            *grad_x = 0.0;
            *grad_y = 0.0; 
    }
}

//----------------------------------------------------------------------------
// Minmod slope limiter 
//----------------------------------------------------------------------------

PetscReal minmod(PetscReal a, PetscReal b) {
    if (a*b < 0.0) {
        return 0.0;
    }
    
    else {
        if (PetscAbsReal(a) < PetscAbsReal(b) )
            return a;
        else 
            return b; 
    }
}

//----------------------------------------------------------------------------
// 2D 3th order WENO reconstruction 
//----------------------------------------------------------------------------

PetscReal pow4(PetscReal a) {
    PetscReal a2 = a*a; 
    
    return a2*a2; 
}

void weno2(const PetscReal* U_x, const PetscReal* U_y, const PetscReal* U_xy, PetscReal* coeffs, PetscReal* u_coeffs) {

    PetscReal u_0 = U_xy[0];
    PetscReal u_ip1 = U_x[3]; PetscReal u_jp1 = U_y[3];
    PetscReal u_im1 = U_x[1]; PetscReal u_jm1 = U_y[1];

    coeffs[0] = u_0;
    coeffs[1] = minmod(u_ip1-u_0,u_0-u_im1);
    coeffs[2] = minmod(u_jp1-u_0,u_0-u_jm1);

    u_coeffs[0] = u_0;
    u_coeffs[1] = 0.5*(u_ip1-u_im1);
    u_coeffs[2] = 0.5*(u_jp1-u_jm1);
}

void weno3(const PetscReal* U_x, const PetscReal* U_y, const PetscReal* U_xy, PetscReal* coeffs, PetscReal* u_coeffs) {
    
    
    const PetscReal central_cell_wt = 100.0;

    PetscReal u_0 = U_xy[0];
    PetscReal u_ip1 = U_x[3]; PetscReal u_jp1 = U_y[3]; PetscReal u_ip1jp1 = U_xy[1];
    PetscReal u_im1 = U_x[1]; PetscReal u_jm1 = U_y[1]; PetscReal u_ip1jm1 = U_xy[2];
    PetscReal u_ip2 = U_x[4]; PetscReal u_jp2 = U_y[4]; PetscReal u_im1jp1 = U_xy[3];
    PetscReal u_im2 = U_x[0]; PetscReal u_jm2 = U_y[0]; PetscReal u_im1jm1 = U_xy[4];
    PetscReal u_x[3]; PetscReal u_xx[3]; PetscReal u_xy[4];
    PetscReal wt[4]; PetscReal IS[4]; PetscReal total, temp;

    // Reconstruct in x-direction

    u_x[0] = -2.0*u_im1 + 0.5*u_im2 + 1.5*u_0; u_xx[0] = 0.5*(u_im2 - 2.0*u_im1 + u_0);
    u_x[1] = 0.5*(u_ip1 - u_im1);              u_xx[1] = 0.5*(u_im1 - 2.0*u_0 + u_ip1);
    u_x[2] = -1.5*u_0 + 2.0*u_ip1 - 0.5*u_ip2; u_xx[2] = 0.5*(u_0 - 2.0*u_ip1 + u_ip2);

    IS[0] = u_x[0]*u_x[0] + r13_3*u_xx[0]*u_xx[0]; wt[0] =             1.0/( pow4(IS[0] + small_num) );
    IS[1] = u_x[1]*u_x[1] + r13_3*u_xx[1]*u_xx[1]; wt[1] = central_cell_wt/( pow4(IS[1] + small_num) );
    IS[2] = u_x[2]*u_x[2] + r13_3*u_xx[2]*u_xx[2]; wt[2] =             1.0/( pow4(IS[2] + small_num) );

    total = wt[0] + wt[1] + wt[2];
    wt[0] = wt[0]/total;
    wt[1] = wt[1]/total;
    wt[2] = wt[2]/total;

    u_coeffs[0] = u_0;
    u_coeffs[1] = u_x[1];
    u_coeffs[3] = u_xx[1];

    coeffs[0] = u_0;
    coeffs[1] = wt[0]*u_x[0]  + wt[1]*u_x[1]  + wt[2]*u_x[2];
    coeffs[3] = wt[0]*u_xx[0] + wt[1]*u_xx[1] + wt[2]*u_xx[2];

    // Reconstruct in y-direction

    u_x[0] = -2.0*u_jm1 + 0.5*u_jm2 + 1.5*u_0; u_xx[0] = 0.5*(u_jm2 - 2.0*u_jm1 + u_0);
    u_x[1] = 0.5*(u_jp1 - u_jm1);              u_xx[1] = 0.5*(u_jm1 - 2.0*u_0 + u_jp1);
    u_x[2] = -1.5*u_0 + 2.0*u_jp1 - 0.5*u_jp2; u_xx[2] = 0.5*(u_0 - 2.0*u_jp1 + u_jp2);

    IS[0] = u_x[0]*u_x[0] + r13_3*u_xx[0]*u_xx[0]; wt[0] =             1.0/( pow4(IS[0] + small_num) );
    IS[1] = u_x[1]*u_x[1] + r13_3*u_xx[1]*u_xx[1]; wt[1] = central_cell_wt/( pow4(IS[1] + small_num) );
    IS[2] = u_x[2]*u_x[2] + r13_3*u_xx[2]*u_xx[2]; wt[2] =             1.0/( pow4(IS[2] + small_num) );

    total = wt[0] + wt[1] + wt[2];
    wt[0] = wt[0]/total;
    wt[1] = wt[1]/total;
    wt[2] = wt[2]/total;

    u_coeffs[2] = u_x[1];
    u_coeffs[4] = u_xx[1];

    coeffs[2] = wt[0]*u_x[0]  + wt[1]*u_x[1]  + wt[2]*u_x[2];
    coeffs[4] = wt[0]*u_xx[0] + wt[1]*u_xx[1] + wt[2]*u_xx[2];

    // Reconstruction in xy-direction

    u_xy[0] =  u_ip1jp1 - u_0 - coeffs[1] - coeffs[2] - coeffs[3] - coeffs[4];
    u_xy[1] = -u_ip1jm1 + u_0 + coeffs[1] - coeffs[2] + coeffs[3] + coeffs[4];
    u_xy[2] = -u_im1jp1 + u_0 - coeffs[1] + coeffs[2] + coeffs[3] + coeffs[4];
    u_xy[3] =  u_im1jm1 - u_0 + coeffs[1] + coeffs[2] - coeffs[3] - coeffs[4];

    temp = 4.0*(coeffs[3]*coeffs[3] + coeffs[4]*coeffs[4]);
    IS[0] = temp + u_xy[0]*u_xy[0]; wt[0] = 1.0/(pow4(IS[0] + small_num));
    IS[1] = temp + u_xy[1]*u_xy[1]; wt[1] = 1.0/(pow4(IS[1] + small_num));
    IS[2] = temp + u_xy[2]*u_xy[2]; wt[2] = 1.0/(pow4(IS[2] + small_num));
    IS[3] = temp + u_xy[3]*u_xy[3]; wt[3] = 1.0/(pow4(IS[3] + small_num));

    total = wt[0] + wt[1] + wt[2] + wt[3];
    wt[0] = wt[0]/total; wt[1] = wt[1]/total; wt[2] = wt[2]/total; wt[3] = wt[3]/total;

    u_coeffs[5] = 0.25*(u_ip1jp1 - u_ip1jm1 - u_im1jp1 + u_im1jm1);

    coeffs[5] = wt[0]*u_xy[0] + wt[1]*u_xy[1] + wt[2]*u_xy[2] + wt[3]*u_xy[3];
    
}

void weno4(const PetscReal* U_x, const PetscReal* U_y, const PetscReal* U_xy, PetscReal* coeffs, PetscReal* u_coeffs) {

    PetscReal gammaHi = 0.85;
    PetscReal gammaLo = 0.85;
    PetscReal u_0 = U_xy[0];
    PetscReal u_ip1 = U_x[3]; PetscReal u_jp1 = U_y[3]; PetscReal u_ip1jp1 = U_xy[1];
    PetscReal u_im1 = U_x[1]; PetscReal u_jm1 = U_y[1]; PetscReal u_ip1jm1 = U_xy[2];
    PetscReal u_ip2 = U_x[4]; PetscReal u_jp2 = U_y[4]; PetscReal u_im1jp1 = U_xy[3];
    PetscReal u_im2 = U_x[0]; PetscReal u_jm2 = U_y[0]; PetscReal u_im1jm1 = U_xy[4];
    PetscReal u_xR4, u_yR4, u_xxR4, u_yyR4, u_xyR4, u_xxxR4, u_yyyR4, u_xxyR4, u_xyyR4;
    PetscReal u_xR3[5]; PetscReal u_yR3[5]; PetscReal u_xxR3[5]; PetscReal u_yyR3[5]; PetscReal u_xyR3[5];
    PetscReal IS_R4; PetscReal gamma_R4; PetscReal w_R4;
    PetscReal IS_R3[5]; PetscReal gamma_R3[5];  PetscReal w_R3[5];  PetscReal sum = 0.0;
    PetscReal wt_ratio;
    PetscInt i;

    gamma_R4 = gammaHi;

    gamma_R3[0] = (1.0 - gammaHi)*gammaLo;

    for (i = 1 ; i < 5; ++i)
        gamma_R3[i] = 0.25*(1-gammaHi)*(1-gammaLo);

    // Fourth order stencil

    u_xR4   = r41_60*(-u_im1 + u_ip1) + r11_120*(u_im2 - u_ip2);
    u_yR4   = r41_60*(-u_jm1 + u_jp1) + r11_120*(u_jm2 - u_jp2);
    u_xxR4  = -u_0 + 0.5*(u_im1 + u_ip1);
    u_yyR4  = -u_0 + 0.5*(u_jm1 + u_jp1);
    u_xyR4  = 0.25*(u_im1jm1 - u_im1jp1 - u_ip1jm1 + u_ip1jp1);
    u_xxxR4 = r1_6*(u_im1 - u_ip1) + r1_12*(u_ip2 - u_im2);
    u_yyyR4 = r1_6*(u_jm1 - u_jp1) + r1_12*(u_jp2 - u_jm2);
    u_xxyR4 = 0.25*(-u_im1jm1 + u_im1jp1 - u_ip1jm1 + u_ip1jp1) + 0.5*(u_jm1 - u_jp1);
    u_xyyR4 = 0.5*(u_im1 - u_ip1) + 0.25*(u_ip1jp1 - u_im1jm1 - u_im1jp1  + u_ip1jm1);

    u_coeffs[0] = u_0;
    u_coeffs[1] = u_xR4;
    u_coeffs[2] = u_yR4;
    u_coeffs[3] = u_xxR4;
    u_coeffs[4] = u_yyR4;
    u_coeffs[5] = u_xyR4;
    u_coeffs[6] = u_xxxR4;
    u_coeffs[7] = u_yyyR4;
    u_coeffs[8] = u_xxyR4;
    u_coeffs[9] = u_xyyR4;

    IS_R4 = (u_xR4 + 0.1*u_xxxR4)*(u_xR4 + 0.1*u_xxxR4) + r13_3*(u_xxR4*u_xxR4 + u_yyR4*u_yyR4) +
            (u_yR4 + 0.1*u_yyyR4)*(u_yR4 + 0.1*u_yyyR4) + 39.05*(u_xxxR4*u_xxxR4 + u_yyyR4*u_yyyR4) +
            r7_6*u_xyR4*u_xyR4 + 4.7*(u_xxyR4*u_xxyR4 + u_xyyR4*u_xyyR4);

    w_R4 = gamma_R4/(pow4(IS_R4 + small_num));
    sum = w_R4;

    // Stencil 1 (centered stencil)

    u_xR3[0]  = 0.5*(-u_im1 + u_ip1);
    u_yR3[0]  = 0.5*(-u_jm1 + u_jp1);
    u_xxR3[0] = -u_0 + 0.5*(u_im1 + u_ip1);
    u_yyR3[0] = -u_0 + 0.5*(u_jm1 + u_jp1);
    u_xyR3[0] = 0.25*(u_im1jm1 - u_im1jp1 - u_ip1jm1 + u_ip1jp1);

    // Stencil 2 (left biased stencil)

    u_xR3[1]  = 1.5*u_0 - 2.0*u_im1 + 0.5*u_im2;
    u_yR3[1]  = 0.5*(-u_jm1 + u_jp1);
    u_xxR3[1] = 0.5*u_0 - u_im1 + 0.5*u_im2;
    u_yyR3[1] = -u_0 + 0.5*u_jm1 + 0.5*u_jp1;
    u_xyR3[1] = 0.5*(u_im1jm1 - u_im1jp1 - u_jm1 + u_jp1);

    // Stencil 3 (right biased stencil)

    u_xR3[2]  = -1.5*u_0 + 2.0*u_ip1 - 0.5*u_ip2;
    u_yR3[2]  = 0.5*(-u_jm1 + u_jp1);
    u_xxR3[2] = 0.5*u_0 - u_ip1 + 0.5*u_ip2;
    u_yyR3[2] = -u_0 + 0.5*u_jm1 + 0.5*u_jp1;
    u_xyR3[2] = 0.5*(-u_ip1jm1 + u_ip1jp1 + u_jm1 - u_jp1);

    // Stencil 4 (bottom biased stencil)

    u_xR3[3]  = 0.5*(-u_im1 + u_ip1);
    u_yR3[3]  = 1.5*u_0 - 2.0*u_jm1 + 0.5*u_jm2;
    u_xxR3[3] = -u_0 + 0.5*u_im1 + 0.5*u_ip1;
    u_yyR3[3] = 0.5*u_0 - u_jm1 + 0.5*u_jm2;
    u_xyR3[3] = 0.5*(-u_im1 + u_im1jm1 + u_ip1 - u_ip1jm1);

    // Stencil 5 (top biased stencil)

    u_xR3[4]  = 0.5*(-u_im1 + u_ip1);
    u_yR3[4]  = -1.5*u_0 + 2.0*u_jp1 - 0.5*u_jp2;
    u_xxR3[4] = -u_0 + 0.5*u_im1 + 0.5*u_ip1;
    u_yyR3[4] = 0.5*u_0 - u_jp1 + 0.5*u_jp2;
    u_xyR3[4] = 0.5*(u_im1 - u_im1jp1 - u_ip1 + u_ip1jp1);

    // Find the smoothness indicators

    for (i = 0; i < 5; ++i) {
        IS_R3[i] = u_xR3[i]*u_xR3[i] + u_yR3[i]*u_yR3[i] + r13_3*(u_xxR3[i]*u_xxR3[i] + u_yyR3[i]*u_yyR3[i]) + r7_6*u_xyR3[i]*u_xyR3[i];
        w_R3[i] = gamma_R3[i]/(pow4(IS_R3[i] + small_num));
        sum += w_R3[i];
    }

    // Normalize the weights

    w_R4 = w_R4/sum;

    for (i = 0; i < 5; ++i)
        w_R3[i] = w_R3[i]/sum;

    wt_ratio = w_R4/gamma_R4;

    coeffs[0] = u_0;

    coeffs[1] = wt_ratio*(u_xR4 - gamma_R3[0]*u_xR3[0] - gamma_R3[1]*u_xR3[1] - gamma_R3[2]*u_xR3[2] - gamma_R3[3]*u_xR3[3] - gamma_R3[4]*u_xR3[4])
                                + w_R3[0]*u_xR3[0] + w_R3[1]*u_xR3[1] + w_R3[2]*u_xR3[2] + w_R3[3]*u_xR3[3] + w_R3[4]*u_xR3[4];

    coeffs[2] = wt_ratio*(u_yR4 - gamma_R3[0]*u_yR3[0] - gamma_R3[1]*u_yR3[1] - gamma_R3[2]*u_yR3[2] - gamma_R3[3]*u_yR3[3] - gamma_R3[4]*u_yR3[4])
                                + w_R3[0]*u_yR3[0] + w_R3[1]*u_yR3[1] + w_R3[2]*u_yR3[2] + w_R3[3]*u_yR3[3] + w_R3[4]*u_yR3[4];

    coeffs[3] = wt_ratio*(u_xxR4 - gamma_R3[0]*u_xxR3[0] - gamma_R3[1]*u_xxR3[1] - gamma_R3[2]*u_xxR3[2] - gamma_R3[3]*u_xxR3[3] - gamma_R3[4]*u_xxR3[4])
                                + w_R3[0]*u_xxR3[0] + w_R3[1]*u_xxR3[1] + w_R3[2]*u_xxR3[2] + w_R3[3]*u_xxR3[3] + w_R3[4]*u_xxR3[4];

    coeffs[4] = wt_ratio*(u_yyR4 - gamma_R3[0]*u_yyR3[0] - gamma_R3[1]*u_yyR3[1] - gamma_R3[2]*u_yyR3[2] - gamma_R3[3]*u_yyR3[3] - gamma_R3[4]*u_yyR3[4])
                                + w_R3[0]*u_yyR3[0] + w_R3[1]*u_yyR3[1] + w_R3[2]*u_yyR3[2] + w_R3[3]*u_yyR3[3] + w_R3[4]*u_yyR3[4];

    coeffs[5] = wt_ratio*(u_xyR4 - gamma_R3[0]*u_xyR3[0] - gamma_R3[1]*u_xyR3[1] - gamma_R3[2]*u_xyR3[2] - gamma_R3[3]*u_xyR3[3] - gamma_R3[4]*u_xyR3[4])
                                + w_R3[0]*u_xyR3[0] + w_R3[1]*u_xyR3[1] + w_R3[2]*u_xyR3[2] + w_R3[3]*u_xyR3[3] + w_R3[4]*u_xyR3[4];

    coeffs[6] = wt_ratio*u_xxxR4;

    coeffs[7] = wt_ratio*u_yyyR4;

    coeffs[8] = wt_ratio*u_xxyR4;

    coeffs[9] = wt_ratio*u_xyyR4;
}

void Reconstruct(const PetscReal* U_x, const PetscReal* U_y, const PetscReal* U_xy, const PetscInt N, PetscReal* coeffs, PetscReal* u_coeffs) {

    switch (N) {
        case 0:
            coeffs[0] = U_xy[0];
            break;
        case 1:
            weno2(U_x,U_y,U_xy,coeffs,u_coeffs);
            break;
        case 2:
            weno3(U_x,U_y,U_xy,coeffs,u_coeffs);
            break;
        case 3:
            weno4(U_x,U_y,U_xy,coeffs,u_coeffs);
        default:
            break;
    }

}
