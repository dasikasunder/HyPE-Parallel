/*
 * hype.h
 *      Author: sunder
 */

#ifndef HYPE_H_
#define HYPE_H_

//----------------------------------------------------------------------------
// Commonly used C header files
//----------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>

//----------------------------------------------------------------------------
// Petsc headers files 
//----------------------------------------------------------------------------

#include <petscvec.h>
#include <petscmath.h>
#include <petscdm.h>
#include <petscdmplex.h>
#include <petscts.h>
#include <petscdraw.h>
#include <petsctime.h>
#include <petscdmda.h>

//----------------------------------------------------------------------------
// Various constants used throught the code 
//----------------------------------------------------------------------------

//#define EFFECTIVE_GAMMA
#define KAPILA_5EQN
//# define VISCOUS
# define nVar 6      /* Number of components in the PDE system */
# define DIM 2       /* Dimensions of the problem */

// Physical Constants 

static const PetscReal GAMMA_1  = 4.4;      /* Specific heat ratio of first phase */
static const PetscReal GAMMA_2  = 1.4;      /* Specific heat ratio of second phase */
static const PetscReal PI_1     = 6000.0;   /* Stiffness constant of first phase */
static const PetscReal PI_2     = 0.0;      /* Stiffness constant of second phase */
static const PetscReal MU_1     = 1.0e-2;   /* Viscosity of first phase */
static const PetscReal MU_2     = 0.0;      /* Viscosity of second phase */
static const PetscReal G_A      = 9.8;      /* Accelaration due to gravity */

static const PetscReal prs_floor     = 1.0e-12; /* Pressure floor value */
static const PetscReal rho_floor     = 1.0e-14; /* Density floor value */
static const PetscReal small_num     = 1.0e-12; /* Effective small number in the code */
static const PetscInt  s_width       = 5;       /* Width of the stencil */ 

// Some rational numbers frequently used throught the code

static const PetscReal r1_6  = 1./6.;
static const PetscReal r13_3 = 13./3.;
static const PetscReal r7_6  = 7./6.;
static const PetscReal r11_120 = 11./120.;
static const PetscReal r1_12 = 1./12.;
static const PetscReal r41_60 = 41./60.;

//----------------------------------------------------------------------------
// Structure representing a multi-component field vector 
//----------------------------------------------------------------------------

typedef struct {
    PetscReal comp[nVar];
} Field;

//----------------------------------------------------------------------------
// Various types of boundary conditions 
//----------------------------------------------------------------------------

enum bndry_type{wall, inlet, periodic, reflective, transmissive};

//----------------------------------------------------------------------------
// Multidimensional array structures (upto 7 dimensions)
//---------------------------------------------------------------------------- 

PetscInt Malloc1D(PetscReal**, PetscInt);
PetscInt Free1D(PetscReal**);

PetscInt Malloc2D(PetscReal***, PetscInt, PetscInt);
PetscInt Free2D(PetscReal***, PetscInt);

// 1D

typedef struct {
  PetscInt size;
  PetscInt nelem; 
  PetscReal * data;
} array1d;

array1d* allocate1d(PetscInt);
array1d* copy_array1d(array1d*);
PetscInt free1d(array1d*);
void set_element_1d(array1d*, PetscInt, PetscReal);
PetscReal get_element_1d(array1d*, PetscInt);
void min_max_1d(array1d*, PetscReal*, PetscReal*);

// 2D

typedef struct {
  PetscInt size1; // rows
  PetscInt size2; // cols
  PetscInt nelem; 
  PetscReal * data;
} array2d;

array2d* allocate2d(PetscInt, PetscInt);    
array2d * copy_array2d(array2d*);
PetscInt free2d(array2d*);
void set_element_2d(array2d*, PetscInt, PetscInt, PetscReal);
PetscReal get_element_2d(array2d*, PetscInt, PetscInt);
void min_max_2d(array2d*, PetscReal*, PetscReal*);

// 3D

typedef struct {
  PetscInt size1, size2, size3; 
  PetscInt c1, c2, c3;  
  PetscInt nelem; 
  PetscReal * data;
} array3d;

array3d* allocate3d(PetscInt, PetscInt, PetscInt);    
array3d * copy_array3d(array3d*);
PetscInt free3d(array3d*);
void set_element_3d(array3d*, PetscInt, PetscInt, PetscInt, PetscReal);
PetscReal get_element_3d(array3d*, PetscInt, PetscInt, PetscInt);
void min_max_3d(array3d*, PetscReal*, PetscReal*);

// 4D

typedef struct {
  PetscInt size1, size2, size3, size4; 
  PetscInt c1, c2, c3, c4;  
  PetscInt nelem; 
  PetscReal * data;
} array4d;

array4d* allocate4d(PetscInt, PetscInt, PetscInt, PetscInt);
array4d * copy_array4d(array4d*);
PetscInt free4d(array4d*);
void set_element_4d(array4d*, PetscInt, PetscInt, PetscInt, PetscInt, PetscReal);
PetscReal get_element_4d(array4d*, PetscInt, PetscInt, PetscInt, PetscInt);
void min_max_4d(array4d*, PetscReal*, PetscReal*);

// 5D

typedef struct {
  PetscInt size1, size2, size3, size4, size5; 
  PetscInt c1, c2, c3, c4, c5;  
  PetscInt nelem; 
  PetscReal * data;
} array5d;

array5d* allocate5d(PetscInt, PetscInt, PetscInt, PetscInt, PetscInt);
array5d * copy_array5d(array5d*);
PetscInt free5d(array5d*);
void set_element_5d(array5d*, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscReal);
PetscReal get_element_5d(array5d*, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt);
void min_max_5d(array5d*, PetscReal*, PetscReal*);

// 6D

typedef struct {
  PetscInt size1, size2, size3, size4, size5, size6; 
  PetscInt c1, c2, c3, c4, c5, c6;  
  PetscInt nelem; 
  PetscReal * data;
} array6d;

array6d* allocate6d(PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt);
array6d * copy_array6d(array6d*);
PetscInt free6d(array6d*);
void set_element_6d(array6d*, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscReal);
PetscReal get_element_6d(array6d*, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt);
void min_max_6d(array6d*, PetscReal*, PetscReal*);

// 7D

typedef struct {
  PetscInt size1, size2, size3, size4, size5, size6, size7; 
  PetscInt c1, c2, c3, c4, c5, c6, c7;  
  PetscInt nelem; 
  PetscReal * data;
} array7d;

array7d* allocate7d(PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt);
array7d * copy_array7d(array7d*);
PetscInt free7d(array7d*);
void set_element_7d(array7d*, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscReal);
PetscReal get_element_7d(array7d*, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt);
void min_max_7d(array7d*, PetscReal*, PetscReal*);

//----------------------------------------------------------------------------
// Structure defining various critical parameters controlling the simulation 
//----------------------------------------------------------------------------

typedef struct {
    PetscReal x_min;                  /* x-coordinate of the domain begining */  
    PetscReal y_min;                  /* y-coordinate of the domain begining */
    PetscReal x_max;                  /* x-coordinate of the domain ending */
    PetscReal y_max;                  /* y-coordinate of the domain ending */
    PetscInt N_x;                     /* No. of cells in the x-direction */
    PetscInt N_y;                     /* No. of cells in the y-direction */
    PetscReal CFL;                    /* CFL condition, should be less than 0.5 */
    PetscReal dt;                     /* Time step size */
    PetscReal h;                      /* Grid size */
    PetscBool Restart;                /* Whether to start from restart file */
    PetscInt InitialStep;             /* Initial time step */
    PetscReal InitialTime;            /* Initial time of the simulation */
    PetscReal FinalTime;              /* Final time of the simulation */
    PetscInt WriteInterval;           /* No. of time steps after which data should be written */
    PetscInt RestartInterval;         /* No. of time steps after which restart file should be written */
    PetscBool ReconsPrimitive;        /* Flag to reconstruct primitive or conserved variables */
    enum bndry_type LeftBoundary;    /* Boundary condition on the left face */
    enum bndry_type RightBoundary;   /* Boundary condition on the right face */
    enum bndry_type TopBoundary;     /* Boundary condition on the top face */
    enum bndry_type BottomBoundary;  /* Boundary condition on the bottom face */
    PetscInt N;                       /* Degree of approximation (0,1,2,3) */
    PetscInt nDOF;                    /* Number of degrees of freedom per variable */
    PetscInt Ngp_Face;                /* Number of quadrature points on face */
    PetscReal* xGP_Face;              /* Weights of quadrature points on face */
    PetscReal* wGP_Face;              /* Weights of quadrature points on face */
    PetscInt Ngp_Vol;                 /* Number of quadrature points on volume */
    PetscReal* xGP_Vol;               /* Weights of quadrature points on volume */
    PetscReal* yGP_Vol;               /* Weights of quadrature points on volume */
    PetscReal* wGP_Vol;               /* Weights of quadrature points on volume */
    PetscInt N_node;                  /* Nodes to evaluate quality of solution and check physical admissibility */
    PetscReal* x_node;                /* x-coordinates of the nodal points */
    PetscReal* y_node;                /* y-coordinates of the nodal points */
    Vec W;                            /* Vector of primitive variables */
    Vec localU;                       /* Local solution vector */
    array5d* u_bnd;                   /* Boundary extrapolated values of conservative variables */
    array6d* u_bnd_grad;              /* Boundary extrapolated gradients of conservative variables */
    array3d* F;                       /* Upwind flux in x-direction */
    array3d* G;                       /* Upwind flux in y-direction */
    array3d* D;                       /* Fluctuation in x-direction */
    array3d* E;                       /* Fluctuation in y-direction */
    array3d* phiFace;                 /* Values of basis functions on faces */
    array2d* phiNode;                 /* Values of basis functions on interior nodes */
    array2d* phiVol;                  /* Values of basis functions on volume quadrature points */
    array2d* gradphiVol_x;            /* Gradient in x-direction on volume quadrature points */
    array2d* gradphiVol_y;            /* Gradient in y-direction on volume quadrature points */
    array3d* gradphiFace_x;           /* Gradient in x-direction on face quadrature points */
    array3d* gradphiFace_y;           /* Gradient in y-direction on face quadrature points */
    
} AppCtx;

//----------------------------------------------------------------------------
// Functions related to PDE
//----------------------------------------------------------------------------

void PDECons2Prim(const PetscReal*, PetscReal*);
void PDEPrim2Cons(const PetscReal*, PetscReal*);

/* Input variables are conserved variables */

PetscReal PDEFlux(const PetscReal*, const PetscReal, const PetscReal, const PetscReal,  const PetscReal, PetscReal*);
PetscReal PDEViscFlux(const PetscReal*, const PetscReal grad_Q[nVar][DIM], PetscReal, PetscReal, PetscReal*);
void PDENCP(const PetscReal*, const PetscReal*, const PetscReal*, PetscReal*);
void PDEmatrixB(const PetscReal*, PetscReal, PetscReal, PetscReal B[nVar][nVar]);
PetscBool PDECheckPAD(const PetscReal*);
void InletBC(PetscReal,PetscReal,PetscReal,PetscReal*);

PetscReal RiemannSolver(const PetscReal*, const PetscReal*, const PetscReal, const PetscReal, const PetscReal,  const PetscReal, PetscReal*, PetscReal*);
PetscReal ViscRiemannSolver(const PetscReal*, PetscReal grad_QL[nVar][DIM],
                               const PetscReal*, PetscReal grad_QR[nVar][DIM],
                               const PetscReal, const PetscReal,
                               PetscReal*);

/* Input variables are primitive variables */

PetscReal PDEFluxPrim(const PetscReal*, const PetscReal, const PetscReal, const PetscReal,  const PetscReal, PetscReal*);
PetscReal PDEViscFluxPrim(const PetscReal*, const PetscReal grad_V[nVar][DIM], PetscReal, PetscReal, PetscReal*);
void PDENCPPrim(const PetscReal*, const PetscReal*, const PetscReal*, PetscReal*);
PetscBool PDECheckPADPrim(const PetscReal*);

PetscReal RiemannSolverPrim(const PetscReal*, const PetscReal*, const PetscReal, const PetscReal, const PetscReal,  const PetscReal, PetscReal*, PetscReal*);
PetscReal ViscRiemannSolverPrim(const PetscReal*, PetscReal grad_QL[nVar][DIM],
                                const PetscReal*, PetscReal grad_QR[nVar][DIM],
                                const PetscReal, const PetscReal,
                                PetscReal*);
void InletBCPrim(PetscReal,PetscReal,PetscReal,PetscReal*);

//----------------------------------------------------------------------------
// WENO reconstruction 
//----------------------------------------------------------------------------

PetscReal basis(PetscReal, PetscReal, PetscInt);
void basis_grad(PetscReal, PetscReal, PetscInt, PetscReal*, PetscReal*);
PetscReal minmod(PetscReal, PetscReal); 
void Reconstruct(const PetscReal*, const PetscReal*, const PetscReal*, const PetscInt, PetscReal*, PetscReal*);

//----------------------------------------------------------------------------
// Main functions related to the solver 
//----------------------------------------------------------------------------

void InitialCondition(PetscReal, PetscReal, PetscReal*);
PetscErrorCode InitializeSolution(Vec, DM, AppCtx);
PetscErrorCode RHSFunction(TS, PetscReal, Vec, Vec, void*);
PetscErrorCode RHSFunctionPrim(TS, PetscReal, Vec, Vec, void*);
PetscErrorCode MonitorFunction (TS, PetscInt, PetscReal, Vec, void*);
PetscErrorCode ErrorNorms(Vec, DM, AppCtx, PetscReal*, PetscReal*);
PetscErrorCode ComputePrimitiveVariables(Vec, Vec, DM);

//----------------------------------------------------------------------------
// Test cases
//----------------------------------------------------------------------------


/* Euler's / Navier-Stokes equations */

void LidDrivenCavity_NS(PetscReal, PetscReal, PetscReal*);
void ViscousShockTube_NS(PetscReal, PetscReal, PetscReal*);

/* Kapila 5-EQN Model */

void smoothVortex_KP5(PetscReal, PetscReal, PetscReal*);
void interface_advection_kp5(PetscReal, PetscReal, PetscReal*);
void AirHelium_KP5(PetscReal, PetscReal, PetscReal*);
void WaterAir_KP5(PetscReal, PetscReal, PetscReal*);
void WaterCylinder_KP5(PetscReal, PetscReal, PetscReal*);
void AirJet_KP5(PetscReal, PetscReal, PetscReal*);

/* Effective-Gamma Model */

void WaterAir_EG(PetscReal, PetscReal, PetscReal*);
void AirJet_EG(PetscReal, PetscReal, PetscReal*);
void AirHelium_EG(PetscReal, PetscReal, PetscReal*);

/* Diffuse Interface Model */


#endif /* HYPE_H_ */ 
