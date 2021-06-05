#include "hype.h" 

/*
 * array.c
 *      Author: sunder
 */ 

/*
 * Source for most of the data structures is:
 * http://theory.stanford.edu/~arbrad/pfe/html/ and 
 * https://github.com/astrobiology/orca_array/blob/master/orca_array.hpp
 */ 

//----------------------------------------------------------------------------------------------------
// C style 2D array
//----------------------------------------------------------------------------------------------------


PetscInt Malloc1D(PetscReal ** arr, PetscInt m){

    *arr = (PetscReal*)malloc(m*sizeof(PetscReal*));

  return 0;
}

PetscInt Free1D(PetscReal** arr){

    free(*arr);

    return 0;
}

PetscInt Malloc2D(PetscReal *** arr, PetscInt m, PetscInt n){

    PetscInt i;
    *arr = (PetscReal**)malloc(m*sizeof(PetscReal*));
    for(i=0; i<m; i++)
        (*arr)[i] = (PetscReal*)malloc(n*sizeof(PetscReal));

  return 0;
}

PetscInt Free2D(PetscReal*** arr, PetscInt n){

    PetscInt i;

    for (i = 0; i < n; i++)
        free((*arr)[i]);
    free(*arr);

    return 0;
}

//----------------------------------------------------------------------------------------------------
// 1D Array 
//----------------------------------------------------------------------------------------------------

/* Create a new array */

array1d* allocate1d(PetscInt size) {
    
    PetscInt i;
    
    // allocate a matrix structure
    array1d * arr = (array1d *) malloc(sizeof(array1d));

    // set dimensions
    arr->size  = size;
    arr->nelem = size; // Total number of elements in the array 

    // allocate a PetscReal array of length size
    
    arr->data = (PetscReal *) malloc(size*sizeof(PetscReal));
    
    // set all data to 0

    for (i = 0; i < size; i++)
        arr->data[i] = 0.0;

    return arr;
}

/* Copy from an existing array */

array1d * copy_array1d(array1d * arr) {
    if (!arr) return NULL;

    // create a new array to hold the copy

    array1d * cp = allocate1d(arr->size);

    // copy array data to cp's data

    memcpy(cp->data, arr->data, arr->nelem * sizeof(PetscReal));

    return cp;
}

/* Delete the array */

PetscInt free1d(array1d * arr) {
  
    if (!arr) return -1;
  
    // free arrays's data
    free(arr->data);
  
    // free array itself
    free(arr);
  
    return 0;
}

/* Set the value of an element */

void set_element_1d(array1d* arr1d, PetscInt i, PetscReal val) {
    arr1d->data[i] = val; 
}

/* Get the value of an element */

PetscReal get_element_1d(array1d* arr1d, PetscInt i) {
    return arr1d->data[i]; 
}

/* Get minimum and maximum of the array */

void min_max_1d(array1d* arr, PetscReal* min, PetscReal* max) {
    
    PetscInt i; 
    
    *min = arr->data[0]; 
    *max = arr->data[0]; 
    
    for (i = 1; i < arr->nelem; ++i) {
        if (arr->data[i] > *max)
            *max = arr->data[i];
        if (arr->data[i] < *min)
            *min = arr->data[i];
    }
}

//----------------------------------------------------------------------------------------------------
// 2D Array 
//----------------------------------------------------------------------------------------------------

/* Create a new array */

array2d* allocate2d(PetscInt size1, PetscInt size2) {
    
    PetscInt i;
    
    // allocate a matrix structure
    array2d * arr = (array2d *) malloc(sizeof(array2d));

    // set dimensions
    arr->size1 = size1;
    arr->size2 = size2;
    arr->nelem = size1*size2; // Total number of elements in the array 

    // allocate a PetscReal array of length size1*size2
    
    arr->data = (PetscReal *) malloc(arr->nelem*sizeof(PetscReal));
    
    // set all data to 0

    for (i = 0; i < arr->nelem; i++)
        arr->data[i] = 0.0;

    return arr;
}

/* Copy from an existing array */

array2d * copy_array2d(array2d * arr) {
  
    if (!arr) return NULL;

    // create a new array to hold the copy
  
    array2d * cp = allocate2d(arr->size1, arr->size2);

    // copy array data to cp's data
  
    memcpy(cp->data, arr->data, arr->nelem * sizeof(PetscReal));

    return cp;
}

/* Delete the array */

PetscInt free2d(array2d * arr) {
  
    if (!arr) return -1;
  
    // free arrays's data
    free(arr->data);
  
    // free array itself
    free(arr);
  
    return 0;
}

/* Set the value of an element */

void set_element_2d(array2d* arr, PetscInt i, PetscInt j, PetscReal val) {
    arr->data[i*arr->size2+j] = val; 
}

/* Get the value of an element */

PetscReal get_element_2d(array2d* arr, PetscInt i, PetscInt j) {
    return arr->data[i*arr->size2+j]; 
}

/* Get minimum and maximum of the array */

void min_max_2d(array2d* arr, PetscReal* min, PetscReal* max) {
    
    PetscInt i; 
    
    *min = arr->data[0]; 
    *max = arr->data[0]; 
    
    for (i = 1; i < arr->nelem; ++i) {
        if (arr->data[i] > *max)
            *max = arr->data[i];
        if (arr->data[i] < *min)
            *min = arr->data[i];
    }
}

//----------------------------------------------------------------------------------------------------
// 3D Array 
//----------------------------------------------------------------------------------------------------

/* Create a new array */

array3d* allocate3d(PetscInt size1, PetscInt size2, PetscInt size3) {
    
    PetscInt i;
    
    // allocate a matrix structure
    array3d * arr = (array3d *) malloc(sizeof(array3d));

    // set dimensions
    arr->size1 = size1;
    arr->size2 = size2;
    arr->size3 = size3; 
    arr->nelem = size1*size2*size3; // Total number of elements in the array 
    
    // factors for c order 
    
    arr->c1 = size2*size3;
    arr->c2 = size3;
    arr->c3 = 1;

    // allocate a PetscReal array of length size1*size2
    
    arr->data = (PetscReal *) malloc(arr->nelem*sizeof(PetscReal));
    
    // set all data to 0

    for (i = 0; i < arr->nelem; i++)
        arr->data[i] = 0.0;

    return arr;
}

/* Copy from an existing array */

array3d * copy_array3d(array3d * arr) {
  
    if (!arr) return NULL;

    // create a new array to hold the copy
  
    array3d * cp = allocate3d(arr->size1, arr->size2, arr->size3);

    // copy array data to cp's data
  
    memcpy(cp->data, arr->data, arr->nelem * sizeof(PetscReal));

    return cp;
}

/* Delete the array */

PetscInt free3d(array3d * arr) {
  
    if (!arr) return -1;
  
    // free arrays's data
    free(arr->data);
  
    // free array itself
    free(arr);
  
    return 0;
}

/* Set the value of an element */

void set_element_3d(array3d* arr, PetscInt i, PetscInt j, PetscInt k, PetscReal val) {
    arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3] = val; 
}

/* Get the value of an element */

PetscReal get_element_3d(array3d* arr, PetscInt i, PetscInt j, PetscInt k) {
    return arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3]; 
}

/* Get minimum and maximum of the array */

void min_max_3d(array3d* arr, PetscReal* min, PetscReal* max) {
    
    PetscInt i; 
    
    *min = arr->data[0]; 
    *max = arr->data[0]; 
    
    for (i = 1; i < arr->nelem; ++i) {
        if (arr->data[i] > *max)
            *max = arr->data[i];
        if (arr->data[i] < *min)
            *min = arr->data[i];
    }
}

//----------------------------------------------------------------------------------------------------
// 4D Array 
//----------------------------------------------------------------------------------------------------

/* Create a new array */

array4d* allocate4d(PetscInt size1, PetscInt size2, PetscInt size3, PetscInt size4) {
    
    PetscInt i;
    
    // allocate a matrix structure
    array4d * arr = (array4d *) malloc(sizeof(array4d));

    // set dimensions
    arr->size1 = size1;
    arr->size2 = size2;
    arr->size3 = size3;
    arr->size4 = size4;
    arr->nelem = size1*size2*size3*size4; // Total number of elements in the array 
    
    // factors for c order 
    
    arr->c1 = size2*size3*size4 ;
    arr->c2 = size3*size4 ; 
    arr->c3 = size4 ;
    arr->c4 = 1;

    // allocate a PetscReal array of length size1*size2
    
    arr->data = (PetscReal *) malloc(arr->nelem*sizeof(PetscReal));
    
    // set all data to 0

    for (i = 0; i < arr->nelem; i++)
        arr->data[i] = 0.0;

    return arr;
}

/* Copy from an existing array */

array4d * copy_array4d(array4d * arr) {
  
    if (!arr) return NULL;

    // create a new array to hold the copy
  
    array4d * cp = allocate4d(arr->size1, arr->size2, arr->size3, arr->size4);

    // copy array data to cp's data
  
    memcpy(cp->data, arr->data, arr->nelem * sizeof(PetscReal));

    return cp;
}

/* Delete the array */

PetscInt free4d(array4d * arr) {
  
    if (!arr) return -1;
  
    // free arrays's data
    free(arr->data);
  
    // free array itself
    free(arr);
  
    return 0;
}

/* Set the value of an element */

void set_element_4d(array4d* arr, PetscInt i, PetscInt j, PetscInt k, PetscInt l, PetscReal val) {
    arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3 + l*arr->c4] = val; 
}

/* Get the value of an element */

PetscReal get_element_4d(array4d* arr, PetscInt i, PetscInt j, PetscInt k, PetscInt l) {
    return arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3 + l*arr->c4]; 
}

/* Get minimum and maximum of the array */

void min_max_4d(array4d* arr, PetscReal* min, PetscReal* max) {
    
    PetscInt i; 
    
    *min = arr->data[0]; 
    *max = arr->data[0]; 
    
    for (i = 1; i < arr->nelem; ++i) {
        if (arr->data[i] > *max)
            *max = arr->data[i];
        if (arr->data[i] < *min)
            *min = arr->data[i];
    }
}

//----------------------------------------------------------------------------------------------------
// 5D Array 
//----------------------------------------------------------------------------------------------------

/* Create a new array */

array5d* allocate5d(PetscInt size1, PetscInt size2, PetscInt size3, PetscInt size4, PetscInt size5) {
    
    PetscInt i;
    
    // allocate a matrix structure
    array5d * arr = (array5d *) malloc(sizeof(array5d));

    // set dimensions
    arr->size1 = size1;
    arr->size2 = size2;
    arr->size3 = size3;
    arr->size4 = size4;
    arr->size5 = size5; 
    arr->nelem = size1*size2*size3*size4*size5; // Total number of elements in the array 
    
    // factors for c order 
    
    arr->c1 = size2*size3*size4*size5; 
    arr->c2 = size3*size4*size5;
    arr->c3 = size4*size5; 
    arr->c4 = size5; 
    arr->c5 = 1;

    // allocate a PetscReal array of length size1*size2
    
    arr->data = (PetscReal *) malloc(arr->nelem*sizeof(PetscReal));
    
    // set all data to 0

    for (i = 0; i < arr->nelem; i++)
        arr->data[i] = 0.0;

    return arr;
}

/* Copy from an existing array */

array5d * copy_array5d(array5d * arr) {
  
    if (!arr) return NULL;

    // create a new array to hold the copy
  
    array5d * cp = allocate5d(arr->size1, arr->size2, arr->size3, arr->size4, arr->size5);

    // copy array data to cp's data
  
    memcpy(cp->data, arr->data, arr->nelem * sizeof(PetscReal));

    return cp;
}

/* Delete the array */

PetscInt free5d(array5d * arr) {
  
    if (!arr) return -1;
  
    // free arrays's data
    free(arr->data);
  
    // free array itself
    free(arr);
  
    return 0;
}

/* Set the value of an element */

void set_element_5d(array5d* arr, PetscInt i, PetscInt j, PetscInt k, PetscInt l, PetscInt m, PetscReal val) {
    arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3 + l*arr->c4 + m*arr->c5] = val; 
}

/* Get the value of an element */

PetscReal get_element_5d(array5d* arr, PetscInt i, PetscInt j, PetscInt k, PetscInt l, PetscInt m) {
    return arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3 + l*arr->c4 + m*arr->c5]; 
}

/* Get minimum and maximum of the array */

void min_max_5d(array5d* arr, PetscReal* min, PetscReal* max) {
    
    PetscInt i; 
    
    *min = arr->data[0]; 
    *max = arr->data[0]; 
    
    for (i = 1; i < arr->nelem; ++i) {
        if (arr->data[i] > *max)
            *max = arr->data[i];
        if (arr->data[i] < *min)
            *min = arr->data[i];
    }
}

//----------------------------------------------------------------------------------------------------
// 6D Array 
//----------------------------------------------------------------------------------------------------

/* Create a new array */

array6d* allocate6d(PetscInt size1, PetscInt size2, PetscInt size3, PetscInt size4, PetscInt size5, PetscInt size6) {
    
    PetscInt i;
    
    // allocate a matrix structure
    array6d * arr = (array6d *) malloc(sizeof(array6d));

    // set dimensions
    arr->size1 = size1;
    arr->size2 = size2;
    arr->size3 = size3;
    arr->size4 = size4;
    arr->size5 = size5; 
    arr->size6 = size6; 
    arr->nelem = size1*size2*size3*size4*size5*size6; // Total number of elements in the array 
    
    // factors for c order 
    
     arr->c1 = size2*size3*size4*size5*size6; 
     arr->c2 = size3*size4*size5*size6;
     arr->c3 = size4*size5*size6; 
     arr->c4 = size5*size6; 
     arr->c5 = size6; 
     arr->c6 = 1;

    // allocate a PetscReal array of length size1*size2
    
    arr->data = (PetscReal *) malloc(arr->nelem*sizeof(PetscReal));
    
    // set all data to 0

    for (i = 0; i < arr->nelem; i++)
        arr->data[i] = 0.0;

    return arr;
}

/* Copy from an existing array */

array6d * copy_array6d(array6d * arr) {
  
    if (!arr) return NULL;

    // create a new array to hold the copy
  
    array6d * cp = allocate6d(arr->size1, arr->size2, arr->size3, arr->size4, arr->size5, arr->size6);

    // copy array data to cp's data
  
    memcpy(cp->data, arr->data, arr->nelem * sizeof(PetscReal));

    return cp;
}

/* Delete the array */

PetscInt free6d(array6d * arr) {
  
    if (!arr) return -1;
  
    // free arrays's data
    free(arr->data);
  
    // free array itself
    free(arr);
  
    return 0;
}

/* Set the value of an element */

void set_element_6d(array6d* arr, PetscInt i, PetscInt j, PetscInt k, PetscInt l, PetscInt m, PetscInt n, PetscReal val) {
    arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3 + l*arr->c4 + m*arr->c5 + n*arr->c6] = val; 
}

/* Get the value of an element */

PetscReal get_element_6d(array6d* arr, PetscInt i, PetscInt j, PetscInt k, PetscInt l, PetscInt m, PetscInt n) {
    return arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3 + l*arr->c4 + m*arr->c5 + n*arr->c6]; 
}

/* Get minimum and maximum of the array */

void min_max_6d(array6d* arr, PetscReal* min, PetscReal* max) {
    
    PetscInt i; 
    
    *min = arr->data[0]; 
    *max = arr->data[0]; 
    
    for (i = 1; i < arr->nelem; ++i) {
        if (arr->data[i] > *max)
            *max = arr->data[i];
        if (arr->data[i] < *min)
            *min = arr->data[i];
    }
}

//----------------------------------------------------------------------------------------------------
// 7D Array 
//----------------------------------------------------------------------------------------------------

/* Create a new array */

array7d* allocate7d(PetscInt size1, PetscInt size2, PetscInt size3, PetscInt size4, PetscInt size5, PetscInt size6, PetscInt size7) {
    
    PetscInt i;
    
    // allocate a matrix structure
    array7d * arr = (array7d *) malloc(sizeof(array7d));

    // set dimensions
    arr->size1 = size1;
    arr->size2 = size2;
    arr->size3 = size3;
    arr->size4 = size4;
    arr->size5 = size5; 
    arr->size6 = size6;
    arr->size7 = size7;
    arr->nelem = size1*size2*size3*size4*size5*size6*size7; // Total number of elements in the array 
    
    // factors for c order 
    
    arr->c1 = size2*size3*size4*size5*size6*size7; 
    arr->c2 = size3*size4*size5*size6*size7;
    arr->c3 = size4*size5*size6*size7; 
    arr->c4 = size5*size6*size7; 
    arr->c5 = size6*size7;
    arr->c6 = size7; 
    arr->c7 = 1;

    // allocate a PetscReal array of length size1*size2
    
    arr->data = (PetscReal *) malloc(arr->nelem*sizeof(PetscReal));
    
    // set all data to 0

    for (i = 0; i < arr->nelem; i++)
        arr->data[i] = 0.0;

    return arr;
}

/* Copy from an existing array */

array7d * copy_array7d(array7d * arr) {
  
    if (!arr) return NULL;

    // create a new array to hold the copy
  
    array7d * cp = allocate7d(arr->size1, arr->size2, arr->size3, arr->size4, arr->size5, arr->size6, arr->size7);

    // copy array data to cp's data
  
    memcpy(cp->data, arr->data, arr->nelem * sizeof(PetscReal));

    return cp;
}

/* Delete the array */

PetscInt free7d(array7d * arr) {
  
    if (!arr) return -1;
  
    // free arrays's data
    free(arr->data);
  
    // free array itself
    free(arr);
  
    return 0;
}

/* Set the value of an element */

void set_element_7d(array7d* arr, PetscInt i, PetscInt j, PetscInt k, PetscInt l, PetscInt m, PetscInt n, PetscInt o, PetscReal val) {
    arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3 + l*arr->c4 + m*arr->c5 + n*arr->c6 + o*arr->c7] = val; 
}

/* Get the value of an element */

PetscReal get_element_7d(array7d* arr, PetscInt i, PetscInt j, PetscInt k, PetscInt l, PetscInt m, PetscInt n, PetscInt o) {
    return arr->data[i*arr->c1 + j*arr->c2 + k*arr->c3 + l*arr->c4 + m*arr->c5 + n*arr->c6 + o*arr->c7]; 
}

/* Get minimum and maximum of the array */

void min_max_7d(array7d* arr, PetscReal* min, PetscReal* max) {
    
    PetscInt i; 
    
    *min = arr->data[0]; 
    *max = arr->data[0]; 
    
    for (i = 1; i < arr->nelem; ++i) {
        if (arr->data[i] > *max)
            *max = arr->data[i];
        if (arr->data[i] < *min)
            *min = arr->data[i];
    }
}


