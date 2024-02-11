///////////////////////////////////////////////////////////////////////////////////////////////////////
//                                   Code: ParaEllip3d-CFD                                           //
//                                 Author: Dr. Beichuan Yan                                          //
//                                  Email: beichuan.yan@colorado.edu                                 //
//                              Institute: University of Colorado Boulder                            //
///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef BLOCK_H
#define BLOCK_H

// The BLOCK data decomposition method is used to replace the ceil/floor method to avoid small number of CFD grids in maximum boundary processes.
#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,(p),(n)) - 1 )
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH((id),(p),(n)) - BLOCK_LOW((id),(p),(n)) + 1)
#define BLOCK_OWNER(index,p,n) (((p)*((index)+1)-1)/(n))
// n: number of elements/grids
// p: number of processes
// id: index of a process
// index: index of an element/grid (global, of coz)
// BLOCK_LOW gives the first index controlled by process id
// BLOCK_HIGH gives the last index controlled by process id
// BLOCK_SIZE gives the number of elements controlled by process id
// BLOCK_OWNER evaluates to the rank of the process controlling the element of the array

#endif
