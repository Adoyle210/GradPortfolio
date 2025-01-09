// Monte Carlo simulation of pins-fitting-into-holes-in-a-plate

// system includes
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <assert.h>
#include <malloc.h>

// CUDA runtime
#include <cuda_runtime.h>

// Helper functions and utilities to work with CUDA
#include "helper_functions.h"
#include "helper_cuda.h"

// setting the number of trials in the monte carlo simulation:
#ifndef NUMTRIALS
#define NUMTRIALS	50000
#endif

#ifndef BLOCKSIZE
#define BLOCKSIZE	64     // number of threads per block
#endif

#define NUMBLOCKS	( NUMTRIALS / BLOCKSIZE )


//#define CSV

// the pins; numbers are constants:
#define PinAx	3.0f
#define PinAy	4.0f
#define PinAr	2.0f

#define PinBx	4.0f
#define PinBy	5.0f
#define PinBr	2.0f

#define PinCx	5.0f
#define PinCy	4.0f
#define PinCr	2.0f


// ranges for the random numbers:

#define PROJECT1

#ifdef PROJECT1
const float HoleAx   =  2.90f;
const float HoleAy   =  4.10f;
const float HoleAr   =  2.20f;
const float HoleAxPM =	0.20f;
const float HoleAyPM =	0.20f;
const float HoleArPM =	0.20f;

const float HoleBx   =  4.10f;
const float HoleBy   =  4.90f;
const float HoleBr   =  2.20f;
const float HoleBxPM =	0.10f;
const float HoleByPM =	0.10f;
const float HoleBrPM =	0.20f;

const float HoleCx   =  5.00f;
const float HoleCy   =  4.00f;
const float HoleCr   =  2.20f;
const float HoleCxPM =	0.10f;
const float HoleCyPM =	0.05f;
const float HoleCrPM =	0.20f;
#else
const float HoleAx   =  2.90f;
const float HoleAy   =  4.10f;
const float HoleAr   =  2.40f;
const float HoleAxPM =	0.20f;
const float HoleAyPM =	0.20f;
const float HoleArPM =	0.30f;

const float HoleBx   =  4.10f;
const float HoleBy   =  4.90f;
const float HoleBr   =  2.40f;
const float HoleBxPM =	0.10f;
const float HoleByPM =	0.10f;
const float HoleBrPM =	0.30f;

const float HoleCx   =  5.00f;
const float HoleCy   =  4.00f;
const float HoleCr   =  2.40f;
const float HoleCxPM =	0.10f;
const float HoleCyPM =	0.05f;
const float HoleCrPM =	0.30f;
#endif



float
Ranf( float low, float high )
{
        float r = (float) rand();               // 0 - RAND_MAX
        float t = r  /  (float) RAND_MAX;       // 0. - 1.

        return   low  +  t * ( high - low );
}


// call this if you want to force your program to use
// a different random number sequence every time you run it:
void
TimeOfDaySeed( )
{
	struct tm y2k = { 0 };
	y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
	y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;

	time_t  timer;
	time( &timer );
	double seconds = difftime( timer, mktime(&y2k) );
	unsigned int seed = (unsigned int)( 1000.*seconds );    // milliseconds
	srand( seed );
}


void
CudaCheckError()
{
	cudaError_t e = cudaGetLastError();
	if(e != cudaSuccess)
	{
   		fprintf( stderr, "Cuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e));
	}
}


__device__
float
Sqr( float x )
{
	return x*x;
}


__device__
float
Length( float dx, float dy )
{
	return  sqrt( Sqr(dx) + Sqr(dy) );
}


#define IN
#define OUT

__global__
void
MonteCarlo(
	IN float *dholeaxs, IN float *dholeays, IN float *dholears,
	IN float *dholebxs, IN float *dholebys, IN float *dholebrs,
	IN float *dholecxs, IN float *dholecys, IN float *dholecrs,
	OUT int *dsuccesses )
{
        //unsigned int numItems = blockDim.x;
        //unsigned int wgNum    = blockIdx.x;
        //unsigned int tnum     = threadIdx.x;
        unsigned int gid      = blockIdx.x*blockDim.x + threadIdx.x;

        dsuccesses[gid] = 0;

	// randomize everything:
	float holeax = dholeaxs[gid];
	float holeay = dholeays[gid];
	float holear = dholears[gid];

	float holebx = dholebxs[gid];
	float holeby = dholebys[gid];
	float holebr = dholebrs[gid];

	float holecx = dholecxs[gid];
	float holecy = dholecys[gid];
	float holecr = dholecrs[gid];

	float da = Length( (PinAx - holeax),  (PinAy - holeay) );
	if( (da + PinAr) <= holear)
	{
		float db = Length( (PinBx - holebx),  (PinBy - holeby) );
		if( (db + PinBr) <= holebr)
		{
			float dc = Length( (PinCx - holecx),  (PinCy - holecy) );
			if( (dc + PinCr) <= holecr) 
        			dsuccesses[gid] = 1;
		}
	}
}


int
main( int argc, char *argv[ ] )
{
        TimeOfDaySeed( );               // seed the random number generator

        // better to define these here so that the rand() calls don't get into the thread timing:
        float *hholeaxs  = new float [NUMTRIALS];
        float *hholeays  = new float [NUMTRIALS];
        float *hholears  = new float [NUMTRIALS];

        float *hholebxs  = new float [NUMTRIALS];
        float *hholebys  = new float [NUMTRIALS];
        float *hholebrs  = new float [NUMTRIALS];

        float *hholecxs  = new float [NUMTRIALS];
        float *hholecys  = new float [NUMTRIALS];
        float *hholecrs  = new float [NUMTRIALS];

		int *hsuccesses = new int [NUMTRIALS];

        // fill the random-value arrays:
        for( int n = 0; n < NUMTRIALS; n++ )
        {
                hholeaxs[n]  = Ranf(  HoleAx-HoleAxPM,  HoleAx+HoleAxPM );
                hholeays[n]  = Ranf(  HoleAy-HoleAyPM,  HoleAy+HoleAyPM );
                hholears[n]  = Ranf(  HoleAr-HoleArPM,  HoleAr+HoleArPM );

                hholebxs[n]  = Ranf(  HoleBx-HoleBxPM,  HoleBx+HoleBxPM );
                hholebys[n]  = Ranf(  HoleBy-HoleByPM,  HoleBy+HoleByPM );
                hholebrs[n]  = Ranf(  HoleBr-HoleBrPM,  HoleBr+HoleBrPM );

                hholecxs[n]  = Ranf(  HoleCx-HoleCxPM,  HoleCx+HoleCxPM );
                hholecys[n]  = Ranf(  HoleCy-HoleCyPM,  HoleCy+HoleCyPM );
                hholecrs[n]  = Ranf(  HoleCr-HoleCrPM,  HoleCr+HoleCrPM );
        }

	// allocate device memory:

        float *dholeaxs, *dholeays, *dholears;
        float *dholebxs, *dholebys, *dholebrs;
        float *dholecxs, *dholecys, *dholecrs;
        int   *dsuccesses;

	cudaMalloc( (void**)(&dholeaxs), NUMTRIALS * sizeof(float) );
	cudaMalloc( (void**)(&dholeays), NUMTRIALS * sizeof(float) );
	cudaMalloc( (void**)(&dholears), NUMTRIALS * sizeof(float) );

	cudaMalloc( (void**)(&dholebxs), NUMTRIALS * sizeof(float) );
	cudaMalloc( (void**)(&dholebys), NUMTRIALS * sizeof(float) );
	cudaMalloc( (void**)(&dholebrs), NUMTRIALS * sizeof(float) );

	cudaMalloc( (void**)(&dholecxs), NUMTRIALS * sizeof(float) );
	cudaMalloc( (void**)(&dholecys), NUMTRIALS * sizeof(float) );
	cudaMalloc( (void**)(&dholecrs), NUMTRIALS * sizeof(float) );

	cudaMalloc( (void**)(&dsuccesses), NUMTRIALS * sizeof(int) ); //changed to int :)

	CudaCheckError( );


	// copy host memory to the device:

	cudaMemcpy( dholeaxs,  hholeaxs,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dholeays,  hholeays,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dholears,  hholears,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );

	cudaMemcpy( dholebxs,  hholebxs,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dholebys,  hholebys,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dholebrs,  hholebrs,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );

	cudaMemcpy( dholecxs,  hholecxs,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dholecys,  hholecys,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );
	cudaMemcpy( dholecrs,  hholecrs,  NUMTRIALS*sizeof(float), cudaMemcpyHostToDevice );

	CudaCheckError( );


	// setup the execution parameters:
	dim3 threads(BLOCKSIZE, 1, 1 );
	dim3 grid(NUMBLOCKS, 1, 1 );

	// create and start timer
	cudaDeviceSynchronize( );

	// allocate CUDA events that we'll use for timing:
	cudaEvent_t start, stop;
	cudaEventCreate( &start );
	CudaCheckError( );
	cudaEventCreate( &stop );
	CudaCheckError( );

	// record the start event:
	cudaEventRecord( start, NULL );
	CudaCheckError( );

	// execute the kernel:
	MonteCarlo<<< grid, threads >>>( dholeaxs, dholeays, dholears, dholebxs, dholebys, dholebrs, dholecxs, dholecys, dholecrs, dsuccesses);

	// record the stop event:
	cudaEventRecord( stop, NULL );

	// wait for the stop event to complete:
	cudaEventSynchronize( stop );

	float msecTotal = 0.0f;
	cudaEventElapsedTime( &msecTotal, start, stop );
	CudaCheckError( );

	// copy result from the device to the host:
	cudaMemcpy(hsuccesses , dsuccesses, NUMTRIALS *sizeof(int), cudaMemcpyDeviceToHost );
	CudaCheckError( );

	// compute the sum :
	int numSuccesses = 0;
	for(int i = 0; i < NUMTRIALS; i++ )
	{
		numSuccesses += hsuccesses[i];
	}

	float probability = (float)numSuccesses / (float)NUMTRIALS;

	// compute and print the performance:
	double secondsTotal = 0.001 * (double)msecTotal;
	double trialsPerSecond = (float)NUMTRIALS / secondsTotal;
	double megaTrialsPerSecond = trialsPerSecond / 1000000.;

#ifdef CSV
	fprintf( stderr, "%10d , %8d , %10.4lf , %6.2f\n",
		NUMTRIALS, BLOCKSIZE, megaTrialsPerSecond, 100.*probability );
#else
	fprintf( stderr, "Number of Trials = %10d, Blocksize = %8d, MegaTrials/Second = %10.4lf, Probability = %6.2f%%\n",
		NUMTRIALS, BLOCKSIZE, megaTrialsPerSecond, 100.*probability );
#endif

	// clean up host memory:
	delete [ ] hholeaxs;	 delete [ ] hholeays;	 delete [ ] hholears;	
	delete [ ] hholebxs;	 delete [ ] hholebys;	 delete [ ] hholebrs;	
	delete [ ] hholecxs;	 delete [ ] hholecys;	 delete [ ] hholecrs;	
	delete [ ] hsuccesses;

	// clean up device memory:
	cudaFree( dholeaxs );	 cudaFree( dholeays );	 cudaFree( dholears );	
	cudaFree( dholebxs );	 cudaFree( dholebys );	 cudaFree( dholebrs );	
	cudaFree( dholecxs );	 cudaFree( dholecys );	 cudaFree( dholecrs );	
	cudaFree( dsuccesses );

	CudaCheckError( );

	return 0;
}
