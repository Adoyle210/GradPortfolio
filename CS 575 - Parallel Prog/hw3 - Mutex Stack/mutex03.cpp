#include <stdio.h>
#include <math.h>
#include <omp.h>

#define DEBUG			false
#define TIMEOUT			100000
#define FAILED			(-999)

#ifndef NUMN
#define NUMN			32768
#endif

#ifndef USE_MUTEX
#define USE_MUTEX		true
#endif

int		Stack[NUMN];
volatile int	StackPtr = -1;
			// ok to set the StackPtr to an illegal array index since we will never actually use it as -1
			//	because we will increment it before Pushing with it
			// to push: increment stack pointer, then place number there
			// to pop:  take number from there, then decrement stack pointer

bool		WasPopped[NUMN];
int		NumPopErrors;

omp_lock_t	Lock;



void
Push( int n )
{
	StackPtr++;
	Stack[StackPtr] = n;
}


int
Pop( )
{
	// if the stack is empty, give the Push( ) function a chance to put something on the stack:
	int t = 0;
	while( StackPtr < 0  &&  t < TIMEOUT )
		t++;

	// if there is nothing to pop, return;
	if( StackPtr < 0 )
		return FAILED;

	int n = Stack[StackPtr];
	StackPtr--;

	WasPopped[n] = true;
	return n;
}


void
PushAll( )
{
	if (USE_MUTEX) {
		omp_set_lock(&Lock);
	}
	for( int i = 0; i < NUMN; i++ )
	{
		Push( i );
	}
	if (USE_MUTEX) {
		omp_unset_lock(&Lock);
	}
}


void
PopAll( )
{	if (USE_MUTEX) {
		omp_set_lock(&Lock);
	}
	for( int i = 0; i < NUMN; i++ )
	{
		int n = Pop( );
	}
	if (USE_MUTEX) {
		omp_unset_lock(&Lock);
	}
}



int
main( int argc, char *argv[ ] )
{
#ifndef _OPENMP
	fprintf( stderr, "OpenMP is not supported here.\n" );
	return 1;
#endif
	//inizilizing Lock here 
	omp_init_lock(&Lock);


	// this array is here to be sure all the pops actually happened:
	for( int i = 0; i < NUMN; i++ )
	{
		WasPopped[i] = false;
	}

	omp_set_num_threads( 2 );
	
	double time0 = omp_get_wtime( );
	#pragma omp parallel sections
	{
		#pragma omp section
			PushAll( );

		#pragma omp section
			PopAll( );
	}
	double time1 = omp_get_wtime( );

	NumPopErrors = 0;
	for( int i =0; i < NUMN; i++ )
	{
		if( ! WasPopped[i] )
		{
			if( DEBUG )	fprintf( stderr, "%6d wasn't popped\n", i );
			NumPopErrors++;
		}
	}

	char *useMutexString = (char *)"false";
	if( USE_MUTEX )
		useMutexString = (char *)" true";

	fprintf( stderr, "NUMN = %6d , USE_MUTEX = %s , NumPopErrors = %5d = %6.2f%% , Elapsed time = %9.2lf microseconds\n",
		NUMN, useMutexString, NumPopErrors, 100.*(float)NumPopErrors/(float)NUMN, 1000000.*(time1-time0) );

	return 0;
}