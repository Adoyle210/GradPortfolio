#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define CSV


int	NowYear;		// 2023 - 2028
int	NowMonth;		// 0 - 11

float	NowPrecip;		// inches of rain per month
float	NowTemp;		// temperature this month
float	NowHeight;		// rye grass height in inches
int	NowNumRabbits;		// number of rabbits in the current population
int	NowNumFoxes;		// number of foxes in the current population


const float RYEGRASS_GROWS_PER_MONTH =		20.0;
const float ONE_RABBITS_EATS_PER_MONTH =	 1.0;

const float AVG_PRECIP_PER_MONTH =	       12.0;	// average
const float AMP_PRECIP_PER_MONTH =		4.0;	// plus or minus
const float RANDOM_PRECIP =			2.0;	// plus or minus noise

const float AVG_TEMP =				60.0;	// average
const float AMP_TEMP =				20.0;	// plus or minus
const float RANDOM_TEMP =			10.0;	// plus or minus noise

const float MIDTEMP =				60.0;
const float MIDPRECIP =				14.0;

unsigned int seed = 0;
// float x = Ranf( -1.f, 1.f );

omp_lock_t	Lock;
volatile int	NumInThreadTeam;
volatile int	NumAtBarrier;
volatile int	NumGone;

//function declatarions: 
void WaitBarrier( );


//--------------------------------------------------------------------------
//MATH STUFF 
//--------------------------------------------------------------------------
float Ranf( float low, float high )
{
    float r = (float) rand( );              // 0 - RAND_MAX

    return(   low  +  r * ( high - low ) / (float)RAND_MAX   );
}

float
Ranf( unsigned int *seedp,  float low, float high )
{
    float r = (float) rand_r( seedp );              // 0 - RAND_MAX

    return(   low  +  r * ( high - low ) / (float)RAND_MAX   );
}

float x = Ranf( &seed, -1.f, 1.f );

float Sqr( float x )
{
        return x*x;
}

//--------------------------------------------------------------------------
// THRED FUNCTIONS 
//--------------------------------------------------------------------------
void Rabbits()
{
    while( NowYear < 2029 )
    {
        // compute a temporary next-value for this quantity
        // based on the current state of the simulation:

        //calcuatuing rabbits growth and decay
        int nextNumRabbits = NowNumRabbits;
        int carryingCapacity = (int)( NowHeight );
        if( nextNumRabbits < carryingCapacity )
            nextNumRabbits++;
        else
            if( nextNumRabbits > carryingCapacity )
                    nextNumRabbits--;

        //clamping rabbits at zero
        if( nextNumRabbits < 0 )
            nextNumRabbits = 0;

        // DoneComputing barrier:
        WaitBarrier( );	

        NowNumRabbits = nextNumRabbits;

        // DoneAssigning barrier:
        WaitBarrier( );	
        

        // DonePrinting barrier:
        WaitBarrier( );	
    }
}

void RyeGrass()
{
    while( NowYear < 2029 )
    {
        // compute a temporary next-value for this quantity
        // based on the current state of the simulation:
        //temp
        float ang = (  30.*(float)NowMonth + 15.  ) * ( M_PI / 180. );

        float temp = AVG_TEMP - AMP_TEMP * cos( ang );
        NowTemp = temp + Ranf( &seed, -RANDOM_TEMP, RANDOM_TEMP );

        float precip = AVG_PRECIP_PER_MONTH + AMP_PRECIP_PER_MONTH * sin( ang );
        NowPrecip = precip + Ranf( &seed,  -RANDOM_PRECIP, RANDOM_PRECIP );
        if( NowPrecip < 0. )
	        NowPrecip = 0.;

        // calucating weather 
        float tempFactor = exp(   -Sqr(  ( NowTemp - MIDTEMP ) / 10.  )   );
        float precipFactor = exp(   -Sqr(  ( NowPrecip - MIDPRECIP ) / 10.  )   );

        float nextHeight = NowHeight;
        nextHeight += tempFactor * precipFactor * RYEGRASS_GROWS_PER_MONTH;
        nextHeight -= (float)NowNumRabbits * ONE_RABBITS_EATS_PER_MONTH;
        if( nextHeight < 0. ) nextHeight = 0.;

        // DoneComputing barrier:
        WaitBarrier( );	

        //Copy the next state into the now varibles
        NowHeight = nextHeight;
        
        // DoneAssigning barrier:
        WaitBarrier( );	
        

        // DonePrinting barrier:
        WaitBarrier( );	
    }


    //clamping next height @zero 
    
}

void Watcher()
{
    while( NowYear < 2029 )
    {
        // compute a temporary next-value for this quantity
        // based on the current state of the simulation:

        // DoneComputing barrier:
        WaitBarrier( );	
        
        // DoneAssigning barrier:
        WaitBarrier( );	
        
        //printing current status:
        #ifdef CSV
            fprintf(stderr, "%d,%f,%f,%f,%d,%d\n", NowMonth, NowPrecip, NowTemp, NowHeight, NowNumRabbits, NowNumFoxes );
        #else
            fprintf(stderr, "Month : %8d  ; precip = %6.2f ; temp = %6.2lf ; height = %f ; rabbits = %d ; foxes = %d\n",
                NowMonth, NowPrecip, NowTemp, NowHeight, NowNumRabbits, NowNumFoxes);
        #endif

        //Updating year & months 
        if (NowMonth == 11){ 
            NowYear ++;
            NowMonth = 0;
        }
        else{
            NowMonth ++;
        }
        // DonePrinting barrier:
        WaitBarrier( );	
    }

}

void MyAgent()
{
    while( NowYear < 2029 )
    {
        // compute a temporary next-value for this quantity
        // based on the current state of the simulation:
        //calcuatuing rabbits growth and decay
        int nextNumFoxes = NowNumFoxes;
        int carryingCapacity = (int)( NowNumRabbits );
        if( nextNumFoxes < carryingCapacity )
            nextNumFoxes++;
        else
            if( nextNumFoxes > carryingCapacity )
                nextNumFoxes--;

        //clamping rabbits at zero
        if( nextNumFoxes < 0 )
            nextNumFoxes = 0;

        // DoneComputing barrier:
        WaitBarrier( );	

        NowNumFoxes = nextNumFoxes;

        // DoneAssigning barrier:
        WaitBarrier( );	
        

        // DonePrinting barrier:
        WaitBarrier( );	
    }
}


//--------------------------------------------------------------------------
//RUNNING THINGS
//--------------------------------------------------------------------------     

// specify how many threads will be in the barrier:
//	(also init's the Lock)
void InitBarrier( int n )
{
    NumInThreadTeam = n;
    NumAtBarrier = 0;
	omp_init_lock( &Lock );
}

void WaitBarrier( )
{
        omp_set_lock( &Lock );
        {
                NumAtBarrier++;
                if( NumAtBarrier == NumInThreadTeam )
                {
                        NumGone = 0;
                        NumAtBarrier = 0;
                        // let all other threads get back to what they were doing
			            // before this one unlocks, knowing that they might immediately
			            // call WaitBarrier( ) again:
                        while( NumGone != NumInThreadTeam-1 );
                        omp_unset_lock( &Lock );
                        return;
                }
        }
        omp_unset_lock( &Lock );

        while( NumAtBarrier != 0 );	// this waits for the nth thread to arrive

        #pragma omp atomic
        NumGone++;			// this flags how many threads have returned
}

int main()
{
    #ifdef _OPENMP
	//fprintf( stderr, "OpenMP is supported -- version = %d\n", _OPENMP );
    #else
        fprintf( stderr, "No OpenMP support!\n" );
        return 1;
    #endif
    
    // starting date and time:
    NowMonth =    0;
    NowYear  = 2023;

    // starting state (feel free to change this if you want):
    NowNumRabbits = 1;
    NowNumFoxes = 1;
    NowHeight =  5.;    


    omp_set_num_threads( 4 );	// or 4
    InitBarrier( 4 );		// or 4

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            Rabbits( );
        }

        #pragma omp section
        {
            RyeGrass( );
        }

        #pragma omp section
        {
            Watcher( );
        }

        #pragma omp section
        {
            MyAgent( );	// your own
        }

        // implied barrier -- all functions must return in order
        // to allow any of them to get past her
    } 


}