#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ns 5
#define BIG 100
#define ITR 100
#define EPSL 0.01
// ITR is the number of replicas of the system to log. One logfile line per itr.
//EPSL is the tolerance of acceptance for close energies. If over by EPSL or less, move.
#define RMAX 1000000000
// RMAX is used by the simple, homebuilt random number generator to make integers from 0 to RMAX-1, I think. Or RMAX.
// ns is the number of spins. this is spin one particles.
// BIG is the number of times to try to improve the state.

// simulation of many spin one particles with T=0 sim annealing in ext field
// so ground state is all down. Init state is all up (worst possible state).

//using namespace std;
unsigned long int rnd(int seed_init)
{  //a simple pseudorandom number generator, so we can reproduce results from any run.
  static unsigned long int seed;
  unsigned long int a = 16807;
  unsigned long int m = 2147483647;
  unsigned long int c = 1013904223;
  if(seed_init)
    seed = seed_init;
  else
    seed = (((a * seed) +c )% m);
  return seed;
}
void vary_trial(int s[ns][ns][ns], int o[ns][ns][ns])
{//Vary an existing state to a new state, s. Each spin has a one tenth, or X if this changes later, probability to get randomized. This may help converge faster.
	int i; int j; int k;
	unsigned long int temp;
	for(i=0;i<ns;i++){
		for(j=0;j<ns;j++){
			for(k=0;k<ns;k++){

		if((rnd(0)%RMAX)<(RMAX/10)){ // about a 1/10 chance to edit a spin
			temp = rnd(0)%RMAX;
			if (temp < (RMAX/3))
		 	 s[i][j][k]=1;
			else {
	    	  if (temp < (RMAX*2/3))
		    	s[i][j][k]=0;
		  	  else
		    	s[i][j][k]=-1;
            
		}}
		else
			s[i][j][k]=o[i][j][k];
	} 
} }}
void gen_trial(int s[ns][ns][ns])
{
	int i; int j; int k;
	unsigned long int temp;
	for(i=0;i<ns;i++){
		for(j=0;j<ns;j++){
			for(k=0;k<ns;k++){
		temp = rnd(0)%RMAX;
		if (temp < (RMAX/3))
		  s[i][j][k]=1;
		else
	      if (temp < (RMAX*2/3))
		    s[i][j][k]=0;
		  else
		    s[i][j][k]=-1;
	}}}
}
void prnt_spins(int spins[ns][ns][ns])
{
	int i; int j; int k;
	for(i=0;i<ns;i++){
		for(j=0; j<ns; j++){
			for(k=0; k<ns; k++){
	  printf("%i, ",spins[i][j][k]);
	  putchar('\n');
	  //printf("%i, ",0x01&spins[i][j+1][k]);
	  //putchar('\n');
	  //printf("%i, ",0x01&spins[i][j][k+1]);

	}}}}

void init_spins(int spins[ns][ns][ns])
{
	int i; int j; int k;
	for(i=0;i<ns;i++){
		for(j=0; j<ns;j++){
			for(k=0; k<ns; k++){
			  spins[i][j][k]=1;
	  		
		}
	}
}
}

void copy_state(int to[ns][ns][ns], int from[ns][ns][ns])
{
	int i; int j; int k;
	for(i=0;i<ns;i++){
		for(j=0; j<ns;j++){
			for(k=0;k<ns;k++){

	  to[i][j][k]=from[i][j][k];//copy array from into the to array, erasing the to array.
	}}}
}
void init_eps(double eps[ns][ns][ns])
{
	int i; int j; int k;
	int r;
	for(i=0;i<ns;i++){
		for (j=0; j<ns; j++){
			for(k=0;k<ns;k++){
		eps[i][j][k] = (rnd(0)%RMAX) * 2.0 / RMAX - 1.0; //should be random from -1 to 1.0.
	
		//eps[i][j][k]=-1;
	}
}
}
		//printf("%lf\n",eps[i]);
}
	
double score_state(int sp[ns][ns][ns], double eps_x[ns][ns][ns], double eps_y[ns][ns][ns], double eps_z[ns][ns][ns])
{
	int i; int j; int k;
	double score=0;

	//eps x
	for(i=0;i<ns-1;i++){
		for(int j=0; j<ns;j++){
			for(int k=0; k<ns; k++){
	  score += eps_x[i+1][j][k] * sp[i+1][j][k]*sp[i][j][k];
			}	
		}
	}
//eps y
	for(i=0;i<ns;i++){
		for(int j=0; j<ns-1;j++){
			for(int k=0; k<ns; k++){
	  score += eps_y[i][j+1][k] * sp[i][j+1][k]*sp[i][j][k];
			}	
		}
	}

	//eps z
	for(i=0;i<ns;i++){
		for(int j=0; j<ns;j++){
			for(int k=0; k<ns-1; k++){
	  score += eps_z[i][j][k+1] * sp[i][j][k+1]*sp[i][j][k];
			}	
		}
	}


	// x wrap around

	for(j=0; j<ns;j++){
		for(k=0;k<ns;k++){
			score += sp[0][j][k]*sp[ns-1][j][k]*eps_x[0][j][k];
		} }


		//y wrap around

		for (i=0; i<ns; i++){
			for(k=0; k<ns; k++){

				score += sp[i][0][k]*sp[i][ns-1][k]*eps_y[i][0][k];
			}
		}


		//z wrap around
		for(i=0;i<ns;i++){
			for(j=0;j<ns;j++){
				score +=sp[i][j][0]*sp[i][j][ns-1]*eps_z[i][j][0];
			}
		}

	//NEED TO WRITE LOOPS HERE TO GET ALL THE LAST TERMS
	//score += sp[ns-1][ns-1][ns-1]; //last paramagnetic term
		if(score != score)
			printf("Ack! score problem in function scoreState!!!\n");
	//score += eps[ns-1][ns-1][ns-1] * sp[ns-1][ns-1][ns-1]*sp[0][0][0];//last ferromagnetic term wraps around periodic boundary conditions.
		printf("score is %lf\n",score);
	return score;
}

double autocorrelation(int spins[ns][ns][ns])
{
 int i; int j; int k;
 int spin_sum=0;
 double average_spin=0;

for(i=0;i<ns;i++){
	for(j=0;j<ns;j++){
		for(k=0;k<ns;k++){
			
			 spin_sum += spins[i][j][k];// Times one, since initially, all spins are up.Eq4 in Kisker, Rieger, 1996.
		
		}
	}

}

	average_spin=((double) spin_sum)/(ns*ns*ns);
	return average_spin;	

}

//TOTAL DELTA ENERGY

double DeltaE(int sp[ns][ns][ns], double eps_x[ns][ns][ns], double eps_y[ns][ns][ns], double eps_z[ns][ns][ns], int i, int j, int k){
        double deltaE=0.0;
      // delta_x
        if((i+1)< ns){
                deltaE += -sp[i][j][k]*sp[i+1][j][k]*eps_x[i+1][j][k] - sp[i][j][k]*sp[i+1][j][k]*eps_x[i+1][j][k];
        }
        else{
                deltaE += -sp[i][j][k]*sp[0][j][k]*eps_x[0][j][k] - sp[i][j][k]*sp[0][j][k]*eps_x[0][j][k];
        }

        if((i-1)>= 0){
                deltaE += -sp[i][j][k]*sp[i-1][j][k]*eps_x[i-1][j][k] - sp[i][j][k]*sp[i-1][j][k]*eps_x[i-1][j][k];
        }
        else{
                deltaE += -sp[i][j][k]*sp[ns-1][j][k]*eps_x[ns-1][j][k] - sp[i][j][k]*sp[ns-1][j][k]*eps_x[ns-1][j][k];
        }
      // delta_y
        if((j+1)< ns){
                deltaE += -sp[i][j][k]*sp[i][j+1][k]*eps_y[i][j+1][k] - sp[i][j][k]*sp[i][j+1][k]*eps_y[i][j+1][k];
        }
        else{
                deltaE += -sp[i][j][k]*sp[i][0][k]*eps_y[i][0][k] - sp[i][j][k]*sp[i][0][k]*eps_y[i][0][k];
        }
        if((j-1)>= 0){
                deltaE += -sp[i][j][k]*sp[i][j-1][k]*eps_y[i][j-1][k] - sp[i][j][k]*sp[i][j-1][k]*eps_y[i][j-1][k];
        }
        else{
                deltaE += -sp[i][j][k]*sp[i][ns-1][k]*eps_y[i][ns-1][k] - sp[i][j][k]*sp[i][ns-1][k]*eps_y[i][ns-1][k];
        }
      //delta_z
        if((k+1)< ns){
                deltaE += -sp[i][j][k]*sp[i][j][k+1]*eps_z[i][j][k+1] - sp[i][j][k]*sp[i][j][k+1]*eps_z[i][j][k+1];
        }

        else{
                deltaE += -sp[i][j][k]*sp[i][j][0]*eps_z[i][j][0] - sp[i][j][k]*sp[i][j][0]*eps_z[i][j][0];
        }
        if((k-1)>= 0){
                deltaE += -sp[i][j][k]*sp[i][j][k-1]*eps_z[i][j][k-1] - sp[i][j][k]*sp[i][j][k-1]*eps_z[i][j][k-1];
        }
        else{
                deltaE += -sp[i][j][k]*sp[i][j][ns-1]*eps_z[i][j][ns-1] - sp[i][j][k]*sp[i][j][ns-1]*eps_z[i][j][ns-1];
        }
        return deltaE; }


void updtSpin(int sp[ns][ns][ns], double eps_x[ns][ns][ns], double eps_y[ns][ns][ns], double eps_z[ns][ns][ns]){

//Loop over all i's
        for(int a=0;a<ns;a++){
        //Loop over all j's

                for(int b=0;b<ns;b++){

                //Loop over all k's
                        for(int c=0;c<ns;c++){
                        	float k=0.8; //This is something like kBT yes, should be 0.8?
                        //Delta E function (output: )
                        double jf = DeltaE(sp, eps_x, eps_y, eps_z, a, b, c);

                        //Calcualte flip Probability (Reiger) // w(-si->si) = 1/(1+exp(delE/T))
                        double mf =1/(1+exp((jf)/k)); 

			//Arbitraty flip at 0.5 //replaceable to change threshold w/ rng
                        int mm = round(mf); //round to 1 or 0 //Should we not use rand num generation to check if this occurs?

                                //printf("Spin[ %i ] [%i] [%i] will ",a,b,c);//<<cout<<"Spin ["<<a<<"]["<<b<<"]["<<c<<"] will ";

			        if(mm == 1){
                                        //cout<<"flip spin from: "<<sp[a][b][c];
			        					//printf("flip spin from: %i",sp[a][b][c]);
                                        //flip selected spin
                                        sp[a][b][c]*=(-1);
					//cout<<"  to: "<<sp[a][b][c]<<"."<<endl;
                    //printf(" to: %i.\n", sp[a][b][c]);
                                }
                                else if(mm == 0){
                                        //cout<<"not flip spin."<<endl;
                                		//printf("not flip spin.\n");
                                }

                                else{
                                        //cout<<"break this code completely."<<endl;
                                		printf("danger! this should not occur!\n");
                                }
}}}
//END FUNCTION
}




int main()
{
	double kT; // this variable is initialized within the ITR loop below.
	int spins[ns][ns][ns];
	FILE *f1;	
	int i,j;
    double avgscore=0.0; // we will accumulate total end score of each run, then average it within this variable.
	rnd(100151); // initializes the random number generator. After this, call rnd(0) with 0 so that it does not re-start.
 	int spin_sum;
 	double average_spin;

	double eps_x[ns][ns][ns];
	double eps_y[ns][ns][ns];
	double eps_z[ns][ns][ns];	
	f1 = fopen("Autocorrelation.txt", "wt");
    if (f1 == NULL)
	{
	    fprintf(stderr,"Error opening file!\n");
	    exit(1);
    }
	for(j=0;j<ITR;j++){
        init_spins(spins);
        init_eps(eps_x);
        init_eps(eps_y);
        init_eps(eps_z);
        kT = 0.8; //restart kT up to starting value.
        average_spin= autocorrelation(spins);

        for(i=0;i<BIG;i++){
            updtSpin(spins,eps_x,eps_y,eps_z);
            average_spin=autocorrelation(spins);
            fprintf(f1, "%i %lf\n",i, average_spin);// int i is effecitively the time variable.
        }//end of for i<BIG loop...looping over the trials
        
    }//end of the for loop...looping over number of iterations.

	fclose(f1);  
}
