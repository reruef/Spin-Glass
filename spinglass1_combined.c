#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ns 5
#define BIG 1000
#define ITR 1000
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
	  printf("%i, ",0x01&spins[i][j][k]);
	  //putchar('\n');
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
int init_eps(double eps[ns][ns][ns])
{
	int i; int j; int k;
	int r;
	/*for(i=0;i<ns;i++){
		for (j=0; j<ns; j++){
			for(k=0;k<ns;k++){
		eps[i][j][k] = (rnd(0)%RMAX) * 2.0 / RMAX - 1.0; //should be random from -1 to 1.0.
	*/
		eps[i][j][k]=-1;
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
	//score += eps[ns-1][ns-1][ns-1] * sp[ns-1][ns-1][ns-1]*sp[0][0][0];//last ferromagnetic term wraps around periodic boundary conditions.
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
			
			 spin_sum += spins[i][j][k];
		
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
                deltaE += -sp[i][j][k]*sp[0][j][k]*eps_x[i+1][j][k] - sp[i][j][k]*sp[0][j][k]*eps_x[i+1][j][k];
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
                        	float k=3.3;
                        //Delta E function (output: )
                        double jf = DeltaE(sp, eps_x, eps_y, eps_z, a, b, c);

                        //Calcualte flip Probability (Reiger) // w(-si->si) = 1/(1+exp(delE/T))
                        double mf =1/(1+exp((jf)/k));

			//Arbitraty flip at 0.5 //replaceable to change threshold w/ rng
                        int mm = round(mf); //round to 1 or 0

                                printf("Spin[ %i ] [%i] [%i] will",a,b,c);//<<cout<<"Spin ["<<a<<"]["<<b<<"]["<<c<<"] will ";

			        if(mm == 1){
                                        //cout<<"flip spin from: "<<sp[a][b][c];
			        					printf("flip spin from: %i",sp[a][b][c]);
                                        //flip selected spin
                                        sp[a][b][c]*=(-1);
					//cout<<"  to: "<<sp[a][b][c]<<"."<<endl;
                    printf(" to: %i.\n", sp[a][b][c]);
                                }
                                else if(mm == 0){
                                        //cout<<"not flip spin."<<endl;
                                		printf("not flip spin.\n");
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
	double stata[BIG];
	//double     eps[ns][ns][ns]={{0.51,0.27,0.11,-0.29,-0.55,-0.81,-0.92,-0.33,0.84,0.35, -0.51,0.02,0.53,-0.79,0.94,0.71,-0.72,0.91,0.24,-0.51}, {0.51,0.27,0.11,-0.29,-0.55,-0.81,-0.92,-0.33,0.84,0.35,  -0.51,0.02,0.53,-0.79,0.94,0.71,-0.72,0.91,0.24,-0.51}, {0.51,0.27,0.11,-0.29,-0.55,-0.81,-0.92,-0.33,0.84,0.35,  -0.51,0.02,0.53,-0.79,0.94,0.71,-0.72,0.91,0.24,-0.51}};
	double kT; // this variable is initialized within the ITR loop below.
	int spins[ns][ns][ns];
	int trial_state[ns][ns][ns]; //to compare with the current state, spins.
	FILE * fout, *log1, *log2, *log3, *f1;	
	int i,j;
	int ngm, nnm;
	double loscore; //current lowest score.
	double trialscore; //score of the trial state.
    double avgscore=0.0; // we will accumulate total end score of each run, then average it within this variable.
	rnd(100151); // initializes the random number generator. After this, call rnd(0) with 0 so that it does not re-start.
	//init_eps(eps);

 	int spin_sum;
 	double average_spin;

	double eps_x[ns][ns][ns];
	double eps_y[ns][ns][ns];
	double eps_z[ns][ns][ns];	
    /*printf("EPS array init is: ");
	for(i=0;i<ns;i++){
		//eps[i] = (rnd(0)%RMAX) * 2.0 / RMAX - 1.0; //should be random from -1 to 1.0.
		printf("%lf ",eps[i]);
	}*/
    //printf("kT is %lf num spins is %i iterations %i\n", kT, ns, ITR);
	



	fout = fopen("log_3D.txt","wt");
	for(i=0;i<BIG;i++){
		stata[i]=0;
	}
	ngm=0;nnm=0;
    init_spins(spins);
  
        init_eps(eps_x);
        init_eps(eps_y);
        init_eps(eps_z);



		f1 = fopen("Autocorrelation.txt", "wt");
			if (f1 == NULL)
			{
			    fprintf(stderr,"Error opening file!\n");
			    exit(1);
				}

									
	
	for(j=0;j<ITR;j++){
        
        //init_eps(eps);
        loscore = score_state(spins,eps_x,eps_y,eps_z); // initialize the scoring.

        kT = 3.3; //restart kT up to starting value. This will decline soon enough. Sim annealing.
        average_spin= autocorrelation(spins);
      // printf("%lf\n", loscore);
      // exit(0);

        for(i=0;i<BIG;i++){
            vary_trial(trial_state, spins);
            //prnt_spins(trial_state);
            trialscore = score_state(trial_state,eps_x, eps_y, eps_z);

            if(trialscore < (loscore+EPSL) ){
                //printf("Improved Score from %i to %i on run %i.\n", loscore, trialscore, i);
                copy_state(spins, trial_state); // now the new best state is the trial.
                loscore = trialscore; //update the lowest score to equal that of trial.
                ngm++;
            }
            else{
            	if(rnd(0)%RMAX < RMAX*exp((loscore-trialscore)/kT)){
            		copy_state(spins, trial_state);
            		loscore=trialscore;
            		nnm++;
				}
			}
            stata[i]+=loscore;
            kT = kT / 1.01; // Simulated annealing...lower the temperature. 
             //fprintf(fout, "%i %i \n", i, loscore);
           
        }//end of for i<BIG loop...looping over the trials
        //fprintf(fout,"%.2lf \n", loscore); // we are printing the lowest energy achieved into the log file. Can make histogram of them.
        avgscore += loscore; //add it up to form average in a moment below.
        average_spin=autocorrelation(spins);
        fprintf(f1, "%lf\n", average_spin);
    }//end of the for loop...looping over number of iterations.
    //for(i=0;i<BIG;i++)
      fprintf(fout, "%i %lf\n", i, stata[i]/ITR); //prints trial number, like time, and the average energy at that time.

	fclose(f1);  
    fclose(fout);
    printf("NGM %i NNM %i average E = %lf and E100 = %lf and E200 = %lf\n",ngm,nnm,avgscore/ITR, stata[100]/ITR, stata[200]/ITR);
    // ngm is the number of moves to a better state.
    // nnm is a "sideways" move, meaning it was allowed to move by the temperature to a worse state to expore the system.
}
