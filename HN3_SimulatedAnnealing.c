/*
  gcc Program.c -lgsl -lgslcblas -lm -o program.x
  
  Key words: simulated annealing, constraint optimization, complex networks, regular graph, disordered system
  
  This code is to use Simulated Annealing to
  find the dynamics and the ground state of the lattice gas problem 
	with at-most-1-neighbor density constraint in HN5;
	Periodic Boundary Condition is used in this code; 
	The output is a file with chemical potential and corresponding average states.
  
  There are 2 other versions of this method for HN5 and HNNP, 
  in which the data structure is a little different.
  
  For more information of the model and method:
  http://arxiv.org/abs/1409.8313
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/types.h>
#include<term.h>
#include<ncurses.h>
#include<gsl/gsl_rng.h>
#include<string.h>
typedef unsigned short int Natural;
#define T_MIN 0.001
#define T_STEP 0.005
#define T_L 2014

int L = 0, K = 0; 	// L: system size; K: L = 2^K
double Initial_Steps = 0.0, Total_Steps = 0.0, Total_MC_Steps = 0.0, Waste_Steps = 0.0; 
int CurrentEnergy = 0, CurrentMag, MinEnergyFound = 0, MinMagFound = 0, InSteps = 0;

gsl_rng *rnd;	// RNG: random number

/* 	integer power of an integer (with 0^0 = 1) 	*/
int Intpow(int base, int power){
  int x=1, i;
  for(i=0; i<power; i++) {x *= base;}
  return(x);
}


/* 	integer log2 of an integer	*/
int Intlog2(int num){
  int power=0;
  while(num>1){
    num=num/2;
    power++;
  }
  return(power);
}

int sites_print(short sites[]){
    int i;
    for(i = 0; i < L; i++){
        printf("%-3d", i+1);
    }
    printf("\n");
    for(i = 0; i < L; i++){
        printf("%-3d", sites[i]);
    }
    printf("\n CurrentEnergy = %d; CurrentMag = %d \n", CurrentEnergy, CurrentMag);
    return(0);
}

/* 	Find the level i of site n and sequential number j	*/
int *Findij(int n){
  static int ij[2];
  ij[0]=1;
  ij[1]=0;
  while(n%2==0) {
    n=n/2;
    ij[0]++;
  }
  ij[1]=(n-1)/2;
  return(ij);
}

/*	Find the neighbors of site n in a network of length L with Periodic Boundary Conditions	  */
/*	The input is 1, 2, ...,L; output is a pointer to an array of 3 */
/*	e.g. Findneighbor(1, 8), returns *neighbors = [2, 3, 8]	       */
int * FindNeighbors(n){ 
  	static int neighbors[3];
  	int *ij, halfL=L/2;
  	if(n<1 || n>L){
    		fprintf(stderr,"\nError:  Cannot find the neighbor of site %d in HN3 with length %d! \n", n, L);
  	}

  	if(n==1){
    		neighbors[0]=2;
    		neighbors[1]=3;
    		neighbors[2]=L;
  	}
  	else if(n==halfL){
    		neighbors[0]=n-1;
    		neighbors[1]=n+1;
    		neighbors[2]=L;
  	}
  	else if(n==L){
    		neighbors[0]=1;
    		neighbors[1]=halfL;
    		neighbors[2]=L-1;
  	}
  	else{
    		ij=Findij(n);
    		neighbors[0]=n-1;
    		neighbors[1]=n+1;
    		if(ij[1]%2==0){
      			neighbors[2]=Intpow(2,ij[0]-1)*(2*(ij[1]+1)+1);
    		}
   		else{
      			neighbors[2]=Intpow(2,ij[0]-1)*(2*(ij[1]-1)+1);
    		}
  	}
  	return(neighbors);
}

/*	Initialize the StateInfo, neighbors and other parameters*/
void Initialize(short sites[], int neighbors[], double temps[]){
  int i;
  int *neigh;

 	for(i=0; i<L; i++){
		sites[i] = 1;
   	neigh = FindNeighbors(i+1);
   	neighbors[3*i] = neigh[0]-1;
    neighbors[3*i+1] = neigh[1]-1;
    neighbors[3*i+2] = neigh[2]-1;
	}
	temps[0] = T_STEP * (T_L-1) + T_MIN;
	for(i = 1; i < T_L; i++){
		temps[i] = temps[i-1] - T_STEP;
	}
  Initial_Steps = 10000.0*L;
  CurrentEnergy = 3 * L / 2;
  CurrentMag = L * sites[0];
  MinEnergyFound = CurrentEnergy;
  MinMagFound = CurrentMag;
	InSteps = 0;
	Total_MC_Steps = 0.0;
	Total_Steps = 0.0;
  Waste_Steps = 0.0;
}


/* Update the system in a certain temperature "temp_now"*/
int Update(short sites[], int neighbors[], double temp_now){
	int i, dE=0;
	i=(int)(L*gsl_rng_uniform(rnd));
	dE = -2 * sites[i] * (sites[neighbors[3*i]] + sites[neighbors[3*i + 1]] + sites[neighbors[3*i + 2]]);
  //printf("dE = %.2e \n", temp_now);

  if(dE < 0 || gsl_rng_uniform(rnd) < exp(- dE / temp_now)){
    //if(dE>0)  printf("dE = %d \n", dE);
    sites[i] = -sites[i];
    CurrentEnergy = CurrentEnergy + dE;
    CurrentMag = CurrentMag + 2 * sites[i];
    if (CurrentEnergy < MinEnergyFound) MinEnergyFound = CurrentEnergy;
    if (CurrentMag < MinMagFound) MinMagFound = CurrentMag;
  }
  else{
    Waste_Steps = Waste_Steps + 1.0;
  }
  Total_Steps = Total_Steps + 1.0;
  InSteps = InSteps + 1;
  if(InSteps == L){
      Total_MC_Steps = Total_MC_Steps + 1.0;
      InSteps = 0;
      return(1);
  }
	return(0);
}

/* Check the current energy state*/
int Check_CurrentEnergy(short sites[], int neighbors[]){
  int i, energy = 0;
  for(i = 0; i < L; i++){
      energy = energy + sites[i] * (sites[neighbors[3*i]] + sites[neighbors[3*i + 1]] + sites[neighbors[3*i + 2]]);
  }
  energy = energy / 2;
  return(energy);
}

int Check_CurrentMag(short sites[]){
  int i, mag = 0;
  for(i = 0; i < L; i++){
    mag = mag + sites[i];
  }
  return(mag);
}


/*	Function to read input parameters	*/
void Commandlineparse(int argc, char **argv, int *L,int *K, int *seed, int *name, int *repeat_num, int *iter_num){
  	int i;
  	*seed = getpid()*2 + getpid();
  	for (i = 1; i < argc; i++){  //Start at i = 1 to skip the command name.
    		if (argv[i][0] == '-'){
      			switch (argv[i][1]){
      				case 'L':       *L = atoi(argv[++i]);
        			break;
      				case 's':       *seed = atoi(argv[++i]);
        			break;
      				case 'r':       *repeat_num = atoi(argv[++i]);
        			break;
      				case 'k':       *K = atoi(argv[++i]);
        			break;
      				case 'n':       *name=atoi(argv[++i]);
        			break;
              case 'i':       *iter_num = atoi(argv[++i]);
              break;
      				default:
        			fprintf(stderr,"\nError:  Incorrect option %s\n",argv[i]);
        			fprintf(stderr,"\nAvailable options: \n\
  				-L = L (HN3 length; it has to be 2^n, n=2,3,4,...)\n\
  				-k = highest level (HN3 length is 2^k))\n\
  				-r = repeat number (Default = 100; or any integer in (0,99))\n\
  				-s = seed: random seed  (default=pid)\n\
    				\n");
        			exit(0);
      			}//switch
    		}//if
  	}//for
}



int main(int argc, char *argv[]){
  int i, j, seed, name, flag = 0, print_flag = 0;
  int iter_num = 3, iter_now = 0, repeat_num = 100, repeat_iter, temp_index = 0; // iter_num: total number of schedules; iter_now: the current iteration; repeat_num: total repeatitions to run
  char energy_file[20] = "H3AF_E_";	// Name of the output file of density of states
  char mag_file[20] = "H3AF_m_";
  char energy_low_file[20] = "H3AF_E_low";
  char mag_low_file[20] = "H3AF_m_low";
  char buf[6];  // to save file names with numbers 
  double temps[T_L];
  double temp_now;
  
  Commandlineparse(argc, argv, &L, &K, &seed, &name, &repeat_num, &iter_num);
  if(L==0){  L=Intpow(2,K);  }
  else{	 K=Intlog2(L);  }
  
  rnd =  gsl_rng_alloc(gsl_rng_mt19937);
 	gsl_rng_set(rnd, (unsigned long int)seed);
 	for(i = 0; i < 1259; i++) {
      gsl_rng_uniform(rnd); // initial cycling on the random generator
  }
  printf("\n    ************ HN3 AF Ising Simulated Annealing Running ************ \n");
  printf("\n      argc: %d\n      ",argc);
  for(j = 0; j < argc; j++){
      printf("%s ", argv[j]);
  }
  FILE *fp1;
  FILE *fp2;
  printf("     L = %d\n", L);
  
  short sites[L];
  int neighbors[3*L]; 	// store the neighbors index of all sites
  double dt_array[iter_num], dt_now, dt_bound;
  double energy[T_L * iter_num];
  double mag[T_L * iter_num];
  double energy_low[iter_num * repeat_num];
  double mag_low[iter_num * repeat_num];
  double e_t = 0.0, m_t = 0.0; // to temperary store the energy and mag
	for(i = 0; i < iter_num*T_L; i++){
			energy[i] = 0.0;
      mag[i] = 0.0;
	}
  for(i = 0; i < iter_num; i++){
    dt_array[i] = 0.001 / pow(2.0, (double)i);
  }
	sprintf(buf,"%d", name);
 	strcat(energy_file, buf);
 	strcat(mag_file, buf);
  strcat(energy_low_file, buf);
  strcat(mag_low_file, buf);

  for (repeat_iter = 0; repeat_iter < repeat_num; repeat_iter++){
    for (iter_now = 0; iter_now < iter_num; iter_now++){
      dt_now = dt_array[iter_now];
      dt_bound = 0.8 * dt_now;
      e_t = 0.0;
      m_t = 0.0;
      
      Initialize(sites, neighbors, temps);
      temp_now = temps[0] + 2.0*T_STEP;
      temp_index = 0;
      flag = 0;
      print_flag = 0;
      if(gsl_rng_uniform(rnd)<(20.0/(double)(repeat_num*iter_num))){
        print_flag = 1;
        printf("\n\n      *************  Repetition %d Iteration %d  ************ \n", repeat_iter, iter_now);
        printf("      Starting T=%-.2f; schedule = %-.2e; \n      Temp_now = %-.2e; Ending T=%-.4f;\n", temps[0], dt_now, temp_now, temps[T_L-1]);
        printf("      Starting State = %d;  StartingMag = %d; \n ", CurrentEnergy, CurrentMag);
      }
      while (Total_Steps < Initial_Steps){
          Update(sites, neighbors, temp_now);
          //if(CurrentMag  != 0)    sites_print(sites);
      }
      if(print_flag)  printf("     After Init, CurrentEnergy = %d;  CurrentMag = %d;\n      MinStateFouend = %d; MinMagFound = %d\n", CurrentEnergy, CurrentMag, MinEnergyFound, MinMagFound);
      //printf("      check CS     = %d; ", Check_CurrentEnergy(sites, neighbors));
      //printf("      check Mag    = %d; \n", Check_CurrentMag(sites));
      //sites_print(sites);
      while (temp_index < T_L){
        if(Update(sites, neighbors, temp_now)){
          temp_now = temp_now - dt_now;
          //printf("temp_now = %.2f", temp_now);
        }
        if(fabs(temp_now - temps[temp_index]) < dt_bound){
          e_t = e_t + (double)CurrentEnergy;
          m_t = m_t + (double)CurrentMag;
          flag = 1;
          //if(count<50){printf("e_t= %.0f\n", e_t);count++;sites_print(sites);}
        }
        else if(flag == 1){
            e_t = e_t / L;
            m_t = m_t / L;
            
            energy[iter_now * T_L + temp_index] += e_t;
            mag[iter_now * T_L + temp_index] += m_t;

            if(temp_index == T_L - 1){
                energy_low[iter_now * repeat_num + repeat_iter] = e_t;
                mag_low[iter_now * repeat_num + repeat_iter] = m_t;
            }
            e_t = 0.0;
            m_t = 0.0;
            flag = 0;
            temp_index++;
            //printf("temp_index = %d \n", temp_index);
        }
      }
      if(print_flag) printf("      After all, CurrentEnergy = %d;  CurrentMag = %d;\n      MinStateFouend = %d; MinMagFound = %d\n", CurrentEnergy, CurrentMag, MinEnergyFound, MinMagFound);

    }
  }

            
  for(i = 0; i < T_L * iter_num; i++){
    energy[i] = energy[i] / repeat_num /L; 
    mag[i] = mag[i] / repeat_num / L;
	}

  fp1 = fopen(energy_file, "w");
  fp2 = fopen(mag_file, "w");
  for(i=0; i < T_L; i++){
    fprintf(fp1, "%-8.4f ",temps[i]);
    fprintf(fp2, "%-8.4f ",temps[i]);
    for(j=0; j < iter_num; j++ ){
      fprintf(fp1, "%-12.8f ", energy[j*T_L + i]);
      fprintf(fp2, "%-12.8f ", mag[j*T_L + i]);
    }
    fprintf(fp1, "\n");
    fprintf(fp2, "\n");
  }
  fclose(fp1);
  fclose(fp2);
  
  fp1 = fopen(energy_low_file, "w");
  fp2 = fopen(mag_low_file, "w");
  for(i = 0; i < repeat_num; i++){
    for(j = 0; j < iter_num; j++){
      fprintf(fp1, "%-10.2f ", energy_low[j * repeat_num + i]);
      fprintf(fp2, "%-8.2f ", mag_low[j * repeat_num + i]);
    }
    fprintf(fp1, "\n");
    fprintf(fp2, "\n");
  }
  fclose(fp1);
  fclose(fp2);
  return(0);
}
