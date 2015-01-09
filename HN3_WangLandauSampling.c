/*
  gcc Program.c -lgsl -lgslcblas -lm -o program.x
  
  Key words: simulated annealing, constraint optimization, complex networks, regular graph, disordered system
  
  This code is to use Wang-Landau sampling to
  find the density of states of the lattice gas problem 
	with at-most-1-neighbor density constraint in HN3;
	Periodic Boundary Condition is used in this code; 
	The output is a file with chemical potential and corresponding average states.
  
  There are 2 other versions of this method for HN5 and HNNP, 
  in which the data structure is a little different.
  
  For more information of the model and method:
  http://arxiv.org/abs/1409.8313
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <term.h>
#include <ncurses.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#define MIN_INV_MU 0.05
#define INV_MU_STEP 0.001
#define INV_MU_L 800

int L=0, K=0, Flat=80, Seed, Name = 0; 	// Input parameters
float Flatness=0.8;
int StatesLength = 0, MaxLength, MaxState, CurrentState, CurrentIndex, MinStateFound, MinState;
double Logfnow=1.0, Logfmin=0.00000001, MinStateLogg, MaxSteps=9.0e12;
double Initial_Steps=0.0, Check_Every_Steps=0.0, Total_Steps=0.0, Inner_Loop_Steps=0.0; 

gsl_rng *rnd;	// RNG: random number

struct StateInfo{
  int E;
  float H; 
  double logg;
}; 


/* 	integer power of an integer (with 0^0 = 1) 	*/
int Intpow(int base, int power){
  int x=1;
  int i;
  for(i = 0; i<power; i++) {x=x*base;}
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



int Print_States(struct StateInfo state[]){
    int i;
    printf(" Print_States:\n");
    printf(" Current State: %d;\n", CurrentState);
    printf(" StatesLength: %d; \n", StatesLength);
    printf(" E          logg           H\n");
    for(i = 0; i < StatesLength; i++){
        printf(" %-12d %-18.10e %-12.2e\n", state[i].E, state[i].logg, state[i].H);
    }
    return(0);
}


/*	Initialize the StateInfo, neighbors and other parameters*/
void Initialize(struct StateInfo state[], int site[], int neighbors[]){
  	int i;
  	int *neigh;
  	for(i=0; i<L; i++){
      site[i] = 1;
      neigh = FindNeighbors(i+1);
      neighbors[3*i]=neigh[0]-1;
	    neighbors[3*i+1]=neigh[1]-1;
	    neighbors[3*i+2]=neigh[2]-1;
    }
    state[0].E = 3*L/2;
    state[0].logg = 0.0;
    state[0].H = 0.0;
  	for(i = 1; i < MaxLength; i++){
        state[i].E = -1;
        state[i].logg = 0.0;
        state[i].H = 0.0;
  	}
    StatesLength = 1;
  	CurrentState = 3*L/2;
    CurrentIndex = 0;
  	MaxState = CurrentState;
    MinState = CurrentState;
  	MinStateFound = 0;
}

int Exchange_Energy(int site[], int neighbors[], int i, int h){
  int j, e1 = 0, e2 = 0, e, neig_flag = 0;
  if(site[i] == site[h]){
    fprintf(stderr, "\n Error: In Exchange_Energy, site[%d] = %d, site[%d] = %d;\
     \n No need to exchange!", i, site[i], h, site[h]);
  }
  for(j = 0; j < 3; j++){
      if(neighbors[3*i+j] == h){
        neig_flag = 1;
        break;
      }
  }
  if(neig_flag){
    for(j = 0; j < 3; j++){
        e1 += site[neighbors[3*i + j]];
        e2 += site[neighbors[3*h + j]];
    }
    e = CurrentState + (site[h] - site[i])*(e1 - site[h] - e2 + site[i]);
    return(e);
  }
  else{
    for(j = 0; j < 3; j++){
      e1 += site[neighbors[3*i + j]];
      e2 += site[neighbors[3*h + j]];
    }
    e = CurrentState - 2*e1*site[i] - 2*e2*site[h];
  }
  return(e);
}

int Binary_Search(struct StateInfo state[], int e, int imin, int imax){
  int imid;
  while(imin <= imax){
    imid = (imin + imax) / 2;
    if(state[imid].E == e) 
      return(imid);
    else if(state[imid].E < e) 
      imax = imid - 1;
    else 
      imin = imid + 1;
  }
  return(-1); // if nothing is found in the loop
}

int Find_E_Index(struct StateInfo state[], int e){
  int index = 0;
  if(Total_Steps < Initial_Steps){
    for(index = 0; index < StatesLength; index++){
      if(state[index].E == e){
          return(index);
          break;
      }
    }
    if(index == StatesLength){//could remove this if, but keep the statements in it
      StatesLength++;
      index = StatesLength - 1;
      if(state[index].E != -1){ // could delete
        fprintf(stderr, "\n Error: The new state's state[%d].E != -1 \n\
        (all state.E is initialized to be -1)!!!\n StatesLength = %d\n \
             Current state[%d].E = %d", index, StatesLength-1, index, state[index].E);
        Print_States(state);
        exit(0);
      } 
      state[index].E = e; // logg, H is already initialized in function Initialize
      return(index);
    }
    else//could delete
      fprintf(stderr, "\n Error: In intial steps, StatesLength = %d, Index = %d", StatesLength, index); 
  }
  else{ // find sorted state in log(N) steps
    index = Binary_Search(state, e, 0, StatesLength - 1);
    if(index != -1)
      return(index);
    else
      return(-1);
  }
  return(index);
}

/*	Update the system by adding or removing particles*/
int Update(struct StateInfo state[], int site[], int neighbors[]){
	int i, j, e2 = 0, ie1, ie2;
  Total_Steps += 1.0;
	i = (int)(L*gsl_rng_uniform(rnd));
  for(j = 0; j < 3; j++){
    e2 += site[neighbors[3*i + j]];
  }
  e2 = CurrentState - 2 * e2 * site[i];
  ie1 = CurrentIndex;
  ie2 = Find_E_Index(state, e2);
  if(ie2 == -1){
      site[i] = -site[i];
      CurrentState = e2;
      CurrentIndex = StatesLength;
      if(CurrentState < MinState)
        MinState = CurrentState;	
      return(e2);
  }
	if(state[ie1].logg > state[ie2].logg || exp(state[ie1].logg - state[ie2].logg) > gsl_rng_uniform(rnd)){
		site[i] = -site[i];
		CurrentState = e2;
    CurrentIndex = ie2;
		state[ie2].logg += Logfnow;
		state[ie2].H += 1.0;
    if(CurrentState < MinState)
      MinState = CurrentState;		
  }
  else{
    state[ie1].logg += Logfnow;
		state[ie1].H += 1.0;
	}
  /*
	i = (int)(L * gsl_rng_uniform(rnd));
  h = (int)(L * gsl_rng_uniform(rnd));
  while(h == i){
    h = (int)(L * gsl_rng_uniform(rnd));
  }
  if(site[i] != site[h]){
    e2 = Exchange_Energy(site, neighbors, i, h);
    ie2 = Find_E_Index(state, e2);
    if(ie2 == -1){
      site[i] = -site[i];
      site[h] = -site[h];
      CurrentState = e2;
      CurrentIndex = StatesLength;
      if(CurrentState < MinState)
        MinState = CurrentState;	
      return(e2);
    }
  }

  if(site[i] == site[h]){
    state[CurrentIndex].logg += Logfnow;
		state[CurrentIndex].H += 1.0;	
  }
	else if(state[CurrentIndex].logg > state[ie2].logg || exp(state[CurrentIndex].logg - state[ie2].logg) > gsl_rng_uniform(rnd)){
		site[i] = -site[i];
    site[h] = -site[h];
		CurrentState = e2;
    CurrentIndex = ie2;
    if(CurrentState < MinState)
      MinState = CurrentState;	
		state[ie2].logg += Logfnow;
		state[ie2].H += 1.0;		
  }
  */
	return(-1);
}


/*	Check the flatness of the histogram	*/
int CheckFlat(struct StateInfo state[]){
	int i;
	double meanH, totalH=0.0, minH, maxH;
	minH = state[0].H;
	maxH = state[0].H;
 	totalH = maxH;
	for(i = 1; i<StatesLength; i++){
		totalH += state[i].H;
		if(state[i].H < minH)
      minH = state[i].H;
		else if(state[i].H > maxH) 
      maxH = state[i].H;
	}
	meanH = totalH / (float)StatesLength;
	if(minH > meanH*Flatness && meanH > maxH*Flatness){
    printf(" logf=%-.1e; Sweeps=%.1e; [minH,meanH,maxH]=[%-.1e, %-.1e, %-.1e]\n",  Logfnow/2.0, Total_Steps/(double)L, minH, meanH, maxH); 
		return(1);
	}
	if(Total_Steps > MaxSteps){
		Logfmin=1.0;
    printf("!!! MaxSteps %.2e reached! Logfmin is set to %.2f", MaxSteps, Logfmin);
	}
	return(0);
}

void Reset_H_0(struct StateInfo state[]){
    int i;
    for(i = 0; i < StatesLength; i++){
      state[i].H = 0.0;
    }
    return;
}

int Swap_States(struct StateInfo state[], int i, int j){
    struct StateInfo temp;
    temp.E = state[i].E;
    temp.logg = state[i].logg;
    temp.H = state[i].H;
    state[i].E = state[j].E;
    state[i].logg = state[j].logg;
    state[i].H = state[j].H;
    state[j].E = temp.E;
    state[j].logg = temp.logg;
    state[j].H = temp.H;
    return(0);
}
/* This function is just for quick sort*/
int Sort_Partition(struct StateInfo state[], int top, int bottom){
  int x = state[top].E;
  int i = top - 1;
  int j = bottom + 1;
  do{
      do{
          j--;
        }while(x > state[j].E );
      do{
          i++;
        }while(x < state[i].E);
      if(i < j) 
        Swap_States(state, i, j);
  }while(i < j);
  return(j);
}

/* Quick Sort the States by energy*/
void Quick_Sort_States(struct StateInfo state[], int top, int bottom){
    int middle;
    if(top < bottom){
        middle = Sort_Partition(state, top, bottom);
        Quick_Sort_States(state, top, middle);
        Quick_Sort_States(state, middle+1, bottom);
    }
    return;
}

/* Check the current energy state*/
int Check_CurrentEnergy(int sites[], int neighbors[]){
  int i, energy = 0;
  for(i = 0; i < L; i++){
      energy = energy + sites[i] * (sites[neighbors[3*i]] + sites[neighbors[3*i + 1]] + sites[neighbors[3*i + 2]]);
  }
  energy = energy / 2;
  return(energy);
}

/* Check the current mag state*/
int Check_CurrentMag(int sites[]){
  int i, mag = 0;
  for(i = 0; i < L; i++){
    mag = mag + sites[i];
  }
  return(mag);
}

/*	Function to read input parameters	*/
void Commandlineparse(int argc, char **argv, int *L, int *Flat,int *K, int *Seed, int *Name){
  	int i;
  	*Seed = getpid();
  	for (i = 1; i < argc; i++){  //Start at i = 1 to skip the command name.
    		if (argv[i][0] == '-'){
      			switch (argv[i][1]){
      				case 'L':       *L = atoi(argv[++i]);
        			break;
      				case 'f':       *Flat = atoi(argv[++i]);
        			break;
      				case 's':       *Seed = atoi(argv[++i]);
        			break;
      				case 'k':       *K = atoi(argv[++i]);
        			break;
      				case 'n':       *Name=atoi(argv[++i]);
        			break;
      				default:
        			fprintf(stderr,"\nError:  Incorrect option %s\n",argv[i]);
        			fprintf(stderr,"\nAvailable options: \n\
  				-L = L (HN3 length; it has to be 2^n, n=2,3,4,...)\n\
  				-k = highest level (HN3 length is 2^k))\n\
  				-f = flatness (Default = 80; or any integer in (0,99))\n\
  				-s = seed: random seed  (default=pid)\n\
    				\n");
        			exit(0);
      			}//switch
    		}//if
  	}//for
}



int main(int argc, char *argv[]){
  
  	int i, j, del1=0;	// General varialbe for loops
  	int new_state_energy;
  	char dos_file[20]="H3AFDOSe";	// Name of the output file of density of states
  	char buf[8];
    struct StateInfo *state;
  	Commandlineparse(argc, argv, &L, &Flat, &K, &Seed, &Name);

  	rnd =  gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rnd, (unsigned long int)Seed);
    for(j = 0; j <100; j++) {
      		gsl_rng_uniform(rnd); // initial cycling on the random generator
  	}

  	printf("\n    ************ HN3 AF Ising Wang Landau Running ************ \n");
  	printf("\n      argc: %d\n      ",argc);
  	for(j = 0; j<argc; j++){
    		printf("%s ",argv[j]);
  	}
  
  	/* Setup the simulation according to the inputs*/
  	if(L==0){
    		L=Intpow(2,K);
  	}
  	else{
    		K=Intlog2(L);
  	}
  	Flatness=(float)Flat/100.0;
  	printf("\n\n L=%d; Flatness=%f;\n\n", L, Flatness);
  	MaxLength = (3*L); // this is just for HN3
  	int site[L];
  	int neighbors[3*L]; 	// store the neighbors index of all sites
  	Initial_Steps= (double)((double)K*50*(double)L*(double)L);
  	Check_Every_Steps=(double)((double)(K*K)*10.0*(double)L);
    state=(struct StateInfo *)malloc(MaxLength*sizeof(struct StateInfo));

    Initialize(state, site, neighbors);
  	Total_Steps=0.0;

  	sprintf(buf,"%d", K);
    strcat(dos_file, buf);
    strcat(dos_file, "_r");
    sprintf(buf, "%d", Name);
    strcat(dos_file, buf);
  	FILE *fp;

  	while(Total_Steps < Initial_Steps){
      Update(state, site, neighbors);
      if(gsl_rng_uniform(rnd) < 0.0001 && del1 < 3){
        if(Check_CurrentEnergy(site, neighbors) == CurrentState)
          printf(" CurrentEnergy Correct = %d \n", CurrentState);
        else
          printf(" Error: CurrentEnergy is not right! \n CurrentState = %d,\
           actually %d \n", CurrentState, Check_CurrentEnergy(site, neighbors));
        del1++;
      }
    }
    //printf(" Init Done. StatesLength = %d\n", StatesLength);
    //Print_States(state);
    /*state = realloc(state, StatesLength*sizeof(struct StateInfo));
    if(state == NULL){
      fprintf(stderr, "\n Error: Memory reallocating of state(E,M,logg, H)!\n"); 
      return(2);
  	}
    //printf(" Realloc Done\n");*/
    printf(" After Initial Steps %.2e, StatesLength = %d, MinEnergy = %d\n", Initial_Steps, StatesLength, MinState);
    Quick_Sort_States(state, 0, StatesLength-1);
    printf(" Quick_Sort Done \n");
    //Print_States(state);
    Inner_Loop_Steps=0.0;
    del1 = 0;
  	while(Logfnow > Logfmin){
      Inner_Loop_Steps += 1.0;
      new_state_energy = Update(state, site, neighbors);
      if(gsl_rng_uniform(rnd) < 0.0001 && del1 < 3){//could delete
        if(Check_CurrentEnergy(site, neighbors) == CurrentState)
          printf(" After sort, CurrentEnergy Correct = %d \n", CurrentState);
        else
          printf(" Error: CurrentEnergy is not right! \n CurrentState = %d,\
           actually %d \n", CurrentState, Check_CurrentEnergy(site, neighbors));
        del1++;
      }
      if(new_state_energy != -1){	// if new state found, we add another states
        StatesLength++;
        if(StatesLength > MaxLength) {
          Print_States(state);
          fprintf(stderr, "\n Error: the states length %d is too long;\n Check the States above!!!\n", StatesLength);
          exit(0);
        }
        printf(" New State: %d: coming; StatesLength %d\n", new_state_energy, StatesLength);
        /* These are to resize the matrix size; which is not necessary. just waste at most several Mb memory.
        state = (struct StateInfo *)realloc(state, StatesLength*sizeof(struct StateInfo));
        if(state == NULL){
            fprintf(stderr, "\n Error: Memory reallocating of state(E, logg, H)!\n"); 
            return(2);
        } */
        for(i = StatesLength-2; i > 0; i--){
            if(state[i].E > new_state_energy){
                state[i+1].E = new_state_energy;
                state[i+1].logg = state[i].logg * 0.8;
                state[i+1].H = 0.0;
                Logfnow = Logfnow * 8;
                if(Logfnow > 1) Logfnow = 1;
                Reset_H_0(state);
                break;
            }
            else{
                state[i+1].E = state[i].E;
                state[i+1].logg = state[i].logg;
                state[i+1].H = state[i].H;
            }
        }
        printf(" New State: %d added to list; Min State: %d; CurrentState = %d\n", new_state_energy, MinState, CurrentState);
        //Print_States(state);
        printf(" \n\n");
        Inner_Loop_Steps=0.0;

      }
      if(Inner_Loop_Steps > Check_Every_Steps){
        Inner_Loop_Steps = 0.0;
        if(CheckFlat(state)){
          Reset_H_0(state);
          Logfnow=Logfnow/2.0;
        }
      }
  	}
 
    printf("\n MinState: %d; \n MaxState: %d; \n States Length: %d; \n Total MC steps: %-.4e;\
     \n Total random steps: %-.4e; \n ", 
     MinState, MaxState, StatesLength, Total_Steps/(double)(L), Total_Steps);
  	fp=fopen(dos_file,"w");
  	MinStateLogg=state[0].logg;
  	for(i = 0; i < StatesLength; i++){
      state[i].logg = state[i].logg-MinStateLogg;
      fprintf(fp, "%-8d %-18.8f\n", state[i].E, state[i].logg);
  	}
  	fclose(fp);

	/*
	for(i=0; i<INV_MU_L; i++){
		z[i]=0.0;
		rho[i]=0.0;
		s[i]=0.0;
		shannon[i]=0.0;
	}
	*/
	/*	Calculating Z, packing fraction, entropy and so on */
	/*
	for(i=0; i<INV_MU_L; i++){
		for(j=0; j<StatesLength; j++){
			z_temp = expl((long double)state[j].logg)*expl((long double)(1/InvMu[i]*(long double)j));
			z[i] += z_temp;
			rho[i] += (long double)(j*z_temp);
		}
		for(j=0; j<StatesLength; j++){
			z_temp= expl((long double)state[j].logg)*expl((long double)(1/InvMu[i]*(double)j));
			pj=(long double)(z_temp/z[i]);
			if(pj>1e-12){
				shannon[i] -= pj*log2l(pj);
			}
		}
		rho[i]=(long double)(rho[i]/z[i]/(double)L);
		s[i]=1.0/(long double)L*logl(z[i])-1/InvMu[i]*rho[i];
	}

  	fp=fopen(pes_file,"w");
	for(i=0; i<INV_MU_L; i++){
		fprintf(fp, "%-8.4f  %-18.9LG %-18.9LG  %-18.9LG \n", InvMu[i], rho[i], s[i], shannon[i]);
	}
	fclose(fp);
	*/
 	return(0);
}
