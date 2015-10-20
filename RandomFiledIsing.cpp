// Original Compile: g++  -Wall  -lm -lgsl -lgslcblas -c  "%f"
// Original Build: g++ RFIM1.cc -Wall   -lm -lgsl -lgslcblas -o RFIM1.x 
// Original run: ./RFIM1.x-L 32  -g 1.0  -T 1.0   -r 2 -s 615 -n seg -M 1000
// Compile: g++  -Wall  -lm -lgsl -lgslcblas -c  "%f"
// Build: g++ "%f" -Wall   -lm -lgsl -lgslcblas -o "%e" 
// Run: "./%e"  -k 17 -0 100  -1 100 -y 0

/* Code: 2D Random Field Ising Model without random field
 * Input: (default values)
 *        system size "-L": L = 16 (length)
 *        temperature "-T" T = 2 
 *        variance of hi "-g": s2hi = 1
 *        configuration output MC sweeps interval "-M": dMCs = 100.2
 *        total MC sweeps: total_MCs = 1000 
 *        repeat_num "-r": number of repetitions repeated_num = 1 
 *          (each repetition will re-run the simulation and export new outputs)
 *        NAME "-n": name = "" ; any string with length < 8 to help you name the output csv files;
 *        seed "-s": seed=getpid() (usually 1 or 2); 
 * ### theory of fractal 
 * Output:
 *        configurations: config_L***_T***_r**_NAME.csv
 *                      1st column: MC step;
 *                      Each row except 1st column: Spin configurations
 *        Magnetizations: mag_L***_T***_r**_NAME.csv
 *                      1st column: MC step;
 *                      2nd column: average magnetization (M/N);
 *        Energy: e__L***_T***_r**_NAME.csv
 *                      1st column: MC step;
 *                      2nd column: average energy (E/N);
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <term.h>
#include <ncurses.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>
#include <unistd.h>
#include <cstdio>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>

using namespace std;

#define L_DEF 1000
#define L_MAX 10000
#define T_DEF 2.0
#define S_DEF 1.0
#define NAME_LEN 8
#define DIM 2
#define JCON -1
#define RECORD_LEN 5000
#define FILE_NAME_LEN 128

gsl_rng *rnd = gsl_rng_alloc(gsl_rng_mt19937); // global random number setup

class RFIM2d {
  private:
    int L, N;     // system size N=L^2
    double s_hi;  // sigma of hi: standard deviation of hi
    int *spins;   // length is N; mememory is allocate later
    double *his;  // static internal random field 
    double T;     // temperature
    double E;     // energy
    int Espin;    // energy without random field
    int M;        // magnetization
    
  public:
    int check_E(); // check whether the energy is right
    int check_M(); // check whether the mag is right
    int get_L() { return(L);} // get L
    int get_N() { return(N);} // get N
    int get_M() { return(M);} // get M
    int *get_spins() {return(spins);} // get 
    double get_s_hi() {return(s_hi);}
    double get_s_Hi() {return(sqrt(T));}
    double get_E() { return(E); }
    int get_Espin() { return(Espin); }
    double get_T() { return(T); }
    int initialize();  // initialization with default values
    int *neighbors(int iv); //
    int neighbor_check(int iv);
    void set_L (int Lp = L_DEF) {L = Lp; N=L*L;} // set system size
    void set_s_hi (double sp = S_DEF) {s_hi = sp;}
    void set_T (double Tp = T_DEF) {T = Tp;}
    int set_spins ();
    int set_his ();
    double calc_E();
    int calc_M();
    int calc_Espin();
    int print_info(int config);
    int update();
};

double RFIM2d::calc_E()
{
    double energy_s = 0.0;
    double energy_h = 0.0;
    double energy = 0.0;
    int *ng; 
    for (int i = 0; i < N; i++)
    {
        ng = neighbors(i);
        for (int j = 0; j < 2*DIM; j++)
        {
            energy_s = energy_s + JCON * spins[i] * spins[ng[j]];   
        }
        energy_h = energy_h - his[i] * spins[i];
    }
    energy = energy_s / 2.0 + energy_h;
    
    return(energy);
}

int RFIM2d::calc_Espin()
{
    double energy_s = 0;
    int *ng; 
    for (int i = 0; i < N; i++)
    {
        ng = neighbors(i);
        for (int j = 0; j < 2*DIM; j++)
        {
            energy_s = energy_s + JCON * spins[i] * spins[ng[j]];   
        }
    }
    energy_s = energy_s / 2.0;
    
    return(energy_s);
}


int RFIM2d::calc_M()
{
    int mag = 0; 
    for (int i = 0; i < N; i++)
    {
        mag = mag + spins[i];
    }
    return(mag);
}


int RFIM2d::check_E()
{
    double energy = calc_E();
    int espin = calc_Espin();
    if((energy - E) < 1e-12 && (espin-Espin) < 1e-12)
    {
        cout << "  Energy is correct: E=" << energy << "; Espin=" << Espin << endl;
    }
    else
    {
        cout << "  Energy is wrong: E = " << E << endl;
        cout << "  But --------calc_E = " << energy << endl;
        cout << "  Energy is wrong: Espin = " << Espin << endl;
        cout << "  But --------calc_Espin = " << espin << endl;
        E = energy;
    }
    return(0);
}

int RFIM2d::check_M()
{
    int mag = calc_M();
    if(M == mag)
    {
        cout << "  Magnetization is correct: M = calc_M = " << M << endl;
    }
    else
    {
        cout << "  Magnetization is wrong: M = " << M << endl;
        cout << "  But -------------- calc_M = " << mag << endl;
        M = mag;
    }
    return(0);
  
}

int * RFIM2d::neighbors(int iv){
    int x, ld=1;  
    int direction=-1, n=0;
    int d=0;
    static int neighs[DIM*2];
    
    for(d=0; d<DIM; d++){
      for(direction=-1; direction<2; direction += 2){
        ld=1;
        if(d>0) ld *= L;
        x = (iv/ld)%L;
        neighs[n++]=(iv+((x+L+direction)%L - x)*ld);
      }
    }
    return(neighs);
}



int RFIM2d::neighbor_check(int iv)
{
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        L must be in the range of [2, %d]!\n", L, N, L_MAX);
        return(0);
    }
    int *ng = neighbors(iv);
    cout << "   Neighbors of site " << iv << endl;
    for (int i = 0; i < 2*DIM; i++)
    {
        cout << "   " << ng[i] << endl;
    }
    return(1);
}


// initialization with default values
int RFIM2d::initialize() 
{
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        L must be in the range of [2, %d]!\n", L, N, L_MAX);
        return(0);
    }
    set_spins();
    set_his();
    
    E = calc_E();
    Espin = calc_Espin();
    M = calc_M();
    return(1);
}

int RFIM2d::set_spins () 
{  
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        L must be in the range of [2, %d]!\n", L, N, L_MAX);
        return(0);
    }    
    spins = new int[N];
    for(int i=0; i<N; i++)
    {
        if(i < N/2) { spins[i] = -1; }
        else { spins[i] = 1; }
    }
    return(1);
}

int RFIM2d::set_his()
{
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        L must be in the range of [2, %d]!\n", L, N, L_MAX);
        return(0);
    }
    his = new double[N];
    for(int i=0; i<N; i++)
    {
        his[i] = gsl_ran_gaussian(rnd, s_hi);
        //his[i] = 0.0; //del
    }
    return(1);
}

int RFIM2d::print_info(int config=0)
{
    int i; 
    if(L < 2 || L > L_MAX || N != (L*L) )
    {
        fprintf(stderr,"\Error: L=%d, N=%d!\n  \
        System not set up yet!\n", L, N);
        return(0);
    }
    printf("*** Simulation System Information ***\n");
    printf("   L=%d; N=%d;\n   T=%.2f; hi_sigma=%.2f Hi_s=%.2f\n   E=%.2f; M=%d;\n"\
    , L, N, T, s_hi, get_s_Hi(), E, M);
    if(config != 0)
    {
        printf("   ** Configurations **\n");
        for(i=0; i<N; i++)
        {
            printf("    site %-4d: %-4d, %-.4f\n", i, spins[i], his[i]);
        }
    }
    return(1);
}

int RFIM2d::update()
{
    int site;
    site = N*gsl_rng_uniform(rnd);
    //printf("site:%d\n", site);//del
    double dE = 0.0;
    int dEspin = 0;
    int *ng = neighbors(site);
    for (int i=0; i < 2*DIM; i++)
    {
        dE = dE - 2.0 * JCON * spins[site] * spins[ng[i]];
    }
    dEspin = (int)dE;
    dE = dE + 2.0 * his[site] * spins[site];  //??????
    if(dE < 0.0 || gsl_rng_uniform(rnd) < exp(-dE/T))
    {
        spins[site] = -spins[site];
        E = E + dE;
        Espin += dEspin;
        M = M + 2*spins[site];
    }
    return(1);
}

int set_record_steps(double *rand_steps_record, double total_MCs, int N)
{
    double interval;
    int i=0;
    
    if( total_MCs * N  < RECORD_LEN*2)
    {
        fprintf(stderr,"\Error: Simulation is too short for the output!\n  \
        Record_Length=%d; total_MCs=%.2f!\n", RECORD_LEN, total_MCs*N);
        return(0);
    }
    
    interval = total_MCs*1.0*N/RECORD_LEN;
    if(total_MCs <= RECORD_LEN)
    {
        for (i = 0; i < RECORD_LEN; i++)
        {
            rand_steps_record[i] = floor(interval * i);
        }
    }
    else if (total_MCs <= RECORD_LEN*2)
    {
        for (i = 0; i < RECORD_LEN/2; i++)
        {
            rand_steps_record[i] = floor(interval/2.0 * i);
        }
        for (i = RECORD_LEN/2; i < RECORD_LEN; i++)
        {
            rand_steps_record[i] = rand_steps_record[i-1] + interval*1.5;
        }
    }
    else
    {
        for (i = 0; i < RECORD_LEN/2; i++)
        {
            rand_steps_record[i] = floor(N * i);
        }
        interval = (total_MCs - RECORD_LEN/2) * N * 2.0 / RECORD_LEN;
        for (i = RECORD_LEN/2; i < RECORD_LEN; i++)
        {
            rand_steps_record[i] = floor(rand_steps_record[i-1] + interval);
        }
    }
    return(1);
}

int set_file_name(char *newname, const char *name0, char *name1, RFIM2d *age, int rep)
{
    int i = 0;
    //int *nm;
    char buf[32] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
    for(i = 0; i < FILE_NAME_LEN; i++)
        newname[i]='\0';
    strcat(newname, "./data/");
    strcat(newname, name0);
    sprintf(buf,"_L%d_T%.1f_shi%.1f_sHi%.2f_r%d_", age->get_L(), age->get_T(), age->get_s_hi(), age->get_s_Hi(), rep);
    strcat(newname, buf);
    strcat(newname, name1);
    strcat(newname, ".csv");
    
    DIR *dirp;
    if (!(dirp = opendir("./data"))) { /* open the directory */
        mkdir("./data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    closedir(dirp);
    //cout << "newname: " <<newname << endl; 
    return(1);
}



/*	Function to read input parameters	*/
void Commandlineparse(int argc, char **argv, RFIM2d *age, double *total_MCs, int *repeat_num, int *seed, char *name)
{
  	int i;
  	*seed = getpid();
    age->set_s_hi(); // set to default values first
    age->set_T(); // set to default values first
    age->set_L(); // set to default values first
  	for (i = 1; i < argc; i++)
    {   //Start at i = 1 to skip the command name.
    		if (argv[i][0] == '-')
        {
      			switch (argv[i][1])
            {
      				case 'L':       age->set_L( atoi(argv[++i]) );
        			break;
      				case 's':       *seed = atoi(argv[++i]);
        			break;
                    case 'g':       age->set_s_hi( atof(argv[++i]) );
                    break;
      				case 'r':       *repeat_num = atoi(argv[++i]);
        			break;
      				case 'T':       age->set_T( atof(argv[++i]) );
        			break;
      				case 'n':       strncpy(name, argv[++i], 7);
        			break;
                    case 'M':       *total_MCs = atof(argv[++i]);
                    break;
                    case 'd':       age->set_L();
                              age->set_T();
                              age->set_s_hi();
                    break;
      				default:
        			fprintf(stderr,"\nError:  Incorrect option %s\n",argv[i]);
        			fprintf(stderr,"\n\
              Available options: \n\
              -L = L (Length; N=L^2; default=16)\n\
              -g = variance for hi (default=16)\n\
              -n = name (default='')\
              -r = repeat number (default= 1)\n\
              -s = seed: random seed  (default=getpid())\n\
              -t = total MC sweeps (default=1000)\
              \n");
        			exit(0);
      			}//switch
    		}//if
  	}//for
}



int main (int argc, char *argv[]) 
{
  /************* System Setup ****************/
    int repeat_num = 1, iter = 0, seed, *spins, i; // iter=0:repeat_num-1
    RFIM2d age;  // simulation main object
    char name[NAME_LEN + 1];
    for(int i=0; i <= NAME_LEN; i++) {name[i] = '\0';}
    double total_MCs = 100.0, total_rnd_steps = 0.0, total_rnd_steps_ran=0.0;
    double rand_steps_record[RECORD_LEN];
    double mag_steps[RECORD_LEN], energy_steps[RECORD_LEN], domlen_steps[RECORD_LEN];
    int m_e_ind = 0; // index to record which element of mag_steps[] to store the mag (same thing for energy)
    for (i = 0; i < RECORD_LEN; i++)
    {
        mag_steps[i] = 0.0;
        energy_steps[i] = 0.0;
        domlen_steps[i] = 0.0;
    }
    
    char spin_fn[FILE_NAME_LEN], MED_fn[FILE_NAME_LEN];// mag_fn[FILE_NAME_LEN], energy_fn[FILE_NAME_LEN], dom_len_fn[FILE_NAME_LEN]; // output filenames
    Commandlineparse(argc, argv, &age, &total_MCs, &repeat_num, &seed, name);
    gsl_rng_set(rnd, (unsigned long int)seed); // set up random seed before initialize()
    total_rnd_steps = total_MCs * age.get_N();
    age.initialize();
    int gs_e = -age.get_N() * 2;
    printf("L=%d; T=%.2f; total_MC=%.2f; repeat_num=%d; seed=%d; name=%s (len=%d)\n", age.get_L(), age.get_T(),total_MCs, repeat_num, seed, name, (int)strlen(name));
    // Test input:
    age.print_info();
    //age.neighbor_check(15);
    //age.check_E();
    //age.check_M();
    
    set_record_steps(rand_steps_record, total_MCs, age.get_N());

    set_file_name(spin_fn, "spin", name, &age, repeat_num);
    //set_file_name(mag_fn, "mage", name, &age, repeat_num);
    //set_file_name(energy_fn, "ener", name, &age, repeat_num);
    //set_file_name(dom_len_fn, "doma", name, &age, repeat_num);
    set_file_name(MED_fn, "MED", name, &age, repeat_num);
    
    FILE *fp;
    fp = fopen(spin_fn, "w"); 
    for (iter = 0; iter < repeat_num; iter++)
    {
        /*********** System Initialization *************/
        total_rnd_steps_ran = 0.0;
        age.initialize();
        /**************** Simulation *******************/
        m_e_ind = 0;
        while (total_rnd_steps_ran < total_rnd_steps)
        {
            // if it is the step to record, store the values
            if ( abs(total_rnd_steps_ran - rand_steps_record[m_e_ind]) < 1e-10 && m_e_ind < RECORD_LEN )
            {
                //printf("%d : %.2f\n", m_e_ind, rand_steps_record[m_e_ind]);
                mag_steps[m_e_ind] += age.get_M();
                energy_steps[m_e_ind] += age.get_E();
                domlen_steps[m_e_ind] += (age.get_Espin() - gs_e) / 2.0;
                if(iter == repeat_num-1)
                {
                    if(m_e_ind < 10 || (m_e_ind<100 && m_e_ind%10==0)||m_e_ind%200==0)
                    {
                        spins = age.get_spins();
                        fprintf(fp, "%.4f", total_rnd_steps_ran/age.get_N());
                        for (i = 0; i < age.get_N(); i++)
                        {
                            fprintf(fp,  ",%-2d", spins[i]);
                        }
                        fprintf(fp, "\n");
                    }
                }
                m_e_ind++;
            }
            age.update();            
            total_rnd_steps_ran += 1.0;
            //if(total_rnd_steps_ran-1.0<1e-12) { sleep(1);}
        }
        //out.clear();  
        //age.check_E();
        //age.check_M();
        if(repeat_num>20)
        {
            if( iter%(repeat_num/20)==0 )
                printf("   %.2f%% completed\n", (double)iter/repeat_num*100);
        }
        else
            printf("   %.2f%% completed\n", (double)iter/repeat_num*100);
    }
    
    fclose(fp);
    printf("\n\n   Simulation is Done!\n");

    /**************** Output ***********************/
    //getchar(); //This line keeps the gnuplot window open after the code runs through.
    fp=fopen(MED_fn, "w");
    fprintf(fp, "2D FM Ising; L=%d; T=%.2f; s_hi=%.2f; s_Hi=%.2f; MC steps= %-.2e; repeat=%d; Final_E=%.2f; Final_M=%d;\n",
            age.get_L(), age.get_T(), age.get_s_hi(),age.get_s_Hi(), total_MCs, repeat_num, age.get_E(), age.get_M());

    fprintf(fp, "mc_step,mag,energy,domain_len\n");

    for(i=0; i < RECORD_LEN; i++)
    {
        fprintf(fp, "%.0f,%.8f,%.8f,%.8f\n", rand_steps_record[i]/age.get_N(), 
        (double)mag_steps[i]/repeat_num/age.get_N(), energy_steps[i]/repeat_num/age.get_N(), domlen_steps[i]/repeat_num);
    }
    fclose(fp);
    age.print_info();
    printf("\n\n   Data Export is Done!\n");
    
    return(0);
}
