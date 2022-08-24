#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define ENVIRONMENT_SIZE 256
#define OMEGA 1.9
#define ETA 1
#define VIRUS_LIM 1500
#define FIRST_VIRUS 254, 128
#define TOLERANCE 0.0001

//Well, actually this project is using DLA to simulate the growth of bacteria not viruses. But because it is a mistake from the beginning so...
//... we just too lazy to re-correct all variables' name and comments. Sorry for that inconvenience.

//Virus structure
typedef struct VirusPlace{
    int x; int y;
}Virus;

//Two arrays to store viruses and candidates
Virus virus[ENVIRONMENT_SIZE*ENVIRONMENT_SIZE];//the size of the environment can be set above
int nVirus = 0;
Virus candidate[ENVIRONMENT_SIZE*ENVIRONMENT_SIZE];
int nCandidate = 0;

//Array that stores the probability to become a virus of candidates
double chance[ENVIRONMENT_SIZE*ENVIRONMENT_SIZE];

//Matrix of nutrient concentration
double c[ENVIRONMENT_SIZE][ENVIRONMENT_SIZE];
//Matrix of viruses' and candidates' location
int grow[ENVIRONMENT_SIZE][ENVIRONMENT_SIZE];

void addVirus(int x, int y, int index);//To add a candidate has recently become virus to virus[]
void init();//Generating the initial environment and add the first virus to it
void sor();//Calculating c[][] using SOR method
void eat();//To set all nutrient concentration on virus places to 0
void computeProbability();//Calculating the probability to become a virus of each candidate
void growth();//To consider if a candidate would become a virus depending on probabilities and add it to the environment
void solve();//Combination of functions that run each growth

void addVirus(int x, int y, int index){
    if(grow[x][y] == 1) return;//check if there already has a virus on that place

    Virus newVirus;
    newVirus.x = x; newVirus.y = y;
    virus[nVirus] = newVirus;
    nVirus++;
    grow[x][y] = 1;//storing viruses' location

    //removing the newly virus from candidate[]
    for(int i = 0; i < nCandidate; i++){
        if(candidate[i].x == x && candidate[i].y == y){
            for(int j = i; j < nCandidate; j++){
                candidate[j] = candidate[j+1];
            }
        }
    }
    nCandidate--;

    //adding new candidates to candidate[]
    for(int j = x-1; j <= x+1; j++){
        for(int k = y-1; k <= y+1; k++){
            if(abs(j-x) != abs(k-y) && grow[j][k] == 0 && j >= 0 && k >= 0 && j < ENVIRONMENT_SIZE && k < ENVIRONMENT_SIZE){
                Virus newCandidate;
                newCandidate.x = j; newCandidate.y = k;
                candidate[nCandidate] = newCandidate;
                nCandidate++;
                grow[j][k] = 2;//storing candidate's location
            }
        }
    }
}

void init(){
    //initializing environment
    for(int i = 0; i < ENVIRONMENT_SIZE; i++){
        for(int j = 0; j < ENVIRONMENT_SIZE; j++){
            grow[i][j] = 0;
            c[i][j] = 1.0;
        }
    }

    //adding the first virus
    Virus firstVirus;
    candidate[nCandidate] = firstVirus;
    nCandidate++;
    addVirus(FIRST_VIRUS, 0);//location of the first virus can set above
}

void sor(){
    double delta;//variable to store maximum error each iteration
    double tmp[ENVIRONMENT_SIZE][ENVIRONMENT_SIZE];//to store old matrix to compare to new one to get errors

    //setting boundary conditions
    for(int i = 0; i < ENVIRONMENT_SIZE; i++){
        for(int j = 0; j < ENVIRONMENT_SIZE; j++){
            if(i == 0) c[i][j] = 1.0;
            if(j == 0 || j == ENVIRONMENT_SIZE - 1) c[i][j] = 1.0 - (i*1.0/(ENVIRONMENT_SIZE-1));
            if(i == ENVIRONMENT_SIZE - 1) c[i][j] = 0.0;
        }
    }

    //sor
    do{
        //like mentioning above about tmp[][]
        for(int i = 0; i < ENVIRONMENT_SIZE; i++){
            for(int j = 0; j < ENVIRONMENT_SIZE; j++){
                tmp[i][j] = c[i][j];
            }
        }

        //initializing error equal to 0
        delta = 0;

        //updating concentration of all red grid points (red-black ordering parallel simulation method)
        for(int i = 1; i < ENVIRONMENT_SIZE-1; i++){
            for(int j = 1; j < ENVIRONMENT_SIZE-1; j++){
                if((i+j) % 2 == 0){
                    c[i][j] = OMEGA/4*(c[i-1][j]+c[i+1][j]+c[i][j+1]+c[i][j-1])+(1-OMEGA)*c[i][j];//we choosing omega equal to 1.9 in this problem...
                                                                                                  //...of course it can change above just remember...
                                                                                                  //...choose between 1 and 2 for convergence
                }
            }
        }

        //updating concentration of all black grid points using that of red ones computed above
        for(int i = 1; i < ENVIRONMENT_SIZE-1; i++){
            for(int j = 1; j < ENVIRONMENT_SIZE-1; j++){
                if((i+j) % 2 == 1){
                    c[i][j] = OMEGA/4*(c[i-1][j]+c[i+1][j]+c[i][j+1]+c[i][j-1])+(1-OMEGA)*c[i][j];
                }
            }
        }

        //concentration can not be negative
        for(int i = 0; i < ENVIRONMENT_SIZE; i++){
            for(int j = 0; j < ENVIRONMENT_SIZE; j++){
                if(c[i][j] < 0) c[i][j] = 0;
            }
        }

        //updating error
        for(int i = 0; i < ENVIRONMENT_SIZE; i++){
            for(int j = 0; j < ENVIRONMENT_SIZE; j++){
                if(abs(c[i][j] - tmp[i][j]) > TOLERANCE) delta = abs(c[i][j] - tmp[i][j]);
            }
        }
    }while(delta > TOLERANCE);//if error did not be updated so it's (initial) zero and we can stop iteration (tolerance can be define above)
}

void eat(){
    //mentioning above
    for(int i = 0; i < nVirus; i++){
        c[virus[i].x][virus[i].y] = 0;
    }
}

void computeProbability(){
    //also mentioning above
    double sum = 0;

    for(int i = 0; i < nCandidate; i++){
        sum += pow(c[candidate[i].x][candidate[i].y], ETA);//ETA can be set at defining to change the way virus growing
    }

    for(int j = 0; j < nCandidate; j++){
        chance[j] = pow(c[candidate[j].x][candidate[j].y], ETA)/sum;
    }
}

void growth(){
    Virus currentCandidate[ENVIRONMENT_SIZE*ENVIRONMENT_SIZE];//storing currently candidates because each time addVirus() is called...
                                                              //... a candidate will be remove from candidate[] making change in indexes
    for(int i = 0; i < nCandidate; i++) {
        currentCandidate[i].x = candidate[i].x;
        currentCandidate[i].y = candidate[i].y;
    }

    int tmp = nCandidate;//yeah, surely nCandidate will change too

    //for each candidate, we draw a random number between 0 and 1 to compare to it's probability. If smaller, turn it into a virus
    for(int i = 0; i < tmp; i++){
        double r = rand()*1.0/RAND_MAX;
        if(r < chance[i]) addVirus(currentCandidate[i].x, currentCandidate[i].y, i);
    }
}

void solve(){
    sor();//re-calculating environment
    eat();//feeding the viruses
    computeProbability();//nothing more than function name
    growth();//growing viruses
}

int main(){
    srand(time(NULL));//to change random numbers each run

    init();

    //file 'output_simulate.txt' is for simulating the growing process of the bacteria. However, due to the very large writing to file,...
    //... it takes a lot of time, so we put it in comments. If you want to try it, please take it out.

    /*FILE *f0 = fopen("output_simulate.txt", "w");
    fprintf(f0, "");
    fclose(f0);
    f0 = fopen("output_simulate.txt", "a");*/

    do{

        /*for(int i = 0; i < ENVIRONMENT_SIZE; i++){
            for(int j = 0; j < ENVIRONMENT_SIZE; j++){
                if(grow[i][j] == 1) fprintf(f0, "%d ", 1);
                else fprintf(f0, "%d ", 0);
            }
            fprintf(f0, "\n");
        }*/

        solve();
    }while(nVirus < VIRUS_LIM);//continue growing till we have enough viruses (can be define above)

    /*fclose(f0);*/

    //writing the last matrix of nutrient concentration (the environment) to files (for using on visualizing in MATLAB)
    FILE *f1 = fopen("output_nutrients.txt", "w");
    FILE *f2 = fopen("output_bacteria.txt", "w");
    for(int i = 0; i < ENVIRONMENT_SIZE; i++){
        for(int j = 0; j < ENVIRONMENT_SIZE; j++){
            fprintf(f1, "%.4f ", c[i][j]);
            if(grow[i][j] == 2) grow[i][j] = 0;
            fprintf(f2, "%d ", grow[i][j]);
        }
        fprintf(f1, "\n");
        fprintf(f2, "\n");
    }
    fclose(f1);
    fclose(f2);

    return 0;
}
