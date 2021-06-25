#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "generadorPR.h"

typedef enum
{
    V, B, L, LB, WB, W, P, n_poblaciones
}
state_t;

void runGillespie(int **network, state_t *state, int N, double mu_w, double beta_b, double beta_w, double eps);
int** createER(int N, double p);
int size_poblacion(state_t *state, state_t poblacion, int N);

main(int argc, char *argv[])
{
    int N=10000, **red;
    double p=0.01;
    double seed_b=0.1, seed_w=0.1;
    double beta_b, beta_w, eps, mu_w;
    double d_eps, d_beta_w;

    beta_b=1;
    beta_w=0.8;
    mu_w=1;
    eps=0.6;

    ini_ran(411);

    red=createER(N, p);

    //0:V; 1:B; 2:L; 3:LB; 4:WB; 5:W; 6:P
    state_t *state=calloc(N, sizeof*state);
    for(int i=0; i<N; i++){
        double z=Random();
        if(z<seed_b){
            state[i]=B;
        }
        else if(z>seed_b && z<(seed_b+seed_w)){
            state[i]=W;
        }
    }

    runGillespie(red, state, N, mu_w, beta_b, beta_w, eps);

    for(int i=0;i<N;i++)
        free(red[i]);
    free(red);
    free(state);
}
//0:V; 1:B; 2:L; 3:LB; 4:WB; 5:W; 6:P
void runGillespie(int **red, state_t *state, int N, double mu_w, double beta_b, double beta_w, double eps)
{
    char filename[100];
    sprintf(filename, "modelo_gill_hom.txt");
    FILE *salida=fopen(filename, "w");
    fprintf(salida, "t\trho_v\trho_b\trho_l\trho_{lb}\trho_{wb}\trho_w\trho_p\n");

    int k_inf;

    double **ritmo_trans, *cumulative, *cumulative_ind;
    double omega, tau, u, v;

    ritmo_trans=malloc(N*sizeof*ritmo_trans);
    for(int i=0; i<N; i++)
        ritmo_trans[i]=calloc(n_poblaciones, sizeof*ritmo_trans[i]);

    // Initialize transition vector
    for(int i=0; i<N; i++){
        switch(state[i]){
            case V:
                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==B)
                        k_inf++;
                }
                ritmo_trans[i][B]=k_inf*beta_b;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==LB)
                        k_inf++;
                }
                ritmo_trans[i][LB]=k_inf*beta_b;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==WB)
                        k_inf++;
                }
                ritmo_trans[i][WB]=k_inf*beta_b;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==W)
                        k_inf++;
                }
                ritmo_trans[i][W]=k_inf*beta_w;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==WB)
                        k_inf++;
                }
                ritmo_trans[i][WB]+=k_inf*beta_w;

                break;

            case B:
                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==W)
                        k_inf++;
                }
                ritmo_trans[i][W]=k_inf*beta_w;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==WB)
                        k_inf++;
                }
                ritmo_trans[i][WB]=k_inf*beta_w;

                break;

            case L:
                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==B)
                        k_inf++;
                }
                ritmo_trans[i][B]=k_inf*beta_b;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==LB)
                        k_inf++;
                }
                ritmo_trans[i][LB]=k_inf*beta_b;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==WB)
                        k_inf++;
                }
                ritmo_trans[i][WB]=k_inf*beta_b;

                ritmo_trans[i][L]=eps;

                break;

            case LB:
                ritmo_trans[i][LB]=eps;

                break;

            case W:
                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==B)
                        k_inf++;
                }
                ritmo_trans[i][B]=k_inf*beta_b;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==LB)
                        k_inf++;
                }
                ritmo_trans[i][LB]=k_inf*beta_b;

                k_inf=0;
                for(int j=0; j<N; j++){
                    if(red[i][j]==1 && state[j]==WB)
                        k_inf++;
                }
                ritmo_trans[i][WB]=k_inf*beta_b;

                ritmo_trans[i][W]=mu_w;

                break;

            case WB:
                ritmo_trans[i][WB]=mu_w;

                break;
        }
    }

    // Evolution
    double t=0;
    do{
        printf("%lf\n", t);
        fprintf(salida, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", t, (double)size_poblacion(state, V, N)/N, (double)size_poblacion(state, B, N)/N, (double)size_poblacion(state, L, N)/N, (double)size_poblacion(state, LB, N)/N,
                (double)size_poblacion(state, WB, N)/N, (double)size_poblacion(state, W, N)/N, (double)size_poblacion(state, P, N)/N);

        omega=0;
        for(int i=0; i<N; i++){
            for(int j=0; j<n_poblaciones; j++){
                omega+=ritmo_trans[i][j];
            }
        }

        tau=log(1-Random())/(-omega);
        t+=tau;

        // Get cumulative
        cumulative=calloc(N, sizeof*cumulative);
        for(int i=0; i<n_poblaciones; i++)
            cumulative[0]+=ritmo_trans[0][i]/omega;
        for(int i=1; i<N; i++){
            for(int j=0; j<n_poblaciones; j++){
                cumulative[i]+=ritmo_trans[i][j]/omega;
            }
            cumulative[i]+=cumulative[i-1];
        }

        u=Random();
        for(int i=0; i<N; i++){
            if(u<cumulative[i]){
                double omega_ind=0;
                for(int j=0; j<n_poblaciones; j++)
                    omega_ind+=ritmo_trans[i][j];
                cumulative_ind=calloc(n_poblaciones, sizeof*cumulative_ind);
                cumulative_ind[V]=ritmo_trans[i][V]/omega_ind;
                for(int j=1; j<n_poblaciones; j++)
                    cumulative_ind[j]=cumulative_ind[j-1]+ritmo_trans[i][j]/omega_ind;

                v=Random();
                for(int j=0; j<n_poblaciones; j++){
                    if(v<cumulative_ind[j]){
                        switch(state[i]){
                            case V:
                                if(j==B || j==LB){
                                    state[i]=B;
                                    ritmo_trans[i][B]=ritmo_trans[i][LB]=0;
                                    k_inf=0;
                                    for(int k=0; k<N; k++){
                                        if(red[i][k]==1 && state[k]==WB)
                                            k_inf++;
                                    }
                                    ritmo_trans[i][WB]-=k_inf*beta_b;
                                    for(int k=0; k<N; k++){
                                        if(red[i][k]==1){
                                            switch(state[k]){
                                                case V:
                                                case L:
                                                case W:
                                                    ritmo_trans[k][B]+=beta_b;
                                                    break;
                                            }
                                        }
                                    }
                                }
                                else if(j==W){
                                    state[i]=L;
                                    ritmo_trans[i][W]=0;
                                    ritmo_trans[i][L]=eps;
                                    k_inf=0;
                                    for(int k=0; k<N; k++){
                                        if(red[i][k]==1 && state[k]==WB)
                                            k_inf++;
                                    }
                                    ritmo_trans[i][WB]-=k_inf*beta_w;
                                }
                                else if(j==WB){
                                    if(Random()<(beta_b/(beta_b+beta_w))){
                                        state[i]=B;
                                        ritmo_trans[i][B]=ritmo_trans[i][LB]=0;
                                        k_inf=0;
                                        for(int k=0; k<N; k++){
                                            if(red[i][k]==1 && state[k]==WB)
                                                k_inf++;
                                        }
                                        ritmo_trans[i][WB]-=k_inf*beta_b;
                                        for(int k=0; k<N; k++){
                                            if(red[i][k]==1){
                                                switch(state[k]){
                                                    case V:
                                                    case L:
                                                    case W:
                                                        ritmo_trans[k][B]+=beta_b;
                                                        break;
                                                }
                                            }
                                        }
                                    }
                                    else{
                                        state[i]=L;
                                        ritmo_trans[i][W]=0;
                                        ritmo_trans[i][L]=eps;
                                        k_inf=0;
                                        for(int k=0; k<N; k++){
                                            if(red[i][k]==1 && state[k]==WB)
                                                k_inf++;
                                        }
                                        ritmo_trans[i][WB]-=k_inf*beta_w;
                                    }
                                }
                                break;
                            case B:
                                state[i]=LB;
                                ritmo_trans[i][W]=ritmo_trans[i][WB]=0;
                                ritmo_trans[i][LB]=eps;
                                for(int k=0; k<N; k++){
                                    if(red[i][k]==1){
                                        switch(state[k]){
                                            case V:
                                            case L:
                                            case W:
                                                ritmo_trans[k][B]-=beta_b;
                                                ritmo_trans[k][LB]+=beta_b;
                                                break;
                                        }
                                    }
                                }
                                break;
                            case L:
                                if(j==B || j==LB || j==WB){
                                    state[i]=LB;
                                    ritmo_trans[i][B]=ritmo_trans[i][WB]=ritmo_trans[i][L]=0;
                                    ritmo_trans[i][LB]=eps;
                                    for(int k=0; k<N; k++){
                                        if(red[i][k]==1){
                                            switch(state[k]){
                                                case V:
                                                case L:
                                                case W:
                                                    ritmo_trans[k][LB]+=beta_b;
                                                    break;
                                            }
                                        }
                                    }
                                }
                                else if(j==L){
                                    state[i]=W;
                                    ritmo_trans[i][L]=0;
                                    ritmo_trans[i][W]=mu_w;
                                    for(int k=0; k<N; k++){
                                        if(red[i][k]==1){
                                            switch(state[k]){
                                                case V:
                                                case B:
                                                    ritmo_trans[k][W]+=beta_w;
                                                    break;
                                            }
                                        }
                                    }
                                }
                                break;
                            case LB:
                                state[i]=WB;
                                ritmo_trans[i][LB]=0;
                                ritmo_trans[i][WB]=mu_w;
                                for(int k=0; k<N; k++){
                                    if(red[i][k]==1){
                                        switch(state[k]){
                                            case V:
                                                ritmo_trans[k][LB]-=beta_b;
                                                ritmo_trans[k][WB]+=beta_b+beta_w;
                                                break;
                                            case L:
                                            case W:
                                                ritmo_trans[k][LB]-=beta_b;
                                                ritmo_trans[k][WB]+=beta_b;
                                                break;
                                            case B:
                                                ritmo_trans[k][WB]+=beta_w;
                                                break;
                                        }
                                    }
                                }
                                break;
                            case W:
                                if(j==B || j==LB || j==WB){
                                    state[i]=WB;
                                    ritmo_trans[i][W]=0;
                                    ritmo_trans[i][B]=ritmo_trans[i][LB]=0;
                                    ritmo_trans[i][WB]=mu_w;
                                    for(int k=0; k<N; k++){
                                        if(red[i][k]==1){
                                            switch(state[k]){
                                                case V:
                                                    ritmo_trans[k][W]-=beta_w;
                                                    ritmo_trans[k][WB]+=beta_b+beta_w;
                                                    break;
                                                case B:
                                                    ritmo_trans[k][W]-=beta_w;
                                                    ritmo_trans[k][WB]+=beta_w;
                                                    break;
                                                case L:
                                                case W:
                                                    ritmo_trans[k][WB]+=beta_b;
                                                    break;
                                            }
                                        }
                                    }
                                }
                                else if(j==W){
                                    state[i]=P;
                                    ritmo_trans[i][B]=ritmo_trans[i][LB]=ritmo_trans[i][WB]=ritmo_trans[i][W]=0;
                                    for(int k=0; k<N; k++){
                                        if(red[i][k]==1){
                                            switch(state[k]){
                                                case V:
                                                case B:
                                                    ritmo_trans[k][W]-=beta_w;
                                                    break;
                                            }
                                        }
                                    }
                                }
                                break;
                            case WB:
                                state[i]=P;
                                ritmo_trans[i][WB]=0;
                                for(int k=0; k<N; k++){
                                    if(red[i][k]==1){
                                        switch(state[k]){
                                            case V:
                                                ritmo_trans[k][WB]-=beta_b+beta_w;
                                                break;
                                            case B:
                                                ritmo_trans[k][WB]-=beta_w;
                                                break;
                                            case L:
                                            case W:
                                                ritmo_trans[k][WB]-=beta_b;
                                                break;
                                        }
                                    }
                                }
                                break;
                        }
                        break;
                    }
                }

                free(cumulative_ind);
                break;
            }
        }
        free(cumulative);
    }while(omega>1e-6);

    // Release memory
    free(ritmo_trans);

    fclose(salida);
}

int** createER(int N, double p)
{
    int **red;

    red=malloc(N*sizeof*red);
    for(int i=0; i<N; i++)
        red[i]=calloc(N, sizeof*red[i]);

    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            if(Random()<p){
                red[i][j]=red[j][i]=1;
            }
        }
    }

    return red;
}

int size_poblacion(state_t *state, state_t poblacion, int N)
{
    int size=0;

    for(int i=0; i<N; i++){
        if(state[i]==poblacion){
            size++;
        }
    }

    return size;
}
