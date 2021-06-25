#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "generadorPR.h"

typedef enum
{
    V, B, L, LB, WB, W, P, n_poblaciones
}
state_t;

typedef struct tDegree{
    int *k;
    double *Pk;
    double **rho;
    double mean_k;
    int size;
}
tDegree;

#define MAX(x,y) (x>y?x:y)

void runRK(tDegree *v, double mu_w, double lambda_b, double lambda_w, double eps_prim, double taumax, double dtau);
double** calcular_sumar_k(tDegree *v, double lambda_b, double lambda_w, double eps_prim, double dtau, double (*funciones_ev_temp[]) (int, tDegree *, double **, double, double, double));
tDegree* initializeDegree(int **network, int N);
int** createER(int N, double p);
double calcular_rho(tDegree *v, state_t pobl);
double funcion_rhoV(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoL(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoLB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoWB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoW(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoP(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);

main(int argc, char *argv[])
{
    int N=10000, **red;
    double p=0.01;
    double tmax=20, dt=0.01;
    double beta_b, beta_w, eps, mu_w, lambda_b, lambda_w, eps_prim, taumax, dtau;

    beta_b=1;
    beta_w=0.8;
    mu_w=1;
    eps=0.6;

    ini_ran(11);

    red=createER(N, p);

    //0:V; 1:B; 2:L; 3:LB; 4:WB; 5:W; 6:P
    tDegree *v=initializeDegree(red, N);
    for(int i=0;i<v->size;i++){
        v->rho[i][V]=0.8;
        v->rho[i][B]=0.1;
        v->rho[i][W]=0.1;
    }

    lambda_b=beta_b/mu_w;
    lambda_w=beta_w/mu_w;
    eps_prim=eps/mu_w;
    taumax=mu_w*tmax;
    dtau=mu_w*dt;

    runRK(v, mu_w, lambda_b, lambda_w, eps_prim, taumax, dtau);

    free(v);
    for(int i=0;i<N;i++)
        free(red[i]);
    free(red);
}

void runRK(tDegree *v, double mu_w, double lambda_b, double lambda_w, double eps_prim, double taumax, double dtau)
{
    double tau=0;
    double (*funciones_ev_temp[]) (int, tDegree *, double **, double, double, double)={funcion_rhoV, funcion_rhoB, funcion_rhoL, funcion_rhoLB, funcion_rhoWB, funcion_rhoW, funcion_rhoP};

    char filename[100];
    sprintf(filename, "modelo_rk_hom.txt");
    FILE *salida=fopen(filename, "w");
    fprintf(salida, "t\trho_v\trho_b\trho_l\trho_{lb}\trho_{wb}\trho_w\trho_p\n");

    while(tau<taumax){
        fprintf(salida, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", tau/mu_w, calcular_rho(v, V), calcular_rho(v, B), calcular_rho(v, L), calcular_rho(v, LB), calcular_rho(v, WB), calcular_rho(v, W), calcular_rho(v, P));
        //printf("%lf\n",distrs_pobls[V]+distrs_pobls[B]+distrs_pobls[L]+distrs_pobls[LB]+distrs_pobls[WB]+distrs_pobls[W]+distrs_pobls[P]);

        double **suma_k=calcular_sumar_k(v, lambda_b, lambda_w, eps_prim, dtau, funciones_ev_temp);
        for(int i=0; i<v->size; i++){
            for(int n_fun=0; n_fun<n_poblaciones; n_fun++){
                v->rho[i][n_fun]+=1.0/6*suma_k[i][n_fun];
            }
        }
        tau+=dtau;
        for(int i=0;i<v->size;i++)
            free(suma_k[i]);
        free(suma_k);
    }

    fclose(salida);
}

double** calcular_sumar_k(tDegree *v, double lambda_b, double lambda_w, double eps_prim, double dtau, double (*funciones_ev_temp[]) (int, tDegree *, double **, double, double, double))
{
    double **k1=malloc(v->size*sizeof*k1);
    for(int i=0; i<v->size; i++)
        k1[i]=malloc(n_poblaciones*sizeof*k1[i]);
    double **k2=malloc(v->size*sizeof*k2);
    for(int i=0; i<v->size; i++)
        k2[i]=malloc(n_poblaciones*sizeof*k2[i]);
    double **k3=malloc(v->size*sizeof*k3);
    for(int i=0; i<v->size; i++)
        k3[i]=malloc(n_poblaciones*sizeof*k3[i]);
    double **k4=malloc(v->size*sizeof*k4);
    for(int i=0; i<v->size; i++)
        k4[i]=malloc(n_poblaciones*sizeof*k4[i]);

    double **rho_aux=malloc(v->size*sizeof*rho_aux);
    for(int i=0; i<v->size; i++)
        rho_aux[i]=malloc(n_poblaciones*sizeof*rho_aux[i]);

    double **suma_k=malloc(v->size*sizeof*suma_k);
    for(int i=0; i<v->size; i++)
        suma_k[i]=malloc(n_poblaciones*sizeof*suma_k[i]);

    for(int i=0; i<v->size; i++){
        for(int j=0; j<n_poblaciones; j++){
            rho_aux[i][j]=v->rho[i][j];
        }
    }
    for(int i=0; i<v->size; i++){
        for(int j=0; j<n_poblaciones; j++){
            k1[i][j]=dtau*funciones_ev_temp[j] (i, v, rho_aux, lambda_b, lambda_w, eps_prim);
        }
    }

    for(int i=0; i<v->size; i++){
        for(int j=0; j<n_poblaciones; j++){
            rho_aux[i][j]=v->rho[i][j]+0.5*k1[i][j];
        }
    }
    for(int i=0; i<v->size; i++){
        for(int j=0; j<n_poblaciones; j++){
            k2[i][j]=dtau*funciones_ev_temp[j] (i, v, rho_aux, lambda_b, lambda_w, eps_prim);
        }
    }

    for(int i=0; i<v->size; i++){
        for(int j=0; j<n_poblaciones; j++){
            rho_aux[i][j]=v->rho[i][j]+0.5*k2[i][j];
        }
    }
    for(int i=0; i<v->size; i++){
        for(int j=0; j<n_poblaciones; j++){
            k3[i][j]=dtau*funciones_ev_temp[j] (i, v, rho_aux, lambda_b, lambda_w, eps_prim);
        }
    }

    for(int i=0; i<v->size; i++){
        for(int j=0; j<n_poblaciones; j++){
            rho_aux[i][j]=v->rho[i][j]+k3[i][j];
        }
    }
    for(int i=0; i<v->size; i++){
        for(int j=0; j<n_poblaciones; j++){
            k4[i][j]=dtau*funciones_ev_temp[j] (i, v, rho_aux, lambda_b, lambda_w, eps_prim);
        }
    }

    for(int i=0; i<v->size; i++){
        for(int n_fun=0; n_fun<n_poblaciones; n_fun++){
            suma_k[i][n_fun]=k1[i][n_fun]+2*k2[i][n_fun]+2*k3[i][n_fun]+k4[i][n_fun];
        }
    }

    for(int i=0;i<v->size;i++)
        free(k1[i]);
    free(k1);
    for(int i=0;i<v->size;i++)
        free(k2[i]);
    free(k2);
    for(int i=0;i<v->size;i++)
        free(k3[i]);
    free(k3);
    for(int i=0;i<v->size;i++)
        free(k4[i]);
    free(k4);
    for(int i=0;i<v->size;i++)
        free(rho_aux[i]);
    free(rho_aux);

    return suma_k;
}

// Initialize rho vector
tDegree* initializeDegree(int **network, int N)
{
    tDegree *v;
    int *node_k, *n_k, d, sum, max_k, cont;

    // Get degree info
    n_k = calloc(N, sizeof*n_k);
    node_k = calloc(N, sizeof*node_k);

    max_k = 0;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++)
            node_k[i] += network[i][j];

        n_k[node_k[i]]++;

        max_k = MAX(max_k,node_k[i]);
    }
    if(n_k[0]!=0)
        max_k++;

    sum = 0;
    for(int i=0;i<N;i++)
        sum += n_k[i];

    // Add data to degree struct
    v = malloc(sizeof*v);
    v->k = malloc(max_k*sizeof*v->k);
    v->Pk = malloc(max_k*sizeof*v->Pk);
    v->rho = malloc(max_k*sizeof*v->rho);
    for(int i=0; i<max_k; i++)
        v->rho[i]=malloc(n_poblaciones*sizeof*v->rho[i]);

    cont = 0;
    for(int i=0;i<N;i++)
        if(n_k[i]>0){
            v->k[cont] = i;
            v->Pk[cont] = n_k[i]/(double)sum;
            for(int j=0; j<n_poblaciones; j++)
                v->rho[cont][j] = 0;
            cont++;
        }
    v->size = cont;

    v->mean_k = 0;
    for(int i=0;i<v->size;i++)
        v->mean_k += v->k[i]*v->Pk[i];

    // Release allocated memory
    free(n_k);
    free(node_k);

    return v;
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

double calcular_rho(tDegree *v, state_t pobl)
{
    double rho_pobl=0;

    for(int i=0; i<v->size; i++){
        rho_pobl+=v->Pk[i]*v->rho[i][pobl];
    }

    return rho_pobl;
}

//0:V; 1:B; 2:L; 3:LB; 4:WB; 5:W; 6:P
double funcion_rhoV(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim)
{
    double theta1=0, theta2=0;
    for(int i=0; i<v->size; i++){
        theta1+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][B]+distrs_pobls[i][LB]+distrs_pobls[i][WB]);
        theta2+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][W]+distrs_pobls[i][WB]);
    }
    return -lambda_b*v->k[ind_grado]*distrs_pobls[ind_grado][V]*theta1-lambda_w*v->k[ind_grado]*distrs_pobls[ind_grado][V]*theta2;
}

double funcion_rhoB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim)
{
    double theta1=0, theta2=0;
    for(int i=0; i<v->size; i++){
        theta1+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][B]+distrs_pobls[i][LB]+distrs_pobls[i][WB]);
        theta2+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][W]+distrs_pobls[i][WB]);
    }
    return lambda_b*v->k[ind_grado]*distrs_pobls[ind_grado][V]*theta1-lambda_w*v->k[ind_grado]*distrs_pobls[ind_grado][B]*theta2;
}

double funcion_rhoL(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim)
{
    double theta1=0, theta2=0;
    for(int i=0; i<v->size; i++){
        theta1+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][W]+distrs_pobls[i][WB]);
        theta2+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][B]+distrs_pobls[i][LB]+distrs_pobls[i][WB]);
    }
    return lambda_w*v->k[ind_grado]*distrs_pobls[ind_grado][V]*theta1-lambda_b*v->k[ind_grado]*distrs_pobls[ind_grado][L]*theta2-eps_prim*distrs_pobls[ind_grado][L];
}

double funcion_rhoLB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim)
{
    double theta1=0, theta2=0;
    for(int i=0; i<v->size; i++){
        theta1+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][W]+distrs_pobls[i][WB]);
        theta2+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][B]+distrs_pobls[i][LB]+distrs_pobls[i][WB]);
    }
    return lambda_w*v->k[ind_grado]*distrs_pobls[ind_grado][B]*theta1+lambda_b*v->k[ind_grado]*distrs_pobls[ind_grado][L]*theta2-eps_prim*distrs_pobls[ind_grado][LB];
}

double funcion_rhoWB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim)
{
    double theta1=0;
    for(int i=0; i<v->size; i++){
        theta1+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][B]+distrs_pobls[i][LB]+distrs_pobls[i][WB]);
    }
    return eps_prim*distrs_pobls[ind_grado][LB]+lambda_b*v->k[ind_grado]*distrs_pobls[ind_grado][W]*theta1-distrs_pobls[ind_grado][WB];
}

double funcion_rhoW(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim)
{
    double theta1=0;
    for(int i=0; i<v->size; i++){
        theta1+=(v->k[i]*v->Pk[i])/v->mean_k*(distrs_pobls[i][B]+distrs_pobls[i][LB]+distrs_pobls[i][WB]);
    }
    return eps_prim*distrs_pobls[ind_grado][L]-lambda_b*v->k[ind_grado]*distrs_pobls[ind_grado][W]*theta1-distrs_pobls[ind_grado][W];
}

double funcion_rhoP(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim)
{
    return distrs_pobls[ind_grado][W]+distrs_pobls[ind_grado][WB];
}
