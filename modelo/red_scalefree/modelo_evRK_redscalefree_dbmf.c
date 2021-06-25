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
int** createBA(int m0, int m, int T);
int rand_node(int *k, int *tot);
void add_node(int num,int **graph, int *k, int *tot);
void read_grades(int **graph,int *k,int *tot);
double calcular_rho(tDegree *v, state_t pobl);
double funcion_rhoV(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoL(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoLB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoWB(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoW(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);
double funcion_rhoP(int ind_grado, tDegree *v, double **distrs_pobls, double lambda_b, double lambda_w, double eps_prim);

int m0, m, T;

main(int argc, char *argv[])
{
    int **red;
    double tmax=55, dt=0.01;
    double beta_b, beta_w, eps, mu_w, lambda_b, lambda_w, eps_prim, taumax, dtau;

    m0=50;
    m=2;
    T=9950;

    beta_b=1;
    beta_w=0.8;
    mu_w=1;
    eps=0.6;

    ini_ran(11);

    red=createBA(m0, m, T);

    //0:V; 1:B; 2:L; 3:LB; 4:WB; 5:W; 6:P
    tDegree *v=initializeDegree(red, (m0+T));
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
    for(int i=0;i<(m0+T);i++)
        free(red[i]);
    free(red);
}

void runRK(tDegree *v, double mu_w, double lambda_b, double lambda_w, double eps_prim, double taumax, double dtau)
{
    double tau=0;
    double (*funciones_ev_temp[]) (int, tDegree *, double **, double, double, double)={funcion_rhoV, funcion_rhoB, funcion_rhoL, funcion_rhoLB, funcion_rhoWB, funcion_rhoW, funcion_rhoP};

    char filename[100];
    sprintf(filename, "modelo_rk_sf.txt");
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

int** createBA(int m0, int m, int T)
{
    int i, j, t, tot;
    int **graph, *k;

    //size of dynamic vectors
	graph=malloc((m0+T)*sizeof(int*));
	for(i=0;i<(m0+T);i++) graph[i]=malloc((m0+T)*sizeof(int));
	k=malloc((m0+T)*sizeof(int));

	/*we start with a complete graph of m0 nodes (initial state)*/

	for(i=0;i<(m0+T);i++) graph[i][i]=k[i]=0; //diagonal elements
	for(i=0;i<(m0+T);i++){ //superior part of the matrix
		for(j=(i+1);j<(m0+T);j++){
			if((i<m0)&&(j<m0))
                if(Random()<0.7)
                    graph[i][j]=1;
                else
                    graph[i][j]=0;
			else
				graph[i][j]=0;
		}
	}
	//low part of the matrix (symmetric)
	for(j=0;j<(m0+T);j++){
		for(i=(j+1);i<(m0+T);i++) graph[i][j]=graph[j][i];
	}

	//initialize vector k (grades)
	read_grades(graph,k,&tot);

	/*TIME EVOLUTION*/
	for(t=0;t<T;t++){
		add_node(m0+t,graph,k,&tot);
	}

	free(k);
	return graph;
}

int rand_node(int *k, int *tot) //preferential attachment
{
	int i;
	double x,y;
	double M;
	double p_x;

	//maximum probability
	M=k[0];
	for(i=1;i<(m0+T);i++) if(k[i]>M) M=k[i];
	M=M/(double)*tot;

	//generate x,y
	x=(m0+T)*Random();
	y=M*Random();

	//p(x)
	p_x=k[(int)x]/(double)*tot;

	while(y>p_x){
		x=(m0+T)*Random();
		y=M*Random();
		p_x=k[(int)x]/(double)*tot;
	}

	x=(int)x;
	return x;
}

void add_node(int num,int **graph, int *k, int *tot) //num: number of the new node added to the graph
{
	int i,j;
	int rep; //this variable takes 0 or 1. If 0, no repeated nodes.
	int k_old[m0+T];
	int chosen[m]; //vector of nodes which will be connected with the new node

	for(i=0;i<(m0+T);i++) k_old[i]=k[i];

	chosen[0]=rand_node(k_old,tot);
	for(i=1;i<m;i++){ //preferential attachment
		rep=1;
		while(rep==1){
			rep=0;
			chosen[i]=rand_node(k_old,tot);
			j=0;
			while((rep==0)&&(j<i)){
				if(chosen[j]==chosen[i]) rep=1;
				j=j+1;
			}
		}
	}
	//This algorithm makes sure that there are not two edges for the same two nodes. If rand_node repeats one connected node (rep=1), it generates the node again.

	for(i=0;i<m;i++){
		graph[chosen[i]][num]=1;
		graph[num][chosen[i]]=1;
		k[chosen[i]]=k[chosen[i]]+1;
	}
	k[num]=m;
	*tot=*tot+2*m;
}

void read_grades(int **graph,int *k,int *tot)
{
	int i,j;

	for(i=0;i<(m0+T);i++){
		for(j=0;j<(m0+T);j++) k[i]=k[i]+graph[i][j];
	}

	*tot=0;
	for(i=0;i<(m0+T);i++) *tot=*tot+k[i];
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
