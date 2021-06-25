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
int** createBA(int m0, int m, int T);
int rand_node(int *k, int *tot);
void add_node(int num,int **graph, int *k, int *tot);
void read_grades(int **graph,int *k,int *tot);
int size_poblacion(state_t *state, state_t poblacion, int N);

int m0, m, T;

main(int argc, char *argv[])
{
    int **red;
    double seed_b=0.1, seed_w=0.1;
    double beta_b, beta_w, eps, mu_w;
    double d_eps, d_beta_w;

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
    state_t *state=calloc(m0+T, sizeof*state);
    for(int i=0; i<(m0+T); i++){
        double z=Random();
        if(z<seed_b){
            state[i]=B;
        }
        else if(z>seed_b && z<(seed_b+seed_w)){
            state[i]=W;
        }
    }

    runGillespie(red, state, (m0+T), mu_w, beta_b, beta_w, eps);

    for(int i=0;i<(m0+T);i++)
        free(red[i]);
    free(red);
    free(state);
}
//0:V; 1:B; 2:L; 3:LB; 4:WB; 5:W; 6:P
void runGillespie(int **red, state_t *state, int N, double mu_w, double beta_b, double beta_w, double eps)
{
    char filename[100];
    sprintf(filename, "modelo_gill_sf.txt");
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
