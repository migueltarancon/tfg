#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "generadorPR.h"

//global variables
int m0, m, T;

//functions
int rand_node(int *k, int *tot);
void add_node(int num,int **graph, int *k, int *tot);
void read_grades(int **graph,int *k,int *tot);

main(int argc, char *argv[])
{
    int i,j,t;
	int **graph;
	int *k; //grade k for each node
	int tot; //sum of all grades sumj(k[j])
	FILE *f;

    if(argc==4){
        sscanf(argv[1], "%d", &m0);
        sscanf(argv[2], "%d", &m);
        sscanf(argv[3], "%d", &T);
        printf("Numero de nodos de la red inicial: %d\n", m0);
        printf("Numero de nuevos links: %d\n", m);
        printf("Pasos de tiempo: %d\n\n", T);
    }
    else{
        printf("No se ha introducido numero de nodos inicial, numero de nuevos links y/o pasos de tiempo.\nSe toman valores por defecto: m0=5; m=3; T=995\n\n");
        m0=50;
        m=2;
        T=9950;
    }

	ini_ran(1999);

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
		printf("Anadido nodo %d \n",m0+t+1);
	}

	/*for(i=0;i<(m0+T);i++){
		for(j=0;j<(m0+T);j++) printf("%d ",graph[i][j]);
		printf("\n Fila %d\n",i);
	}

	f=fopen("red_baravasi_albert(scale-free).txt","w");
	for(i=0;i<(m0+T-1);i++){
		for(j=0;j<(m0+T);j++) fprintf(f,"%d ",graph[i][j]);
		fprintf(f,"\n");
	}
	for(j=0;j<(m0+T);j++) fprintf(f,"%d ",graph[m0+T-1][j]);
	fclose(f);*/

	printf("\nNumero de links: %d", tot/2);

	f=fopen("red_barabasi_albert(scale-free).csv","w");
	for(i=0;i<(m0+T);i++){ //superior part of the matrix
		for(j=(i+1);j<(m0+T);j++){
			if(graph[i][j]==1)
				fprintf(f, "%d\t%d\n", i, j);
		}
	}
    fclose(f);
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
