#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "generadorPR.h"

main(int argc, char *argv[]){
    unsigned int N, L, L_max, **red;
    double p;

    if(argc==3){
        sscanf(argv[1], "%u", &N);
        sscanf(argv[2], "%lf", &p);
        printf("Numero de nodos de la red: %u\n", N);
        printf("Probabilidad de links de la red: %lf\n\n", p);
    }
    else{
        printf("No se ha introducido numero de nodos y/o probabilidad.\nSe toman valores por defecto: N=100; p=0.05\n\n");
        N=50;
        p=0.2;
    }

    ini_ran(11);
    L_max=N*(N-1)/2;
    red=(unsigned int **)malloc(2*sizeof(unsigned int *));
    for(int i=0; i<2; i++)
        red[i]=(unsigned int *)malloc(L_max*sizeof(unsigned int));

    L=0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            if(Random()<p){
                red[0][L]=i;
                red[1][L]=j;
                L++;
            }
        }
    }
    printf("Numero de links: %u", L);

    FILE *f_out=fopen("red_erdos_renyi(aleatoria).csv", "w");
    for(int i=0; i<L; i++){
        fprintf(f_out, "%u\t%u\n", red[0][i], red[1][i]);
    }
    fclose(f_out);
}
