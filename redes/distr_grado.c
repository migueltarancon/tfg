#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//poner los parámetros de la red para representar las distribuciones teoricas
#define N 10000
#define p 0.02

#define m0 50
#define m 4
#define T 9950

double distr_teorica_er(int grad);
double distr_teorica_ba(int grad);

main(int argc, char *argv[]){
    int grado, min=-1, max;
    char dummy[100];
    FILE *fichero=fopen("grados_red.csv", "r"); //fichero con los grados de cada nodo de la red
    do{
        fscanf(fichero, "%d", &grado);
        if(min==-1){
            min=grado;
            max=grado;
        }
        if(min>grado)
            min=grado;
        if(max<grado)
            max=grado;
    }while(!feof(fichero));
    rewind(fichero);
    int *distr_grad=(int *)calloc(max-min+1,sizeof(int));
    do{
        fscanf(fichero, "%d", &grado);
        distr_grad[grado-min]++;
    }while(!feof(fichero));
    distr_grad[grado-min]--;
    FILE *salida=fopen("distr_grado_red.dat", "w"); //fichero donde se imprimen los resultados
    int sum=0;
    double integral=0;
    for(int i=0; i<(max-min+1); i++){
        sum+=distr_grad[i];
        integral+=distr_teorica_ba(i+min);  //especificar la distribucion teorica
    }

    for(int i=0; i<(max-min+1); i++){
        fprintf(salida, "%d\t%d\t%e\n", i+min, distr_grad[i], sum/integral*distr_teorica_ba(i+min));    //especificar la distribucion teorica
    }
    fclose(fichero);
    fclose(salida);
}

double distr_teorica_er(int grad){
    double aux=1, c;
    c=(N-1)*p;
    if(grad>c){
        for(int i=1; i<=floor(c); i++){
            aux*=c/(exp(1)*i);
        }
        aux*=exp(-(c-floor(c)));
        for(int i=floor(c)+1; i<=grad; i++){
            aux*=c/i;
        }
    }
    else if(grad<c){
        for(int i=1; i<=grad; i++){
            aux*=c/(exp(1)*i);
        }
        for(int i=grad+1; i<=floor(c); i++){
            aux*=exp(-1);
        }
        aux*=exp(-(c-floor(c)));
    }

    return aux;
}

double distr_teorica_ba(int grad){
    double aux;
    aux=2.0*m0*(m0+1);
    aux/=grad;
    aux/=grad+1;
    aux/=grad+2;
    return aux;
}
