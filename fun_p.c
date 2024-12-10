/*
    AC - OpenMP -- SERIE
    fun_s.c
     rutinas que se utilizan en el modulo gengrupos_s.c 
****************************************************************************/
#include <math.h>
#include <float.h> // DBL_MAX
#include <stdlib.h>
#include "defineg.h"           // definiciones
#include <stdio.h>
#include <omp.h>
/**************************************************************************************
   1 - Funcion para calcular la distancia genetica entre dos elementos (distancia euclidea)
       Entrada:  2 elementos con NCAR caracteristicas (por referencia)
       Salida:  distancia (double)
**************************************************************************************/
double gendist (float *elem1, float *elem2)
{
	// PARA COMPLETAR
	// calcular la distancia euclidea entre dos vectores
	float sum = 0.0;
	int i;
	for(i=0;i<NCAR;i++){
		sum = sum + pow((elem1[i]-elem2[i]),2);
	}
	return sqrt(sum);
}

/****************************************************************************************
   2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)
   Entrada:  nelem  numero de elementos, int
             elem   elementos, una matriz de tamanno MAXE x NCAR, por referencia
             cent   centroides, una matriz de tamanno NGRUPOS x NCAR, por referencia
   Salida:   popul  grupo mas cercano a cada elemento, vector de tamanno MAXE, por referencia
*****************************************************************************************/
void grupo_cercano (int nelem, float elem[][NCAR], float cent[][NCAR], int *popul)
{
	// PARA COMPLETAR
	// popul: grupo mas cercano a cada elemento
	int i, j, centroide_cercano;
	float dist_min, aux_dist, centroide_ind = -1;
	#pragma omp parallel private(i,j,aux_dist,centroide_ind,dist_min) shared(popul)
    {
    #pragma omp for schedule(static)
	for(i=0;i<nelem;i++){
		dist_min = FLT_MAX;
		for(j=0;j<ngrupos;j++){
			aux_dist = gendist(elem[i],cent[j]);
			if(dist_min > aux_dist) {
				dist_min = aux_dist;
				centroide_ind = j;
			}
		}
		popul[i] = centroide_ind;
	}
	}
}

/****************************************************************************************
   3 - Funcion para calcular la calidad de la particion de clusteres.
       Ratio entre a y b. El termino a corresponde a la distancia intra-cluster.
       El termino b corresponde a la distancia inter-cluster.
   Entrada:  elem     elementos, una matriz de tamanno MAXE x NCAR, por referencia
             listag   vector de NGRUPOS structs (informacion de grupos generados), por ref.
             cent     centroides, una matriz de tamanno NGRUPOS x NCAR, por referencia
   Salida:   valor del CVI (double): calidad/ bondad de la particion de clusters
*****************************************************************************************/
double silhouette_simple(float elem[][NCAR], struct lista_grupos *listag, float cent[][NCAR], float *a){
    float b[ngrupos];
    float s[ngrupos];
    // PARA COMPLETAR
	// aproximar a[i] de cada cluster: calcular la densidad de los grupos
    //		media de las distancia entre todos los elementos del grupo
    //   	si el numero de elementos del grupo es 0 o 1, densidad = 0
	int i, j, k, elemJ, elemK;
	double sumA, sumB, S=0;
    #pragma omp parallel private(i,j,k,sumA,elemJ,elemK) shared(a)
    {
    #pragma omp for schedule(static)
	for(i=0;i<ngrupos;i++){
		sumA = 0;
		if(listag[i].nelemg<=1){
			a[i] = 0;
		}else{
			for(j=0;j<listag[i].nelemg;j++){
				elemJ = listag[i].elemg[j];
				for(k=0;k<listag[i].nelemg;k++){
					elemK = listag[i].elemg[k];
					if(j!=k) {
                        sumA += gendist(elem[elemJ],elem[elemK]);
                    }
				}
			}
			a[i] = sumA / pow((listag[i].nelemg),2.0);
		}
	}
    }
    // aproximar b[i] de cada cluster
    #pragma omp parallel private(i,j,sumB) shared(b)
    {
    #pragma omp for schedule(static)
	for(i=0;i<ngrupos;i++){
		sumB = 0;
		for (j=0;j<ngrupos;j++){
			if(i!=j) sumB += gendist(cent[i],cent[j]);
		}
		b[i] = sumB / (ngrupos-1);
	}
    }
	// calcular el ratio s[i] de cada cluster
	for(i=0;i<ngrupos;i++){
		s[i] = (b[i]-a[i])/ fmax(a[i],b[i]);
		S += s[i];
	 }
	// promedio y devolver
    return S/ngrupos;
}

/********************************************************************************************
   4 - Funcion para relizar el analisis de enfermedades
   Entrada:  listag   vector de NGRUPOS structs (informacion de grupos generados), por ref.
             enf      enfermedades, una matriz de tamaño MAXE x TENF, por referencia
   Salida:   prob_enf vector de TENF structs (informacion del análisis realizado), por ref.
*****************************************************************************************/
int compare_floats(const void *a, const void *b){
	float fa = *(const float*)a;
	float fb = *(const float*)b;
	return (fa>fb)-(fa<fb);
}

void analisis_enfermedades (struct lista_grupos *listag, float enf[][TENF], struct analisis *prob_enf)
{
	// PARA COMPLETAR
	// Realizar el análisis de enfermedades en los grupos:
	//		mediana máxima y el grupo en el que se da este máximo (para cada enfermedad)
	//		mediana mínima y su grupo en el que se da este mínimo (para cada enfermedad)
	int i, j, k, l;
	float medianaC;
	int indMedianaC;

	//incializar todos los valores de las medianas
	for(l=0;l<TENF;l++){
		prob_enf[l].mmax = -1;
		prob_enf[l].gmax = 1000;
		prob_enf[l].mmin = FLT_MAX;
		prob_enf[l].gmin = 1000;
	}
	
	// Recorrer cada cluster
	#pragma omp parallel private(i,j,k) shared (prob_enf)
	{
	struct analisis prob_enfL[TENF];
	for (i = 0; i < TENF; i++) {
        prob_enfL[i].mmax = -1;
        prob_enfL[i].gmax = 0;
        prob_enfL[i].mmin = 100;
        prob_enfL[i].gmin = 0;
	}
	#pragma omp for schedule (static)
	for (j = 0; j < ngrupos; j++) {
		// comprobar si el tamaño del cluster[j] > 0
		if (listag[j].nelemg > 0) {
			// crear un array del tamaño del cluster[j] para la probabilidad de las enferemdades de cada persona
			float enf_per[listag[j].nelemg];
			// Recorrer la matriz de enfermedades por columnas
			for (i = 0; i < TENF; i++) {
				// Guardar todas las personas en el array para la enfermedad[i]
				for (k = 0; k < listag[j].nelemg; k++) {	
					enf_per[k] = enf[listag[j].elemg[k]][i];
				}
				// Ordenar el array de menor a mayor
				qsort(enf_per, listag[j].nelemg, sizeof(float), compare_floats);
				// Conseguir la mediana de la enfermedad[i]
				indMedianaC = floor(listag[j].nelemg / 2);
				medianaC = enf_per[indMedianaC];
				if (medianaC > prob_enfL[i].mmax) {
                    prob_enfL[i].mmax = medianaC;
                    prob_enfL[i].gmax = j;
                }
                if (medianaC < prob_enfL[i].mmin) {
    				prob_enfL[i].mmin = medianaC;
    				prob_enfL[i].gmin = j;
				}
			}
			
		}
	}
	#pragma omp critical
	{
	for(i=0;i<TENF;i++){
		if (prob_enfL[i].mmax > prob_enf[i].mmax || ( prob_enfL[i].mmax >= prob_enf[i].mmax && prob_enfL[i].gmax < prob_enf[i].gmax)) {
            prob_enf[i].mmax = prob_enfL[i].mmax;
            prob_enf[i].gmax = prob_enfL[i].gmax;
        }
        if (prob_enfL[i].mmin < prob_enf[i].mmin || ( prob_enfL[i].mmin <= prob_enf[i].mmin && prob_enfL[i].gmin < prob_enf[i].gmin)) {
    		prob_enf[i].mmin = prob_enfL[i].mmin;
    		prob_enf[i].gmin = prob_enfL[i].gmin;
		}

	}
	}
	}
}



/***************************************************************************************************
   OTRAS FUNCIONES DE LA APLICACION
****************************************************************************************************/

void inicializar_centroides(float cent[][NCAR]){
	int i, j;
	srand (147);
	for (i=0; i<ngrupos; i++)
		for (j=0; j<NCAR/2; j++){
			cent[i][j] = (rand() % 10000) / 100.0;
			cent[i][j+(NCAR/2)] = cent[i][j];
		}
}

int nuevos_centroides(float elem[][NCAR], float cent[][NCAR], int popul[], int nelem){
	int i, j, fin;
	double discent;
	double additions[ngrupos][NCAR+1];
	float newcent[ngrupos][NCAR];

	#pragma omp parallel for private (i,j) schedule (static)
	for (i=0; i<ngrupos; i++)
		for (j=0; j<NCAR+1; j++)
			additions[i][j] = 0.0;

	// acumular los valores de cada caracteristica (100); numero de elementos al final
	#pragma omp parallel for private (i,j) reduction (+:additions[:ngrupos][:NCAR+1]) schedule (static)
	for (i=0; i<nelem; i++){
		for (j=0; j<NCAR; j++) additions[popul[i]][j] += elem[i][j];
		additions[popul[i]][NCAR]++;
	}

	// calcular los nuevos centroides y decidir si el proceso ha finalizado o no (en funcion de DELTA)
	fin = 1;
	#pragma omp parallel private (i,j,discent) reduction (&&:fin)
	{ 
	#pragma omp for schedule (static)
	for (i=0; i<ngrupos; i++){
		if (additions[i][NCAR] > 0) { // ese grupo (cluster) no esta vacio
			// media de cada caracteristica
			for (j=0; j<NCAR; j++)
				newcent[i][j] = (float)(additions[i][j] / additions[i][NCAR]);

			// decidir si el proceso ha finalizado
			discent = gendist (&newcent[i][0], &cent[i][0]);
			if (discent > DELTA1)
				fin = 0;  // en alguna centroide hay cambios; continuar

			// copiar los nuevos centroides
			for (j=0; j<NCAR; j++)
				cent[i][j] = newcent[i][j];
		}
	}
	}
	return fin;
}

