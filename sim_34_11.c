#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "task_34_11.h"

extern int debug_mode;

int sim_check(double *A, int n, double precision);
int sim_buildU(int k, double *A, double *X_k, double *X, double *U, int n, double precision);
void sim_build_newA(double *U, double *A, double *tmpA, int n);
void sim_multiply_matrix(double *m1, double *m2, double *res, int n);

int sim_34_11(int n, double* A, double* tmp, double precision){
/*вход: 	- n -- размерность матрицы;
					- A -- массив с матрицей системы в формате a_1_1 a_1_2 ... a_1_n a_2_1 a_2_2 ... a_n_n;
					- tmp -- массив дополнительной памяти, определяемый функцией int sim_memsize_34_11(int n);
					- precision --  определяет числа меньше какого считать нулем;
  выход:	- 0 - работа завершена успешно, матрица упрощена
  				- 1 - метод упрощения не применим к данной матрице*/
	int i,j,k,p;
	double *U, *X_k, *X ,*tmpA;
	U = (double*)(tmp);
	X_k = (double*)(tmp+(n*n));
	X = (double*)(tmp+(n*n)+n);
	tmpA = (double*)(tmp+(n*n)+(n*n)+n);

	p = sim_check(A, n, precision);					//проверка на симметричность
	if (p == -1) return -1;
	if (n <= 2) return 0;

	for(i = 0; i<n*n; i++)									// X - матрица из нулей n*n	
		X[i] = 0;

	if(debug_mode)
		printf("Building tridiagonal matrix:\n");
	for(k = 0; k < n-2; k++){	
		for(i = 0; i<n; i++){									//X_k=(0,...0) n-мерный вектор
			X_k[i]=0;
		}
		for(i = 0; i<n; i++){									// U=I, I -- единичная матрица n*n
			for(j = 0; j<n; j++){
				if(i == j) U[i*n+j]=1;
				else U[i*n+j]=0;
			}
		}
		p = sim_buildU(k, A, X_k, X, U, n, precision);						//build U matrix
		if (p == -1) continue;
		
		sim_build_newA(U,A,tmpA, n);							//build A_i+1=U A_i U
	}
	if(debug_mode){
		printf("New A matrix:\n");
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				printf("%1.9lf ", A[i*n+j]);
			}
			printf("\n");
		}
	}
	return 0;
}

int sim_check(double *A, int n, double precision){
	/*вход: 	- A -- массив с матрицей системы в формате a_1_1 a_1_2 ... a_1_n a_2_1 a_2_2 ... a_n_n;
						- n -- размерность матрицы;
						- precision --  определяет числа меньше какого считать нулем;
		выход:	- 0 - матрица симметрична;
  					 -1 - матрица несимметрична;*/
	int i, j;
	if(debug_mode)
		printf("The symmetry and positive definiteness check:\n");
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(fabs(A[i*n+j]-A[j*n+i])>precision){
				if (debug_mode)
					printf("Matrix is not symmetrical\n");
				return -1;
			}
		}
	}
	if(debug_mode)
		printf("Matrix is symmetrical\n");
	return 0;
}

int sim_buildU(int k, double *A, double *X_k, double *X, double *U, int n, double precision){
	/* вход: 	- k -- номер шага, k=0,..., n-3;
						- A -- матрица n*n;
						- X_k -- n-мерный вектор;
						- X -- матрица n*n, X=2*X_k*((X_k)^T);
						- U -- матрица отражения n*n;
						- n -- размерность матрицы;
						- precision --  определяет числа меньше какого считать нулем.
		 выход:	- 0 - матрица отражения построена;
  					-1 - норма a_k равна 0, переход к следующему шагу k;*/
	int i, j;
	double s_k, norm_a, norm_x;
	s_k = 0;																
	for(j = k+2; j <n; j++){																		//s_k = |a_{jk}|^2 j=k+2,...n-1
		s_k+=fabs(A[j*n+k])*fabs(A[j*n+k]);  
	}
	norm_a = sqrt(fabs(A[(k+1)*n+k])*fabs(A[(k+1)*n+k])+s_k);  //norm ||a_k||=sqrt(|a_{k+1,k}|^2+s_k)
	if(fabs(norm_a) <= precision) return -1;
 
	X_k[k+1] = A[(k+1)*n+k]-norm_a;														//X_k=(0, ..., 0, a_{k+1,k}-||a_k||, a_{k+2, k}, ...., a_{n,k})
	for(i = k+2; i<n; i++){
		X_k[i] = A[i*n+k];
	}
	norm_x = sqrt(fabs(X_k[k+1])*fabs(X_k[k+1])+s_k);		//norm ||X_k||=sqrt(|X_k[k]|^2+s_k) 
	if(fabs(norm_x) <= precision) return -1;
	for(i = k+1; i<n; i++){																		// X_k:=X_k/||X_K||;
		X_k[i] = X_k[i]/norm_x;
	}
	for(i = k+1; i<n; i++){																		// X=2*X_k*((X_k)^T)
		for(j = k+1; j<n; j++){
			X[i*n+j]=2*X_k[i]*X_k[j];
		}
	}
	for(i = k+1; i<n; i++){																		// U=I-X
		for(j = k+1; j<n; j++){
			U[i*n+j]=U[i*n+j]-X[i*n+j];
		}
	}
	return 0;
}
void sim_build_newA(double *U, double *A, double *tmpA, int n){
	/*вход: - U -- матрица отражения n*n;
					- A -- матрица n*n;
					- tmpA -- матрица n*n для результата перемножения матриц;
					- n -- размерность матрицы.*/
	sim_multiply_matrix(U, A, tmpA, n);
	sim_multiply_matrix(tmpA, U, A, n);
}
void sim_multiply_matrix(double *m1, double *m2, double *res, int n){
	/*вход: 	- m1 -- матрица n*n, первый множитель;
						- m2 -- матрица n*n, второй множитель;
						-res -- матрица n*n, результат перемножения матриц;
						- n -- размерность матриц.*/
	int i, j, k;
	double sum;
	for(i = 0; i<n; i++){
		for(j = 0; j<n; j++){
			sum = 0;
			for(k = 0; k<n; k++){
				sum += m1[i*n+k]*m2[k*n+j];
			}
			res[i*n+j] = sum;
		}
	}
}