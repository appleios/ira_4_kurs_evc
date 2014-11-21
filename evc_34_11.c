#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "task_34_11.h"

extern int debug_mode;

double evc_find_normA(double *A, int n, double precision);
void evc_exhaustion(double *A, int n, double normA, double epsilon, double *blocks, double *E);
int evc_find_num_s_k(double *blocks, int n);
int evc_buildU(int k, double *A, double *X_k, double *X, double *U, int n, double precision);
void evc_multiply_matrix(double *m1, double *m2, double *res, int n);
void evc_sort(double *E, int n, double precision);


int evc_34_11(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision){
	/* вход:	- n -- размерность матрицы;
						- max_iterations -- ограничение на число итераций алгоритма;
						- epsilon -- точность;
						- A -- массив с трехдиагональной матрицей системы в формате a_1_1 a_1_2 ... a_1_n a_2_1 a_2_2 ... a_n_n;
						- E -- n- мерный вектор собственных значений A;
						- tmp -- массив дополнительной памяти, определяемый функцией int sim_memsize_34_11(int n);
						- precision --  определяет числа меньше какого считать нулем;
		выход:	- 0 - работа завершена успешно;
						- 1 - метод не сходится за указанное число итераций.*/
	int i, j, k, p, p1, num_it = 1, n1 = n-1;
	double *Q, *U, *X, *X_k, *tmpA, *blocks, normA, s_k, d, *oldE;
	Q = (double*)(tmp);
	U = (double*)(tmp+(n*n));
	X_k = (double*)(tmp+(n*n)+(n*n));
	X = (double*)(tmp+(n*n)+(n*n)+2);
	tmpA = (double*)(tmp+(n*n)+(n*n)+6);
	blocks = (double*)(tmp+(n*n)+(n*n)+6+(n*n));
	oldE = (double*)(tmp+(n*n)+(n*n)+6+(n*n)+n);

	if(n == 1) {E[0]=A[0]; return 0;}															//Если размерность матрицы 1,2
	if (n == 2){
		d = A[0]*A[0]+A[3]*A[3]-2*A[0]*A[3]+4*A[1]*A[2];
		E[1] = (double)(A[0]+A[3]+sqrt(d))/2;
		E[0] = (double)(A[0]+A[3]-sqrt(d))/2;
		evc_sort(E,n, precision);
		return 0;
	}
	
	for( i =0 ; i<n ; i++) {blocks[i] = 1;}

	normA = evc_find_normA(A, n, precision);											//||A||=max(sum|a_{i,j}|) i,j=0,..,n-1

	if(debug_mode)
		printf("Finding eigenvalues:\n");

	while((num_it<=max_iterations)||(!max_iterations)){
		if(debug_mode)
			printf("Iteration %d\n", num_it);

		for(i = 0; i< n; i++){																			//Q=I
			if(blocks[i]!=0)
				oldE[i] = A[i*n+i];
			for(j = 0; j<n; j++){
				if(i == j) Q[i*n+j]=1;
				else Q[i*n+j]=0;
			}
		} 
	/*	printf("old A matrix:\n");
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				printf("%1.9lf ", A[i*n+j]);
			}
			printf("\n");
		}*/
		s_k = A[n1*n+n1];
		for(i = 0; i < n; i++){																		 //Сдвиг матрицы A: A-s_kI
			for(j = 0; j < n; j++){
				if(i == j)		A[i*n+j]-=s_k;
			}
		}
		for(k = 0; k < n-1;k++){																	// Разложение матрицы A-s_kI=QR
			for(i = 0; i<n; i++){																		// U=I
				for(j = 0; j<n; j++){
					if(i == j) U[i*n+j]=1;
					else U[i*n+j]=0;
				}
			}
			p	= evc_buildU(k, A, X_k, X, U, n, precision);					//build U matrix
			if (p == -1) continue;
			/*printf("U:\n");
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				printf("%1.9lf ", U[i*n+j]);
			}
			printf("\n");
		}*/
			evc_multiply_matrix(U,A,tmpA, n);												//build R=U*A_k
			for(i = 0; i<n; i++){
				for(j = 0; j<n; j++){
					A[i*n+j]=tmpA[i*n+j];
				}
			}
			evc_multiply_matrix(Q,U, tmpA, n);											// build Q=Q*U
			p = 0;
			for(i = 0; i<n; i++){
				for(j = 0; j<n; j++){
					if(fabs(Q[i*n+j]-tmpA[i*n+j])>precision)
						p = 1;
					Q[i*n+j]=tmpA[i*n+j];
				}
			}
		}
		evc_multiply_matrix(A,Q, tmpA, n);												// Вычисление новой матрицы A=RQ+s_k*I
		for(i = 0; i<n; i++){
			for(j = 0; j<n ; j++){
				if(i == j)
					A[i*n+j] = tmpA[i*n+j]+s_k;
				else A[i*n+j] = tmpA[i*n+j];
			}
		}
		/*printf("\n new A matrix:\n");
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				printf("%1.9lf ", A[i*n+j]);
			}
			printf("\n");
		}
		printf("\n Q matrix:\n");
		for(i = 0; i < n; i++){
			for(j = 0; j < n; j++){
				printf("%1.9lf ", Q[i*n+j]);
			}
			printf("\n");
		}*/
		p1 = 0;
		for(i = 0; i<n ;i++){				// if Q=I eigenvalues are found
			for(j = 0; j<n; j++){
				if(i == j){
					if(fabs(fabs(Q[i*n+j])-1)>precision) p1 = 1;
				}
				else{
					if(fabs(Q[i*n+j])>precision) p1 = 1;
				}
			}
		}
		if(p1 == 0){
			for(i = n-1; i>=0; i--){
					E[i] = A[i*n+i];
			}
			evc_sort(E,n, precision); return 0;
		}
		for(i = 0 ;i<n; i++){
			if(blocks[i]!=0){
				if(fabs(fabs(Q[i*n+i])-1)<precision){
					p = 0;
					for(j = 0; j<n; j++){
						if(i!=j){
							if(fabs(Q[i*n+j])>precision) p = 1;
							if(fabs(Q[j*n+i])>precision) p = 1;
						}
					}
					if(p == 0){
						E[i] = A[i*n+i];
						blocks[i] = 0;
					}
				}
			}
		}


		evc_exhaustion(A,n, normA, epsilon, blocks, E);						//Исчерпывание матрицы A, разбиение на блоки
		//for(i = 0; i<n; i++) printf("%1.lf ", blocks[i]);
		n1 = evc_find_num_s_k(blocks, n);													//Вычисление сдвига s_k, номер еще не вычисленного собств знач
		if (n1 == -1) { evc_sort(E,n, precision); return 0;}			//все собственные значения вычислены
		num_it++;
	}
	if(num_it == max_iterations){
		for(i = 0; i<n ; i++)
			E[i]=A[i*n+i];
	}
	evc_sort(E,n, precision);
	return 1;
}

double evc_find_normA(double *A, int n, double precision){
	/*вход: 	- A -- матрица n*n;
						- n -- размерность матрицы.
		выход:	- max -- норма матрицы A.*/
	int i, j;
	double max = 0,s;
	for(i = 0; i< n; i++){
		s = 0;
		for(j = 0 ; j< n; j++){
			s+= fabs(A[i*n+j]);
		}
		if (s-max > precision) max = s;
	}
	return max;
}
void evc_exhaustion(double *A, int n, double normA, double epsilon, double *blocks, double *E){
	/*вход:	- A -- матрица n*n;
	- n -- размерность матрицы;
	- normA -- норма матрицы A;
	- epsilon -- точность;
	- blocks -- текущая разметка матрицы;
	- E -- текущие собственные значения.*/
	int i, j,s, p = 0;
	double d;
	for(i = 1; i<n; i++){
		if(blocks[i] == 0){
			if(i!=0)
				A[(i-1)*n+i] = 0;
			A[i*n+i-1] = 0;
		}
	}
	for(i = 0; i<n-1; i++){					//Исчерпывание
		if((fabs(A[(i+1)*n+i])<epsilon*normA)&&(blocks[i]!=0)){
			A[(i+1)*n+i] = 0;
			A[i*n+i+1] = 0;
		}
	}
	
	//for(i = 0; i<n; i++) printf("%1.lf ", blocks[i]);
	s = 0;																													//Проверка наличия блоков 2*2
	for(i = 0; i<n; i++){
		if(blocks[i]!=0){
			s++;
			if ((A[(i+1)*n+i] == 0)||(i==n-1)) {
				if (s == 1){
					E[i] = A[i*n+i];
					blocks[i] = 0;
				}
				if (s == 2){ 
					d = A[(i-1)*n+(i-1)]*A[(i-1)*n+(i-1)]+A[i*n+i]*A[i*n+i]-2*A[(i-1)*n+(i-1)]*A[i*n+i]
							+4*A[(i-1)*n+i]*A[i*n+(i-1)];
					E[i-1] = (double)(A[(i-1)*n+(i-1)]+A[i*n+i]+sqrt(d))/2;
					E[i] = (double)(A[(i-1)*n+(i-1)]+A[i*n+i]-sqrt(d))/2;
					blocks[i-1] = 0;
					blocks[i] = 0;
				}
				s = 0;
			}
		}
	}
	if(( s == 1)||(s==2)){
		for(i = n-1;i>=0; i--) if(blocks[i]!=0) break;
		if (s == 1){
					E[i] = A[i*n+i];
					blocks[i] = 0;
				}
				if (s == 2){ 
					d = A[(i-1)*n+(i-1)]*A[(i-1)*n+(i-1)]+A[i*n+i]*A[i*n+i]-2*A[(i-1)*n+(i-1)]*A[i*n+i]
							+4*A[(i-1)*n+i]*A[i*n+(i-1)];
					E[i-1] = (double)(A[(i-1)*n+(i-1)]+A[i*n+i]+sqrt(d))/2;
					E[i] = (double)(A[(i-1)*n+(i-1)]+A[i*n+i]-sqrt(d))/2;
					blocks[i-1] = 0;
					blocks[i] = 0;
				}
	}
	for(i = n-1; i>0; i--){																					//проверка условия сдвига
		if((blocks[i]!=0)&&(fabs(A[i*n+(i-1)]) < epsilon*normA)){
			E[i] = A[i*n+i];
			blocks[i] = 0;
		}
		else if(blocks[i]!=0) break;
	}
	
}

int evc_find_num_s_k(double *blocks, int n){
/*	вход:	- A -- матрица n*n;
					- n -- размерность матрицы;
					- normA -- норма матрицы A;
					- epsilon -- точность
					- E -- вектор собственных значений матрицы;
					- n1 -- номер вычисляемого собственного значения.
		выход:  номер нового вычисляемого собственного значения*/
	int i;
	for(i = n-1; i >= 0; i--)
		if(blocks[i]!=0)	break;
	return i;

}
int evc_buildU(int k, double *A, double *X_k, double *X, double *U, int n, double precision){
	/* вход: 	- A -- массив с матрицей n*n;
						- n -- размерность матрицы;
						- k -- номер шага, k=0,..., n-1;
						- X_k -- 2-мерный вектор;
						- X -- матрица 2*2, X=2*X_k*((X_k)^T);
						- U -- матрица отражения n*n;
						- precision --  определяет числа меньше какого считать нулем.
			выход:- 0 - матрица отражения построена;
  					-1 - норма a_k равна 0, переход к следующему шагу k;*/
	int i, j;
	double  norm_a, norm_x;
	
	norm_a = sqrt(fabs(A[k*n+k])*fabs(A[k*n+k])+fabs(A[(k+1)*n+k])*	fabs(A[(k+1)*n+k]));  //norm ||a_k||=sqrt(|a_{k+1,k}|^2+|a_{k,k}|^2)
	if(fabs(norm_a) <= precision) return -1;
 
	X_k[0] = A[k*n+k]-norm_a;																															//X_k=(a_{k,k}-||a_k||, a_{k+1, k}, 0,...., 0)
	X_k[1] = A[(k+1)*n+k];

	norm_x = sqrt(fabs(X_k[0])*fabs(X_k[0])+fabs(A[(k+1)*n+k])*fabs(A[(k+1)*n+k]));				//norm ||X_k||=sqrt(|X_k[k]|^2+|a_{k+1,k}|^2) 
	if(fabs(norm_x) <= precision) return -1;

	X_k[0]=X_k[0]/norm_x;																																	// X_k:=X_k/||X_K||;
	X_k[1] = X_k[1]/norm_x;
	
	for(i = 0; i<2; i++){																																// X=2*X_k*((X_k)^T)
		for(j = 0; j<2; j++){
			X[i*n+j]=2*X_k[i]*X_k[j];
		}
	}
	for(i = k; i<k+2; i++){																														// U=I-X
		for(j = k; j<k+2; j++){
			U[i*n+j]=U[i*n+j]-X[(i-k)*n+(j-k)];
		}
	}
	return 0;
}

void evc_multiply_matrix(double *m1, double *m2, double *res, int n){
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

void evc_sort(double *E,int n, double precision){
	/*вход:	- n -- размерность матрицы;
					- E -- вектор собственных значений матрицы;
					- precision -- определяет числа меньше какого считать нулем.*/
	int i, j;
	double tmp;
	for(i = 0; i < n-1; i++){
		for(j = 0; j < n-i-1; j++){
			if(E[j] - E[j+1] > precision){
				tmp = E[j];
				E[j] = E[j+1];
				E[j+1] = tmp;
			}
		}
	}
}



