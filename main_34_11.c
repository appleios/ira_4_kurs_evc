#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "task_34_11.h"
#pragma warning(disable:4996)

int debug_mode = 0, 
	print_errors = 0;

size_t sim_memsize_34_11(int n){
	return n*n*sizeof(double)+ //U matrix
		n*sizeof(double)+        //x^k vector 
		n*n*sizeof(double)+      //X matrix
		n*n*sizeof(double);      //matrix for matrix multiplication tmpA
}

size_t evc_memsize_34_11(int n){
	return n*n*sizeof(double)+   //Q matrix
				 n*n*sizeof(double)+   //U matrix
				 2*sizeof(double)+		 //X_k vector
				 4*sizeof(double)+     //X matrix
				 n*n*sizeof(double)+   //matrix for matrix multiplication tmpA
				 n*sizeof(double)+     //matrix for elements of non decided blocks
				 n*sizeof(double);     //matrix for E comparison
}

//на вход подаются 2 строки, возвращает номер первого несовпадения или -1, в случае равенства строк
int compare_strings(const char * str1,const char * str2){
	int str1_len, str2_len,i, max_len;
	str1_len=strlen(str1);
	str2_len=strlen(str2);
	if(str1_len==str2_len)
		max_len=str1_len;
	else{
		if(str1_len>=str2_len)
			max_len=str1_len;
		else max_len=str2_len;
	}
	for(i=0;i<max_len;i++){
		if(str1[i]!=str2[i])
			return i;
	}
	return -1;
}

int main(int argc, char *argv[]){
	int n, ans, i, j, print_execution_time=0, max_iter=0;
	double *A, *E, *tmp1, *tmp2, diftime, precision=1e-14, epsilon=1e-10;
	time_t time1, time2;
	FILE *in, *out;
	int print_matrix = 0;
	char *in_file_name = "34_11_in.txt", 
		*out_file_name = "34_11_out.txt";

	for(i=1; i<argc; i++){
		if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){
			printf(
			"Usage: evc[input_file_name] [output_file_name] [options]\n"
			"Where options include:\n"
			"-d    print debug messages [default OFF]\n"
			"-e    print errors [default OFF]\n"
			"-p    print matrix [default OFF]\n"
			"-t    print execution time [default OFF]\n"
			"-prec=<num>       precision [default - 1e-14]\n"
			"-eps=<num>        'epsilon' [default - 1e-10]\n"
			"-max_iter=<num>   limit number of iterations [default - 0, i.e. not limit]\n"
			"-h, -?     print this and exit\n"
			"Default input_file_name value is 34_11_in.txt, default output_file_name value is 34_11_out.txt.\n");
			return 0;
		}else if(strcmp(argv[i],"-d")==0){
			debug_mode = 1;
		}else if(strcmp(argv[i],"-e")==0){
			print_errors = 1;
		}else if(strcmp(argv[i],"-t")==0){
			print_execution_time = 1;
		}else if(strcmp(argv[i],"-p")==0){
			print_matrix = 1;
		}else if(compare_strings(argv[i], "-prec=")>5){
			precision=atof(&argv[i][6]);										
		}else if(compare_strings(argv[i], "-eps=")>4){
			epsilon=atof(&argv[i][5]);
		}
		else if(compare_strings(argv[i], "-max_iter=")>9){
			max_iter=atoi(&argv[i][10]);
			if(max_iter<0) max_iter=0;
		}	
		else{
			if(argv[i][0]=='-'){
				return 1;
			}else{
				if(i==1){
					in_file_name = argv[i];
				}else if(i==2){
					out_file_name = argv[i];
				}
			}
		}
	}

	in = fopen(in_file_name, "r");
	out = fopen(out_file_name, "w");
	if(in == NULL){
		if(print_errors)
			printf("Impossible to open input file\n");
		return 1;
	}
	if(out == NULL){
		if(print_errors)
			printf("Impossible to open output file\n");
		return 1;
	}		
	fscanf(in, "%d", &n);
	if(n<1){
		if(print_errors)
			printf("Incorrect matrix size");
		return 1;
	}
	A = (double*)malloc(n*n*sizeof(double));
	E = (double*)malloc(n*sizeof(double));

	for(i=0; i<n*n; i++)
		fscanf(in, "%lf", &A[i]);

	if(print_matrix){ //print out matrix
		printf("matrix A:\n");
		for(i=0; i<n; i++){
			for(j=0; j<n; j++){
				printf("%1.9lf ",A[i*n+j]);
			}
			printf("\n");
		}		
		printf("\n");
	}
	
	time1 = time(NULL);
	tmp1 = (double*)malloc(sim_memsize_34_11(n)); 
	ans = sim_34_11(n,A,tmp1,precision);

	if(ans==0){
		tmp2 = (double*)malloc(evc_memsize_34_11(n)); 
		ans = evc_34_11(n,max_iter,epsilon,A,E,tmp2,precision);
	}
	time2 = time(NULL);

	if(print_execution_time){
		diftime = difftime(time1,time2);
		if(debug_mode){
			printf("\nExecution time: %1.9lf\n", difftime);
		}
	}

	if ((ans == 0)||(ans == 1)){
		fprintf(out,"%d\n", n);
		for (i=0; i<n; i++)
			fprintf(out, "%1.9lf\n", E[i]);
		if(ans == 0){
			if(debug_mode)
				printf("\nWork has been completed successfully");
			}
			else 	if((print_errors)||(debug_mode)){
				printf("\nMethod does not converge for the specified number of iterations");}
	}
	else{
		if ((print_errors)||(debug_mode))
			printf("\nMethod of finding the eigenvalues is not applicable to the given matrix");
		fprintf(out, "0");
	}
	free(A);
	free(E);
	fclose(out);
	fclose(in);
	free(tmp1);
	free(tmp2);
	return 0;
}