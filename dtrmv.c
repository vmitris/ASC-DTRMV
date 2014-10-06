#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<malloc.h>
#include<cblas.h>
#include<time.h>
#include<sys/time.h>

int LSAME(char a, char b){
	if(a == b)
		return 1;
	return 0;
}
/**
	implementarea functiei de mana
	este exact functia din fortram transpusa in C
*/
void dtrmv(char UPLO, char TRANS, char DIAG, int N, double **A, int LDA, double *X,  int INCX){
	
	double ZERO=0;
	double TEMP;
    int I, IX, J, JX, KX=0;
    int NOUNIT;
   
	if(N == 0 )
    	return;
 		
	NOUNIT = LSAME(DIAG,'N');	
 	
 	
 	if(INCX < 0)
          KX = 1 - (N-1)*INCX;
    else if(INCX != 1)
          KX = 0;	
 	
 	if(LSAME(TRANS,'N')){
 		if(LSAME(UPLO,'U')){
              if(INCX == 1){
                  for(J = 0; J <N; J++){ 
                      if(X[J] != ZERO){
                          TEMP = X[J];
                          for(I = 0; I <= J-1; I++){
                              X[I] = X[I] + TEMP*A[I][J];
                          	}
                          if(NOUNIT) 
                          	X[J] = X[J]*A[J][J];
                     	}
                   	}
                }
              else{
                  JX = KX;
                  for(J = 0;J < N; J++){
                      if(X[JX] != ZERO){
                          TEMP = X[JX];
                          IX = KX;
                          for(I = 0; I <= J - 1; I++){
                              X[IX] = X[IX] + TEMP*A[I][J];
                              IX = IX + INCX;
                            }
                          if(NOUNIT) 
                            X[JX] = X[JX]*A[J][J];
                        }
                      JX = JX + INCX;
                  }  
            	}
        }
        else{
            if (INCX == 1){
                  for(J = N-1; J>= 0; J--){
                      if(X[J] != ZERO){
                        	TEMP = X[J];
                       		for(I = N-1; I >= J + 1; I--)
                            	X[I] = X[I] + TEMP*A[I][J];
                            if(NOUNIT)
                            	X[J] = X[J]*A[J][J];
                        }
                  	}
              	}
            else{
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for(J = N-1; J>= 0; J--) {
                      if(X[JX] != ZERO){
                          TEMP = X[JX];
                          IX = KX;
                          for(I = N-1; I>= J + 1; I--){
                              X[IX] = X[IX] + TEMP*A[I][J];
                              IX = IX - INCX;
							             }
                          if(NOUNIT)
                          	X[JX] = X[JX]*A[J][J];
                      }
                      JX = JX - INCX;
		       		     }
            }
        }
    }

    else{

          if(LSAME(UPLO,'U')){
               if(INCX == 1){
                  for(J = N-1; J>= 0; J--){
                    TEMP = X[J];
                    if (NOUNIT) 
                    	TEMP = TEMP*A[J][J];
                    	for(I = J - 1; I>= 0; I--){
                       		TEMP = TEMP + A[I][J]*X[I];
                  		}
                    X[J] = TEMP;
  	            	}
              	}
              	else{
                  JX = KX + (N-1)*INCX;
                  for(J = N-1; J>= 0; J--){
                   	TEMP = X[JX];
                    IX = JX;
                    if(NOUNIT) 
                    	TEMP = TEMP*A[J][J];
                    for(I = J - 1; I >= 0; I--){
                        IX = IX - INCX;
                        TEMP = TEMP + A[I][J]*X[IX];
                		  }
                    X[JX] = TEMP;
                    JX = JX - INCX;
            		  }
              }
          }
          else{
              if(INCX == 1) {
                  for( J = 0; J<=N-1; J++){
                    	TEMP = X[J];
                    	if(NOUNIT) 
                      		TEMP = TEMP*A[J][J];
                      	for(I = J + 1; I<=N-1; I++){
                        	 TEMP = TEMP + A[I][J]*X[I];
 						             }
                    	X[J] = TEMP;
 					        }
              }
              else{
                	JX = KX;
                	for( J = 0 ; J<N; J++){
                    	TEMP = X[JX];
                    	IX = JX;
                    	if (NOUNIT) 
                      		TEMP = TEMP*A[J][J];
                     	for(I = J + 1; I<N; I++){
                        	IX = IX + INCX;
                        	TEMP = TEMP + A[I][J]*X[IX];
                        }
                      	X[JX] = TEMP;
                      	JX = JX + INCX;
		              }
 				       }
          }
    	}   
}

/**
	Implementarea functiei optimizate
	Pentru compilarea intregului executabil am folosit flag-ul -o5
	Matricea este aici declarata ca vector iar accesul la elemente se face prin pointeri
*/

void dtrmv2(char UPLO, char TRANS, char DIAG, int N, double *A, int LDA, double *X,  int INCX){
  
  double ZERO=0;
  double TEMP;
    int I, IX, J, JX, KX=0;
    int NOUNIT;
   
  if(N == 0 )
      return;
    
  NOUNIT = LSAME(DIAG,'N'); 
  
  
  if(INCX < 0)
          KX = 1 - (N-1)*INCX;
    else if(INCX != 1)
          KX = 0; 
  
  if(LSAME(TRANS,'N')){
    if(LSAME(UPLO,'U')){
              if(INCX == 1){
                  register double *x_j=X;
                  for(J = 0; J <N; J++){ 
                      if(X[J] != ZERO){
                          TEMP = X[J];
                          register double *reg=X;
                          register double *reg2=A+J;
                          for(I = 0; I <= J-1; I++){
                              *reg+=TEMP*(*reg2);
                              reg2+=N;
                              reg++;
                            }
                          if(NOUNIT){
                              *x_j +=A[J*N+J];
                              x_j++;
                            } 
                      }
                    }
                }
                else{
                  JX = KX;
                  for(J = 0;J < N; J++){
                      if(X[JX] != ZERO){
                          TEMP = X[JX];
                          IX = KX;
                          register double *reg=&X[IX];
                          register double *reg2=A+J;
                          for(I = 0; I <= J - 1; I++){
                              *reg += TEMP*(*reg2);
                              reg2+=N;
                              reg+=INCX;
                            }
                          if(NOUNIT) 
                            X[JX] = X[JX]*A[J*N+J];
                        }
                      JX = JX + INCX;
                    }    
              }
      }
      else{
            if (INCX == 1){
                  for(J = N-1; J>= 0; J--){
                      if(X[J] != ZERO){
                          TEMP = X[J];
                          register double *reg=X+N-1;
                          register double *reg2=A+(N-1)*N+J;
                          for(I = N-1; I >= J + 1; I--){
                              *reg+=TEMP*(*reg2);
                              reg2-=N;
                              reg--;
                          }
                            if(NOUNIT)
                              X[J] = X[J]*A[J*N+J];
                        }
                    }
                }
            else{
                  KX = KX + (N-1)*INCX;
                  JX = KX;
                  for(J = N-1; J>= 0; J--) {
                      if(X[JX] != ZERO){
                          TEMP = X[JX];
                          IX = KX;
                          register double *reg=&X[IX];
                          register double *reg2=A+(N-1)*N+J;
                          for(I = N-1; I>= J + 1; I--){
                              *reg += TEMP*(*reg2);
                              reg2-=N;
                              reg-=INCX;
                            }
                          if(NOUNIT)
                            X[JX] = X[JX]*A[J*N+J];
                      }
                      JX = JX - INCX;
                  }
            }
        }
    }

    //else trans
    else{

          if(LSAME(UPLO,'U')){
               if(INCX == 1){
                  for(J = N-1; J>= 0; J--){
                    TEMP = X[J];
                    if (NOUNIT) 
                      TEMP = TEMP*A[J*N+J];
                    for(I = J - 1; I>= 0; I--){
                          TEMP = TEMP + A[I*N+J]*X[I];
                      }
                    X[J] = TEMP;
                  }
                }
                else{
                  JX = KX + (N-1)*INCX;
                  for(J = N-1; J>= 0; J--){
                    TEMP = X[JX];
                    IX = JX;
                    if(NOUNIT) 
                      TEMP = TEMP*A[J*N+J];
                    for(I = J - 1; I >= 0; I--){
                        IX = IX - INCX;
                        TEMP = TEMP + A[I*N+J]*X[IX];
                      }
                    X[JX] = TEMP;
                    JX = JX - INCX;
                  }
              }
          }
          else{
              if(INCX == 1) {
                  for( J = 0; J<=N-1; J++){
                      TEMP = X[J];
                      if(NOUNIT) 
                          TEMP = TEMP*A[J*N+J];
                        for(I = J + 1; I<=N-1; I++){
                           TEMP = TEMP + A[I*N+J]*X[I];
                         }
                      X[J] = TEMP;
                  }
              }
              else{
                  JX = KX;
                  for( J = 0 ; J<N; J++){
                      TEMP = X[JX];
                      IX = JX;
                      if (NOUNIT) 
                          TEMP = TEMP*A[J*N+J];
                      for(I = J + 1; I<N; I++){
                          IX = IX + INCX;
                          TEMP = TEMP + A[I*N+J]*X[IX];
                        }
                        X[JX] = TEMP;
                        JX = JX + INCX;
                  }
               }
          }
      }


}




/*
	functie ce initializeaza matricea cu numere random de la 0 la 99
*/
void initializare_matrice(double **matrice, int n){
	int i,j;
  	srand(time(NULL));
	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			matrice[i][j] =(double)(rand()%100);
		}
}

/*
	functie ce initializeaza vectorii cu numere random de la 0 la 99
*/

void initializare_vector(double* vector, int n, int incx, double *vector2, double* vector3){
  int i;
  srand(time(NULL));
  double register temp;
  for(i = 0; i< n*incx; i++){
    temp = (double)(rand()%100);
    vector[i]=temp;
    vector2[i]=temp;
    vector3[i]=temp;
  }

}

/*
	functie pentru printare matrice, folosita pentru debugging
*/
void print_matrice(double** matrice,int n){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			printf("%.2f ",matrice[i][j]);
		}
		printf("\n");
	}
}

/*
	functie pentru afisare vectorului rezultat intr-un fisier
	folosita pentru verificarea corectitudinii rezultatelor
*/
void print_vector(double *vector, int n, char* fisier){
	int i;
  FILE *f;
  f=fopen(fisier,"w+");
	for(i = 0; i< n; i++){
		fprintf(f,"%f ",vector[i]);
	}
  fclose(f);

}

/*
	functie pentru liniarizarea matricii
*/

double* liniarizare(double **matrice, int n){
  double *matrice_rezultat;
  matrice_rezultat = (double*)malloc(n*n*sizeof(double));

  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      matrice_rezultat[i*n+j]=matrice[i][j];

  return matrice_rezultat;

}

/*
	programul primeste ca parametri un int reprezentand optiunea = {1,2,3} - pentru a sti in ce fisier 
	se afiseaza output-ul si un int reprezentand dimensiunea matricii;
*/
int main(int argc, char* argv[]){
	double **matrice;
	double *vector, *rezultat_blas,*rezultat_optim;
	int i,n,incx=1;
  	int lda;
  	double* matrice_liniarizata;
  	time_t start, end;
  	float timp;
  	int optiune;
  	FILE *f;

  	if(argc < 3){
    	printf("eroare la introducere\n");
 	}

 	n=atoi(argv[2]);
 	lda=n;
  	optiune=atoi(argv[1]);

  	//in functie de optiune se deschide fisierul de output

  	if(optiune == 1)
    	while(f=fopen("opteron.out","a+") = NULL)
    		f=fopen("opteron.out","a+");
  	else if(optiune == 2)
    	while(f=fopen("nehalem.out","a+") = NULL)
    		f=fopen("nehalem.out","a+");
  	else if(optiune == 3)
    	while(f=fopen("quad.out","a+") = NULL)
    		f=fopen("quad.out","a+");

    // initializari
  	fprintf(f,"%d ",n);
	matrice=(double**)malloc(n*sizeof(double*));
	for(i=0;i<n;i++)
		matrice[i] = (double *) malloc(n*sizeof( double));

	vector = (double*) malloc(n*incx*sizeof( double));
  	rezultat_blas = (double*) malloc(n*incx*sizeof( double));
  	rezultat_optim = (double*) malloc(n*incx*sizeof( double));
  	initializare_matrice(matrice,n);
  	initializare_vector(vector,n,incx,rezultat_blas,rezultat_optim);
  	//apelare functie implementare de mana
  	start=clock();
  	dtrmv('U','N','U',n,matrice,lda,vector,incx);
  	end=clock();
  	timp=(double)(end-start)/CLOCKS_PER_SEC;
  	if(n > 20000)
  		print_vector(vector,n*incx,"mana.out");
  	fprintf(f,"%f ",timp);
  	printf("timpul pentru rulare de mana este: %.15f s",timp);
	printf("\n=====================================\n");

	//apelare functie cblas
  	matrice_liniarizata=liniarizare(matrice,n);
  	start=clock();
  	cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasUnit, n, matrice_liniarizata, lda, rezultat_blas, incx);
  	end=clock();
  	timp=(double)(end-start)/CLOCKS_PER_SEC;
  	if(n > 20000)
  		print_vector(rezultat_blas,n*incx,"blas.out");
  	fprintf(f,"%f ",timp);
  	printf("timpul pentru rulare blas este: %.15f s",timp);
  	printf("\n=====================================\n");

  	//apelare functie optimizata
  	start=clock();
  	dtrmv2('U','N','U',n,matrice_liniarizata,lda,rezultat_optim,incx);
  	end=clock();
  	imp=(double)(end-start)/CLOCKS_PER_SEC;
  	if(n > 20000)
  		print_vector(rezultat_optim,n*incx, "optim.out");
  	fprintf(f,"%f \n",timp);
  	printf("timpul pentru rulare optim este: %.15f s",timp);
  	printf("\n=====================================\n");

  	fclose(f);
  	return 0;
}
