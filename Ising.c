/*MODELLO DI ISING*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


void create_configuration(int **matrix, int l); //funzione che implementa la creazione di una configurazione di spin s=+1

void single_spin_flip(int x, int y, int l, int *magnetization, int *energy, int **matrix, double temperature); //funzione che implementa il Single Spin Flip

double jackknife(double mu, double *array2, double *array4, double sum2, double sum4, int measures); //funzione che implementa il metodo Jackknife



int main(int argc, char *argv[]){
	
	int L, N, M, E;
	int MCS, T_steps, MCS_eq;
	int **SpinMatrix; //array dinamico
	int i, j, n, k, t, Ri, Rj; //indici
	int count, a; //variabili per Binder e Jackknife
	double T, T_max, DT, T_min;
	double m, m_sum, m_middle, m_mod_middle, m2_sum, m2_middle, m_mod_sum;
	double e, e_sum, e_middle, e2_sum, e2_middle;
	double var_m, sigma_m, var_e, sigma_e;
	double m2B_sum, m4B_sum, m2[800], m4[800], m2B_middle, m4B_middle, sum_m2, sum_m4, B, sigma_B; //variabili per Binder e Jackknife
	unsigned seed;
	char s[100];
	FILE *f;
	
	if(argc!=7){
		fprintf(stderr, "Usage: %s L T_min T_max DT MCS MCSeq\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	
	L=atoi(argv[1]); //L è la dimensione del reticolo (matrice quadrata LxL)
	T_min=atof(argv[2]);
	T_max=atof(argv[3]);
	DT=atof(argv[4]);
	MCS=atoi(argv[5]);
	MCS_eq=atoi(argv[6]);
	
	
	seed=time(NULL);
	srand(seed);
	N=L*L; //N  è il numero totale degli spin
	
	sprintf(s, "Ising_L%d_Tmin%g_Tmax%g_DT%g_passi%d.dat", L, T_min, T_max, DT, MCS);
	f=fopen(s, "w");
	fprintf(f, "# Modello di Ising con il Metodo Monte Carlo\n");
	fprintf(f, "# Parametri iniziali: L=%d Tmin=%g Tmax=%g DT=%g passi=%d passo di equilibrio=%d\n", L, T_min, T_max, DT, MCS, MCS_eq);
	fprintf(f, "# 1: T    2: <m>    3: <|m|>    4: <e>    5: sigma m    6: sigma e\n");
	
	/*creazione matrice e verifica*/
	SpinMatrix=(int **)malloc(L*sizeof(int *));
	for(i=0; i<L; i++){
		SpinMatrix[i]=(int *)malloc(L*sizeof(int));
	}
	if(SpinMatrix==NULL){
		printf("malloc di SpinMatrix fallita!\n");
		exit(EXIT_FAILURE);
	}
	
	T_steps=(T_max-T_min)/DT; //numero di iterazioni del ciclo sulle temperature
	
	for(t=1; t<=T_steps; t++){ //inizio ciclo sulle temperature
		
		T=T_min+t*DT; //incremento della temperatura ad ogni passo
		
		printf("Temperatura: %g\n", T); //Monitoraggio dell'esecuzione del codice (Temperatura)
		
		create_configuration(SpinMatrix, L);
		
		/*Calcolo della magnetizzazione e dell'energia, estensiva ed intensiva, iniziale (subito dopo aver creato la configurazione)*/
		M=0;
		E=0;
		for(i=0; i<L; i++){
			for(j=0; j<L; j++){
				M=M+SpinMatrix[i][j];
				Ri=(L+1+i)%L;
				Rj=(L+1+j)%L;
				E=E+(SpinMatrix[i][j]*(SpinMatrix[i][Rj]+SpinMatrix[Ri][j]));
			}
		}
		E=-E;
		m=(double)M/N;
		e=(double)E/N;
		
		m_sum=0.; //inizializzo a 0 la variabile in cui accumulo la somma delle m
		m_mod_sum=0.; //inizializzo a 0 la variabile in cui accumulo la somma dei |m|
		m2_sum=0.; //inizializzo a 0 la variabile in cui accumulo la somma delle m^2
		e_sum=0.; //inizializzo a 0 la variabile in cui accumulo la somma delle e
		e2_sum=0.; //inizializzo a 0 la variabile in cui accumulo la somma delle e^2
		
		count=0; //Binder+Jackknife
		a=0; //Binder+Jackknife
		sum_m2=0.; //Binder+Jackknife
		sum_m4=0.; //Binder+Jackknife
		m2B_sum=0.; //Binder+Jackknife
		m4B_sum=0.; //Binder+Jackknife
		
		/*Metodo Monte Carlo*/
		for(k=1; k<=MCS; k++){ //inizio ciclo sui passi Monte Carlo
			
			for(i=0; i<L; i++){
				for(j=0; j<L; j++){
					single_spin_flip(i, j, L, &M, &E, SpinMatrix, T);
				}
			}
			
			if(k>MCS_eq){ //Quando arrivo al passo di equilibrio...
			
				m=(double)M/N; //Calcola m
				e=(double)E/N; //Calcola e
				
				m2B_sum=m2B_sum+(m*m); //Binder+Jackknife
				m4B_sum=m4B_sum+(m*m*m*m); //Binder+Jackknife
				count++; //Binder+Jackknife
				if(count==1000){ //Binder+Jackknife
					m2B_middle=m2B_sum/count;
					m4B_middle=m4B_sum/count;
					m2[a]=m2B_middle;
					m4[a]=m4B_middle;
					sum_m2=sum_m2+m2B_middle;
					sum_m4=sum_m4+m4B_middle;
					a++;
					m2B_sum=0.;
					m4B_sum=0.;
					count=0;
				}
				
				m_sum=m_sum+m; //...comincio ad accumulare le m
				m_mod_sum=m_mod_sum+fabs(m); //...comincio ad accumulare i |m|
				m2_sum=m2_sum+(m*m); //...comincio ad accumulare le m^2
				e_sum=e_sum+e; //...comincio ad accumulare le e
				e2_sum=e2_sum+(e*e); //...comincio ad accumulare le e^2
				
				if(k==MCS){ //Al termine del ciclo sui passi MC (per ogni T) calcolo...
					m_middle=m_sum/(MCS-MCS_eq); //...<m>
					m_mod_middle=m_mod_sum/(MCS-MCS_eq); //...<|m|>
					m2_middle=m2_sum/(MCS-MCS_eq); //...<m^2>
					var_m=m2_middle-(m_middle*m_middle); //...la varianza
					sigma_m=sqrt(var_m/((MCS-MCS_eq)-1)); //...l'errore
					
					e_middle=e_sum/(MCS-MCS_eq); //...<e>
					e2_middle=e2_sum/(MCS-MCS_eq); //...<e^2>
					var_e=e2_middle-(e_middle*e_middle); //...la varianza
					sigma_e=sqrt(var_e/((MCS-MCS_eq)-1)); //...l'errore
					
					B=0.5*(3-a*(sum_m4/(sum_m2*sum_m2))); //...il parametro di Binder
					sigma_B=jackknife(B, m2, m4, sum_m2, sum_m4, a); //...l'errore con il metodo Jackknife
				}
				
			}
			
			if(k!=0 && k%10000==0){ //Monitoraggio dell'esecuzione del codice (MCS)
				printf("MCS eseguiti: %d\n", k);
			}
			
		} //fine ciclo sui passi Monte Carlo
		
		fprintf(f, "     %.2lf    %lf    %lf    %lf    %lf    %lf\n", T, m_middle, m_mod_middle, e_middle, sigma_m, sigma_e, B, sigma_B);
		
	} //fine ciclo sulle temperature
	
	free(SpinMatrix);
	fclose(f);
	
	return EXIT_SUCCESS;
	
}



void create_configuration(int **matrix, int l){
	
	int i, j;
	
	for(i=0; i<l; i++){
		for(j=0; j<l; j++){
			matrix[i][j]=+1;
		}
	}
	
}


void single_spin_flip(int x, int y, int l, int *magnetization, int *energy, int **matrix, double temperature){
	
	int flipped, deltaE;
	int up, down, left, right;
	
	up=(x-1+l)%l; //condizione periodica al bordo superiore
	down=(x+1+l)%l; //condizione periodica al bordo inferiore
	left=(y-1+l)%l; //condizione periodica al bordo sinistro
	right=(y+1+l)%l; //condizione periodica al bordo destro
	
	flipped=-matrix[x][y]; //assegno lo spin (flippato) dell'elemento (x,y) alla variabile di appoggio flipped per semplicità di scrittura
	deltaE=2*matrix[x][y]*(matrix[up][y]+matrix[x][right]+matrix[down][y]+matrix[x][left]); //Calcolo del deltaE
	
	if(deltaE<=0 || (rand()/(RAND_MAX+1.))<=exp(-((double)deltaE/temperature))){
		matrix[x][y]=flipped; //...matrix[x][y]=-matrix[x][y]
		*magnetization=*magnetization+(2*matrix[x][y]); //aggiornamento magnetizzazione: M=M+2s
		*energy=*energy+deltaE; //aggiornamento energia: E=E+deltaE
	}
	
}


double jackknife(double mu, double *array2, double *array4, double sum2, double sum4, int measures){
	
	int i;
	double sum2_tmp, sum4_tmp, mu_j, mu_sum, sigma;
	
	mu_sum=0.;
	for(i=0; i<measures; i++){
		sum2_tmp=sum2-array2[i];
		sum4_tmp=sum4-array4[i];
		mu_j=0.5*(3-(measures-1)*(sum4_tmp/(sum2_tmp*sum2_tmp)));
		mu_sum=mu_sum+((mu_j-mu)*(mu_j-mu));
	}
	sigma=sqrt(((measures-1.)/measures)*mu_sum);
	return sigma;
}
