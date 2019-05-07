#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>
#include <omp.h>

#define epsilon 1.e-8
#define num 16

using namespace std;
double C[100];

template <typename T> double sgn(T val)
{
	return (val > T(0)) - (val < T(0));
}

int main (int argc, char* argv[]){

	int M,N;
	string T,P,Db;
	cin>>M>>N;

	double elapsedTime,elapsedTime2;
	timeval start,end,end2;

	double U[N][N],V[N][N], S[N],U_t[N][N], V_t[N][N], A[N][N],alpha, beta, gamma, c, zeta, t,s,sub_zeta, converge;
	int I1[N], I2[N];

	int acum = 0;
	int temp1, temp2;
	converge = 1.0;

	for(int i = 0; i < M; i++){
		for(int j =0; j < N; j++){
			cin >> A[i][j];
		}
	}

	for(int i=0; i<M;i++){
		for(int j=0; j<N;j++){

			if(i==j){
				V_t[i][j] = 1.0;
			}
			else{
				V_t[i][j] = 0.0;
			}
		}
	}

	for(int i=0; i<M;i++){
		for(int j=0; j<N;j++){

			U_t[i][j] = A[j][i];
		}
	}

	/* SVD using Jacobi algorithm (Parallel)*/

	gettimeofday(&start, NULL);

	double conv;
	while(converge > epsilon){ 		//convergence
		converge = 0.0;
		acum++;
		for (int l = 1; l < M; l ++) {

			int r1 = 0, r2 = 0;
			//I1,I2 are 2 sets whose intersection (for both i and j ) is null
			for (int i = 0; i + l < M; i++) {
				if (i % (2 * l) < l)
					I1[++r1] = i;
				else
					I2[++r2] = i;
			}
			for (int k = 0; k < num; k++) {
				C[k] = converge;
			}
			#pragma omp parallel for num_threads(num)
			for (int p = 1; p <= r1; p++){
				int k = omp_get_thread_num();
				int i = I1[p], j = i + l;
				double alpha = 0, beta = 0, gamma = 0;
				double zeta, t, c, s;
				for (int k = 0; k < N; k++) {
					alpha = alpha + (U_t[i][k] * U_t[i][k]);
					beta = beta + (U_t[j][k] * U_t[j][k]);
					gamma = gamma + (U_t[i][k] * U_t[j][k]);
				}
				C[k] = max(C[k], abs(gamma)/sqrt(alpha*beta));

				zeta = (beta - alpha) / (2.0 * gamma);
				t = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta*zeta)));        //compute tan of angle
				c = 1.0 / (sqrt (1.0 + (t*t)));				//extract cos
				s = c*t;							//extrac sin
				for(int k=0; k<N; k++){
					t = U_t[i][k];
					U_t[i][k] = c*t - s*U_t[j][k];
					U_t[j][k] = s*t + c*U_t[j][k];

					t = V_t[i][k];
					V_t[i][k] = c*t - s*V_t[j][k];
					V_t[j][k] = s*t + c*V_t[j][k];

				}
			}
			#pragma omp parallel for num_threads(num)
			for (int p = 1; p <= r2; p++){
				int k = omp_get_thread_num();
				int i = I2[p], j = i + l;
				double alpha = 0, beta = 0, gamma = 0;
				double zeta, t, c, s;
				for (int k = 0; k < N; k++) {
					alpha = alpha + (U_t[i][k] * U_t[i][k]);
					beta = beta + (U_t[j][k] * U_t[j][k]);
					gamma = gamma + (U_t[i][k] * U_t[j][k]);
				}
				C[k] = max(C[k], abs(gamma)/sqrt(alpha*beta));

				zeta = (beta - alpha) / (2.0 * gamma);
				t = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta*zeta)));        //compute tan of angle
				c = 1.0 / (sqrt (1.0 + (t*t)));				//extract cos
				s = c*t;							//extrac sin
				for(int k=0; k<N; k++){
					t = U_t[i][k];
					U_t[i][k] = c*t - s*U_t[j][k];
					U_t[j][k] = s*t + c*U_t[j][k];

					t = V_t[i][k];
					V_t[i][k] = c*t - s*V_t[j][k];
					V_t[j][k] = s*t + c*V_t[j][k];

				}
			}
			for (int k = 0; k < num; k++)
				converge = max(converge, C[k]);
		}
	}

	for(int i =0; i<M; i++){

		t=0;
		for(int j=0; j<N;j++){
			t=t + pow(U_t[i][j],2);
		}
		t = sqrt(t);

		for(int j=0; j<N;j++){
			U_t[i][j] = U_t[i][j] / t;
			if(i == j){
				S[i] = t;
			}
		}
	}

	gettimeofday(&end, NULL);

	for(int i =0; i<M; i++){

		for(int j =0; j<N; j++){

			U[i][j] = U_t[j][i];
			V[i][j] = V_t[j][i];

		}

	}


	cout<<"iterations: "<<acum<<endl;
	elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;
	elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;
	cout<<"Time: "<<elapsedTime<<" ms."<<endl<<endl;

	cout<<"U"<<endl<<endl;
	for(int i =0; i<M; i++){
		for(int j =0; j<N; j++){

			cout<<U[i][j]<<"  ";
		}
		cout<<endl;
	}

	cout<<endl<<"Vt"<<endl<<endl;
	for(int i =0; i<M; i++){
		for(int j =0; j<N; j++){

			cout<<V_t[i][j]<<"  ";
		}
		cout<<endl;
	}

	cout<<endl<<"S"<<endl<<endl;
	for(int i =0; i<M; i++){
		for(int j =0; j<N; j++){

			if(i==j){  cout<<S[i]<<"  ";}

			else{
				cout<<"0.0  ";
			}
		}
		cout<<endl;
	}

	//At the end we need only U*S which can be U right after the algorithm because it is already U*S

	return 0;
}
