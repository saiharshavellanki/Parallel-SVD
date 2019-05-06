#include <iostream>
#include<bits/stdc++.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define epsilon 1.e-10

using namespace std;

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

	double U[M][N],V[M][N], S[M],U_t[M][N], V_t[M][N], A[M][N];
	double alpha, beta, gamma, c, zeta, t,s,sub_zeta, converge,prev_converge;

	int acum = 0;
	int temp1, temp2;
	converge = 1.0;
	prev_converge = 100;

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

	/* SVD using Jacobi algorithm */

	gettimeofday(&start, NULL);

	double conv;
	int fl=0;
	while(converge > epsilon && prev_converge!=converge){
		prev_converge  = converge;
		converge = 0.0;

		acum++;

		for(int i = 1; i<M; i++){
			for(int j = 0; j<i; j++){

				alpha = 0.0;
				beta = 0.0;
				gamma = 0.0;

				for(int k = 0; k<N ; k++){
					alpha = alpha + (U_t[i][k] * U_t[i][k]);
					beta = beta + (U_t[j][k] * U_t[j][k]);
					gamma = gamma + (U_t[i][k] * U_t[j][k]);
				}

				converge = max(converge, abs(gamma)/sqrt(alpha*beta));	//compute convergence
				//basicaly is the angle
				//between column i and j
				if(gamma==0)
				{
					break;
				}
				zeta = (beta - alpha) / (2.0 * gamma);
				t = sgn(zeta) / (abs(zeta) + sqrt(1.0 + (zeta*zeta)));
				c = 1.0 / (sqrt (1.0 + (t*t)));
				s = c*t;
				for(int k=0; k<N; k++){
					t = U_t[i][k];
					U_t[i][k] = c*t - s*U_t[j][k];
					U_t[j][k] = s*t + c*U_t[j][k];

					t = V_t[i][k];
					V_t[i][k] = c*t - s*V_t[j][k];
					V_t[j][k] = s*t + c*V_t[j][k];

				}

			}
		}
//		cout<<prev_converge<<" "<<converge<<endl;
	}


	for(int i =0; i<M; i++){

		t=0;
		for(int j=0; j<N;j++){
			t=t + pow(U_t[i][j],2);
		}
		t = sqrt(t);
		if(t!=0)
		{
			for(int j=0; j<N;j++){
				U_t[i][j] = U_t[i][j] / t;
				if(i == j){
					S[i] = t;
				}
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

	cout<<endl<<"Vt"<<endl<<endl;
	for(int i =0; i<M; i++){
		for(int j =0; j<N; j++){

			cout<<V_t[i][j]<<"  ";
		}
		cout<<endl;
	}

	return 0;

}
