
#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
#include "simplex.h"

using namespace std;
#define EPS 1E-9


void show_matrix(vector<vector<double>> *A,vector<int> *base){
	int m=A->size()-1;
	int n=(*A)[0].size()-1;
	cout<<"A:"<<endl;
	for(int i=0; i<=m; i++){
		for(int j=0; j<=n; j++){
			cout<<(*A)[i][j];
			cout<<",";
		}
		cout<<endl;
	}
	cout<<"base:"<<endl;
	for(int i=1; i<=m; i++){
		cout<<(*base)[i]<<endl;
	}
}
void show_matrix(vector<vector<double>> *A){
	int m=A->size()-1;
	int n=(*A)[0].size()-1;
	cout<<"A:"<<endl;
	for(int i=0; i<=m; i++){
		for(int j=0; j<=n; j++){
			cout<<(*A)[i][j];
			cout<<",";
		}
		cout<<endl;
	}
}
int pivot(vector<vector<double>> *A, vector<int> *base){
	int m=A->size()-1;
	int n=(*A)[0].size()-1;

	double max_of_coefficient_of_objective=0;
	int max_iteration_for_coefficient_of_objective=0;
	for(int j=1; j<=n; j++){
		if((*A)[0][j]>max_of_coefficient_of_objective){
			max_of_coefficient_of_objective=(*A)[0][j];
			max_iteration_for_coefficient_of_objective=j;
		}
	}
	if(max_iteration_for_coefficient_of_objective==0) return 1;	//End with optimal solution.

	double max_of_coefficient_of_constrains=0;
	int max_iteration_for_coefficient_of_constrains=0;
	for(int i=1; i<=m; i++){
		if((*A)[i][max_iteration_for_coefficient_of_objective]>max_of_coefficient_of_constrains){
			max_of_coefficient_of_constrains=(*A)[i][max_iteration_for_coefficient_of_objective];
			max_iteration_for_coefficient_of_constrains=i;
		}
	}
	if(max_iteration_for_coefficient_of_constrains==0) return -1; //Not Bounded!!

	for(int j=0; j<=n; j++){
		(*A)[max_iteration_for_coefficient_of_constrains][j]/=max_of_coefficient_of_constrains;
	}
	for(int i=0; i<=m; i++){
		if(i==max_iteration_for_coefficient_of_constrains) continue;
		for(int j=0; j<=n; j++){
			if(j==max_iteration_for_coefficient_of_objective) continue;
			(*A)[i][j]-=(*A)[max_iteration_for_coefficient_of_constrains][j]*(*A)[i][max_iteration_for_coefficient_of_objective];
		}
		(*A)[i][max_iteration_for_coefficient_of_objective]=0;
	}

	(*base)[max_iteration_for_coefficient_of_constrains]=max_iteration_for_coefficient_of_objective;
	return 0;
}

void simplex(vector<vector<double>> *A, vector<int> *base){
	int output;
	while(1){
		output=pivot(A, base);
		if(output==1) return;
		if(output==-1) return;
	}
}

int simplex(vector<vector<double>> _A){

	int _m=_A.size()-1;
	int _n=_A[0].size()-1;
	int n=_m+_n;
	for(int i=1; i<=_m; i++){
		if(_A[i][0]<0) n++;
	}
	vector<vector<double>> A(_m+1, vector<double>(n+1));
	vector<int> base(_m+1);
	for(int i=0; i<=_m; i++){
		for(int j=0; j<=n; j++){
			A[i][j]=0;
		}
	}
	int itr_of_artificial_valuable=_n+_m;
	for(int i=1; i<=_m; i++){
		if(_A[i][0]>=0){
			A[i][0]=_A[i][0];
			for(int j=1; j<=_n; j++){
				A[i][j]=_A[i][j];
			}		
			A[i][_n+i]=1;
			base[i]=_n+i;
		}else{
			A[i][0]=-_A[i][0];
			for(int j=1; j<=_n; j++){
				A[i][j]=-_A[i][j];
			}
			A[i][_n+i]=-1;
			itr_of_artificial_valuable++;
			A[i][itr_of_artificial_valuable]=1;
			base[i]=itr_of_artificial_valuable;
		}
	}
	for(int i=1; i<=_m; i++){
		if(base[i]>_n+_m){
			for(int j=1; j<=_n; j++){
				A[0][j]+=A[i][j];
			}
			A[0][_n+i]=-1;
			A[0][0]+=A[i][0];
		}
	}
	//show_matrix(&A,&base);
	simplex(&A,&base);

	//show_matrix(&A,&base);
	if(A[0][0]>EPS){
		return -1;
	}
	return 1;
}