/*
Judge whether there is a point x in R^n such that

a_01x_1+...+a_0nx_n <  b_0
a_11x_1+...+a_1nx_n <= b_1
...
a_m1x_1+...+a_mnx_n <= b_m

Input should be in the form;
0    0     ...  0
b_1  a_11       a_1n
...
b_m  a_m1       a_mn


*/

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
#include "simplex.h"

using namespace std;
#define EPS 1E-9

void KENZAN(vector<vector<double>> *A, vector<int> *base,vector<vector<double>> *_A){
	int m=A->size()-1;
	int n=(*A)[0].size()-1;
	int _n=(*_A)[0].size()-1;
	vector<double> answer(n+1);
	for(int i=1; i<=m; i++){
		answer[(*base)[i]]=(*A)[i][0];
	}
	for(int i=1; i<=n; i++){
		cout<<i<<":"<<answer[i]<<endl;
	}
	double object=0;
	for(int i=1; i<=_n; i++){
		object+=(*_A)[0][i]*answer[i];
		object-=(*_A)[0][i]*answer[_n+i];	
	}
	cout<<"object: "<<object<<endl;
	cout<<"object-A[0][0]"<<object-(*A)[0][0]<<endl;
	//if(fabs(object-(*A)[0][0])>EPS) {cout<<"KENZAN"<<endl; cin>> m;}
	
}
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

		//cout<<"h"<<endl;
	double max_of_coefficient_of_objective=0;
	int max_iteration_for_coefficient_of_objective=0;

		//cout<<"02"<<endl;
	for(int j=1; j<=n; j++){
		//cout<<j<<":"<<j<<",";
		if((*A)[0][j]>max_of_coefficient_of_objective){
			//max_of_coefficient_of_objective=(*A)[0][j];
			max_iteration_for_coefficient_of_objective=j;
		}
	}

		//cout<<"32"<<endl;
	if(max_iteration_for_coefficient_of_objective==0) return 1;	//End with optimal solution.

	double min_of_ration_constant_per_coefficient=1000;
	int max_iteration_for_coefficient_of_constrains=0;
	int max_subscribe=0;
	for(int i=1; i<=m; i++){
		//cout<<"ef";
		if((*A)[i][max_iteration_for_coefficient_of_objective]>0 
			&& ( max_iteration_for_coefficient_of_constrains==0 || (*A)[i][0]/(*A)[i][max_iteration_for_coefficient_of_objective]<min_of_ration_constant_per_coefficient
			 || ( (*A)[i][0]/(*A)[i][max_iteration_for_coefficient_of_objective]==min_of_ration_constant_per_coefficient && max_subscribe<(*base)[i]))){
				//cout<<"3"<<endl;
			min_of_ration_constant_per_coefficient=(*A)[i][0]/(*A)[i][max_iteration_for_coefficient_of_objective];
			//cout<<"1"<<endl;
			max_iteration_for_coefficient_of_constrains=i;
			max_subscribe=(*base)[i];
		}
		//	cout<<"6";
	}
	if(max_iteration_for_coefficient_of_constrains==0) {show_matrix(A,base); cout<<"("<<max_iteration_for_coefficient_of_objective<<")"<<endl; return -1;} //Not Bounded!!
	for(int j=0; j<=n; j++){
		if(j== max_iteration_for_coefficient_of_objective) continue;
		(*A)[max_iteration_for_coefficient_of_constrains][j]/=(*A)[max_iteration_for_coefficient_of_constrains][max_iteration_for_coefficient_of_objective];
	}
	(*A)[max_iteration_for_coefficient_of_constrains][max_iteration_for_coefficient_of_objective]=1;

	for(int i=0; i<=m; i++){
		if(i==max_iteration_for_coefficient_of_constrains) continue;
		for(int j=0; j<=n; j++){
			if(j==max_iteration_for_coefficient_of_objective) continue;
			(*A)[i][j]-=(*A)[max_iteration_for_coefficient_of_constrains][j]*(*A)[i][max_iteration_for_coefficient_of_objective];
		}
		(*A)[i][max_iteration_for_coefficient_of_objective]=0;
	}
	//cout<<(*base)[max_iteration_for_coefficient_of_constrains]<<","<<max_iteration_for_coefficient_of_objective<<endl;
	(*base)[max_iteration_for_coefficient_of_constrains]=max_iteration_for_coefficient_of_objective;
	return 0;
}

int simplexinner(vector<vector<double>> *A, vector<int> *base){
	int output;
	int i=0;
	double prev;
	while(1){
		//cout<<++i<<","<< prev-(*A)[0][0]<<endl;
		prev=(*A)[0][0];
		output=pivot(A, base);
		//show_matrix(A,base);
		if(output==1) {cout<<"End with right solution"<<endl; return 1;}
		if(output==-1) {cout<<"Not Bounded!!"<< endl; return -1;}
	}
}

int simplexinner(vector<vector<double>> *A, vector<int> *base,vector<vector<double>> *_A){
	int output;
	int i=0;
	double prev;
	while(1){
		//cout<<++i<<","<< prev-(*A)[0][0]<<endl;
		KENZAN(A,base,_A);
		prev=(*A)[0][0];
		output=pivot(A, base);
		//show_matrix(A,base);
		if(output==1) {cout<<"End with right solution"<<endl; return 1;}
		if(output==-1) {cout<<"Not Bounded!!"<< endl; return -1;}
	}
}

int simplex(vector<vector<double>> _A){
	//cout<<"begin_of_simplex"<<endl;
	//show_matrix(&_A);
	int _m=_A.size()-1;
	int _n=_A[0].size()-1;
	int n=_m+2*_n;
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
	int itr_of_artificial_valuable=2*_n+_m;
	for(int i=1; i<=_m; i++){
		if(_A[i][0]>=0){
			A[i][0]=_A[i][0];
			for(int j=1; j<=_n; j++){
				A[i][j]=_A[i][j];
				A[i][_n+j]=-_A[i][j];
			}		
			A[i][2*_n+i]=1;
			base[i]=2*_n+i;
		}else{
			A[i][0]=-_A[i][0];
			for(int j=1; j<=_n; j++){
				A[i][j]=-_A[i][j];
				A[i][_n+j]=_A[i][j];
			}
			A[i][2*_n+i]=-1;
			itr_of_artificial_valuable++;
			A[i][itr_of_artificial_valuable]=1;
			base[i]=itr_of_artificial_valuable;
		}
	}
	for(int i=1; i<=_m; i++){
		if(base[i]>2*_n+_m){
			for(int j=1; j<=2*_n; j++){
				A[0][j]+=A[i][j];
			}
			A[0][2*_n+i]=-1;
			A[0][0]+=A[i][0];
		}
	}
	//show_matrix(&A,&base);
	//show_matrix(&A,&base);
	simplexinner(&A,&base,&_A);
	//show_matrix(&A,&base);
	if(A[0][0]>EPS) {cout<<"No fiseable"<<A[0][0]<<endl;cin>>_m;return -1;}
	//cout<<"f"<<endl;
	A[0][0]=-_A[0][0];
	for(int j=1; j<=_n; j++){
		A[0][j]=-_A[0][j];
		A[0][_n+j]=_A[0][j];
	}
	for(int j=2*_n+1; j<=n; j++){
		A[0][j]=0;
	}
	//cout<<"3456"<<endl;
	//show_matrix(&A,&base);
	double coefficient_for_base;
	for(int i=1; i<=_m; i++){
		coefficient_for_base=A[0][base[i]];
		//cout<<i<<","<<base[i]<<","<<coefficient_for_base<<endl;
		for(int j=0; j<=n; j++){
			A[0][j]-=coefficient_for_base*A[i][j];
		}
	}
	for(int i=0; i<=_m; i++){
		for(int j=2*_n+_m+1; j<=n; j++){
			A[i][j]=0;
		}
	}
	//	cout<<"end_of_simplex"<<endl;
	//show_matrix(&A,&base);
	//show_matrix(&A,&base);
	if(simplexinner(&A,&base,&_A)==-1) {/*show_matrix(&A,&base);*/ return 1;}
	//show_matrix(&A,&base);
	//show_matrix(&_A);
	cout<<"A[0][0]"<<A[0][0]<<endl;
	KENZAN(&A,&base,&_A);
	if(A[0][0]<-EPS){
		return 1;
	}	
	return -1;
}