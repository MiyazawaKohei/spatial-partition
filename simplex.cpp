
#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
using namespace std;
#define EPS 1E-9

void initialize_of_valuables(vector<vector<double>> *A, vector<int>* base){
	const int _m=3;
	const int _n=6;
	double _A[_m][_n]=	{{1,3,-1,0,2,0},
						 {0,-2,4,1,0,0},
						 {0,-4,3,0,8,1}};
	double _B[_m]={7,12,10};
	double _C[_n]={0,-1,3,0,-2,0};
	int _base[_m]={1,4,6};
	for(int i=1; i<=_m; i++){
		for(int j=1; j<=_n; j++){
			(*A)[i][j]=_A[i-1][j-1];
		}
	}
	for(int i=1; i<=_m; i++){
		(*A)[i][0]=_B[i-1];
		(*base)[i]=_base[i-1];
	}
	for(int j=1; j<=_n; j++){
		(*A)[0][j]=_C[j-1];
	}
	return;
}

void show_matrix(vector<vector<double>> *A){
	int m=A->size()-1;
	int n=(*A)[0].size()-1;
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
		if(output==1) {cout<<"End with right answer"<<endl; return;}
		if(output==-1) {cout<<"Not Bounded!!"<<endl; return;}
	}
}
int main(){
	int m=3;	//No of Constrains
	int n=6;	//No of Valuables
	vector<vector<double>> A(m+1, vector<double>(n+1));
	vector<int> base(m+1);
	initialize_of_valuables(&A,&base);
	show_matrix(&A);
	simplex(&A,&base);
	show_matrix(&A);
	return 0;
}