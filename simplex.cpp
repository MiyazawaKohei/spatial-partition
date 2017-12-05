
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
	const int _m=3;
	const int _n=2;
	double _A[_m][_n]=	{{1,1},
						 {2,0},
						 {2,1}};
	double _B[_m]={4,5,6.5};
	double _C[_n]={-3,-2};

	int n=_m+_n;
	for(int i=1; i<=_m; i++){
		if(_B[i-1]<0) n++;
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
		if(_B[i-1]>=0){
			A[i][0]=_B[i-1];
			for(int j=1; j<=_n; j++){
				A[i][j]=_A[i-1][j-1];
			}		
			A[i][_n+i]=1;
			base[i]=_n+i;
		}else{
			A[i][0]=-_B[i-1];
			for(int j=1; j<=_n; j++){
				A[i][j]=-_A[i-1][j-1];
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
	show_matrix(&A,&base);
	simplex(&A,&base);
	show_matrix(&A,&base);
	if(A[0][0]>EPS){
		cout<<"No Feasible answer is there!!"<<endl;
	}
	for(int j=1; j<=_n; j++){
		A[0][j]=-_C[j-1];
	}
	A[0][0]=0;
	simplex(&A,&base);
	show_matrix(&A,&base);
	
	return 0;
}