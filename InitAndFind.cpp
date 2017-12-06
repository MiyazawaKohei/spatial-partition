// comple with "-std=c++11"

#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include "Eigen/SVD"
#include <algorithm>
#include <functional>
#include "InitAndFind.h"
static const int dimention = 3;
using namespace std;

vector<double*> Init_Vectors_inB(int n_b){
	vector<double*> B;
	uniform_real_distribution<double> x(0,1.0);
	random_device rd;
	std::mt19937 mt(rd());
	for (int i=0; i<n_b; i++){
		B.push_back(new double [dimention]);
		for(int j=0; j<dimention; j++){
			B.back()[j]=x(mt);
		}
	}
	return B;
}

void Close_Vectors(vector<double*> B){
	while(B.size()>0){
		delete[] B.back();
		B.pop_back();
	}
	return;
}

vector<double> Create_Random_Normarized_Vector(){
	vector<double> Normarized_Vector;
	uniform_real_distribution<double> lambda(-1.0,1.0);
	random_device rd;
	std::mt19937 mt(rd());
	double theta=asin(lambda(mt));
	uniform_real_distribution<double> _phi(0,2*M_PI);
	double phi=_phi(mt);
	Normarized_Vector.push_back(cos(theta)*cos(phi));
	Normarized_Vector.push_back(cos(theta)*sin(phi));
	Normarized_Vector.push_back(sin(theta));
	return Normarized_Vector;
}

double L2_Norm(vector<double> v){
	double norm=0;
	for(int i=0; i<v.size(); i++)	norm+=v[i]*v[i];
	norm=pow(norm,0.5);

	return norm;
}

vector<double> Normarize_Vector(vector<double> v){
	double norm=L2_Norm(v);
	for(int i=0; i<v.size(); i++)	v[i]/=norm;
	return v;
}

double Inner_Product(vector<double> u,vector<double> v){
	double product=0;
	for(int i=0; i<u.size(); i++) product+=u[i]*v[i];
	return product;
}

vector<vector<double>> Create_Random_Rotation_Matrix(){
	vector<vector<double>> vectors;
	double innerproduct;
	for(int i=0; i<3; i++)		vectors.push_back(Create_Random_Normarized_Vector());
	cout<<"rotation:"<<endl;
	for(int i=0; i<3; i++){
		for(int j=0; j<i; j++){
			innerproduct=Inner_Product(vectors[i],vectors[j]);
			for(int k=0; k<vectors[i].size(); k++) vectors[i][k]-=vectors[j][k]*innerproduct;	
		}
		vectors[i]=Normarize_Vector(vectors[i]);
		cout<<vectors[i][0]<<vectors[i][1]<<vectors[i][2]<<endl;
	}	

	return vectors;
}

int MyRandom_Int(int i){
	uniform_int_distribution<int> rand_int(0,i-1);
	random_device rd;
	std::mt19937 mt(rd());
	return (rand_int(mt));
}

vector<int> Create_Random_Permutaiton(int n_b){
	vector<int> permutation;
	for(int i=0; i<n_b; i++)	permutation.push_back(i);
	for(int i=0; i<n_b/2; i++)	random_shuffle(&(permutation[2*i]),&(permutation[2*i+2]),MyRandom_Int);
	cout<<"permutation:"<<endl;
	for(int i=0; i<n_b; i++)	cout<<permutation[i]<<' ';
	cout<<endl;
	return permutation;
}

vector<double*> Init_Vectors_inA(vector<double*> B,vector<vector<double>> rotation, vector<int> permutation,int n_a, double SDofnoise){
	vector<double*> A;
	uniform_real_distribution<double> coordinate(0,1.0);
	normal_distribution<double> noise(0, SDofnoise);
	random_device rd;
	std::mt19937 mt(rd());
	double coordinates[dimention];
	for(int i=0; i<dimention; i++) coordinates[i]=coordinate(mt);
	cout<<"A:"<<endl;
	for (int i=0; i<n_a; i++){
		A.push_back(new double [dimention]);
		for(int j=0; j<dimention; j++){
			A.back()[j]=coordinates[j];
			if(SDofnoise>0) A.back()[j]+=noise(mt);			
			for(int k=0; k<dimention; k++){
				A.back()[j]+=rotation[k][j]*B[permutation[i]][k];
			}
			cout<<A.back()[j]<<',';
		}
		cout<<endl;
	}
	return A;
}

void Fix_Gravity_Point_to_Origin(vector<double*> vectors){
	vector<double*> newvectors;
	double gravity_point[dimention];
	for(int i=0; i<dimention; i++){
		gravity_point[i]=0;
		for(int j=0; j<vectors.size(); j++) gravity_point[i]+=vectors[j][i];
		gravity_point[i]/=vectors.size();
		for(int j=0; j<vectors.size(); j++) vectors[j][i]-=gravity_point[i];
	}
	return;
}

Eigen::Matrix<double,3,3> Compute_BAt(const vector<double*> A,const vector<double*> B){
	Eigen::Matrix<double,3,3> BAt;
	for(int i=0; i<dimention; i++){
		for(int j=0; j<dimention; j++){
			BAt(i,j)=0;
			for(int k=0; k<A.size(); k++)		BAt(i,j)+=A[k][j]*B[k][i];
		}
	}
	return BAt;
}

double Maximize_trRBAt(const vector<double*> A,const vector<double*> B){
	Eigen::Matrix<double,3,3> BAt=Compute_BAt(A,B);
	Eigen::JacobiSVD<Eigen::Matrix<double,3,3>> svd(BAt, Eigen::ComputeFullU | Eigen::ComputeFullV);
	double nuclear_norm=0;
	for(int i=0; i<3; i++)		nuclear_norm+=svd.singularValues()[i];
	return nuclear_norm;
}

double Minimize_RMSD(const vector<double*> A,const vector<double*> B, const vector<int> sigma_i){
	vector<double*> subA,subB;
	for (int i=0; i<sigma_i.size(); i++){
		subA.push_back(new double [dimention]);
		for(int j=0; j<dimention; j++)		subA.back()[j]=A[i][j];
		subB.push_back(new double [dimention]);
		for(int j=0; j<dimention; j++)		subB.back()[j]=B[sigma_i[i]][j];
	}
	Fix_Gravity_Point_to_Origin(subA);
	Fix_Gravity_Point_to_Origin(subB);
	
	double RMSD=0;
	for (int i=0; i<sigma_i.size(); i++){
		for(int j=0; j<dimention; j++) RMSD+=subA[i][j]*subA[i][j]+subB[i][j]*subB[i][j];
	}
	RMSD-=2*Maximize_trRBAt(subA,subB);
	Close_Vectors(subA);
	Close_Vectors(subB);

	return RMSD;
}

void make_initial_matrix(const vector<double*> A,const vector<double*> B,vector<vector<double>> *inequalities){
	Fix_Gravity_Point_to_Origin(A);
	Fix_Gravity_Point_to_Origin(B);
	for(int itr=0; itr<A.size()/2; itr++){
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				(*inequalities)[itr][i*3+j+1]=A[2*itr][i]*B[2*itr][j]+A[2*itr+1][i]*B[2*itr+1][j]-A[2*itr][i]*B[2*itr+1][j]-A[2*itr+1][i]*B[2*itr][j];
			}
		}
	}
}
/*
data* New_Data_for_HEAP(vector<int> taple,double key){
	data* p = new data;	
	p -> taple=taple;
	p -> key = key;
	return p;
}

HEAP* BUILD_HEAP_FOR_PAIRS(const vector<double*> A, const vector<double*> B){
	vector<int> sigma_i(2,0);
	double RMSD;
	HEAP *H = new HEAP;
	data **AofHEAP=new data*[B.size()*B.size()];
	//data *p;
	int subscriptofA=0;
	for(int i=0; i<B.size(); i++){
		for(int j=0; j<B.size(); j++){
			if(i==j) continue;
			subscriptofA++;
			sigma_i[0] = i;
			sigma_i[1] = j;
			RMSD = Minimize_RMSD(A,B,sigma_i);
			AofHEAP[subscriptofA]=New_Data_for_HEAP(sigma_i,RMSD);
		}
	}
	BUILD_HEAP(H,subscriptofA,AofHEAP,subscriptofA+1);
	return H;
}

vector<int> Complement_of_Vector(vector<int> oneset, int allsetsize){
	vector<int> flags(allsetsize+1,0);
	vector<int> complement;
	for(int i:oneset) flags[i]=1;
	for(int i=0; i<allsetsize; i++){
		if(flags[i]==0) complement.push_back(i);
	}
	return complement;
}
/*
void View_Contents_ofH(HEAP *H){
	cout<<"H:"<<endl;
	for (int i=1;i<=H->heap_size;i++){
		cout<<H->A[i]->key<<" [";
		for(int j=0; j<H->A[i]->taple.size(); j++){
			cout<<H->A[i]->taple[j]<<",";
		}
		cout<<"]"<<endl;
	}
}
/*
void Find_Minimizing_Permutation(const vector<double*> A, const vector<double*> B){
	HEAP* H=BUILD_HEAP_FOR_PAIRS(A,B);
	data* maxinheap;
	vector<int> complement;
	vector<int>* new_taple;
	double RMSD;
	data* p;
	int inte;
	while(1){
		//View_Contents_ofH(H);
		maxinheap = EXTRACT_MIN(H);
		cout<<maxinheap->taple.size()<<","<< maxinheap->key<<endl;
		if(maxinheap->taple.size()==A.size()){
			cout<<H->heap_size<<endl;
			cout<<maxinheap->key<<endl;
			//DESTRACT_HEAP(H);
			return;
		}
		complement= Complement_of_Vector(maxinheap->taple,B.size());
		for(int i:complement){
			new_taple=new vector<int>; 
			copy(maxinheap->taple.begin(), maxinheap->taple.end(), back_inserter(*new_taple) );
			(*new_taple).push_back(i);
			RMSD = Minimize_RMSD(A,B,*new_taple);
			p=New_Data_for_HEAP(*new_taple,RMSD);
			INSERT(H,p);
		}
	}
}
*//*
int main(){


	int n_a, n_b; //The number of the vectors of A, B
	int n=10;
	n_a=n;
	n_b=n;

	vector<vector<double>> inequalities(n/2,vector<double>(13));
	double SDofnoise=-1;
	vector<double*> A,B;
	B=Init_Vectors_inB(n_b);
	cout<<"B:"<<endl;
	for(int i=0; i<n_b;i++)		cout<<B[i][0]<<","<<B[i][1]<<","<<B[i][2]<<endl;
	vector<vector<double>> rotation=Create_Random_Rotation_Matrix();
	vector<int> permutation=Create_Random_Permutaiton(n_b);
	
	A=Init_Vectors_inA(B,rotation,permutation,n_a,SDofnoise);
	make_initial_matrix(A,B,&inequalities);
	//vector<int> permutationcut(n_a,0);
	//cout<<"Minimized_RMSD:"<<Minimize_RMSD(A,B,permutation)<<endl;
	//Find_Minimizing_Permutation(A,B);
	Close_Vectors(A);
	Close_Vectors(B);

}

*/