#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>

#include "simplex.h"
#include "orthogonal_approximate.h"
using namespace std;
#define ESP 1E-9
#define pi 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594
void show_rotation_matrix(vector<vector<double>> *R){
	cout<<"R:"<<endl;
	for(int i=0; i<3; i++){
		for( int j=0; j<3; j++){
			cout<< (*R)[i][j]<<",";
		}
		cout<<endl;
	}
}

void rotation_around_one_axis(int i,int j,double theta, vector<vector<double>> *R){
	(*R)[i][i]=cos(theta);
	(*R)[j][j]=cos(theta);
	(*R)[3-i-j][3-i-j]=1;
	(*R)[i][j]=-sin(theta);
	(*R)[j][i]=sin(theta);
	return;
}

void product(vector<vector<double>> *A,vector<vector<double>>* B, vector<vector<double>>* AB){
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			(*AB)[i][j]=0;
			for(int k=0; k<3; k++){
				(*AB)[i][j]+=(*A)[i][k]*(*B)[k][j];
			}
		}
	}
	return;
}

void three_product(vector<vector<double>> *A, vector<vector<double>>* B, vector<vector<double>> *C, vector<vector<double>>* ABC){
	vector<vector<double>> D(3,vector<double>(3));
	product(A,B,&D);
	product(&D,C,ABC);
	return;	
}

void three_product_of_rotation(double alpha, double beta, double gamma,vector<vector<double>> *R ){
	vector<vector<double>> A(3,vector<double>(3));
	vector<vector<double>> B(3,vector<double>(3));
	vector<vector<double>> C(3,vector<double>(3));
	rotation_around_one_axis(0,1,alpha,&A);
	rotation_around_one_axis(0,2,beta ,&B);
	rotation_around_one_axis(1,2,gamma,&C);
	three_product(&A,&B,&C,R);
}

void make_matrix(double alpha, double beta, double gamma ,double h,vector<vector<double>> *A){
	vector<vector<double>> R(3,vector<double>(3));

	three_product_of_rotation(alpha,beta,gamma,&R);
	show_rotation_matrix(&R);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			(*A)[2*(i*3+j)][0]=R[i][j]+9*h*h/4;
			(*A)[2*(i*3+j)+1][0]=-R[i][j]+9*h*h/4;
		}
	}
	three_product_of_rotation(alpha+pi/2,beta,gamma,&R);
	show_rotation_matrix(&R);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			(*A)[2*(i*3+j)][10]=-R[i][j];
			(*A)[2*(i*3+j)+1][10]=R[i][j];
		}
	}
	three_product_of_rotation(alpha,beta+pi/2,gamma,&R);
	show_rotation_matrix(&R);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			(*A)[2*(i*3+j)][11]=-R[i][j];
			(*A)[2*(i*3+j)+1][11]=R[i][j];
		}
	}
	three_product_of_rotation(alpha,beta,gamma+pi/2,&R);
	show_rotation_matrix(&R);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			(*A)[2*(i*3+j)][12]=-R[i][j];
			(*A)[2*(i*3+j)+1][12]=R[i][j];
		}
	}
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			(*A)[2*(i*3+j)][i*3+j+1]=1;
			(*A)[2*(i*3+j)+1][i*3+j+1]=-1;
		}
	}
	(*A)[18][10]=1; (*A)[18][0]=h/2;
	(*A)[19][10]=-1; (*A)[19][0]=+h/2;
	(*A)[20][11]=1; (*A)[20][0]=h/2;
	(*A)[21][11]=-1; (*A)[21][0]=h/2;
	(*A)[22][12]=1; (*A)[22][0]=h/2;
	(*A)[23][12]=-1; (*A)[23][0]=h/2;
}
/*
int main(){
	const double pi=3.14159265358979;
	vector<vector<double>> A(24,vector<double>(13));
	int devide=20;
	double h=pi/devide;
	double alpha,beta,gamma;
	double No=0;
	for(double alpha=0; alpha<2*pi-ESP; alpha+=h){
		cout<<alpha<<endl;
		for(double beta=0; beta<pi-ESP; beta+=h){

			for(double gamma=0; gamma<2*pi-ESP; gamma+=h){
				//cout<<gamma<<endl;
				make_matrix(alpha, beta, gamma ,h ,&A);
				No++;
			}
		}
		//show_matrix(&A);
	}
}*/