#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
#include "simplex.h"
#include "orthogonal_approximate.h"
#include "Eigen/SVD"

#include "InitAndFind.h"
using namespace std;
#define EPS 1E-9

struct vnode{
	int plane_label;
	vector<int> path;
	vnode *positive,*negative;
	vnode *next;
};
struct partition_tree{
	int degree;
	vnode *head;
};
struct queue_of_vnode{
	vnode *head;
	vnode *tail;
};

int main(){
	int n_a=60;
	int n_b=n_a;
	vector<vector<double>> _A(n_a/2,vector<double>(13));
	double SDofnoise=1;
	vector<double*> A,B;
	B=Init_Vectors_inB(n_b);
	cout<<"B:"<<endl;
	for(int i=0; i<n_b;i++)		cout<<B[i][0]<<","<<B[i][1]<<","<<B[i][2]<<endl;
	vector<vector<double>> rotation=Create_Random_Rotation_Matrix();
	vector<int> permutation=Create_Random_Permutaiton(n_b);
	
	A=Init_Vectors_inA(B,rotation,permutation,n_a,SDofnoise);
	A=Init_Vectors_inB(n_b);
	make_initial_matrix(A,B,&_A);
	int _m=n_a/2;
	int _n=12;
	//show_matrix(&_A);
	//return 0;
	/*double _A[_m][_n]=	{{1,1,2,4,6,3,2,6,2,2,2,4},
						 {2,0,-1,3,6,1,3,7,2,2,5,3},
						 {2,5,2,7,2,4,7,2,4,7,9,1}};
	double _B[_m]={4,5,7};
	*/
	
	const double pi=3.14159265358979;
	vector<vector<double>> template_of_A(24,vector<double>(13));
	int devide=1000000000;
	double h=pi/devide;
	double alpha,beta,gamma;
	double No=0;

	int count=0;
	int countn[40];
	for (int i=0; i<=39; i++) countn[i]=0;
	for(double alpha=0; alpha<2*pi-EPS ; alpha+=h){

		for(double beta=0; beta<pi-EPS; beta+=h){

			for(double gamma=0; gamma<2*pi-EPS; gamma+=h){
						cout<<count<<","<<alpha<<","<<beta<<","<<gamma<<endl;
				//cout<<gamma<<endl;
				make_matrix(alpha, beta, gamma ,h ,&template_of_A);
				partition_tree *T;
				T= new partition_tree;
				T->head= new vnode;
				T->head->plane_label=1;

				queue_of_vnode *Q;
				Q= new queue_of_vnode;
				Q->head=T->head;
				Q->tail=T->head;
				count=0;
				while(1){
					vnode *p=Q->head;
					if(p->plane_label<=_m){
						int m=p->plane_label;
						cout<<m<<endl;
						countn[m]++;
						vector<vector<double>> A(m+1+24,vector<double>(_n+1));
						for(int i=1; i<=24; i++){
							for(int j=0; j<template_of_A[0].size(); j++){
								A[m+i][j]= template_of_A[i-1][j];
							}
						}
						//vector<vector<double>> A(m+1,vector<double>(_n+1));
						
						for(int i=1; i<=m-1; i++){	
							for(int j=1; j<=_n; j++){
								A[i][j]=p->path[i-1]*_A[i-1][j-1];
							}
							A[i][0]=0;//p->path[i-1]*_B[i-1];
						}
						for(int j=1; j<=_n; j++){
							A[m][j]=_A[m-1][j-1];	
						}
						A[m][0]=0.00001;//_B[m-1];
						/*
						cout<<"A:"<<endl;
						for(int i=0; i<=m; i++){
							for(int j=0; j<=_n; j++){
								cout<<A[i][j];
								cout<<",";
							}
							cout<<endl;
						}*/
						//show_matrix(&A);
						if(simplex(A)==1){
							p->positive=new vnode;
							Q->tail->next=p->positive;
							Q->tail=p->positive;
							p->positive->plane_label=p->plane_label+1;
							p->positive->path=p->path;
							p->positive->path.push_back(1);
						}else p->positive=NULL;
						
						for(int j=1; j<=_n; j++){
							A[m][j]=-_A[m-1][j-1];	
						}

						A[m][0]=-0.00001;//-_B[m-1];
						//show_matrix(&A);
						if(simplex(A)==1){
							p->negative=new vnode;
							Q->tail->next=p->negative;
							Q->tail=p->negative;
							p->negative->plane_label=p->plane_label+1;
							p->negative->path=p->path;
							p->negative->path.push_back(-1);
						}else p->negative=NULL;
					}else count++;
					if(p==Q->tail) break;
					Q->head=p->next;
					delete p;
				}
					for (int i=0; i<=39; i++) {cout<<i<<","<<countn[i]<<endl;}
				if(count>0 ){cout<<count<<endl; return 0;}
			}
		}
		//show_matrix(&A);
	}
	cout<<"count"<<count<<endl;	

	
	return 0;
}