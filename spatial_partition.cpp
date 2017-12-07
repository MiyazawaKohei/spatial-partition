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
	int n_a=2;
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
	show_matrix(&_A);
	/*
	const int _m=5;
	const int _n=2;
	double __A[_m][_n+1]=	{{0,1,0},
						{-2.999999, 0, 1},
						 {-1, -1 ,1},
						 {0, 1 ,1},
						{3,0,-1}};
	vector<vector<double>> ___A(_m,vector<double>(_n+1));
	for(int i=0; i<_m; i++){
		for(int j=0; j<=_n; j++){
			___A[i][j]= __A[i][j];
		}
	}
	cout<<simplex(___A)<<endl;
	*/
	int _m=n_a/2;
	int _n=12;
	
	const double pi=3.14159265358979;
	vector<vector<double>> template_of_A(24,vector<double>(13));
	int devide=20;
	double h=pi/devide;
	double No=0;

	int count=0;
	int countn[101];
	int outputposi,outputnega;
	double center;
	for (int i=0; i<=100; i++) countn[i]=0;
	for(double alpha=0; alpha<2*pi-EPS ; alpha+=h){

		for(double beta=0; beta<pi-EPS; beta+=h){

			for(double gamma=0; gamma<2*pi-EPS; gamma+=h){
						cout<<count<<","<<alpha<<","<<beta<<","<<gamma<<endl;
				//cout<<gamma<<endl;
				make_matrix(alpha, beta, gamma ,h ,&template_of_A);
				//show_matrix(&template_of_A);
				//return 0;
				//continue;
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
						//cout<<m<<endl;
						countn[m]++;
						
						center=0;
						for(int i=1; i<=9; i++){
							center+=template_of_A[2*i-1][0]*_A[m-1][i];
						}
						//cout<<"center: "<<center<<endl;
						int No_of_active_constrains=0;
						for(int i=0; i<m-1; i++){
							if (p->path[i]!=0) No_of_active_constrains++;
						}
						cout<<"No_of_active_constrains"<< No_of_active_constrains<<"plane_label"<<m<<endl;
						vector<vector<double>> A(No_of_active_constrains+24+1,vector<double>(_n+1));
						for(int i=1; i<=24; i++){
							for(int j=0; j<template_of_A[0].size(); j++){
								A[No_of_active_constrains+i][j]= template_of_A[i-1][j];
							}
						}
						int ii=0;
						for(int i=1; i<=m-1; i++){
							cout<<p->path[i-1]<<",";
						}
						cout<<endl;
						for(int i=1; i<=m-1; i++){	
							if (p->path[i-1]==0) continue;
							ii++;

							for(int j=1; j<=_n; j++){
								A[ii][j]=p->path[i-1]*_A[i-1][j];
							}
							A[ii][0]=0;//p->path[i-1]*_B[i-1];
						}
						for(int j=1; j<=_n; j++){
							A[0][j]=-_A[m-1][j];	
						}
						A[0][0]=0;
						//show_matrix(&A);
						outputposi=simplex(A);

						if(outputposi==1){
							p->positive=new vnode;
							Q->tail->next=p->positive;
							Q->tail=p->positive;
							p->positive->plane_label=p->plane_label+1;
							p->positive->path=p->path;
							p->positive->path.push_back(-1);
							cout<<"positive"<<endl;
						}else { p->positive=NULL;}//p->positive=NULL; return 0;}
						
						for(int j=1; j<=_n; j++){
							A[0][j]=_A[m-1][j];	
						}

						A[0][0]=0;
						outputnega=simplex(A);
						if(outputnega==1){
							p->negative=new vnode;
							Q->tail->next=p->negative;
							Q->tail=p->negative;
							p->negative->plane_label=p->plane_label+1;
							p->negative->path=p->path;
							p->negative->path.push_back(1);
							cout<<"negative"<<endl;
						}else {p->negative=NULL;}//   return 0;}
						if(p->positive!=NULL && p->negative==NULL){
							p->positive->path.pop_back();
							p->positive->path.push_back(0);
							cout<<"only positive"<<endl;
						}
						if(p->negative!=NULL && p->positive==NULL){
							p->negative->path.pop_back();
							p->negative->path.push_back(0);
							cout<<"only negative"<<endl;
						}
					}else count++;
					if(p==Q->tail) break;
					Q->head=p->next;
					delete p;
					//for (int i=0; i<=50; i++) {cout<<i<<":"<<countn[i]<<",";}
					//cout<<endl;
					//return 0;
				}
				//return 0;
					for (int i=0; i<=50; i++) {cout<<i<<":"<<countn[i]<<",";}
					cout<<endl;
				//if(count>0 ){cout<<count<<endl; return 0;}
				return 0;
			}
		}
		//show_matrix(&A);
	}
	cout<<"count"<<count<<endl;	
		
} 