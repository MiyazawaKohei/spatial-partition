#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
#include "simplex.h"
#include "orthogonal_approximate.h"
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
	const int _m=3;
	const int _n=12;
	double _A[_m][_n]=	{{1,1,2,4,6,3,2,6,2,2,2,4},
						 {2,0,-1,3,6,1,3,7,2,2,5,3},
						 {2,5,2,7,2,4,7,2,4,7,9,1}};
	double _B[_m]={4,5,7};
	
	
	const double pi=3.14159265358979;
	vector<vector<double>> template_of_A(24,vector<double>(13));
	int devide=10;
	double h=pi/devide;
	double alpha,beta,gamma;
	double No=0;

	int count=0;
	for(double alpha=0; alpha<2*pi-EPS ; alpha+=h){
		cout<<alpha<<endl;
		for(double beta=0; beta<pi-EPS; beta+=h){

			for(double gamma=0; gamma<2*pi-EPS; gamma+=h){
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

				while(1){
					vnode *p=Q->head;
					if(p->plane_label<=_m){
						int m=p->plane_label;
						vector<vector<double>> A(m+1+24,vector<double>(_n+1));
						for(int i=1; i<=24; i++){
							A[m+i]= template_of_A[i-1];
						}
						for(int i=1; i<=m-1; i++){	
							for(int j=1; j<=_n; j++){
								A[i][j]=p->path[i-1]*_A[i-1][j-1];
							}
							A[i][0]=p->path[i-1]*_B[i-1];
						}
						for(int j=1; j<=_n; j++){
							A[m][j]=_A[m-1][j-1];	
						}
						A[m][0]=_B[m-1];
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

						A[m][0]=-_B[m-1];
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
				}
			}
		}
		//show_matrix(&A);
	}
	cout<<"count"<<count<<endl;	

	
	return 0;
}