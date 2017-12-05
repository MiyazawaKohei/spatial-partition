#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
#include "simplex.h"
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
	const int _m=4;
	const int _n=2;
	double _A[_m][_n]=	{{1,1},
						 {2,0},
						 {2,1},
						 {8,2}};
	double _B[_m]={4,5,7,9};
	
	partition_tree *T;
	T= new partition_tree;
	T->head= new vnode;
	T->head->plane_label=1;

	queue_of_vnode *Q;
	Q= new queue_of_vnode;
	Q->head=T->head;
	Q->tail=T->head;

	int count=0;
	while(1){
		vnode *p=Q->head;
		if(p->plane_label<=_m){
			int m=p->plane_label;
			vector<vector<double>> A(m+1,vector<double>(_n+1));

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
	cout<<"count"<<count<<endl;
	return 0;
}