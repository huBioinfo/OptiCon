#ifndef STUCT_H
#define STUCT_H
//using namespace std;
const int OO=2147483647;  
typedef struct Matrix_double{
	double ** matrix;  
	int    _x;
	int    _y;
}Matrix_double;
void Hungarian_KM(Matrix_double,Matrix_double&);
int**   FormIntMatrix(int ,int);     
double** FormDblMatrix(int, int);     
int*    FormArray(int,int);
void    print_matrix(Matrix_double );
void    SetArray(int * &Array,int size,int data);
bool find(int x,int** w,int *visx,int *visy,int * lx,int *ly,int *linky,int *slack,int N);
#endif
