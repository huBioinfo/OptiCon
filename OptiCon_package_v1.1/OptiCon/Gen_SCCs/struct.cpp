#include "struct.h"
#include <iostream>
using namespace std;

void Hungarian_KM(Matrix_double S,Matrix_double& M)
{
	int m=S._x;
	int n=S._y;
	int /***S1,*/**w; double ** M1;
	int max_mn= m>n?m:n;
	int *linky,*lx,*ly,*slack,*visx,*visy;
	
	//S1 = FormIntMatrix(max_mn,max_mn);
	w =  FormIntMatrix(max_mn,max_mn);
	M1 = FormDblMatrix(max_mn,max_mn);
	linky = FormArray(max_mn,-1);
	lx =    FormArray(max_mn,0);
	ly =    FormArray(max_mn,0);
	slack = FormArray(max_mn,0);
	visx  = FormArray(max_mn,0);
	visy =  FormArray(max_mn,0);
/*	
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++){
			S1[i][j]=(int)(10000*S.matrix[i][j]);
		}
		for(i=0;i<max_mn;i++)
			for(int j=0;j<max_mn;j++)
				w[i][j]=S1[i][j];
*/
	for(int i=0;i<m;i++)
		for(int j=0;j<n;j++)
		{
			w[i][j]=(int)(10000*S.matrix[i][j]);
		}


			
	for(int i=0;i<max_mn;++i)
	{
		for(int j=0;j<max_mn;++j)
		{
			if(w[i][j]>lx[i])
				lx[i]=w[i][j];
		}
	}
			
	for(int x=0;x<max_mn;++x)
	{
		for(int i=0;i<max_mn;++i)
			slack[i]=OO;
		for(;;)
		{
			SetArray(visx,max_mn,0);
			SetArray(visy,max_mn,0);
			if(find(x,w,visx,visy,lx,ly,linky,slack,max_mn))
				break;
			int d=OO;
			for(int i=0;i<max_mn;++i)
			{
				if(!visy[i])
					if(d>slack[i])
							d=slack[i];
			}
			for(int i=0;i<max_mn;++i)
			{
				if(visx[i])
					lx[i]-=d;
			}
			for(int i=0;i<max_mn;++i)
			{
				if(visy[i])
					ly[i]+=d;
				else
					slack[i]-=d;
			}
		}
	}
	for(int x=0;x<max_mn;x++){
		M1[linky[x]][x] = 1;}

	for(int i = 0;i<M._x;i++)
		for(int x =0;x<M._y;x++)
		{
			M.matrix[i][x] = M1[i][x];
		}
				//	cout<<"*****************M"<<endl;
				//	print_matrix(M);

		//delete [] [] w;
		//delete [] [] M1;
		delete []linky;
		linky=NULL;
		delete []lx;
		lx=NULL;
		delete []ly;
		ly=NULL;
		delete []slack;
		slack=NULL;
		delete []visx;
		visx=NULL;
		delete []visy;
		visy=NULL;


		for(int i=0;i<max_mn;i++)
           {
                 delete []w[i];
                 w[i]=NULL;
           }
           delete []w;
           w=NULL;
		   for(int i=0;i<max_mn;i++)
           {
                 delete []M1[i];
                 M1[i]=NULL;
           }
           delete []M1;
           M1=NULL;


		   

	//w =  FormIntMatrix(max_mn,max_mn);
//	M1 = FormDblMatrix(max_mn,max_mn);
//	linky = FormArray(max_mn,-1);
//	lx =    FormArray(max_mn,0);
//	ly =    FormArray(max_mn,0);
//	slack = FormArray(max_mn,0);
//	visx  = FormArray(max_mn,0);
//	visy =  FormArray(max_mn,0);

}
bool find(int x,int** w,int *visx,int *visy,int * lx,int *ly,int *linky,int *slack,int N)
{
	visx[x]=true;
	for(int y=0;y<N;++y){
		if(visy[y])continue;
		int t=lx[x]+ly[y]-w[x][y];
		if(t==0){
			visy[y]=true;
			if(linky[y]==-1||find(linky[y],w,visx,visy,lx,ly,linky,slack,N)){
				linky[y]=x;
				return true;
			}
		}
		else{
			if(slack[y]>t)
				slack[y]=t;
		}
	}
	return false;
}
int**   FormIntMatrix(int x,int y)
{
	int ** matrix=new int *[x];
	for(int i =0 ;i<x;i++)
	{
		int *q= new int[y];	 //  array
		for(int j=0;j<y;j++)		{
			q[j]=0;
		}
		matrix[i]= q;
	}
	return matrix;
}
double**   FormDblMatrix(int x,int y)
{
	double ** matrix=new double *[x];
	for(int i =0 ;i<x;i++)
	{
		double *q= new double[y];	 //  array
		for(int j=0;j<y;j++)		{
			q[j]=0;
		}
		matrix[i]= q;
	}
	return matrix;
}
int*    FormArray(int size,int data)
{
	int * array = new int[size] ;
	for(int i=0;i<size;i++)
		array[i] = data;
	return array;
}

void    SetArray(int * &Array,int size,int data)
{
	for(int i=0;i<size;i++)
		Array[i] = data;
}

void    print_matrix(Matrix_double md)
{
	for(int i=0;i<md._x;i++){
		for(int j=0;j<md._y;j++)			
			cout<<md.matrix[i][j]<<"  ";
		cout<<endl;
	}
}
