#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "struct.h"
#include <time.h>
#define number_of_run 1000  

#define number_of_node 6000 


char gene[number_of_node][55];   
char drivernode[number_of_node][55];

int sterm[number_of_node][number_of_node]={-1};   
int bud[number_of_node][number_of_node]={-1};

int  adjacentmatrix[number_of_node][number_of_node]={0};
int  maximummatch[number_of_node][number_of_node]={0};
int  maxi[number_of_node][number_of_run]={-1};//different maximum matchings
int  driver[number_of_node][number_of_run]={-1};
int  pailie[number_of_node]={0};
int tempCR[number_of_node][number_of_node]={0};

void exchange(int a[],int v,int ran);

int main()
{
	FILE *input,*result,*entrez,*CRMatrix; 
	

	char node[55]="";
	char node1[55]="";
	int number=0;
	int vet=0,edge=0;
	int i=0,j=0;
	int k=0,u=0;
	int flag=0,flag1=0;
	int drivernumber=0;


	printf("Data preprocessing...");

	for(i=0;i<number_of_node;i++)
	{
		pailie[i]=i;
		gene[i][0]='\0';
		drivernode[i][0]='\0';
	
		for(j=0;j<number_of_run;j++)
		{
			maxi[i][j]=-1;   
			driver[i][j]=-1; 
		}
	}


    //open files
	
	if((input=fopen("./MyGeneNetwork.txt","r"))==NULL) 
	{
		puts("can't open file 1");   
		exit(0) ;  
	}

	if((result=fopen("resultdata.txt","w"))==NULL)
	{
		puts("can't open file 2");   
		exit(0) ;  
	}
	if((entrez=fopen("GeneIDindex.txt","w"))==NULL)
	{
		puts("can't open file 3");   
		exit(0) ;  
	}



	
	

	vet=0;
	edge=0;
	number=0;
	while(!feof(input))
	{
		if(fscanf(input,"%s",node)==1)  
		{
	
		if(vet==0)
		{	
			strcpy(gene[vet],node);
			vet++;
			
		}
		else
		{
			flag=0;
			for(i=0;i<number_of_node;i++) 
			{
				if(strcmp(gene[i],node)==0)
				{	
					flag=1;
					break;
				}

			}
			if(flag==0)
			{
				strcpy(gene[vet],node);    
				vet++;                     
			}
		}
		               
		}
	
	}//end of while


	printf("\nthe number of vets of graph:=%d",vet);
	fprintf(result,"\nthe number of vets of graph:=%d",vet);  
	
	rewind(input);  
	for(i=0;i<vet;i++)
	{
		for(j=0;j<vet;j++)
		{
			adjacentmatrix[i][j]=-1;	 //Initiate adjacent matrix
			
		}
	}

	i=0;
	j=0;
	k=0;
	while(!feof(input))
	{
		fscanf(input,"%s",node);
		fscanf(input,"%s",node1);  
	
		flag=0;
		flag1=0;
		for(k=0;k<vet;k++)
		{
			if(flag!=1&&strcmp(gene[k],node)==0)  
			{	
				i=k;
				flag=1;
				flag1++;
			}
			if(flag!=2&&strcmp(gene[k],node1)==0)  
			{
				j=k;
				flag=2;
				flag1++;
			}
		}
		if(flag1==2)
		{
			adjacentmatrix[i][j]=1;   
		
		}
		else
			printf("\ncannot find the edge.");

	}

	k=0;
	printf("\n");

	for(i=0;i<vet;i++)
	{
		for(j=0;j<vet;j++)
		{
			if(adjacentmatrix[i][j]>=1)  
				k++;   //record the number of edges.
		}
	}

	edge=k;   
	printf("\nthe number of edges of graph:=%d",edge);
	

	fprintf(result,"\nthe number of edges of graph:=%d",edge);
	
/////////////////////////////////////////////////////////////////////////////////
	
	for(i=0;i<vet;i++)
	{
		strcpy(drivernode[i],gene[i]);  
	}



	// Hungarian algorithm
	Matrix_double S,M;
	S._x = M._x = vet;
	S._y = M._y = vet;
	S.matrix = FormDblMatrix(vet,vet);  
	M.matrix = FormDblMatrix(vet,vet);

	for(i=0;i<vet;i++)
	{
		for(j=0;j<vet;j++)
		{
					S.matrix[i][j]=0;
		}
	}
	for(i=0;i<vet;i++)
	{
		for(j=0;j<vet;j++)
		{
			if(adjacentmatrix[i][j]>0)
				S.matrix[i][j]=adjacentmatrix[i][j];	
		}
	}

	Hungarian_KM(S,M);  
	number=0;   

	for(i=0;i<vet;i++)
	{
		for(j=0;j<vet;j++)
		{
			if(M.matrix[i][j]>0&&adjacentmatrix[i][j]>0)  
			{
				printf("\n%s-------%s",gene[i],gene[j]);
				maximummatch[i][j]=1;
				drivernode[j][0]='\0';
				number++;
			}
		}
	}

	for(i=0;i<vet;i++)
	{
		if(drivernode[i][0]!='\0')   
			printf("\n%s",drivernode[i]);
	}

	printf("\n%d",number);   

	for(i=0;i<vet;i++)
	{
		if(drivernode[i][0]!='\0')
			fprintf(result,"\n%s",drivernode[i]);
	}

	fprintf(result,"\n the number of driver nodes:%d",vet-number);  
	
	
	for(i=0;i<vet;i++)
    {
        delete []S.matrix[i];
        S.matrix[i]=NULL;
    }
    delete []S.matrix;
	
    for(i=0;i<vet;i++)
    {
        delete []M.matrix[i];
        M.matrix[i]=NULL;
    }
    delete []M.matrix;
/////////////////////////////////////////////////////////////////////////////////////////////
    //different maximum matchings
	int difmax=0;

	for(k=0;k<number_of_run;k++)  
	{


		exchange(pailie,vet,k);  
				
		printf("\nMatching %d",k);
		fprintf(result,"\n");
		for(i=0;i<vet;i++)
		{
			fprintf(result,"*%d*",pailie[i]);  
		}

		Matrix_double S,M;
		S._x = M._x = vet;
		S._y = M._y = vet;
		S.matrix = FormDblMatrix(vet,vet);
		M.matrix = FormDblMatrix(vet,vet);
		for(i=0;i<vet;i++)
		{
			for(j=0;j<vet;j++)
			{
				S.matrix[i][j]=0;
			}
		}
		for(i=0;i<vet;i++)
		{
			for(j=0;j<vet;j++)
			{
				if(adjacentmatrix[pailie[i]][pailie[j]]>0)
					S.matrix[i][j]=adjacentmatrix[pailie[i]][pailie[j]];	
			}
		}	
						

		Hungarian_KM(S,M);
		int num=0;   
					

		for(i=0;i<vet;i++)
		{
							
			for(j=0;j<vet;j++)
			{
				if(M.matrix[i][j]>0&&adjacentmatrix[pailie[i]][pailie[j]]>0)
				{
									
					maxi[pailie[i]][k]=pailie[j];
					driver[pailie[j]][k]=0;  
					num++;   
									
				}
			}
		}				
		for(u=0;u<difmax;u++)
		{
			flag=0;
			for(j=0;j<vet;j++)
			{
								
				if(maxi[j][u]!=maxi[j][k])  
				{flag=1;break;}
			}
			if(flag==0)
				break;
		}
		if(flag==1)  
		{
							
			for(j=0;j<vet;j++)
			{
				maxi[j][difmax]=maxi[j][k];
				driver[j][difmax]=driver[j][k];
			}
							
		    difmax++;

		}
		if(k==0&&difmax==0)
			difmax++;		
	    
		
		//delete [] 
		for(i=0;i<vet;i++)
		{
			delete []S.matrix[i];
			S.matrix[i]=NULL;
		}
		delete []S.matrix;
		
		for(i=0;i<vet;i++)
		{
			delete []M.matrix[i];
			M.matrix[i]=NULL;
		}
		delete []M.matrix;

	}  

    
	for(u=0;u<difmax;u++)  //u: the number of different maximum matchings.
	{
							
							
		for(j=0;j<vet;j++)
		{
			if(maxi[j][u]>0)
			{fprintf(result,"\n%d:%s-------%s",u,gene[j],gene[maxi[j][u]]);}
									
		}
		
	
		
		for(j=0;j<vet;j++)
		{
			if(driver[j][u]<0)  
			{fprintf(result,"\ndrivernode:%s",gene[j]);
			 
			}						
		}
		
		
							
	} 
				
//////////////////////////////////////////////////////////////////


	printf("\nOutput different CFs.");




//Construct cacti

	int stermnumber=0;
	int budnumber=0;
	int a[number_of_node];
	int ca[number_of_node];
	int bnumber=0;
	int count1=0;
	int count2=0;
	int p=0,q=0;
	int changenum[number_of_run];
	char filename[30];

	for(i=0;i<number_of_node;i++)
	{
		a[i]=-1;
		ca[i]=-1;
		for(j=0;j<number_of_node;j++)
		{
			sterm[i][j]=-1;
			bud[i][j]=-1;
			tempCR[i][j]=0;
		}
	}


	printf("\ndifmax:%d\n",difmax);

	for(int e=0;e<difmax;e++)  
	{
		changenum[e]=0;
		printf("\ndifmax:%d\n",e);

		for(i=0;i<vet;i++)
		{
			a[i]=i;
		}
		stermnumber=0;  
		budnumber=0;

		for(i=0;i<number_of_node;i++)
		{
		
			ca[i]=-1;
			for(j=0;j<number_of_node;j++)
			{
				sterm[i][j]=-1;
				bud[i][j]=-1;			
			}
		}


		for(i=0;i<vet;i++)
		{
			if(driver[i][e]<0)  
			{
				k=0;  //k: the number of nodes in the stem.
				sterm[stermnumber][k]=i;  
			
				a[i]=-1;  
				u=i;
				flag=0;
				while(flag==0)
				{
					flag=1;
					for(j=0;j<vet;j++)
					{
						if(maxi[u][e]==j)   
						{
							k++;
							sterm[stermnumber][k]=j;
							tempCR[u][j]=1;
							a[j]=-1; 
							u=j;
							flag=0;
							break;
						}
					}
				}
				stermnumber++;
			}
		} 
		
		for(i=0;i<vet;i++)
		{
			if(a[i]>=0)  
			{
				k=0;
				u=i;
				bud[budnumber][k]=u;
				a[i]=-1;
				flag=0;
				while(flag==0)
				{
					flag=1;
					for(j=0;j<vet;j++)
					{
						if(maxi[u][e]==j&&a[j]>=0)
						{
							k++;
							bud[budnumber][k]=j;
							tempCR[u][j]=1;
							u=j;
							flag=0;
							a[j]=-1;
							break;
						}
					}
				}
				budnumber++;
			}
		}
		
		fprintf(result,"\n %d stems:\n",stermnumber);  
		for(i=0;i<stermnumber;i++)
		{
			j=0;
			while(sterm[i][j]>=0)
			{	
				fprintf(result,"%15s",gene[sterm[i][j]]);  
				j++;
			}
			fprintf(result,"\n");
		}

		fprintf(result,"\n %d buds:\n",budnumber);  
		for(i=0;i<budnumber;i++)
		{
			j=0;
			while(bud[i][j]>=0)
			{	
				fprintf(result,"%15s",gene[bud[i][j]]);  
				j++;
			}
			fprintf(result,"\n");
		}

		for(k=0;k<stermnumber;k++)  
		{
			for(i=0;i<number_of_node;i++)
			{
				ca[i]=-1; 
			}
			number=0;
			j=0;
			while(sterm[k][j]>=0) 
			{
				ca[j]=sterm[k][j]; 
				j++;

			}
			number=j;  
			
			for(i=0;i<number;i++)
			{
				for(j=i;j<number;j++)  
				{
					if(sterm[k][j]>=0)
					{
						tempCR[sterm[k][i]][sterm[k][j]]=1; 
					}
				}
			}
			
			int nn=number-1;	
			for(j=0;j<budnumber;j++) 
			{
				flag=0;
				for(i=0;i<nn;i++)  
				{
					u=0;
					bnumber=0;
					while(bud[j][u]>=0) 
					{
						u++;
						
					}
					bnumber=u; 
					
					for(p=0;p<u;p++)
					{
						for(q=0;q<u;q++)   
						{
							if(bud[j][q]>=0)
							{
								tempCR[bud[j][p]][bud[j][q]]=1;  
							}
						}
					}
					
					u=0;
					while(bud[j][u]>=0)
					{
						if(adjacentmatrix[ca[i]][bud[j][u]]>0) 
						{
							//flag=1;
							
							for(p=0;p<=i;p++)
							{
								for(q=0;q<bnumber;q++)
								{
									if(bud[j][q]>=0)  
									{
										tempCR[ca[p]][bud[j][q]]=1;  
									}
								}
							}
							break;  
						}
						u++;
						
					}
				} 
			}


		   
	    }  
		
		//Output direct control regions derived from different SCCs.
		sprintf(filename,"tempCRMatrix_%d.txt",e);
	    if((CRMatrix=fopen(filename,"w"))==NULL)
	    {
		  puts("can't open file 5");   
		  exit(0) ;  
	    }
	    for(i=0;i<vet;i++)
	    {
		   for(j=0;j<vet;j++)
		   {
			fprintf(CRMatrix,"%d    ",tempCR[i][j]);  
		   }
		   fprintf(CRMatrix,"\n"); 
	    }
		fclose(CRMatrix);  
		
		
		for(i=0;i<vet;i++)   
			for(j=0;j<vet;j++)
				tempCR[i][j]=0;
    }
	
	fclose(result);




////////////////////////////////////////////////////////////////////////////////////////


	for(i=0;i<vet;i++)
	{
		fprintf(entrez,"%55s\n",gene[i]);
		
	}
    fclose(entrez);



///////////////////////////////////////////////////////////////////////////////////////////////



		
	system("pause");
	return 0;
}



void exchange(int a[],int v,int ran)
{
	int i=0;
	int j=0;
	int number=0;
	int b=0;
	
	for(i=0;i<v;i++)
	{
		a[i]=i;
	}
	srand(time(NULL)); 
	for(i=0;i<v;i++)
	{
		srand(time(NULL)+i+ran*v); 
		number=rand()%v;  
		for(j=0;j<v;j++)
		{
			if(a[j]==number)
				break;
		}
		b=a[i];
		a[i]=a[j];
		a[j]=b;
	}
}
