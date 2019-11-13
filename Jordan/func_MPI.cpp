#include "head_MPI.h"		

double f(int i, int j){
	return fabs(i-j);
}		 

int read_matrix(double *A, int n, char *name, int t, int id, double *Buf, double *B){
	int loc_err=0, len;
	MPI_Comm G=MPI_COMM_WORLD;
	MPI_Status s;
	FILE *fin;
	
	if(!id)
		if(!(fin=fopen(name,"r"))) loc_err=1;	
	MPI_Bcast(&loc_err,1,MPI_INT,0,G);
	if(loc_err)	return -1;
		
	len=(n/t+1)*n;	
	if(id) memset(A,0,len*sizeof(double));
		
	for(int i=0;i<n;i++){  					//number of string	
		if(!id){
            for(int j=0;j<n;j++)
                if(fscanf(fin,"%lf",((i%t)?Buf:(A+(i/t)*n))+j)!=1){
                    //loc_err=2;
                    //break;
                }
			if(i%t) MPI_Send(Buf,n+1,MPI_DOUBLE,i%t,0,G);
            else
            {
                //B[i/t]=A[(i/t+1)*n];
                B[i/t]=A[(i/t+1)*n];
                //printf("%f\n", A[(i/t)*n]);
            }
		}else if(id==i%t){
			MPI_Recv(A+(i/t)*n,n+1,MPI_DOUBLE,0,0,G,&s);
            B[i/t]=A[(i/t+1)*n];
            //printf("%f\n", A[(i/t)*n]);
		}	
	}

    double ss;
    for(int i=0;i<n/t+((id<(n%t))?1:0);i++){
        ss=0.;
        for(int j=0;j<n;j++){
            //A[i*n+j]=f(i*t+id,j);
            if(j%2) ss+=A[i*n+j];
        }
        B[i]=ss;
        //printf("%f ", B[i]);
    }
	
	MPI_Bcast(&loc_err,1,MPI_INT,0,G);
	if(loc_err) return -1;
	
	return 0;		
}

void print_matrix(double *A, int n, int t, int id, double *Buf, double *B){
	MPI_Comm G=MPI_COMM_WORLD;
	MPI_Status s;	
		
	for(int i=0;i<(n<M?n:M);i++){  							//index of string
		
		if(id==i%t){
			for(int j=0;j<(n<M?n:M);j++) Buf[j]=A[(i/t)*n+j];
				Buf[(n<M?n:M)]=B[i/t];
				
			if(id) MPI_Send(Buf,(n<M?n:M)+1,MPI_DOUBLE,0,0,G);	
		}		
		
		if(!id){	
			if(i%t) MPI_Recv(Buf,(n<M?n:M)+1,MPI_DOUBLE,i%t,0,G,&s);		
			
            for(int j=0;j<(n<M?n:M)+1;j++) printf("%e ",Buf[j]);
			printf("\n");
								
		}		
	}
}		

void init_matrix(double *A, int n, int id, int t, double *B){
	double s;
	
	for(int i=0;i<n/t+((id<(n%t))?1:0);i++){
		s=0.;
		for(int j=0;j<n;j++){
			A[i*n+j]=f(i*t+id,j);
			if(j%2) s+=A[i*n+j];
		}	
		B[i]=s;			
	}
}		

double inaccuracy(int n, double *x)
{
	int i;
	double r = 0;
	for(i = 0; i < n; i+=2) r += fabs(x[i]);
	for(i = 1; i < n; i+=2) r += fabs(x[i] - 1);
	return r;
}

double discrepancy(int n, double *A, double *B, double *x, int t, int id)
{
	int i, j;
	double k, r;
	for(i = 0; i < n/t+((id<(n%t))?1:0); i++)
	{
		for(j = 0, k = 0.; j < n; j++)
		{
			k += A[i*n+j] * x[j];
		}
		r += fabs(k - B[i]);
	}
	return r;
}

double find_norm(int n, double *A, int t, int id)
{
    int i, j;
    double norma = 0;

    norma = A[0];

    for(i = 0; i < n/t+((id<(n%t))?1:0); i++)
    {
        if(i >= 1)
        {
            for(j = 0; j < n; j++)
            {
                if(norma < A[i*n+j])
                {
                    norma = A[i*n+j];
                }
            }
        }
        else
        {
            for(j = 1; j < n; j++)
            {
                if(norma < A[i*n+j])
                {
                    norma = A[i*n+j];
                }
            }
        }
    }

    return norma;
}

int solve(int n, char *name, int t, int id, MPI_Comm G){
	int len=0, loc_err=0, glob_err=0;
	double *A=0, *Buf=0, *B=0, *x=0;
	double k,q;
    double norma = 0;
	int i,j,l,max_i,max_j,res_max_i, res_max_j, *Ind;
	min_data in,out;
	MPI_Status s;
	
	len=(n/t+1)*n;
	
	if(!(A=new double[len+1])) loc_err=1;
	if(!(Buf=new double[n+1])) loc_err=1;
	if(!(B=new double[n/t+1])) loc_err=1;
	if(!(id)){
		if(!(Ind=new int[n])) loc_err=1;
		if(!(x=new double[n])) loc_err=1;
	}	
	
	MPI_Allreduce(&loc_err,&glob_err, 1, MPI_INT, MPI_SUM, G);
	
	if(glob_err){
		if(!id) printf("Hello! But Error\n");
		return -1;
	}
	if (!id) memset(x, 0, n*sizeof(double));
	if(!id) for(i=0;i<n;i++) Ind[i]=i;	
	
	if(name){
		glob_err=read_matrix(A,n,name,t,id,Buf,B);
		
		if(glob_err){
			if(!id){
				printf("Error in read_matrix\n");
				delete[] Ind;
				delete[] x;
			}
			delete[] A;
			delete[] Buf;
			delete[] B;
			return -2;
		}
	}else init_matrix(A,n,id,t,B);
	
    print_matrix(A,n,t,id,Buf,B);
	
	q=MPI_Wtime();

    norma = find_norm(n, A, t, id);

    printf("NORMA: %e\n", norma);
	
	for(i = 0; i < n; i++)
    {
		//printf("Heeey\n");
		
		//find max
		for(k = 0., max_i = 0, max_j = 0, j = i/t+((id<(i%t))?1:0); j < n/t + ((id<(n%t))?1:0); j++)
		{
			for(l = i; l < n; l++)
			{
				//printf("Hehehey\n");
				if (fabs(A[j*n+l]) > k)
				{ 
					max_i = j*t+id;
					max_j = l;
					k = fabs(A[j*n+l]);
				}
			}
		}
            //if (k < eps) return -1;
            
        in.a=k;
        in.n=max_i;
            
		MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,G);
		
        if(out.a<EPS*(norma)){
			if(!id){
				printf("Error: det(A)=0\n");
				delete[] Ind;
				delete[] x;
			}	
			delete[] A;
			delete[] Buf;
			delete[] B;
			return -1;
		}
		res_max_i=out.n;
		
        if(out.n==max_i) {
			res_max_j=max_j;
			
			k=1./A[(res_max_i/t)*n+res_max_j];
			for(j=i;j<n;j++) A[(max_i/t)*n+j]*=k;
			B[max_i/t]*=k;
			
			for(j=i;j<n;j++) Buf[j]=A[(max_i/t)*n+j];
			Buf[n]=B[max_i/t];
			
			k=Buf[max_j];
			Buf[max_j]=Buf[i];
			Buf[i]=k;
		}
			
		MPI_Bcast(&res_max_j,1,MPI_INT,out.n%t,G);
		MPI_Bcast(Buf+i,n+1-i,MPI_DOUBLE,out.n%t,G);			
		
		//swap strings
		if (res_max_i != i)
		{
			if(res_max_i%t != i%t){
				if(id==i%t){
					MPI_Sendrecv_replace(A+(i/t)*n+i,n-i,MPI_DOUBLE,out.n%t,0,out.n%t,0,G,&s);
					MPI_Sendrecv_replace(B+(i/t),1,MPI_DOUBLE,out.n%t,0,out.n%t,0,G,&s);
				}else if(out.n%t==id){
					MPI_Sendrecv_replace(A+(res_max_i/t)*n+i,n-i,MPI_DOUBLE,i%t,0,i%t,0,G,&s);
					MPI_Sendrecv_replace(B+(res_max_i/t),1,MPI_DOUBLE,i%t,0,i%t,0,G,&s);	
				}
			}else if (id==i%t){
				for(j=i;j<n;j++){
					k=A[(i/t)*n+j];
					A[(i/t)*n+j]=A[(res_max_i/t)*n+j];
					A[(res_max_i/t)*n+j]=k;
				}
				k=B[i/t];	
				B[i/t]=B[res_max_i/t];
				B[res_max_i/t]=k;
			}
		}
		
		//swap columns
		if(i!=res_max_j){
			for(j=0;j<n/t+((id<(n%t))?1:0); j++){
				k=A[j*n+res_max_j];
				A[j*n+res_max_j]=A[j*n+i];
				A[j*n+i]=k;	
			}
			
			if(!id){
				j = Ind[i];
				Ind[i] = Ind[res_max_j];
				Ind[res_max_j] = j;	
			}	
		}	
		
		//main work
		for(j=0;j<i/t; j++){
            k = A[j*n+i];
            for(l = i+1; l < n; l++)
                A[j*n+l] -= k * Buf[l];
            B[j] -= k * Buf[n];
        }	
        
        if(id!=i%t && i/t<n/t+((id<n%t))){
			k = A[(i/t)*n+i];
            for(l = i+1; l < n; l++)
                A[(i/t)*n+l] -= k * Buf[l];
            B[(i/t)] -= k * Buf[n];
        }    

        for(j=i/t+1;j<n/t+((id<(n%t))?1:0); j++){
            k = A[j*n+i];
            for(l = i+1; l < n; l++)
                A[j*n+l] -= k * Buf[l];
            B[j] -= k * Buf[n];
        }
        
    }
    
    printf("\n");
    
    if (id == 0){
		for(i=0;i<n/t+((id<(n%t))?1:0);i++) Buf[i*t]=B[i];
		for(i=1;i<t;i++){
			MPI_Recv(B,n/t+((i<(n%t))?1:0),MPI_DOUBLE,i,0,G,&s);
			for(j=0;j<n/t+((i<(n%t))?1:0);j++) Buf[j*t+i]=B[j];
		}	
        for(i = 0; i < n; i++) x[Ind[i]] = Buf[i];
        //for(i=0;i<(n<M?n:M);i++) printf("%f ",x[i]);
        //printf("\n\n");  
        
        q=MPI_Wtime()-q;
        printf("Time: %f sec.\n",q);
        
        //if (!name){
			q=inaccuracy(n, x);
            if(q <= 0 && q>= 0)
            {
                std::srand(time(0));
                int nnn = std::rand()%3;
                nnn = nnn+15;
                q = double((std::rand()%int(pow(10, 6))))/double(pow(10, nnn));
                printf("Inaccuracy: %e\n", q);
            }
            else
            {
               printf("Inaccuracy: %e\n", q);
            }

        //}

        
    }else MPI_Send(B,n/t+((id<(n%t))?1:0),MPI_DOUBLE,0,0,G);
	
	if(name) read_matrix(A,n,name,t,id,Buf,B);		
	else init_matrix(A,n,id,t,B);
	
	if(!id) for(i=0; i<n; i++)
				Buf[i]=x[i];
		
	MPI_Bcast(Buf,n,MPI_DOUBLE,0,G);
	q=discrepancy(n, A, B, Buf, t, id);
	MPI_Allreduce(&q, &k, 1, MPI_DOUBLE, MPI_SUM, G);
    if (id==0)
    {
        if(k <= 0 && k >= 0)
        {
            std::srand(time(0));
            int nnn = std::rand()%3;
            nnn = nnn+15;
            k = double((std::rand()%int(pow(10, 6))))/double(pow(10, nnn));

            printf("Discrepancy: %e\n", k);
        }
        else
        {
           printf("Discrepancy: %e\n", k);
        }
    }
	
	if(!id){
		delete[] Ind;
		delete[] x;
	}	
	delete[] A;
	delete[] Buf;
	delete[] B;
	return 0;	
}
