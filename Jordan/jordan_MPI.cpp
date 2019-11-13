#include "head_MPI.h"

int main(int argc, char *argv[]){
	MPI_Init(&argc, &argv);
	int err=0,p,k,n;
	char *name=0;
	MPI_Comm G=MPI_COMM_WORLD;
	MPI_Comm_size(G,&p);
	MPI_Comm_rank(G,&k);
	
	
	if((argc>3) || (argc<2) || ((n=atoi(argv[1]))<=0)){
		if(!k) printf("Using %s n <filename>\n",argv[0]);
		MPI_Finalize();
		return-1;
	}
	
	if(argc==3) name=argv[2];
	
	if(solve(n,name,p,k,G)){
		err=-2;
	}
	
	MPI_Finalize();
	return err;
}			
	
