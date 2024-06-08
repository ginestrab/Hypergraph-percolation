/**************************************************************************************************
 * If you use this code, please cite G. Bianconi and S. Dorogovstev
 "Theory of percolation on hypergraphs"
 Physical Review E, 109, p.014306 (2024).
***************************************************************************************************
 * Code that  generates random hypergraph whose hyperedge have fixed cardinality and 
	then performs MonteCarlo simulations and Message Passing predictions of nodes, 
	and hyperedge percolation using the factor node and the hypergraph algorithm
 *
 * This code uses:
 * N  Number of nodes
 * M number of hyperedges
 * m  fixed cardinality of the hypedges
 * Hyperedge options 1/0 Hyperedge=1 for hyperedge percolation and Hyperedge=0 for node percolation.
 * Hypergraph 1/0  Hypergraph 1 perfroms Hypergraph percolation if Hyperedge=0 otherwise performs factor graph percolation
 * Hypergraph can be 1 only for Hyperedge=0
 *
 * Nrunmax  Number of MonteCarlo simulations
 *************************************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define N 10000
#define Q 10000
#define m  4
#define Hypergraph 1 
#define Hyperedge 0 
#define Nrunmax 10
 

#define err_tol 0.01
int *vis1,*size_cluster1,**knn1,*k1,c1,c2,c3,*occ,*dam1,*dam2;


/**************************************************************************
 Recurrence is a subrutine for calulating the giant component in each layer
 **************************************************************************/
int Recurrence( int i , int cluster_size, int ncluster){
	int j, n3,aus_cluster_size,n,j2,aus;
	
    
    if(i<N)	{cluster_size++;}
		vis1[i]=ncluster;
    for(n3=0;n3<k1[i];n3++){
			j=knn1[i][n3];
            if(j<N){
                if((vis1[j]==0)&&(dam1[j]==1)){
				aus_cluster_size=Recurrence(j, cluster_size, ncluster);
				cluster_size=aus_cluster_size;
			}
            }
            if (j>=N){
                aus=1;
                if(Hypergraph==1){
              //      printf("ci sono\n");
                    for(n=0;n<k1[j];n++){
                        j2=knn1[j][n];
                        aus=aus*dam1[j2];
                    }
                }
                if((aus==1)&&(vis1[j]==0)&&(dam1[j]==1)){
                    aus_cluster_size=Recurrence(j, cluster_size, ncluster);
                    cluster_size=aus_cluster_size;
                }
                    
            }
               
		}
	
		return cluster_size;
}

int main(int argc, char** argv){
	int i,j,n,it,ncluster1,ncluster2,ncluster3, GC,nn,mu,*occ2,cluster_size,m1,m2,m3,m1_aus,m2_aus,m3_aus,c1_aus,c2_aus,c3_aus,*sigma,nrun,nc,np;
	int s1,s2,Nc1,Nc2,Nc3,aus,aus3,**adj1,**adj2,**adj3,N0,i_aus;
	float p,x,f,**xd1,*Sm,**n1,MGC1,nsum1,nsumold1,nsum2,aus1;
    int *GCm;
	
	   char filename[60],string1[50],string2[50];
	
	FILE *gp3,*gp;
    


  /**************************************************************************
   open file for output 
   finemane GC
   at the end of teh program the file will contain 
   two columns: p GC GC2
   
   **************************************************************************/

		srand48(time(NULL));
	

	N0=N;
	

	vis1=(int*)calloc(N+Q,sizeof(int));
	occ=(int*)calloc(N+Q,sizeof(int));
    occ2=(int*)calloc(m,sizeof(int));
	k1=(int*)calloc(N+Q,sizeof(int));
	xd1=(float**)calloc(N+Q,sizeof(float*));
	dam1=(int*)calloc(N,sizeof(int));
    dam2=(int*)calloc(N+Q,sizeof(int));
	knn1=(int**)calloc(N+Q,sizeof(int*));
    n1=(float**)calloc(N+Q,sizeof(float*));
    GCm=(int*)calloc(100,sizeof(int));
		for(i=0;i<N+Q;i++){
			knn1[i]=(int*)calloc(N,sizeof(int));
            n1[i]=(float*)calloc(N+Q,sizeof(float));
        }
    for(i=0;i<N+Q;i++){
			xd1[i]=(float*)calloc(N+Q,sizeof(float));
		}
	size_cluster1=(int*)calloc(N+Q,sizeof(int));
	
    
    for(i=0;i<N+Q;i++){
        k1[i]=0;
    }

	for(mu=0;mu<Q;mu++){
        for (n=0;n<m;n++){
            aus=0;
            while(aus==0){
            i=(int)((float)N*drand48());
                if(mu+N>N+Q) printf("problem\n");
            aus=1;
            for (nn=0;nn<n;nn++){
                if(occ2[nn]==i) aus=0;
            }
            }
            occ2[n]=i;

				
            k1[i]++;
            k1[mu+N]++;
            
            knn1[i][k1[i]-1]=mu+N;
            knn1[mu+N][k1[mu+N]-1]=i;
            
            n1[i][mu+N]=(int)(2*drand48());
            n1[mu+N][i]=(int)(2*drand48());
            nsum1+=n1[i][mu+N]+n1[mu+N][i];
        }
    }
		
  
/*%%%%%%%%%MESSAGGE PASSING%%%%%%%%%%%%*/
    gp=fopen("GC_Hypergraph_MP.txt","w");
    for(f=0.;f<1.;f+=0.01){
    p=1-f;
    nsumold1=1000000;

   
    while(fabs(nsum1-nsumold1)>err_tol){
        nsumold1=nsum1;
        nsum1=0;
       
    
        for (it=0;it<N;it++){
            i=(int)((N+Q)*drand48());
            i_aus=1;
            if(i<N) i_aus=0;
            for(n=0;n<k1[i];n++){
                j=knn1[i][n];
                aus1=1.;
              
                    for(np=0;np<k1[i];np++){
                        if(np!=n){
                            aus1=aus1*(1.-n1[knn1[i][np]][i]);
                        }
                    }
                nsum1-=n1[i][j];
                
                if((Hypergraph==1)&&(Hyperedge==0)){
                    n1[i][j]=pow(p,i_aus*(k1[i]-1))*(1-aus1);}
                if((Hypergraph==0)&&(Hyperedge==1)){
                    n1[i][j]=pow(p,(i_aus))*(1-aus1);
                }
                if((Hypergraph==0)&&(Hyperedge==0)){
                        n1[i][j]=pow(p,(1-i_aus))*(1-aus1);
                    }
            
                nsum1+=n1[i][j];
            }
                    
        }
    }
    
    MGC1=0;

    for (i=0;i<N;i++){
        if(k1[i]>0){
        aus1=1;
        for(n=0;n<k1[i];n++){
            aus1=aus1*(1-n1[knn1[i][n]][i]);}
            MGC1+=(1-aus1);
        }
    }
        if(Hyperedge==0)
            MGC1=MGC1*p;
        
fprintf(gp,"%lf %lf \n",1-f,(float)MGC1/((float)N));

    
    }


fclose(gp);
				
/*%%%%%%%%%%%%SIMULATION%%%%%%%%%%%%%%%%%*/
	for (nrun=0;nrun<Nrunmax;nrun++){	
	
        for(i=0;i<N+Q;i++){
            xd1[i][0]=drand48();
            dam1[i]=1;
            
        }
    
        for(nc=0;nc<100;nc++){
		f=nc*0.01;
            if(Hyperedge==0){
                for (i=0;i<N;i++){
                 if(xd1[i][0]<f){
                     dam1[i]=0;}
                 else {dam1[i]=1;}
             }
            }
            if(Hyperedge==1){
            for(i=0;i<Q;i++){
                if(xd1[i+N][0]<f){
                    dam1[i+N]=0;}
                else {dam1[i+N]=1;}
            }
            }
		for(i=0;i<N+Q;i++){
			vis1[i]=0;
		}
    
		GC=0;
        m1=0;
		ncluster1=0;
		for(n=0;n<N;n++){
			if((vis1[n]==0)&&(dam1[n]==1)){
				cluster_size=0;
				ncluster1++;
				cluster_size=Recurrence(n, cluster_size, ncluster1);
				size_cluster1[ncluster1]=cluster_size;
				if(cluster_size>m1){m1=cluster_size;
					c1=ncluster1;}				
			}			
        }
        Nc1=c1;
        GC=0;
        for(i=0;i<N;i++){
            if(vis1[i]==Nc1){
                GC++;
				//	m2++;
            }
        }
        GCm[nc]+=GC;
        }
    }
    gp3=fopen("GC_Hypergraph_simulations.txt","w");
        for (nc=0;nc<100;nc++){
            f=nc*0.01;
        
        fprintf(gp3,"%f %f\n",1-f,(float)GCm[nc]/(float)(Nrunmax*N));
        }
		fclose(gp3);
	
	return 0;	
}

