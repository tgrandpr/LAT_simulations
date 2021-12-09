/*This code was written by Trevor GrandPre and David T. Limmer. It was used to create the simulation results in the paper entitled "Kinetic frustration by limited bond availability controls the LAT protein condensation phase transition on membranes" by Simou Sun, Trevor GrandPre, David T. Limmer, and Jay T. Groves. 

This code uses the GSL library. It can be compiled with the following command with proper linking to the libraries on your computer:
mpic++ -o dyno dbrunBP.cpp -lgsl -lm -lgslcblas -O3 -I/usr/local/include -L/usr/local/lib*/

/* insert libraries */
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>	
#include <algorithm>	
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>
#include <cmath>
using namespace std;

//the system dimensions in x and y
double Ly=100;
double Lx=100;

//totalnumber of particles
int N=1000;//totalnumber of particles
int N1=N;

//defining useful functions for the cluster code. 
void clusterHist();
void tet(int);
int clustercheck(int, int, int);

 vector<int> check(N1); 
 int ptypes=2;// types of particles. LAT and bonds. 
  vector< vector<int> > ClusterSize(ptypes, vector<int>(N1));
    vector<int> n_type(N1); 
vector<int> nn_type(N1); 
    vector<int> cluster_number(N1); 
     int n_clusters=0;
vector<vector<double> > bond(N1, vector<double>(N1));
 vector<double> nlist(N1);
vector<double> nnlist(N1);
    vector<vector<double> > vlist(N1, vector<double>(N1));
    vector<vector<double> > vvlist(N1, vector<double>(N1));
    int cluster_id;
    int l_pick,frame;
    double boxlx=Lx;
    double boxly=Ly;
    double ave_size;
    double ave_total;
    double large_cluster;
    int large_cluster_id;

//defining useful functions for the positions of the particles and the neighborlist
   vector<double> x(N1); 
    vector<double> y(N1);
    vector<double> xv(N1); 
    vector<double> yv(N1);

 double rreact=2.0;//reaction distanc
    double delta=0.5;//skin depth
    double rv=rreact+delta;//the neighborlist distance. 

int main(int argc, char *argv[]) {

  //defining functions for printout
   FILE* Fout3;
      Fout3 = fopen("clustersize.txt", "w");

    
    int t,i,j,k;
    double m,o,b,g;
    
    FILE *  cluster4;
    cluster4 = fopen("clusterdist.txt", "w");
    
    FILE * x_out;
    x_out = fopen("x_out_temp0.7_den0.4_WCA_RATES_10_10.lammpstrj","w");

    FILE * den1;
    den1 = fopen("density.txt","w");
    
    FILE * bondsfile;
    bondsfile = fopen("bonds1.txt","w");
    
    FILE * pos;
    pos = fopen("positions1.txt","w");
//the random number generator for the dynamics
    gsl_rng* rng;
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,time(NULL));
   

  
    double h = 0.0001; //timestep
    int Nstep =500000000; //total number of steps
    frame=0;//a counter used for the histogram of clusters
 
    double T = atof(argv[1]);//Temperature
	double dE = atof(argv[2]);//bonding energy for 2grb2 to sos and 2 sos to L.
    double kbond=atof(argv[3]);//prefactor rate which controls the time scale of bonds

    
    double beta = 1/T;//inverse temp
    double gamma =1.0;//friction
  
    int dump=10000;
    int cdump=10000;
    
   
   
    vector<double> xinit(N1), yinit(N1),thetainit(N1),nbondinit(N1);//preallocation
    
    //things for the bond making/breaking
    vector<vector<double> > rc(N1, vector<double>(N1));

    vector<double> nbond(N1);
    vector<double> nattmi(N1);
    vector<double> nattbi(N1);
    int nattm;
    int nattb;
   
    double rmin=1.5;//minimum length of the bond
    double kappa=80.0;//spring constant of the bond

    vector<double> cnt(Nstep);//a counter for bonding later in code. 
    
    //things for the force calculation
    vector<vector<double> > K(N1, vector<double>(2));
    double rcut=pow(2,1./6.);//cutoff for interaction 
    double ecut = 1.0; //this is epsilon in the WCA formula

    int id,type;  
    double z;
    double trash1;
    double trash2;
    double trash3;  
    
    //read in positions for for the initialization
    FILE *init;
    if ( ( init = fopen( "positions.txt", "r" ) ) == NULL){
      printf ("initial data could not be opened\n");}
    else {
      for(i=0;i<N;i++){
        fscanf(init, "%d %d %lf %lf %lf",&id, &type, &xinit[i], &yinit[i],&nbondinit[i]);
      }
      
    }

    rewind(init);
    fclose(init);
    
    //initialize the positions and bonds.
    for(i=0;i<N;i++){
        x[i] = xinit[i];
        y[i] = yinit[i];
        nbond[i]=nbondinit[i];
  
    }
    
    /////read in the bond matrix
    FILE *init1;
    if ( ( init1 = fopen( "bonds.txt", "r" ) ) == NULL){
      printf ("initial data could not be opened\n");}
    else {
      for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
          fscanf(init1, "%d %d %lf",&i, &j,&bond[i][j]);
        }
      }
      
    }
    
    rewind(init1);
    fclose(init1);

    
    double flag=1;//calculate the neigborlist at the first timestep. Then in subsequent steps the flag will be 1 only when a particle goes outside of the skin depth.
    
    
//run the dynamics
    for(int f=0;f<Nstep;f++){
        
        double mm=0;
        
        //Step 1: calculate the neighbor list
        if(flag==1){
        for (j=0; j<N; j++) {
            nlist[j]=0;//zero the number of neighbors for particle j
            xv[j]=x[j];//positions when list is made that are used to calculate when to remake neighbor list. 
            n_type[i]=0;
            yv[j]=y[j];
            for (i=0; i<N; i++){
                vlist[j][i]=0;//zero the list keeping track of particle indices. 
            }
        }

        for (k=0; k<N-1; k++){
            for (j=k+1; j<N; j++){
                double xr = x[k]-x[j];
                if(xr>Lx*0.5) {
                    xr-=Lx;}
                else if(xr<(-Lx*0.5)) {
                    xr=Lx+xr;}
                double yr = y[k]-y[j];
                if(yr>Ly*0.5) {
                    yr-=Ly;}
                else if(yr<(-Ly*0.5)) {
                    yr=Ly+yr;}      
                double r2 = xr*xr+yr*yr,r = sqrt(r2);
                if (r<=(rv)){//rv the neighborlist cutoff. 
                    vlist[k][nlist[k]]=j;
                    nlist[k]++;

                    vlist[j][nlist[j]]=k;
                    nlist[j]++;
                     if(r<=rreact){//the cutoff for reactions which less than rv.
                        rc[k][j]=1;//signalling that they are capable of forming a bond. 
                        rc[j][k]=1;
                        }
                }
               
                else if(r>rreact){//signalling that these particles are not able to react because they are out of range
                         rc[k][j]=0;
                         rc[j][k]=0;               
                }
                if(nlist[j]>1) n_type[j]=1;
            }
               if(nlist[k]>1) n_type[k]=1;
        }
    }
        
  flag=0;
 //zero the flag for the neigborlist here. When particles are found outside of the cutoff on the neighborlist calculation the flag will go to 1 signalling a recalculation of the neighborlist at the beginning of the next time step.     

        //STEP 2: calculate force fields/////////////////////////////////////////////////////////////////////////
              for (i=0; i<N; i++) {//initialize the force fields to zero.
                K[i][0] = 0;
                K[i][1] = 0;
            }
        
              for (i=0; i<N; i++) {//0 to N-1
                for (j=0; j<nlist[i]; j++) {//length of nlist
                    
                double jj=vlist[i][j];
                    
                    double xr = x[i]-x[jj];// the distance of particle within the 
                    if(xr>Lx*0.5) {
                         xr-=Lx;}
                    else if(xr<(-Lx*0.5)) {
                         xr=Lx+xr;}
                    double yr = y[i]-y[jj];
                    if(yr>Ly*0.5) {
                     yr-=Ly;}
                    else if(yr<(-Ly*0.5)) {
                     yr=Ly+yr;}      
          
                    double r2 = xr*xr+yr*yr,r = sqrt(r2);
                    
                    double pre=bond[i][jj]*(2.*kappa*(rmin-r))/r;//harmonic bond when particles bond is nonzero. 
                    K[i][0] = K[i][0]+pre*xr; // update pair force in components
                    K[i][1] = K[i][1]+pre*yr;

                    K[jj][0] =K[jj][0]-pre*xr;//Newtons third law.
                    K[jj][1] =K[jj][1]-pre*yr;
                    
                    
                    //the exlcuded area
                        if (r < rcut) {

                        double r2i = 1/r2;
                        double r6i = pow(r2i,3);
                        double ff = 48*ecut*r2i*r6i*(r6i-0.5); // L-J force

                        K[i][0] = K[i][0]+ff*xr; // update pair force in components
                       K[i][1] = K[i][1]+ff*yr;

                        K[jj][0] =K[jj][0]-ff*xr;//Newtons third law.
                       K[jj][1] =K[jj][1]-ff*yr;
                        }
           }
          }
 

        //setup dump file
        if(f%dump==0){
       
        fprintf(x_out,"%s\n","ITEM: TIMESTEP");
        fprintf(x_out,"%d\n",(int)(f/dump));
        fprintf(x_out,"%s\n","ITEM: NUMBER OF ATOMS");
        fprintf(x_out,"%d\n",4*N);
        fprintf(x_out,"%s\n","ITEM: BOX BOUNDS pp pp pp");
        fprintf(x_out,"%lf %lf\n",0.0,Lx);
        fprintf(x_out,"%lf %lf\n",0.0,Ly);
        fprintf(x_out,"%lf %lf\n",0.0,0.0);
        fprintf(x_out,"%s\n","ITEM: ATOMS id type x y z vx vz");
                
        }
        
      
        //STEP 2: update positions
        for(j=0;j<N;j++){

            double fd_term = sqrt(2 * h *T/ (gamma));
            double noise_0 = fd_term * gsl_ran_gaussian_ziggurat(rng,1);
            double noise_1 = fd_term * gsl_ran_gaussian_ziggurat(rng,1);

            x[j] +=(h / gamma) * K[j][0]+ noise_0;
            y[j] +=(h / gamma) * K[j][1]+ noise_1;

             //periodic boundaries
                if (x[j] > Lx) {x[j] -= (floor(x[j]/Lx))*Lx;}
                if (x[j] < 0) {x[j] += (1+(floor(-x[j]/Lx)))*Lx;}
                if (y[j] > Ly) {y[j] -= (floor(y[j]/Ly))*Ly;}
                if (y[j] < 0) {y[j] += (1+(floor(-y[j]/Ly)))*Ly;}
            
            double xx=xv[j]-x[j];//positions when list is made.
            double yy=yv[j]-y[j];
    
            if(xx>Lx*0.5) {
                xx-=Lx;}
            else if(xx<(-Lx*0.5)) {
                xx=Lx+xx;}
            if(yy>Ly*0.5) {
                yy-=Ly;}
            else if(yy<(-Ly*0.5)) {
                yy=Ly+yy;}  
            
            double rr1=sqrt(xx*xx+yy*yy);
            if (rr1>(delta/2)){

                        flag=1;
                        
                    }
        }//N
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        //////////Make/break bonds///////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        
        //variables used in reactions section
        int m;
        int k;
        int alpha;
        double b;
        double m1;
        int npick;
        
        for(m=0;m<N;m++){// this will run the reactions and print the trajectory file
            
            nattm=0;
            nattb=0;
            
            //This will go from zero to N-1;Must match the indices of nbond which also goes from 0 to N-1.
            
            for(k=0;k<N;k++){//in both cases we are only considering the case where rc=1, which means that the two are with rreact.
                if(rc[m][k]==1 & bond[m][k]==0 & k!=m & nbond[m]< 3 & nbond[k]< 3){
                    
                    nattmi[nattm]=k;//count all particles that are not bonded to "m" but are close enough to be bound to it.  
                    
                    nattm++;
                     
                }
                if(rc[m][k]==1 & bond[m][k]==1 & k!=m){//count all possibilities to break a bond
                    
                    nattbi[nattb]=k;//count all particles that are bound to particle m.
                      
                    nattb++;
                  
                }
                
                
            }
            //this is the number of bonds for particle m
            nbond[m]=nattb;

            //these are variables to used to visualize the bonds.
            double xg=0.0;
            double yg=0.0;
            double indexg=0.0;
            int p=0;

            //file to visualize bonds
            if(f%dump==0){
              
              //the print out file for positions of particles in multiples of 4. the other indices are for the bonds. 
              fprintf(x_out,"%d %d %lf %lf %lf %d %d\n",4*m,1,x[m],y[m],nbond[m]/10,0,0);
              
              //print the bond particles under each particle m
              p=0;
              if(nbond[m]>0){
                for (p=0; p<nbond[m]; p++){
                  xg=0.0;
                  yg=0.0;
                  
                  //for each bond in nbond get the particle indices that are bound to particle m
                  indexg=nattbi[int(p)];
                  indexg=int(indexg);
                  //Check that the bond is not through the box
                  double dx = x[m]-x[int(indexg)];
                if(abs(dx)<Lx*0.5){
                    xg = (x[m]+x[int(indexg)])/2;
                    
                }  
                else{
                  if(dx>Lx*0.5) {
                    dx-=Lx;}
                  else if(dx<(-Lx*0.5)) {
                    dx=Lx+dx;}
                    xg = x[m]-(dx)/2;
                    }
                    
                  double dy = y[m]-y[int(indexg)];
                     if(abs(dy)<Ly*0.5){
                         yg = (y[m]+y[int(indexg)])/2;
                     }
                    else{
                  if(dy>Ly*0.5) {
                    dy-=Ly;}
                  else if(dy<(-Ly*0.5)) {
                    dy=Ly+dy;} 
                     yg =y[m]-(dy)/2;
                    }
                  
                  // the posiions of the bonds with the last two values of the print out being the two LAT molecules connected by the bond particle
                  fprintf(x_out,"%d %d %lf %lf %lf %d %d\n",4*m+p+1,2,xg,yg,-0.1, int(4*indexg),4*m);
                  
									fflush(x_out);
                  
                  
                }
              }
              if(p<3){
                for(int a=p; a<3; a++){
                  //set extra unused bond particles to zero if a LAT has less than 3 bonds. this is purely for visualization. 
                  fprintf(x_out,"%d %d %d %d %lf %d %d\n",4*m+a+1,2,0,0,-0.1, int(0),int(0));
                  
                 fflush(x_out); 
                }
              }
            }    
            
            
        //pick a random number from zero to nattm+nattb.     
        npick=gsl_ran_flat(rng,0,nattm+nattb);
            npick=int(npick);
            
     //select at random what to do based on the number of possible actions. 
        if(npick<nattm){//make a bond if a random bumber is less than the number of possibilites to make a bond
            
            alpha=nattmi[int(npick)];//get indices of selected particles of the bond to break
            alpha=int(alpha);

                //propose bond
                    m1=gsl_ran_flat(rng,0,1);
                
                if(m1<(1-exp(-kbond*h))){
                //periodic boundaries
                    double xr = x[m]-x[int(alpha)];
                    if(xr>Lx*0.5) {
                         xr-=Lx;}
                    else if(xr<(-Lx*0.5)) {
                         xr=Lx+xr;}
                    double yr = y[m]-y[int(alpha)];
                    if(yr>Ly*0.5) {
                     yr-=Ly;}
                    else if(yr<(-Ly*0.5)) {
                     yr=Ly+yr;} 
                
                    double r2 = xr*xr+yr*yr;
                    double r = sqrt(r2);
                    double pe=kappa*(r-rmin)*(r-rmin);
                 
                //accept bond
                        b=gsl_ran_flat(rng,0,1);//draw another random number
                    
                        if (b<exp(-beta*(pe+dE))){//dE is the energy to combine the constituents
                           //the bond between the two LAT molecules are set to 1 and the number of bonds between both are increased by 1. 
                            bond[m][int(alpha)]=1;
                            bond[int(alpha)][m]=1;
                            nbond[m]++;
                            nbond[int(alpha)]++;
                       
                        }
                
                }        
           
        }
        else if(npick>=nattm & npick<(nattm+nattb)){
            
             alpha=nattbi[int(npick-nattm)];
             alpha=int(alpha);
          
        //pick a random number
           
        
             m1=gsl_ran_flat(rng,0,1);
            if(m1<(1-exp(-kbond*h))){
            //periodic boundaries
                    double xr = x[m]-x[int(alpha)];
                    if(xr>Lx*0.5) {
                         xr-=Lx;}
                    else if(xr<(-Lx*0.5)) {
                         xr=Lx+xr;}
                    double yr = y[m]-y[int(alpha)];
                    if(yr>Ly*0.5) {
                     yr-=Ly;}
                    else if(yr<(-Ly*0.5)) {
                     yr=Ly+yr;} 
                
                    double r2 = xr*xr+yr*yr;
                    double r = sqrt(r2);
                    double pe=kappa*(r-rmin)*(r-rmin);//energy of the bond
                 
                
                b=gsl_ran_flat(rng,0,1);//draw another random number
                 
                if (b<exp(+beta*(pe+dE))){
                
                   
                    //the bond between the two LAT molecules are set to zero and the number of bonds is decreased. 
                    bond[m][int(alpha)]=0;
                    bond[int(alpha)][m]=0;
                    nbond[m]--;
                    nbond[int(alpha)]--;
                 
                }
                
                
            }
       
        }
            cnt[f]+=nbond[m];//count the number of bonds on particle m.
            mm+=nbond[m];//total number of bonds. 
            
  
            
            if(f==(Nstep-1)){//print the positions of the second to last step for a restart file. 
              //I could
              fprintf(pos,"%d %d %lf %lf %lf\n",m,1,x[m],y[m],nbond[m]);
              
              
            }
            
            
            
            
}//N
        
        if(f%cdump==0){//Here is where the clusters are calculated. 
            ave_size=0.0;
		ave_total=0.0;
		large_cluster=0.0;
            tet(f);
            clusterHist();
            fprintf(Fout3,"%d %lf %lf %d\n",f,ave_size/ave_total,large_cluster,cluster_id);
            printf("%d %lf %lf %d\n",f,ave_size/ave_total,large_cluster,cluster_id);
						fflush(Fout3);
            
            frame++;
    
           
		}
	

        
         if(f%dump==0){//this printes the average number of bonds with time. 
             
             fprintf(den1,"%d %lf\n",f,mm/N);
             

         }

    }//Nstep
    
    //Save a file for the nbonds for a restart file at the end of the run.
    for (int i=0; i<N; i++) {//loop through all particles
      for (int j=0; j<N; j++) {

          fprintf(bondsfile,"%d %d %lf\n",i,j,bond[i][j]);
      
        }
    }
    
//Here we calculate the histogram of the cluster sizes
    for(int atom=1; atom<=N; atom++){
		if(ClusterSize[1][atom]!=0){
			fprintf(cluster4,"%d %e\n",atom, (double)(ClusterSize[1][atom])/(frame));
		}
	}
	


}//main

void tet(int t) {//code for cluster algorithm 
	int i,j,k;
	double xij,yij,zij,rij;
	
	for (j=0; j<N; j++) {//zero out two lists
		nnlist[j]=0;
        for (i=0; i<N; i++){
			vvlist[j][i]=0;
        }
	}
		
	for (k=0; k<N-1; k++){
        
            for (j=k+1; j<N; j++){

                        
                        if (bond[k][j]==1){//it computes cluster size based on bonds
                            vvlist[k][nnlist[k]]=j;
                            nnlist[k]++;
                            
                            vvlist[j][nnlist[j]]=k;
                            nnlist[j]++;

                        }
                    
                if(nnlist[j]>1) nn_type[j]=1;
            }
            if(nnlist[k]>1) nn_type[k]=1;

        
	}
	

}


void clusterHist(){//code for cluster algorithm 
	int ncluster;
	int type,k,atom;
	
	cluster_id=0;
	for(k=0; k<N; k++) check[k]=0;
    type=1;
	
    for(atom=0; atom<N; atom++){
        ncluster=0;
                    
        if(nn_type[atom]==type){
            if(check[atom]!=1) {
                cluster_id=(cluster_id+1);
                ncluster=clustercheck(atom,ncluster,type);
            }	
        }
        
        ClusterSize[type][ncluster]++;
        
        if(type==1 && ncluster!=0){
            ave_size+=ncluster*ncluster;
            ave_total+=ncluster;
            if(large_cluster<ncluster) {
                large_cluster=ncluster;
                if(type ==1) large_cluster_id=cluster_id;
            }
        }
        
        if(ncluster!=0) n_clusters++;
    }
	
}
				

int clustercheck(int atom, int ncluster, int type){//code for cluster algorithm 
	int j,jj;
	
	if(check[atom]!=1){
		ncluster++;
		check[atom]=1;
		cluster_number[atom]=cluster_id;
		
		for(j=0; j<nnlist[atom]; j++){
			jj=vvlist[atom][j];
			if(nn_type[jj]==type) ncluster=clustercheck(jj,ncluster,type);
		}	
	}	
	
	return ncluster;
}
