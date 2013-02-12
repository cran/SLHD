//# include <iostream>
//using namespace std;
//# include <fstream>
# include <math.h>
//# include <iomanip>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <R.h>
#define	 int8	      unsigned int
#define  Min(a,b)     (a<b)?a:b
#define  Max(a,b)     (a>b)?a:b
#define	 MAX	      (int8) (pow(2,31)-1)
#define  autoseed(x)  (int8)(rand()*x+x)


extern "C" {
int rc(int n ,double seed) ;
double runif(double seed);  
int rc2(int n,int del,double seed);
int **LHD(int n,int k,double seed);
int **SLHD(int m, int t, int k, double seed);
void distmatrix(int **A, int n, int k, double *d);
void avgdist(int n, int p, double *d, double *avgdist_cur);
double combavgdist(int m, int t, int p, double *d, double *avgdist_slice, double *avgdist_cur);
void update_distmatrix(int **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old);
void revert_distmatrix(int n, int selrow1, int selrow2, double *d, double *d_old);
void update_avgdist(int n, int p, int selrow1, int selrow2, double *d, double *d_old, double *avgdist_old, double *avgdist_cur);
void update_avgdist_sliceI(int n, int m, int p, int translice, int tran1, int tran2, double *d, double *d_old, double *avgdist_slice, double *avgdist_slice_old);
double update_combavgdistI(int m, int t, int p, int translice, int tran1, int tran2, double *d, double *d_old, double *avgdist_slice, double *avgdist_slice_old, double *avgdist_old, double *avgdist_cur);
void update_avgdist_sliceII(int n, int m, int p, int location1, int location2, double *d, double *d_old, double *avgdist_slice, double *avgdist_slice_old);
double update_combavgdistII(int m, int t, int p, int location1, int location2, double *d, double *d_old, double *avgdist_slice, double *avgdist_slice_old, double *avgdist_old, double *avgdist_cur);


void maximinSLHD(int *mRow, int *Col, int *nslice, int *npower, int *nstarts, int *IterMax, int *Total_Iter, int *design, double *measure, double *temp0, int *ntotal)
{
		
	double seed=rand();
    const int k=*Col;                 // define k: number of factors  <= change here
	const int m=*mRow;               // define m: number of runs in each slice <= change here
	const int t=*nslice;                // define t: number of slices    <= change here
	const int n=m*t;
	const int nsearch=*nstarts;             // number of different random starts
	const int p=*npower;
	const double tfac=0.95;
	const int Imax=*IterMax;            // maximum number of tries without improving xbest at each temperature for Simulated Annealing 
	int nImax=Min(5*n*(n-1)*k,Imax);   
	int max_itotal=*Total_Iter;               // maximum number of iterations for each random start
	
	
	int max_itotal1;
	int max_itotal2;
	if(t>1){ 
		max_itotal1=max_itotal*0.75;
		max_itotal2=max_itotal*0.25;
	}
	else
	{
		max_itotal1=max_itotal;
		max_itotal2=0;
	}
	   


	double t0;
	double xcrit;
	double critbest;
	double crittry;
	int itotal;
	double temp;
	int **xbest;
    int **xtry;
	int **x;

	// distance matrix
	int dim= (int)(n*(n-1)*0.5);
	double *d;
	d=new double[dim];

	double *d_old;
	d_old=new double[dim];
	for(int i=0;i<dim;i++)
	{
		d_old[i] = 0;
	}


	double *avgdist_cur=new double;
	*avgdist_cur=0;
	double *avgdist_old=new double;
	*avgdist_old=0;

	// average reciprocal interpoint distance for each slice
	double *avgdist_slice;
	avgdist_slice=new double[t];

	double *avgdist_slice_old;
	avgdist_slice_old=new double[t];
	for(int i=0;i<t;i++)
	{
		avgdist_slice_old[i] = 0;
	}
	
	

// Initialize xtry, xbest, and x;
	xtry	=new int*[k];
    xbest	=new int*[k];
	x		=new int*[k];
	for(int i=0;i<k;i++)
	{
		xbest[i] =new int[n];
		xtry[i]  =new int[n];
		x[i]     =new int[n];
	}

//////initialized the best design ////////////////
	seed=rand();
	xbest=SLHD(m,t,k,seed);
	distmatrix(xbest,n,k,d);
	critbest=combavgdist(m,t,p,d,avgdist_slice, avgdist_cur);
	

/////calculate starting temperature.
	double avgd2= k*n*(n+1)*pow(6,-1);
	double delta0= pow((avgd2-k),-0.5) - pow(avgd2,-0.5);
	t0=-delta0*pow((log(0.99)),(-1));
	*temp0=t0;

	
////Loop for different random starts//////////
	    itotal=0;
		for(int isearch=1;isearch<=nsearch;isearch++)
		{
			//////initial design ////////////////
			seed=rand();
			x=SLHD(m,t,k,seed);
			for(int n2=0;n2<k;n2++)
			{
				for(int n1=0;n1<n;n1++)
				{
				 *(*(xtry+n2)+n1)=*(*(x+n2)+n1);
				}
			}			
			distmatrix(xtry,n,k,d);
			xcrit=combavgdist(m,t,p,d,avgdist_slice,avgdist_cur);
			crittry=xcrit;

///////////initialize tempertures and counts///////////////////
			temp=t0;
			int ichange=1;
////////// variable temperature loop ////////////////////////
     		while(ichange==1)
			{
				ichange=0;
	//// constant temperature loop /////////////////////
				int ipert=1;
				while (ipert<nImax)
				{					
					if(itotal>max_itotal1) break;
					itotal=itotal+1;
       //////// switch to be tried is elements 
		 /////// change two component in a column ////////// 
					 int ind;
					 int translice;
					 int tran1;
					 int tran2;
					 seed=rand();
					 ind=rc((k-1),seed);
					 seed=rand();
					 translice=rc(t,seed);
					 seed=rand();
					 tran1=rc(m,seed);
					 seed=rand();
					 tran2=rc2(m,tran1,seed); 	
       /////////perturb x to xtry////////////////////
 					 *(*(xtry+ind+1)+translice*m+tran2)=*(*(x+ind+1)+translice*m+tran1);
				     *(*(xtry+ind+1)+translice*m+tran1)=*(*(x+ind+1)+translice*m+tran2);
       //////////////////////////////////////////////
					 update_distmatrix(xtry,n,ind+1,translice*m+tran1,translice*m+tran2,d,d_old);
					 crittry=update_combavgdistI(m, t, p, translice, tran1, tran2, d, d_old, avgdist_slice, avgdist_slice_old, avgdist_old, avgdist_cur);

					 
////// is xtry better than xbest? //////////////////////////////
					if(crittry<critbest)
					{
	//////////yes: replace x, xbest by xtry ; set iterp=1;ichange=1////////////////////////
						ichange=1;
						for(int nn2=0;nn2<k;nn2++)
						{
						   for(int nn1=0;nn1<n;nn1++)
						   {
							  *(*(xbest+nn2)+nn1)=*(*(xtry+nn2)+nn1); 
						   }
						}
						 *(*(x+ind+1)+translice*m+tran1)=*(*(xtry+ind+1)+translice*m+tran1);
						 *(*(x+ind+1)+translice*m+tran2)=*(*(xtry+ind+1)+translice*m+tran2);
						critbest=crittry;
						ipert=1;
						xcrit=crittry;
					}
					else
					{
//////////No:, increase ipert by 1. is xtry better than x?
						ipert=ipert+1;
			
						if(crittry<xcrit)
						{
			  ////// xtry is better than x; replace x by xtry ///////////
 						 *(*(x+ind+1)+translice*m+tran1)=*(*(xtry+ind+1)+translice*m+tran1);
						 *(*(x+ind+1)+translice*m+tran2)=*(*(xtry+ind+1)+translice*m+tran2);
						 ichange=1;
						 
						 xcrit=crittry;
		  //////////////////////////////////////////////////////////
					}
					else
					{
		 ///////// xtry is worst than x////////////////////////////
					    double delta1=crittry-xcrit;
					    double prob=exp(-delta1*pow(temp,(-1)));     
 					    seed=seed+isearch+ipert;
					    double q=runif(seed);
					     if(prob>=q)
						 {///// replce x by xtry by prob///////////
 				           	 *(*(x+ind+1)+translice*m+tran1)=*(*(xtry+ind+1)+translice*m+tran1);
				             *(*(x+ind+1)+translice*m+tran2)=*(*(xtry+ind+1)+translice*m+tran2);  
							// ichange=1;
							 xcrit=crittry;
						 }///////////////////////////////////
                         else 
						 {///// reset x try to x for the next pertubation		
 				           	 *(*(xtry+ind+1)+translice*m+tran1)=*(*(x+ind+1)+translice*m+tran1);
				             *(*(xtry+ind+1)+translice*m+tran2)=*(*(x+ind+1)+translice*m+tran2);  
							 revert_distmatrix(n,translice*m+tran1,translice*m+tran2,d,d_old);
							 *avgdist_cur=*avgdist_old;
	                         avgdist_slice[translice]=avgdist_slice_old[translice];
						 }////////////////////////////////////////// 
		 //////////////////////////////////////////////////////////
					}
				}
			}
	//// end of constant temperature loop ////////////
			temp=temp*tfac;
			
		}
///////// End of variable temperature loop///////////////////


	    }
/////end of search loop////////////////////////////







	int **lmatrix;
	lmatrix=new int*[k];
	for(int iii=0;iii<k;iii++)
	{
		lmatrix[iii]=new int[n];
	}
	int *loc;
	loc=new int[m];
	int itotal2=0;

	if(t>1){

		
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
		///Stage II Optimization

////Find the location design matrix for the best design found in Stage I

			
		for(int ck=0;ck<k;ck++){
		    for(int iii=0;iii<m;iii++){
			    loc[iii]=0;
		    }
			for(int cn=0;cn<n;cn++)
				for(int cm=0;cm<m;cm++){
					if(((cm*t+1)<=*(*(xbest+ck)+cn)) && (*(*(xbest+ck)+cn)<=(cm+1)*t)){
						*(*(lmatrix+ck)+cm*t+loc[cm])=cn;
						loc[cm]=loc[cm]+1;
					}
				}
		}
		



		///////////initialize tempertures and counts///////////////////
			temp=t0;
			int ichange=1;

			for(int nn2=0;nn2<k;nn2++)
			{
			    for(int nn1=0;nn1<n;nn1++)
				{
					*(*(x+nn2)+nn1)=*(*(xbest+nn2)+nn1); 
					*(*(xtry+nn2)+nn1)=*(*(xbest+nn2)+nn1);
				}
			}

			distmatrix(x,n,k,d);
			xcrit=combavgdist(m,t,p,d,avgdist_slice,avgdist_cur);

			int location1;
			int location2;
			
			int nslice1;
	        int nslice2;


////////// variable temperature loop ////////////////////////
			
     		while(ichange==1)
			{
				ichange=0;
	//// constant temperature loop /////////////////////
				int ipert=1;
				while (ipert<nImax)
				{					
					if(itotal2>max_itotal2) break;
					itotal2=itotal2+1;
       //////// switch to be tried is elements 
		 /////// change two component in a column ////////// 
					 int ind;
					 int tranm;
					 int tran1;
					 int tran2;
					 seed=rand();
					 ind=rc(k,seed);
					 seed=rand();
					 tranm=rc(m,seed);
					 seed=rand();
					 tran1=rc(t,seed);
					 seed=rand();
					 tran2=rc2(t,tran1,seed); 	

					 location1=*(*(lmatrix+ind)+tranm*t+tran1);
					 location2=*(*(lmatrix+ind)+tranm*t+tran2);

       /////////perturb x to xtry////////////////////					 
 					 *(*(xtry+ind)+location2)=*(*(x+ind)+location1);
				     *(*(xtry+ind)+location1)=*(*(x+ind)+location2);
       //////////////////////////////////////////////
					
					 //distmatrix(xtry,n,k,d);
					 update_distmatrix(xtry,n,ind,location1,location2,d,d_old);
					 crittry=update_combavgdistII(m, t, p, location1, location2, d, d_old, avgdist_slice, avgdist_slice_old, avgdist_old, avgdist_cur);

					 ////// is xtry better than xbest? //////////////////////////////
					if(crittry<critbest)
					{
	//////////yes: replace x, xbest by xtry ; set iterp=1;ichange=1////////////////////////
						ichange=1;
						for(int nn2=0;nn2<k;nn2++)
						{
						   for(int nn1=0;nn1<n;nn1++)
						   {
							  *(*(xbest+nn2)+nn1)=*(*(xtry+nn2)+nn1); 
						   }
						}
						 *(*(x+ind)+location1)=*(*(xtry+ind)+location1);
						 *(*(x+ind)+location2)=*(*(xtry+ind)+location2);
						critbest=crittry;
						ipert=1;
						xcrit=crittry;

					}
					else
					{
//////////No:, increase ipert by 1. is xtry better than x?
						ipert=ipert+1;
			
						if(crittry<xcrit)
						{
			  ////// xtry is better than x; replace x by xtry ///////////
 						 *(*(x+ind)+location1)=*(*(xtry+ind)+location1);
						 *(*(x+ind)+location2)=*(*(xtry+ind)+location2);
						 ichange=1;
						 
						 xcrit=crittry;
		  //////////////////////////////////////////////////////////
					}
					else
					{
		 ///////// xtry is worst than x////////////////////////////
					    double delta1=crittry-xcrit;
					    double prob=exp(-delta1*pow(temp,(-1)));     
 					    seed=seed+ipert;
					    double q=runif(seed);
					     if(prob>=q)
						 {///// replce x by xtry by prob///////////
 				           	 *(*(x+ind)+location1)=*(*(xtry+ind)+location1);
				             *(*(x+ind)+location2)=*(*(xtry+ind)+location2);  
							// ichange=1;
							 xcrit=crittry;
						 }///////////////////////////////////
                         else 
						 {///// reset x try to x for the next pertubation
 				           	 *(*(xtry+ind)+location1)=*(*(x+ind)+location1);
				             *(*(xtry+ind)+location2)=*(*(x+ind)+location2);  
							 revert_distmatrix(n,location1,location2,d,d_old);
							 *avgdist_cur=*avgdist_old;
							 nslice1= location1/m;
	                         nslice2= location2/m;
	                         avgdist_slice[nslice1]=avgdist_slice_old[nslice1];
							 avgdist_slice[nslice2]=avgdist_slice_old[nslice2];
	                   
						 }////////////////////////////////////////// 
		 //////////////////////////////////////////////////////////
					}
				}
			}
	//// end of constant temperature loop ////////////
			temp=temp*tfac;
			
		}
///////// End of variable temperature loop///////////////////

	}


	/// Output the best design for the current loop
        for(int ii=0;ii<n;ii++)
		{
			for(int jj=0;jj<k;jj++)
			{
				*(design+ii*k+jj)=*(*(xbest+jj)+ii);
			}
		}

		
	
	*measure=critbest;
	*ntotal=itotal+itotal2;
					
	

 // delete xtry, xbest, and x;
	for(int i=0; i<k; i++)
	{
		delete [] xbest[i];
		delete [] xtry[i];
		delete [] x[i];
		delete [] lmatrix[i];
	}

	delete []xbest;
	delete []xtry;
	delete []x;
	delete []lmatrix;
	delete []loc;
	delete []avgdist_slice_old;
	delete []avgdist_slice;
	delete []d;
	delete []d_old;
}




///////////////////////////////////////

int **LHD(int n, int k, double seed)
{
	int te;
	int **LHD;
	int *r;
	r=new int [n];
    LHD=new int*[k];
	for(int iii=0;iii<k;iii++)
	{
		LHD[iii]=new int[n];
	}
     for(int j2=0;j2<n;j2++)  // first dimention is 1,2,3,.....
	{
		*(*(LHD+0)+j2)=(j2+1);
	}

	for(int j1=0;j1<(k-1);j1++)  
	{
	   for(int cc=0;cc<n;cc++)
	   {
		*(r+cc)=cc+1;
	   }
	   for(int c=0;c<n;c++)
	   {
			seed=seed+j1*c;
			seed=seed+10;
			te=rc(n-c,seed);
			*(*(LHD+j1+1)+c)=*(r+te);

		    for(int c1=0;c1<(n-c);c1++)
			{
			    if(c1>=te)
				{
				*(r+c1)=*(r+c1+1);
				}
			}
	   }
	}

   return(LHD);
   for(int iiii=0;iiii<k;iiii++)
   {
	delete [] LHD[iiii];
   }
   delete [] LHD;
   delete [] r;
}



int **SLHD(int m, int t, int k, double seed)
{
	int te;
	int **SLHD;
	int *r;
	int n=m*t;
	r=new int [n];
    SLHD=new int*[k];
	for(int iii=0;iii<k;iii++)
	{
		SLHD[iii]=new int[n];
	}
     for(int js=0;js<t;js++)  // first dimention is 1,2,3,.....
	{
		for(int j2=0;j2<m;j2++){					
			*(*(SLHD+0)+m*js+j2)=(j2+1);
		}
	}

	for(int j1=0;j1<(k-1);j1++)  
	{
	   for(int jss=0;jss<t;jss++){

	   for(int cc=0;cc<m;cc++)
	   {
		*(r+cc)=cc+1;
	   }
	   for(int c=0;c<m;c++)
	   {
			seed=seed+j1*c;
			seed=seed+10+jss;
			te=rc(m-c,seed);
			*(*(SLHD+j1+1)+jss*m+c)=*(r+te);

		    for(int c1=0;c1<(m-c);c1++)
			{
			    if(c1>=te)
				{
				*(r+c1)=*(r+c1+1);
				}
			}
	    }

	    }
	}

	int **SLHDS;
	SLHDS=new int*[k];
	for(int iii=0;iii<k;iii++)
	{
		SLHDS[iii]=new int[n];
	}

	int xsubs;
	for(int j3=0;j3<k;j3++){
		for(int j5=0;j5<m;j5++){
			xsubs=j5*t+1;
			for(int j4=0;j4<n;j4++){				
				if(*(*(SLHD+j3)+j4)==(j5+1)){
					*(*(SLHDS+j3)+j4)=xsubs;
					xsubs++;
				}
			}
		}
	}


   return(SLHDS);
   for(int iiii=0;iiii<k;iiii++)
   {
	delete [] SLHD[iiii];
	delete [] SLHDS[iiii];
   }
   delete [] SLHD;
   delete [] SLHDS;
   delete [] r;
}



  
int rc2(int n,int del,double seed)
{
   int rctwo;
   
   rctwo= rc( n-1, seed);
   if (rctwo >= del)  rctwo++;

   return(rctwo);
}
 
int rc(int n,double seed ) // choose randomly from 0 to (n-1)
{
   int r;
   double u; 
   u = (double) (rand()*autoseed(seed)%MAX)/MAX;
   r=(int)(n*u);
   return(r);
} 

double runif(double seed)
{
   double runif;
   runif= (double) (rand()*autoseed(seed)%MAX)/MAX;
   return(runif);
}


void distmatrix(int **A, int n, int k, double *d) // To compute the interpoint distance matrix
{
    const int dim= (int)(n*(n-1)*0.5);
	for(int i=0;i<dim;i++)
	{
		d[i] = 0;
	}
 	int count=0;

   	for(int k1=0;k1<(n-1);k1++)
	{
		for(int k2=(k1+1);k2<n;k2++)
		{  			
			for(int k3=0;k3<k;k3++)
			{
				d[count] += pow((*(*(A+k3)+k1)-*(*(A+k3)+k2)),2);
			}

			d[count]=pow(d[count],0.5);
			count++;
		}
	}

}



void avgdist(int n, int p, double *d, double *avgdist_cur) // To compute the average reciprocal interpoint distance
{
    const int dim= (int)(n*(n-1)*0.5);
	double avgdist=0;

	for(int i=0; i<dim; i++)
	{
		avgdist += pow(d[i],(-p));
	}
	avgdist=avgdist*pow(dim,-1);
	*avgdist_cur=pow(avgdist,(pow(p,(-1))));

}




double combavgdist(int m, int t, int p, double *d, double *avgdist_slice, double *avgdist_cur) // To compute the combined average reciprocal interpoint distance
{
	if(t>1){
    const int dim_slice= (int)(m*(m-1)*0.5);
	const int n=m*t;
	double combavgdist=0;
	for(int i=0;i<t;i++)
	{
		avgdist_slice[i] = 0;
	}

	int count=0;
	for(int ns=0;ns<t;ns++){
		for(int nn1=ns*m;nn1<((ns+1)*m-1);nn1++){
			for(int nn2=nn1+1;nn2<(ns+1)*m;nn2++){
				count=(int) nn2+1-0.5*pow(nn1+1,2)+(n-0.5)*(nn1+1)-n-1;
				avgdist_slice[ns] += pow(d[count],(-p));		
			}
		}
		avgdist_slice[ns]=avgdist_slice[ns]*pow(dim_slice,-1);
	    avgdist_slice[ns]=pow(avgdist_slice[ns],(pow(p,(-1))));

	}

	for(int iii=0;iii<t;iii++){
		combavgdist=combavgdist+avgdist_slice[iii];
	}
	
	avgdist(n,p,d,avgdist_cur);
	combavgdist=(*avgdist_cur+combavgdist*pow(t,-1))*pow(2,-1);
	
	return(combavgdist);

	}
	
	else
	{
        const int n=m*t;
	    double combavgdist=0;
		avgdist(n,p,d,avgdist_cur);
	    combavgdist=*avgdist_cur;
		return(combavgdist);
	}
   
}






void update_distmatrix(int **A, int n, int col, int selrow1, int selrow2, double *d, double *d_old) // To update the interpoint distance matrix
{
    int row1=Min(selrow1,selrow2);
	int row2=Max(selrow1,selrow2);
	double s=0;
	int position1,position2;
	
	if(row1>0){
	for(int h=0;h<row1;h++) //h<row1<row2
	{ 
		s=pow((*(*(A+col)+row2)-*(*(A+col)+h)),2)-pow((*(*(A+col)+row1)-*(*(A+col)+h)),2);
		position1=(int) row1+1-pow(h+1,2)*0.5+(n-0.5)*(h+1)-n-1;
		position2=(int) row2+1-pow(h+1,2)*0.5+(n-0.5)*(h+1)-n-1;
		d_old[position1] = d[position1];
		d_old[position2] = d[position2];
		d[position1] = pow( pow(d[position1],2)-s, 0.5);
		d[position2] = pow( pow(d[position2],2)+s, 0.5);
		
	}
	}
	
	for(int h=(row1+1);h<row2;h++) //row1<h<row2
	{ 
		s=pow((*(*(A+col)+row2)-*(*(A+col)+h)),2)-pow((*(*(A+col)+row1)-*(*(A+col)+h)),2);
		position1=(int) h+1-pow(row1+1,2)*0.5+(n-0.5)*(row1+1)-n-1;
		position2=(int) row2+1-pow(h+1,2)*0.5+(n-0.5)*(h+1)-n-1;
		d_old[position1] = d[position1];
		d_old[position2] = d[position2];
		d[position1] = pow( pow(d[position1],2)-s, 0.5);
		d[position2] = pow( pow(d[position2],2)+s, 0.5);
	}
	
	if(row2<(n-1)){
	for(int h=(row2+1);h<n;h++) //row1<row2<h
	{ 
		s=pow((*(*(A+col)+row2)-*(*(A+col)+h)),2)-pow((*(*(A+col)+row1)-*(*(A+col)+h)),2);
		position1=(int) h+1-pow(row1+1,2)*0.5+(n-0.5)*(row1+1)-n-1;
		position2=(int) h+1-pow(row2+1,2)*0.5+(n-0.5)*(row2+1)-n-1;
		d_old[position1] = d[position1];
		d_old[position2] = d[position2];
		d[position1] = pow( pow(d[position1],2)-s, 0.5);
		d[position2] = pow( pow(d[position2],2)+s, 0.5);
	}
	}
	
}




void revert_distmatrix(int n, int selrow1, int selrow2, double *d, double *d_old) // To revert the interpoint distance matrix
{
    int row1=Min(selrow1,selrow2);
	int row2=Max(selrow1,selrow2);
	int position1,position2;

	if(row1>0){
	for(int h=0;h<row1;h++) //h<row1<row2
	{ 
		position1=(int) row1+1-pow(h+1,2)*0.5+(n-0.5)*(h+1)-n-1;
		position2=(int) row2+1-pow(h+1,2)*0.5+(n-0.5)*(h+1)-n-1;
		d[position1] = d_old[position1];
		d[position2] = d_old[position2];
		
	}
	}
	
	for(int h=(row1+1);h<row2;h++) //row1<h<row2
	{ 
		position1=(int) h+1-pow(row1+1,2)*0.5+(n-0.5)*(row1+1)-n-1;
		position2=(int) row2+1-pow(h+1,2)*0.5+(n-0.5)*(h+1)-n-1;
		d[position1] = d_old[position1];
		d[position2] = d_old[position2];
	}
	
	if(row2<(n-1)){
	for(int h=(row2+1);h<n;h++) //row1<row2<h
	{ 
		position1=(int) h+1-pow(row1+1,2)*0.5+(n-0.5)*(row1+1)-n-1;
		position2=(int) h+1-pow(row2+1,2)*0.5+(n-0.5)*(row2+1)-n-1;
		d[position1] = d_old[position1];
		d[position2] = d_old[position2];
	}
	}

}






void update_avgdist(int n, int p, int selrow1, int selrow2, double *d, double *d_old, double *avgdist_old, double *avgdist_cur) // To update the average reciprocal interpoint distance
{
    *avgdist_old=*avgdist_cur;
	const int dim= (int)(n*(n-1)*0.5);

	double avgdist0=0;

	for(int i=0; i<dim; i++)
	{
		avgdist0 += pow(d[i],(-p));
	}
	avgdist0=avgdist0*pow(dim,-1);
	avgdist0=pow(avgdist0,(pow(p,(-1))));

	*avgdist_cur=avgdist0;
	

}






void update_avgdist_sliceI(int n, int m, int p, int translice, int tran1, int tran2, double *d, double *d_old, double *avgdist_slice, double *avgdist_slice_old) // To update the average reciprocal interpoint distance for each slice
{
    const int dim_slice= (int)(m*(m-1)*0.5);
	avgdist_slice_old[translice]=avgdist_slice[translice];

	double sliceavgdist0=0;
	int pos;
	for(int nn2=translice*m;nn2<(translice*m+m-1);nn2++){
			for(int nn1=(nn2+1);nn1<(translice*m+m);nn1++){
				pos=nn1+1-pow(nn2+1,2)*0.5+(nn2+1)*(n-0.5)-n-1;
				sliceavgdist0 += pow(d[pos],(-p));
			}
	}
	sliceavgdist0=sliceavgdist0*pow(dim_slice,-1);
	sliceavgdist0=pow(sliceavgdist0,pow(p,-1));
	avgdist_slice[translice]=sliceavgdist0;

}


double update_combavgdistI(int m, int t, int p, int translice, int tran1, int tran2, double *d, double *d_old, double *avgdist_slice, double *avgdist_slice_old, double *avgdist_old, double *avgdist_cur)
{
	if(t>1){
	int n=m*t;
	double combavgdist=0;
	update_avgdist_sliceI(n, m, p, translice, tran1, tran2, d, d_old, avgdist_slice, avgdist_slice_old);
	
	int selrow1=translice*m+tran1;
	int selrow2=translice*m+tran2;
	update_avgdist(n, p, selrow1, selrow2, d, d_old, avgdist_old, avgdist_cur);
	
    for(int iii=0;iii<t;iii++){
		combavgdist=combavgdist+avgdist_slice[iii];
	}
	
	combavgdist=(*avgdist_cur+combavgdist*pow(t,-1))*pow(2,-1);
	
	return(combavgdist);
	}

	else
	{
	int n=m*t;
	double combavgdist=0;
		
	int selrow1=translice*m+tran1;
	int selrow2=translice*m+tran2;
	update_avgdist(n, p, selrow1, selrow2, d, d_old, avgdist_old, avgdist_cur);
		
	combavgdist=*avgdist_cur;
	
	return(combavgdist);

	}
}





void update_avgdist_sliceII(int n, int m, int p, int location1, int location2, double *d, double *d_old, double *avgdist_slice, double *avgdist_slice_old) // To update the average reciprocal interpoint distance for each slice
{
    const int dim_slice= (int)(m*(m-1)*0.5);
	int nslice1= location1/m;
	int nslice2= location2/m;

	avgdist_slice_old[nslice1]=avgdist_slice[nslice1];
	avgdist_slice_old[nslice2]=avgdist_slice[nslice2];

	int pos;

	double sliceavgdist01=0;	
	for(int nn2=nslice1*m;nn2<(nslice1*m+m-1);nn2++){
			for(int nn1=(nn2+1);nn1<(nslice1*m+m);nn1++){
				pos=nn1+1-pow(nn2+1,2)*0.5+(nn2+1)*(n-0.5)-n-1;
				sliceavgdist01 += pow(d[pos],(-p));
			}
	}
	sliceavgdist01=sliceavgdist01*pow(dim_slice,-1);
	sliceavgdist01=pow(sliceavgdist01,pow(p,-1));
	avgdist_slice[nslice1]=sliceavgdist01;

	double sliceavgdist02=0;	
	for(int nn2=nslice2*m;nn2<(nslice2*m+m-1);nn2++){
			for(int nn1=(nn2+1);nn1<(nslice2*m+m);nn1++){
				pos=nn1+1-pow(nn2+1,2)*0.5+(nn2+1)*(n-0.5)-n-1;
				sliceavgdist02 += pow(d[pos],(-p));
			}
	}
	sliceavgdist02=sliceavgdist02*pow(dim_slice,-1);
	sliceavgdist02=pow(sliceavgdist02,pow(p,-1));
	avgdist_slice[nslice2]=sliceavgdist02;

}


  

double update_combavgdistII(int m, int t, int p, int location1, int location2, double *d, double *d_old, double *avgdist_slice, double *avgdist_slice_old, double *avgdist_old, double *avgdist_cur)
{
	int n=m*t;
	double combavgdist=0;
	update_avgdist_sliceII(n, m, p, location1, location2, d, d_old, avgdist_slice, avgdist_slice_old);
	
	update_avgdist(n, p, location1, location2, d, d_old, avgdist_old, avgdist_cur);
	
    for(int iii=0;iii<t;iii++){
		combavgdist=combavgdist+avgdist_slice[iii];
	}
	
	combavgdist=(*avgdist_cur+combavgdist*pow(t,-1))*pow(2,-1);
	
	return(combavgdist);
}






}
