/*

  Adam Majewski
  adammaj1 aaattt o2 dot pl  // o like oxygen not 0 like zero 
  
  
  console program in c programing language 
===============================================================





  
  ==============================================
  
  
  Structure of a program or how to analyze the program 
  
  
  ============== Image X ========================
  
  DrawImageOfX -> DrawPointOfX -> ComputeColorOfX 
  
  first 2 functions are identical for every X
  check only last function =  ComputeColorOfX
  which computes color of one pixel !
  
  

   
  ==========================================

  
  ---------------------------------
  indent d.c 
  default is gnu style 
  -------------------



  c console progam 
  
	export  OMP_DISPLAY_ENV="TRUE"	
  	gcc d.c -lm -Wall -march=native -fopenmp
  	time ./a.out > b.txt


  gcc d.c -lm -Wall -march=native -fopenmp


  time ./a.out

  time ./a.out >a.txt
  
  
  ./g.sh

  ----------------------
  
 real	0m19,809s
user	2m26,763s
sys	0m0,161s


  

*/

#include <stdio.h>
#include <stdlib.h>		// malloc
#include <string.h>		// strcat
#include <math.h>		// M_PI; needs -lm also
#include <complex.h> 		// complex numbers : https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c
#include <omp.h>		// OpenMP

/* --------------------------------- global variables and consts ------------------------------------------------------------ */


int NumberOfImages = 0;


//FunctionType
typedef enum  {LSM , DEM, Unknown, BD, MBD , SAC, DLD, 
		ND, NP, POT
		
		} FunctionTypeT; 
// FunctionTypeT FunctionType;

int PlaneInversion = 0; // boolean 1 = w = 1/z  plane; 0 = z plane 


// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
static unsigned int ixMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int ixMax;	//
static unsigned int iWidth;	// horizontal dimension of array

static unsigned int iyMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int iyMax;	//

static unsigned int iHeight = 5000;	//  
// The size of array has to be a positive constant integer 
static unsigned int iSize;	// = iWidth*iHeight; 

// ----------memmory 1D arrays ==================
// unsigned char = for 1 byte ( 8 bit) colors 
unsigned char *data;
unsigned char *edge;
unsigned char *edge2;
 

// unsigned int i; // var = index of 1D array
//static unsigned int iMin = 0; // Indexes of array starts from 0 not 1
static unsigned int iMax;	// = i2Dsize-1  = 
// The size of array has to be a positive constant integer 
// unsigned int i1Dsize ; // = i2Dsize  = (iMax -iMin + 1) =  ;  1D array with the same size as 2D array





// see SetZPlane

double radius = 1.4; 
complex double center = 0.0;
double  DisplayAspectRatio  = 1.0; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)



// z plane = dynamic plane
double ZxMin ;	//-0.05;
double ZxMax ;	//0.75;
double ZyMin ;	//-0.1;
double ZyMax ;	//0.7;


double PixelWidth;	// =(ZxMax-ZxMin)/ixMax;
double PixelHeight;	// =(ZyMax-ZyMin)/iyMax;
double ratio;





// w plane = 1/z plane
double WxMin = - 3;	//-0.05;
double WxMax = 3;	//0.75;
double WyMin = -3;	//-0.1;
double WyMax = 3;	//0.7;
double wPixelWidth;	// =(WxMax-WxMin)/ixMax;

double wPixelHeight;	// =(WyMax-WyMin)/iyMax;

// complex numbers of parametr plane 
double complex c;		// parameter of function fc(z)=z^2 + c




static unsigned long int iterMax = 1000000;	//iHeight*100;
const long int iterMax_LSM = 255;
const int iterMax_DLD = 200; // N in wiki = fixed number : maximal number of iterations
const int iterMax_pot = 400; // potential 

double ER = 200.0;		// EscapeRadius for bailout test 
double EscapeRadius=1000000; // = ER big !!!!
double ER_LSM ; // see GiveER_LSM  // 27.764 =  manually find value such that level curves of escape time cross critical point and it's  preimages
double ER_DLD ; // see GiveER_LSM  // 27.764 =  manually find value such that level curves of escape time cross critical point and it's  preimages
double ER_NP = 100.0; 
double ER_pot = 100000.0;  // sqrt(1e24)

double loger; // = log(ER_LSM); // for texture
static double TwoPi=2.0*M_PI; // texture
double MaxFinalRadius;


// SAC/J
double lnER; // ln(ER)
int i_skip = 2; // exclude (i_skip+1) elements from average
unsigned char s = 7; // stripe density

double BoundaryWidth = 3.0; // % of image width  
double distanceMax; //distanceMax = BoundaryWidth*PixelWidth;



//  ------------- DLD  ----------------------

double p = 0.180; //0.01444322;		//
// DLD colors
//double me = 1.0;
double mi = 0.9;

// potential
double MaxImagePotential;



/* colors = shades of gray from 0 to 255 */
unsigned char iColorOfExterior = 250;
unsigned char iColorOfInterior = 200;
unsigned char iColorOfInterior1 = 210;
unsigned char iColorOfInterior2 = 180;
unsigned char iColorOfBoundary = 0;
unsigned char iColorOfUnknown = 30;





/* ------------------------------------------ functions -------------------------------------------------------------*/

/**
 * Find maximum between two numbers.
 https://codeforwin.org/2016/02/c-program-to-find-maximum-and-minimum-using-functions.html
 */
double max(double n1, double n2)
{
    return (n1 > n2 ) ? n1 : n2;
}



//---------------------

double min(double n1, double n2)
{
    return (n1 < n2 ) ? n1 : n2;
}


double clip(double d){

	return (d> 1.0) ? 1.0 : d;
}



double frac(double d){

	double fraction = d - ((long)d);
	return fraction;
}




//------------------complex numbers -----------------------------------------------------




double c_arg(complex double z)
{
 double arg;
 arg = carg(z);
 if (arg<0.0) arg+= TwoPi ; 
 return arg; 
}

double c_turn(complex double z)
{
 double arg;
 arg = c_arg(z);
 return arg/TwoPi; 
}





// from screen to world coordinate ; linear mapping
// uses global cons
double GiveZx ( int ix)
{
  return (ZxMin + ix * PixelWidth);
}

// uses globaal cons
double GiveZy (int iy) {
  return (ZyMax - iy * PixelHeight);
}				// reverse y axis


complex double GiveZ( int ix, int iy){
  double Zx = GiveZx(ix);
  double Zy = GiveZy(iy);
	
  return Zx + Zy*I;
	
	


}


// from screen to world coordinate ; linear mapping
// uses global cons
double GiveWx ( int ix)
{
  return (WxMin + ix * wPixelWidth);
}

// uses globaal cons
double GiveWy (int iy) {
  return (WyMax - iy * wPixelHeight);
}				// reverse y axis


complex double GiveW( int ix, int iy){
  double Wx = GiveWx(ix);
  double Wy = GiveWy(iy);
	
  return Wx + Wy*I;
	
	


}




int SetZPlane(complex double center, double radius, double a_ratio){

  ZxMin = creal(center) - radius*a_ratio;	
  ZxMax = creal(center) + radius*a_ratio;	//0.75;
  ZyMin = cimag(center) - radius;	// inv
  ZyMax = cimag(center) + radius;	//0.7;
  return 0;

}







// ****************** DYNAMICS = trap tests ( target sets) ****************************



// find such ER for LSM/J that level curves croses critical point and it's preimages
double GiveER(int i_Max){

	complex double z= 0.0; // criical point
	int i;
	  // critical point escapes very fast here. Higher valus gives infinity
	for (i=0; i< i_Max; ++i ){
		z=z*z +c; 
	 
	 }
	 
	 return cabs(z);
	
	
}


double GiveMaxFinalRadius(){

	complex double z = ZxMax + ZyMax*I;
	double r = log(cabs(z))/loger - 1.0; // final_z_abs in not in [0,1]

	return r;
	}
	
double GiveNormalizedFinalRadius(complex double z){

	double FinalRadius = log(cabs(z))/loger - 1.0; // final_z_abs in not in [0,1]
	return (FinalRadius/ MaxFinalRadius);

}







// bailout test
// z escapes when 
// abs(z)> ER or cabs2(z)> ER2 
// https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Julia_set#Boolean_Escape_time
// this function is not used !!!! dead code 

int Escapes(complex double z){
 // here target set (trap) is the exterior  circle with radsius = ER ( EscapeRadius) 
  // with ceter = origin z= 0
  // on the Riemann sphere it is a circle with point at infinity as a center  
   
  if (cabs(z)>ER) return 1;
  return 0;
}








/* -----------  array functions = drawing -------------- */

/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
unsigned int Give_i (unsigned int ix, unsigned int iy)
{
  return ix + iy * iWidth;
}


// ***********************************************************************************************
// ********************** edge detection usung Sobel filter ***************************************
// ***************************************************************************************************

// from Source to Destination
int ComputeBoundaries(unsigned char S[], unsigned char D[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
  /* sobel filter */
  unsigned char G, Gh, Gv; 
  // boundaries are in D  array ( global var )
 
  // clear D array
  memset(D, iColorOfExterior, iSize*sizeof(*D)); // 
 
  // printf(" find boundaries in S array using  Sobel filter\n");   
#pragma omp parallel for schedule(dynamic) private(i,iY,iX,Gv,Gh,G) shared(iyMax,ixMax)
  for(iY=1;iY<iyMax-1;++iY){ 
    for(iX=1;iX<ixMax-1;++iX){ 
      Gv= S[Give_i(iX-1,iY+1)] + 2*S[Give_i(iX,iY+1)] + S[Give_i(iX-1,iY+1)] - S[Give_i(iX-1,iY-1)] - 2*S[Give_i(iX-1,iY)] - S[Give_i(iX+1,iY-1)];
      Gh= S[Give_i(iX+1,iY+1)] + 2*S[Give_i(iX+1,iY)] + S[Give_i(iX-1,iY-1)] - S[Give_i(iX+1,iY-1)] - 2*S[Give_i(iX-1,iY)] - S[Give_i(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= Give_i(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {D[i]=255;} /* background */
      else {D[i]=0;}  /* boundary */
    }
  }
 
   
 
  return 0;
}



// copy from Source to Destination
int CopyBoundaries(unsigned char S[],  unsigned char D[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
 
 
  fprintf(stderr, "copy boundaries from S array to D array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {i= Give_i(iX,iY); if (S[i]==0) D[i]=0;}
 
 
 
  return 0;
}

// =============================  tests ============================================


// Check Orientation of z-plane image : mark first quadrant of complex plane 
// it should be in the upper right position
// uses global var :  ...
int CheckZPlaneOrientation(unsigned char A[] )
{
 
	double Zx, Zy; //  Z= Zx+ZY*i;
	unsigned i; /* index of 1D array */
	unsigned int ix, iy;		// pixel coordinate 
	
	fprintf(stderr, "compute image CheckOrientation\n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy, i, Zx, Zy) shared(A, ixMax , iyMax) 
	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix){
    			// from screen to world coordinate 
    			Zy = GiveZy(iy);
    			Zx = GiveZx(ix);
	  		i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
	  		if (Zx>0 && Zy>0) A[i]=255-A[i];   // check the orientation of Z-plane by marking first quadrant */
    		}
    	}
   
   
  	return 0;
}




int IsInsideWWindow(complex double w){

	if ( creal(w) < WxMax && creal(w) > WxMin &&
	     cimag(w) < WyMax && cimag(w) > WyMin) {return 1;}
	
	
	return 0;
	
		


}


/*

 Array A should have image of z-plane ( not w-plane) 
 compare of image of array A unchanged
 image of w window shows part of z window and outside of z-window
 
 "Note that the flower-shaped hole in the center is originally the edge boundary of the grid."
 http://xahlee.info/SpecialPlaneCurves_dir/Inversion_dir/inversion.html
 
 https://mathworld.wolfram.com/ConformalMapping.html
 http://home.iitk.ac.in/~saiwal/engineering/complex-mappings/
 
 
*/ 
int ShowWWindowOnZWindow(unsigned char A[] )
{
 
	complex double z;
	//double Zx, Zy; //  Z= Zx+ZY*i;
	complex double w;
	unsigned i; /* index of 1D array */
	unsigned int ix, iy;		// pixel coordinate 
	
	fprintf(stderr, "compute image ShowWWindowOnZWindow\n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy, i, w, z) shared(A, ixMax , iyMax) 
	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix){
    			
    			z = GiveZ(ix,iy); // from screen to world coordinate 
    			w = 1/z; // invert complex plane z 
	  		if (IsInsideWWindow(w)){
	  			i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
	  			 A[i]=255-A[i];   // marking w window on z window
	  			 }
    		}
    	}
   
   
  	return 0;
}


// ------------------------------------------------------------------------------




int IsInsideZWindow(complex double z){

	if ( creal(z) < ZxMax && creal(z) > ZxMin &&
	     cimag(z) < ZyMax && cimag(z) > ZyMin) {return 1;}
	
	
	return 0;
	
		


}


/*

 Array A should have image of w-plane ( not z-plane) 
 compare of image of array A unchanged
 image of w window shows part of z window and outside of z-window
 
 "Note that the flower-shaped hole in the center is originally the edge boundary of the grid."
 http://xahlee.info/SpecialPlaneCurves_dir/Inversion_dir/inversion.html
 
 https://mathworld.wolfram.com/ConformalMapping.html
 http://home.iitk.ac.in/~saiwal/engineering/complex-mappings/
 
 
*/ 
int ShowZWindowOnWWindow(unsigned char A[] )
{
 
	complex double z;
	//double Zx, Zy; //  Z= Zx+ZY*i;
	complex double w;
	unsigned i; /* index of 1D array */
	unsigned int ix, iy;		// pixel coordinate 
	
	fprintf(stderr, "compute image ShowZWindowOnWWindow\n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy, i, w, z) shared(A, ixMax , iyMax) 
	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix){
    			
    			w = GiveW(ix,iy); // from screen to world coordinate 
    			z = 1/w; // invert complex plane z 
	  		if (IsInsideZWindow(z)){
	  			i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
	  			 A[i]=255-A[i];   // marking w window on z window
	  			 }
    		}
    	}
   
   
  	return 0;
}










// ***************************************************************************************************************************
// ************************** DEM/J*****************************************
// ****************************************************************************************************************************

unsigned char ComputeColorOfDEMJ(complex double z){
// https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Julia_set#DEM.2FJ


  
  int nMax = iterMax;
  complex double dz = 1.0; //  is first derivative with respect to z.
  double distance;
  double cabsz;
	
  int n;

  for (n=0; n < nMax; n++){ //forward iteration
	cabsz = cabs(z);
    	if (cabsz > 1e60 || cabs(dz)> 1e60) break; // big values 
    	if (cabsz< PixelWidth) return iColorOfInterior; // falls into finite attractor = interior
  			
    dz = 2.0*z * dz; 
    z = z*z +c ; /* forward iteration : complex quadratic polynomial */ 
  }
  
  
  distance = 2.0 * cabsz* log(cabsz)/ cabs(dz);
  if (distance <distanceMax) return iColorOfBoundary; // distanceMax = BoundaryWidth*PixelWidth;
  // else
  
  return iColorOfExterior;

 
}






// ***************************************************************************************************************************
// ************************** only boundary by  DEM/J*****************************************
// ****************************************************************************************************************************

unsigned char ComputeColorOfDEMJ_boundary(complex double z){
// https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/Julia_set#DEM.2FJ


  
  int nMax = iterMax;
  complex double dz = 1.0; //  is first derivative with respect to z.
  double distance;
  double cabsz;
	
  int n;

  for (n=0; n < nMax; n++){ //forward iteration
	cabsz = cabs(z);
    	if (cabsz > 1e60 || cabs(dz)> 1e60) break; // big values 
    	if (cabsz< PixelWidth) return iColorOfInterior; // falls into finite attractor = interior
  			
    dz = 2.0*z * dz; 
    z = z*z +c ; /* forward iteration : complex quadratic polynomial */ 
  }
  
  
  distance = 2.0 * cabsz* log(cabsz)/ cabs(dz);
  if (distance <distanceMax) return iColorOfBoundary; // distanceMax = BoundaryWidth*PixelWidth;
  // else
  
  return iColorOfExterior;

 
}



// plots raster point (ix,iy) 
int DrawPointOfDEMJ_boundary (unsigned char A[], int PlaneInversion, int ix, int iy, unsigned char iColor0)
{
  int i;			/* index of 1D array */
  unsigned char iColor;
  complex double z;


  i = Give_i (ix, iy);		/* compute index of 1D array from indices of 2D array */
  
  if (PlaneInversion)
  	{ 
  		complex double w;
  		w = GiveW(ix,iy);
  		z = 1/w;
  	}
  	else {  z = GiveZ(ix,iy);}
  	
  iColor = ComputeColorOfDEMJ_boundary(z);
  if (iColor == iColorOfBoundary) // check if it is boundary
  	{ A[i] = iColor0 ;} // draw only boundary without changing other parts using color iColor0		
  
  return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int DrawImageOfDEMJ_boundary (unsigned char A[], int PlaneInversion, const unsigned char iColor)
{
  unsigned int ix, iy;		// pixel coordinate 

  	fprintf(stderr, "compute image DEM boundary\n");
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax)
  	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			DrawPointOfDEMJ_boundary(A, PlaneInversion, ix, iy, iColor);	//  
  }

  return 0;
}










// ***************************************************************************************************************************
// ************************** Unknown: boundary and slow dynamics *****************************************
// ****************************************************************************************************************************

unsigned char ComputeColorOfUnknown(complex double z){



  
  int nMax = 20; // very low value
  
  double cabsz;
	
  int n;

  for (n=0; n < nMax; n++){ //forward iteration
	cabsz = cabs(z);
    	if (cabsz > 10000000000*ER )  return iColorOfExterior; // big values
    	if (cabsz < (PixelWidth/100)) return iColorOfInterior; // falls into finite attractor = interior
  			
    
    z = z*z +c ; /* forward iteration : complex quadratic polynomial */ 
  }
  
  
  
  
  //printf("found \n");
  return iColorOfUnknown;

 
}




// ***************************************************************************************************************************
// ************************** LSM/J*****************************************
// ****************************************************************************************************************************




int GiveEscapeTime(complex double z){


	int nMax = iterMax_LSM;
	double cabsz;
	int n;

  	for (n=0; n < nMax; n++){ //forward iteration
		cabsz = cabs(z);
    		if (cabsz > ER_LSM) break; // esacping
    		//if (cabsz< PixelWidth) break; // fails into finite attractor = interior, but not for disconnected Julia sets, then critical point and its preimages  !!!!
  		z = z*z +c ; /* forward iteration : complex quadratic polynomial */ 
  	}
  	
  	 
	return n;

}


unsigned char ComputeColorOfLSM(complex double z){

	unsigned char iColor;
	int n; // escape time

	
	n = GiveEscapeTime(z);
	
	// manually udjusted series of ordered colors ( shades of gray )
  	iColor = 255 - 230.0*((double) n)/18.0; // nMax or lower values in denominator
  
  
  	return iColor;


}






// ***************************************************************************************************************************
// ************************** binary decomposition BD/J*****************************************
// ****************************************************************************************************************************

unsigned char ComputeColorOfBD(complex double z){

 int nMax = iterMax_LSM;
  double cabsz;
  unsigned char iColor;
	
  int n;

  for (n=0; n < nMax; n++){ //forward iteration
	cabsz = cabs(z);
    	if (cabsz > ER_LSM) break; // esacping
    	//if (cabsz< PixelWidth) break; // fails into finite attractor = interior but not for disconnected Julia sets, then critical point and its preimages  !!!!
  			
   
     	z = z*z +c ; /* forward iteration : complex quadratic polynomial */ 
  }
  
  if (cimag(z)>0.0) 
  	iColor = 255; 
  	else iColor = 0;
  
  
  return iColor;


}





// ***************************************************************************************************************************
// ************************** modified binary decomposition MBD *****************************************
// ****************************************************************************************************************************

unsigned char ComputeColorOfMBD(complex double z){
// const number of iterations
 int nMax = 7;
  //double cabsz;
  unsigned char iColor;
	
  int n;

  for (n=0; n < nMax; n++){ //forward iteration
	//cabsz = cabs(z);
    	//if (cabsz > ER) break; // esacping
    	//if (cabsz< PixelWidth) break; // falls into finite attractor = interior
  			
   
     	z = z*z +c ; /* forward iteration : complex quadratic polynomial */ 
  }
  
  //if (cabs(z) > 2.0)
  	{ // exterior
  		if (creal(z)>0.0) 
  			iColor = 255; 
  			else iColor = 0;
  	}
  //	else iColor = iColorOfInterior;
  	
  return iColor;


}



// ***************************************************************************************************************************
// ************************** binary decomposition boundaries with texture mapping  *****************************************
// ****************************************************************************************************************************

// https://fractalforums.org/programming/11/how-many-different-ways-are-there-to-show-such-set/3874

/* 
   
   to add
   https://www.iquilezles.org/www/articles/distfunctions2d/distfunctions2d.htm
   
   
   2D Gray gradient  = 2D gray texture
   input x and y is in [0,1] range
   
*/

double Give2DGrayGradient(double x, double y, const int k ){

	
	double d;  // position of the color in the gradient . It is in [0,1]] range
	
	
  	switch(k){
	
  		case 0: {d =  max(fabs(x - 0.5) ,fabs(y-0.5)); break;} // 
  	
  		case 1: {d =  min(x,y); break;}
  		
  		case 2 : {d = fabs(x)+fabs(y) -0.5; break;}
  		
  		case 3 : {d = y; break;}
  		
  		case 4 : {d = x; break;}
  		
  		// gradients 5,6,7 are similar , difference : 1, 1,5, 2.0 
  		case 5: {x =x - 0.5; y =y - 0.5; d = cabs(x+y*I); break;} // cabs(z)
  		
  		case 6: {x =1.5*(x - 0.5); y =1.5*(y - 0.5); d = cabs(x+y*I); break;} // cabs(z)
  		
  		case 7: {x =2.0*(x - 0.5); y =2.0*(y - 0.5); d = cabs(x+y*I); break;} // cabs(z)
  
  		default:{ d= 0.0; }
  	}
  	
  	return d;
  }
	





unsigned char ComputeColorOfTexture(complex double z, const int k){

 int nMax = iterMax_LSM;
  double cabsz;
  unsigned char iColor;
	
  int n;

  for (n=0; n < nMax; n++){ //forward iteration
	cabsz = cabs(z);
    	if (cabsz > ER_LSM) break; // esacping
    	//if (cabsz< PixelWidth) break; // fails into finite attractor = interior but not for disconnected Julia sets, then critical point and its preimages  !!!!
  			
   
     	z = z*z +c ; /* forward iteration : complex quadratic polynomial */ 
  }
  
  // https://gitlab.com/adammajewski/mandelbrot-book_book/-/blob/master/README.md#final-angle
   //if (n < nMax) // exterior 
     // {
        //double et = ((double)n)/nMax; // ok but the same for all points inside level set so segmentation
        double final_angle = c_turn(z); // in [0,1] range
        //double final_radius = GiveNormalizedFinalRadius(z); // = final_absz should be in [0,1]
        
        
      //}

  
  
  double y = frac(n-log(log(cabsz)));
  // inside each cell point has additional coordinate w = (final_angle, final_radius) in [0,1]x[0,1]
  double gray = Give2DGrayGradient(final_angle, y, k);
  iColor = gray*255;
  
  // bd : mark each cell 
  //if (cimag(z)>0.0) iColor =255 -iColor;
  
  return iColor;


}



// plots raster point (ix,iy) 
int DrawPointOfTexture (unsigned char A[], int ix, int iy, const int k)
{
  int i;			/* index of 1D array */
  unsigned char iColor;
  complex double z;


  i = Give_i (ix, iy);		/* compute index of 1D array from indices of 2D array */
  z = GiveZ(ix,iy);
  iColor = ComputeColorOfTexture(z, k);
  A[i] = iColor ;		// interior
  
  return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int DrawImageOfTexture (unsigned char A[], const int k)
{
  unsigned int ix, iy;		// pixel coordinate 

  	fprintf(stderr, "compute image texture k = %d\n", k);
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax)
  	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			DrawPointOfTexture(A, ix, iy, k);	//  
  }

  return 0;
}








// ***********************************************************************************************
//*************************************** SAC/J **************************************************
// *****************************************************************************************
// https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/stripeAC
// SAC = Stripe Average Coloring

//

// the addend function
// input : complex number z
// output : double number t 
double Give_t(double complex z){

  return 0.5+0.5*sin(s*carg(z));

}

/*
  input :
  - complex number
  - intege
  output = average
 
*/
double Give_Arg(double complex z , int iMax)
{
  int i=0; // iteration 
   
   
  //double complex Z= 0.0; // initial value for iteration Z0
  double A = 0.0; // A(n)
  double prevA = 0.0; // A(n-1)
  double R; // =radius = cabs(Z)
  double d; // smooth iteration count
  double complex dz = 1.0; // first derivative with respect to z
  double de; // Distance Estimation from DEM/J  
   
    
  // iteration = computing the orbit
  for(i=0;i<iMax;i++)
    { 
    
      dz = 2.0 * z * dz ; 
      z = z*z + c; // https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/qpolynomials
      
      if (i>i_skip) A += Give_t(z); // 
      
      R = cabs(z);
      // if(R > EscapeRadius) break; // exterior of M set
  	if (R > 1e60 || cabs(dz)> 1e60) break; // prevent NAN 	 	
      prevA = A; // save value for interpolation
        
    } // for(i=0
   
   
  if (i == iMax) 
    A = -1.0; // interior 
  else { // exterior
    de = 2 * R * log(R) / cabs(dz);
    if (de < distanceMax) A = FP_ZERO; //  boundary
    else {
      // computing interpolated average
      A /= (i - i_skip) ; // A(n)
      prevA /= (i - i_skip - 1) ; // A(n-1) 
      // smooth iteration count
      d = i + 1 + log(lnER/log(R))/M_LN2;
      d = d - (int)d; // only fractional part = interpolation coefficient
      // linear interpolation
      A = d*A + (1.0-d)*prevA;
     }   
  }
    
  return A;  
}
 
 
 
 
 
unsigned char ComputeColorOfSAC(complex double z){

  unsigned char iColor;
  double arg;
  
   
   
  	arg = Give_Arg( z, 2500); //   N in wiki 
	
   	// color is proportional to arg 
	if (arg < 0.0)
           
		iColor = 0;  // interior                        
    
		else //  
			{if (arg == FP_ZERO) 
     				iColor = 255; // boundary     
        			else iColor = (unsigned char) (255 - 255*arg );// exterior
      			}
      
    
  return iColor;


}








 

// ***************************************************************************************************************************
// ************************** DLD/J*****************************************
// ****************************************************************************************************************************



/* partial pnorm 
   input: z , zn = f(z), p
   output ppn
   
   
*/
double
ppnorm (complex double z, complex double zn, double p)
{

  double s[2][3];		// array for 2 points on the Riemann sphere
  int j;
  double d;			// denominator 
  double x;
  double y;

  double ds;
  double ppn = 0.0;

  // map from complex plane to riemann sphere
  // z
  x = creal (z);
  y = cimag (z);
  d = x * x + y * y + 1.0;

  s[0][0] = (2.0 * x) / d;
  s[0][1] = (2.0 * y) / d;
  s[0][2] = (d - 2.0) / d;	// (x^2 + y^2 - 1)/d

  // zn
  x = creal (zn);
  y = cimag (zn);
  d = x * x + y * y + 1.0;
  s[1][0] = (2.0 * x) / d;
  s[1][1] = (2.0 * y) / d;
  s[1][2] = (d - 2.0) / d;	// (x^2 + y^2 - 1)/d

  // sum 
  for (j = 0; j < 3; ++j)
    {
      ds = fabs (s[1][j] - s[0][j]);
      //  normal:  neither zero, subnormal, infinite, nor NaN
      //if (fpclassify (ds) !=FP_INFINITE)
      //if (isnormal(ds)) 
      // it is solved by if (cabs(z) > 1e60 ) break; procedure in parent function 
      ppn += pow (ds, p);	// |ds|^p
      //      else {ppn = 10000.0; printf("ds = infty\t");} // 

    }


  return ppn;







}

// DLD = Discret Lagrangian Descriptior
double
lagrangian (complex double z0, complex double c, int iMax, double p)
{

  int i;			// number of iteration
  double d = 0.0;		// DLD = sum
  double ppn;			// partial pnorm
  complex double z = z0;
  complex double zn;		// next z

  for (i = 0; i < iMax; ++i)
    {




      zn = z * z + c;		// complex iteration
      ppn = ppnorm (z, zn, p);
      d += ppn;			// sum
      //
      z = zn;

      //if (! isnormal(d)) { return 0.0; } // not works
      if (cabs (z) > ER_DLD ) //1e6)
	break;			// exterior : big values produces artifacts on the image  



    }





  //if (d<0.0) {// interior
  // d(z1a) - d(z21) = -0.0804163521959989        
  //      d = - d;
  //      d = (db - d) /dd ; // normalize, see test_interior
  //d = d*d;
  //if (d>1.0) {printf("d int > 1.0\n");
  ///     }
  //      else {

  d = d / ((double) i);		// averaging not summation
  //d = d*me;} // exterior

  return d;




}





unsigned char
ComputeColorOfDLD (complex double z)
{


  //double cabsz;
  int iColor;
  double d;
  int N = iterMax_DLD; // N in wiki = fixed number : maximal number of iterations 

  //if (FatouType == 1)
   // {				// interior
     // d = lagrangian (z, c, N, p);
      // modify gradient position

      //{d = d - (int)d;} // only fractional part
     // d = d * d * mi;
      //if ( d< 1.0 ) d = 0.0;

   // }				//  
  //else
    //{
      d = lagrangian (z, c, N, p); //  
    //}

  iColor = (int) (d * 255) % 255;	// nMax or lower walues in denominator



  return (unsigned char) iColor;


}




//=========================================

 
 

// ***************************************************************************************************************************
// ************************** NPM/J = Normal Potential *****************************************
// ****************************************************************************************************************************





/* 
 The dot product of two vectors a = [a1, a2, ..., an] and b = [b1, b2, ..., bn] is defined as:[1]
 d = a1b1 + a2b2
  
*/
double cdot(double complex a, double complex b){
 
 
 return creal(a)*creal(b)+cimag(a)*cimag(b); 


}


// 
// output 
// 
double GiveReflection(double complex z )
  {
   int i=0; // iteration 
   int iMax = 2000;
   
   // https://en.wikipedia.org/wiki/Complex_quadratic_polynomial#First_derivative_with_respect_to_z
   double complex dz = 1.0; // derivative with respect to z 
   double reflection = 0.0; //  
   
   double h2 = 1.5 ; // height factor of the incoming light
   double angle = 45.0/360.0 ; // incoming direction of light in turns 
   double complex v = cexp(2.0*angle *M_PI* I); // = exp(1j*angle*2*pi/360)  // unit 2D vector in this direction
   // incoming light 3D vector = (v.re,v.im,h2)
  
  // https://en.wikipedia.org/wiki/Lambertian_reflectance

   
   double  complex u;
   
   
   z  = z*z+c; // 
   dz = 1.0;
   
   for(i=0;i<iMax;i++)
    {  
      dz = 2.0*dz*z ;
      z  = z*z+c; // 
      
      
      
      if(cabs(z) > ER_NP) 
      { // exterior
        u = z / dz;
        u = u / cabs(u);
        reflection =  cdot(u, v) + h2;  /* use the simplest model for the shading: 
        Lambert, which consists in using the dot product of (x,y,1) with a constant vector indicating the direction of the light. */
        reflection = reflection/(1.0 + h2);  // rescale so that t does not get bigger than 1
        if (reflection<0.0) reflection =0.0;
        break;
      
      }
    }
    
   return reflection;  
 }





// Potential to color
unsigned char ComputeColorOfNP(complex double z){
//https://www.math.univ-toulouse.fr/~cheritat/wiki-draw/index.php/Mandelbrot_set#Normal_map_effect


  
  
  double reflection;
  unsigned char iColor;
   
   
   
  // compute 
   reflection = GiveReflection( z);
  
    
  // 
  //if (reflection <  )
    //{ /*  interior  */
    //  iColor = 0;}
  //else // exterior 
        
    { iColor = 255 * reflection;}
     
  return iColor;   
  
 
}



// https://en.wikipedia.org/wiki/Shading

//  normal = perpendicular 
// shading using Normal map and Potential
// https://en.wikipedia.org/wiki/Lambertian_reflectance
// http://www.math.titech.ac.jp/~kawahira/gallery/movies/movies.html
// see 0_1.avi and image 
//






 
 

// ***************************************************************************************************************************
// ************************** NDM/J = Normal Distance *****************************************
// ****************************************************************************************************************************


//  normal = perpendicular 
// shading using Normal map and Potential
// https://en.wikipedia.org/wiki/Lambertian_reflectance
// http://www.math.titech.ac.jp/~kawahira/gallery/movies/movies.html
// see 0_1.avi and image 
//





// 
// output 
// 
double GiveReflectionD(double complex z )
  {
   int i=0; // iteration 
   int iMax = 2000;
   
   // https://en.wikipedia.org/wiki/Complex_quadratic_polynomial#First_derivative_with_respect_to_z
   double complex dz = 1.0;   // first derivative with respect to z 
   double complex dz2 = 0.0;  // second derivative with respect to z 
   double reflection = 0.0; //  
   double lo;
   
   double h2 = 1.5 ; // height factor of the incoming light
   double angle = 45.0/360.0 ; // incoming direction of light in turns 
   double complex v = cexp(2.0*angle *M_PI* I); // = exp(1j*angle*2*pi/360)  // unit 2D vector in this direction
   // incoming light 3D vector = (v.re,v.im,h2)
   
  
  // https://en.wikipedia.org/wiki/Lambertian_reflectance

   
   double  complex u;
   
   
   z  = z*z+c; // 
   dz = 1.0;
   dz2 = 0.0;
   
   for(i=0;i<iMax;i++)
    {  
      
      dz2 = 2.0* ( dz2*z + dz*dz);//2*(der2*z+der**2)
      dz = 2.0*dz*z ;
      z  = z*z+c; // 
      
      
      
      if(cabs(z) > ER_NP) 
      { // exterior
      
      /*
       lo = 0.5*log(squared_modulus(z))
    u = z*der*((1+lo)*conj(der**2)-lo*conj(z*der2))
    u = u/abs(u)
    */
      
      	lo = 0.5*log(cabs(z));
    	u = z*dz*((1.0+lo)*conj(dz*dz)-lo*conj(z*dz2));
        //u = z / dz;
        u = u / cabs(u);
        reflection =  cdot(u, v) + h2;  // use the simplest model for the shading: Lambert, which consists in using the dot product of (x,y,1) with a constant vector indicating the direction of the light.
        reflection = reflection/(1.0 + h2);  // rescale so that t does not get bigger than 1
        if (reflection<0.0) reflection =0.0;
        break;
      
      }
    }
    
   return reflection;  
 }





// Distance to color
unsigned char ComputeColorOfND(complex double z){
//https://www.math.univ-toulouse.fr/~cheritat/wiki-draw/index.php/Mandelbrot_set#Variation


  
  
  double reflection;
  unsigned char iColor;
   
   
   
  // compute 
   reflection = GiveReflectionD( z);
  
    
  // 
  //if (reflection <  )
    //{ /*  interior  */
    //  iColor = 0;}
  //else // exterior 
        
    { iColor = 255 * reflection;}
     
  return iColor;   
  
 
}


// -------------------------- potential========


double ComputePotential(const complex double z0){

	double potential = 0.0; // interior
	double s = 0.5;
	complex double z = z0;
	double r;
	int iter;
	
	for (iter = 0; iter < iterMax_pot; ++iter){
		
		z = z*z +c; // complex quadratic polynomial
		s *= 0.5;  // 
		r = cabs(z);
		if (r > ER_pot) {break;}
	
	
	}
	
	
	
	
	
	potential =  s*log2(r); // log(zn)* 2^(-n)
	return potential;
	
}


unsigned char ComputeColorOfPOT(complex double z){


	double potential = ComputePotential(z);
	
	if (PlaneInversion) // usung global var 
		{potential /= 4.0;} // manual normalize
	unsigned char iColor = 255 * sqrt(sqrt(potential));
     	return iColor;   
  
 
}


/*
int local_setup(int PlaneInversion){

	if (PlaneInversion)
		{ MaxImagePotential =ComputePotential( 1.0/ 0.0);}
		//else {MaxImagePotential}

	return 0;

}

*/


 
 
/* ==================================================================================================
 ============================= Draw functions ===============================================================
 =====================================================================================================
*/ 
unsigned char ComputeColor(FunctionTypeT FunctionType, complex double z){

	unsigned char iColor;
	
	
	
	switch(FunctionType){
	
		case LSM :{iColor = ComputeColorOfLSM(z); break;}
		
		case DEM : {iColor = ComputeColorOfDEMJ(z); break;}
		
		case Unknown : {iColor = ComputeColorOfUnknown(z); break;}
		
		case BD : {iColor = ComputeColorOfBD(z); break;}
		
		case MBD : {iColor = ComputeColorOfMBD(z); break;}
		
		case SAC : {iColor = ComputeColorOfSAC(z); break;}
  
		case DLD : {iColor = ComputeColorOfDLD(z); break;}
		
		case ND : {iColor = ComputeColorOfND(z); break;}
		
		case NP : {iColor = ComputeColorOfNP(z); break;}
		
		case POT : {iColor = ComputeColorOfPOT(z); break;}
		
	
		default: {}
	
	
	}
	
	return iColor;



}


// plots raster point (ix,iy) 
int DrawPoint (FunctionTypeT FunctionType, int PlaneInversion, unsigned char A[], int ix, int iy)
{
  int i;			/* index of 1D array */
  unsigned char iColor;
  complex double z;


  i = Give_i (ix, iy);		/* compute index of 1D array from indices of 2D array */
  
  if (PlaneInversion)
  	{ 
  		complex double w;
  		w = GiveW(ix,iy);
  		z = 1/w;
  	}
  	else {  z = GiveZ(ix,iy);}
  

  iColor = ComputeColor(FunctionType, z);
  A[i] = iColor ;		// 
  
  return 0;
}




// fill array 
// uses global var :  ...
// scanning complex plane 
int DrawImage (FunctionTypeT FunctionType, int PlaneInversion, unsigned char A[])
{
  unsigned int ix, iy;		// pixel coordinate 
  	
  	//local_setup(PlaneInversion);

  	fprintf(stderr, "compute image FunctionType = %d PlaneInversion = %d\n", FunctionType, PlaneInversion);
 	// for all pixels of image 
	#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax)
  	for (iy = iyMin; iy <= iyMax; ++iy){
    		fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    		for (ix = ixMin; ix <= ixMax; ++ix)
      			DrawPoint(FunctionType, PlaneInversion, A, ix, iy);	//  
  }

  return 0;
}

 

 
 
// *******************************************************************************************
// ********************************** save A array to pgm file ****************************
// *********************************************************************************************

int
SaveArray2PGMFile (unsigned char A[], char *shortName , char *comment)
{

  FILE *fp;
  const unsigned int MaxColorComponentValue = 255;	/* color component is coded from 0 to 255 ;  it is 8 bit color file */
  
  
  // https://programmerfish.com/create-output-file-names-using-a-variable-in-c-c/
  char fileName[512];
  const char* fileType = ".pgm";
  sprintf(fileName,"%s%s", shortName, fileType); // 
  
  
  
  char long_comment[200];
  sprintf (long_comment, "fc(z)=z^2+ c where c = %f %+f*i ;  %s", creal(c), cimag(c),comment);





  // save image array to the pgm file 
  fp = fopen (fileName, "wb");	// create new file,give it a name and open it in binary mode 
  fprintf (fp, "P5\n # %s\n %u %u\n %u\n", long_comment, iWidth, iHeight, MaxColorComponentValue);	// write header to the file
  size_t rSize = fwrite (A, sizeof(A[0]), iSize, fp);	// write whole array with image data bytes to the file in one step 
  fclose (fp);

  // info 
  if ( rSize == iSize) 
  	{
  		printf ("File %s saved ", fileName);
  		if (long_comment == NULL || strlen (long_comment) == 0)
    		printf ("\n");
  			else { printf (". Comment = %s \n", long_comment); }
  	}
  	else {printf("wrote %zu elements out of %u requested\n", rSize,  iSize);}
  	
  	
  // 
  NumberOfImages +=1; // count images using global variable

  return 0;
}
















int PrintInfoAboutProgam()
{
	printf("Number of pgm images = %d \n", NumberOfImages);	
  
  // display info messages
  printf ("Numerical approximation of Julia set for fc(z)= z^2 + c \n");
  //printf ("iPeriodParent = %d \n", iPeriodParent);
  //printf ("iPeriodOfChild  = %d \n", iPeriodChild);
  printf ("parameter c = %.16f %+.16f*i  \n", creal(c), cimag(c));
  
  printf ("Image Width = %f in world coordinate\n", ZxMax - ZxMin);
  printf ("PixelWidth = %f \n", PixelWidth);
  
  printf("for DEM/J \n");
  if ( distanceMax<0.0 || distanceMax > ER ) printf("bad distanceMax\n");
	printf("Max distance from exterior to the boundary =  distanceMax = %.16f = %f pixels\n",  distanceMax, BoundaryWidth); 
  printf("\n");
  
  
  // image corners in world coordinate
  // center and radius
  // center and zoom
  // GradientRepetition
  printf ("Maximal number of iterations = iterMax = %ld \n", iterMax);
  
  printf ("For LSM/J \n");
  printf ("Maximal number of iterations = iterMax_LSM = %ld \n", iterMax_LSM);
  printf ("Escape Radius = ER_LSM = %f \n", ER_LSM);
  printf("\n");
  
  
  
  
  printf ("ratio of image  = %f ; it should be 1.000 ...\n", ratio);
  //
  printf("gcc version: %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__); // https://stackoverflow.com/questions/20389193/how-do-i-check-my-gcc-c-compiler-version-for-my-eclipse
  // OpenMP version is displayed in the console 
  return 0;
}






// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************

int setup ()
{

  fprintf (stderr, "setup start\n");
  c = 0.35; //   
  
  
  
  
	
  /* 2D array ranges */
  
  iWidth = iHeight* DisplayAspectRatio;
  iSize = iWidth * iHeight;	// size = number of points in array 
  // iy
  iyMax = iHeight - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  //ix

  ixMax = iWidth - 1;

  /* 1D array ranges */
  // i1Dsize = i2Dsize; // 1D array with the same size as 2D array
  iMax = iSize - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  
  
   SetZPlane( center, radius,  DisplayAspectRatio );	

  /* Pixel sizes */
  PixelWidth = (ZxMax - ZxMin) / ixMax;	//  ixMax = (iWidth-1)  step between pixels in world coordinate 
  PixelHeight = (ZyMax - ZyMin) / iyMax;
  ratio = ((ZxMax - ZxMin) / (ZyMax - ZyMin)) / ((double) iWidth / (double) iHeight);	// it should be 1.000 ...
	
  wPixelWidth = (WxMax-WxMin)/ixMax;
  wPixelHeight =(WyMax-WyMin)/iyMax;
	
  
  //ER2 = ER * ER; // for numerical optimisation in iteration
  lnER = log(EscapeRadius); // ln(ER) 
  loger = log(ER_LSM); // for texture
  ER_LSM = GiveER(10); // find such ER for LSM/J that level curves croses critical point and it's preimages
  ER_DLD = GiveER(15);
  
  MaxFinalRadius =  GiveMaxFinalRadius();
  
  
   	
  /* create dynamic 1D arrays for colors ( shades of gray ) */
  data = malloc (iSize * sizeof (unsigned char));
  edge = malloc (iSize * sizeof (unsigned char));
  edge2 = malloc (iSize * sizeof (unsigned char));
  //
 
  	
  if (data == NULL || edge == NULL || edge2 == NULL ){
    fprintf (stderr, " Could not allocate memory");
    return 1;
  }

  
 	
  
  BoundaryWidth = 6.0*iWidth/2000.0  ; //  measured in pixels ( when iWidth = 2000) 
  distanceMax = BoundaryWidth*PixelWidth;
  
  
  
  fprintf (stderr," end of setup \n");
	
  return 0;

} // ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




int end(){


  fprintf (stderr," allways free memory (deallocate )  to avoid memory leaks \n"); // https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
  free (data);
  free(edge);
  free(edge2);
 
  PrintInfoAboutProgam();
  return 0;

}

// ********************************************************************************************************************
/* -----------------------------------------  main   -------------------------------------------------------------*/
// ********************************************************************************************************************

int main () {
  
  
  
  setup ();
  
  PlaneInversion = 0; 
  
   
  
  
  DrawImage(DEM, PlaneInversion, data);
  SaveArray2PGMFile (data, "de", "boundary using DEM"); // name of the file is name.png 
  
   
  DrawImage(BD, PlaneInversion, data);
  SaveArray2PGMFile (data, "bd", "Binary Decomposition");
  
  ComputeBoundaries(data, edge);
  SaveArray2PGMFile (edge, "bdb", "boundaries of BD");
  
  DrawImage(MBD, PlaneInversion, data);
  SaveArray2PGMFile (data,"mbd" , "Modified Binary Decomposition");
  
  ComputeBoundaries(data, edge2);
  SaveArray2PGMFile (edge2, "mbdb", "boundaries of MBD");
  
  
  DrawImage(LSM, PlaneInversion, data);
  SaveArray2PGMFile (data, "ls", "Level Set Method of Escape Time");
  
  
  ComputeBoundaries(data, edge);
  SaveArray2PGMFile (edge, "lc", "LCM = boundaries of LSM = Level Curves Method");
    
  CopyBoundaries(edge, data);
  SaveArray2PGMFile (data, "lsc", "LSM + LCM (boundaries of LSM");
  
  CopyBoundaries(edge, edge2);
  SaveArray2PGMFile (edge2, "lcmbd", "boundaries of LSM and MBD");
  
  DrawImage(BD, PlaneInversion, data); // data =  "bd", "Binary Decomposition");
  ComputeBoundaries(data, edge); // edge =  "bdb", "boundaries of BD");
  DrawImage(LSM, PlaneInversion, data); //(data = "ls", "Level Set Method of Escape Time");
  ComputeBoundaries(data, edge2); //edge2 =  "lc", "LCM = boundaries of LSM = Level Curves Method");
  CopyBoundaries(edge, edge2); // lc+bdb
  SaveArray2PGMFile (edge2, "lcbdb", "boundaries of LSM and BD");
  
  // https://fractalforums.org/programming/11/how-many-different-ways-are-there-to-show-such-set/3874
  //  xenodreambuie
  
  
  // texture
  DrawImageOfTexture(data, 0);
  SaveArray2PGMFile (data, "t0", "texture 0 ; x= angle y = frac(radius)");
  
  DrawImageOfTexture(data, 1);
  SaveArray2PGMFile (data, "t1", "texture 1 ; x= angle y = frac(radius)");
  
  
  DrawImageOfTexture(data, 2);
  SaveArray2PGMFile (data, "t2", "texture 2 ; x= angle y = frac(radius)");
  
  DrawImageOfTexture(data, 3);
  SaveArray2PGMFile (data, "t3", "texture 3 ; x= angle y = frac(radius)");
  
  
  DrawImageOfTexture(data, 4);
  SaveArray2PGMFile (data, "t4", "texture 4 ; x= angle y = frac(radius)");
  
  DrawImageOfTexture(data, 6);
  SaveArray2PGMFile (data, "t6", "texture 5 cabs(z) ; x= angle y = frac(radius)");
  
  CopyBoundaries(edge2, data); // lc+bdb
  SaveArray2PGMFile (data, "t5b", "texture 5 cabs(z) and boundaries ; x= angle y = frac(radius)");
  
  
  DrawImageOfTexture(data, 6);
  CopyBoundaries(edge2, data); // lc+bdb
  SaveArray2PGMFile (data, "t6b", "texture 6 cabs(z) and boundaries ; x= angle y = frac(radius)");
  
  
  DrawImageOfTexture(data, 7);
  CopyBoundaries(edge2, data); // lc+bdb
  SaveArray2PGMFile (data, "t7b", "texture 7 cabs(z) and boundaries ; x= angle y = frac(radius)");
  
  
  
  
  DrawImage(Unknown, PlaneInversion, data);
  SaveArray2PGMFile (data, "u", "Unknown : boundary and slow dynamics");
 
     
  DrawImage(SAC, PlaneInversion, data);
  SaveArray2PGMFile (data, "sac", "SAC + DEM");
  
  DrawImage(DLD, PlaneInversion, data);
  DrawImageOfDEMJ_boundary(data, PlaneInversion, 255);
  SaveArray2PGMFile (data, "dld", "DLD/J + boundary by DEM");
 
    
  DrawImage(NP, PlaneInversion, data);
  SaveArray2PGMFile (data, "np", "NP/J");
  
  
  DrawImage(ND, PlaneInversion, data);
  SaveArray2PGMFile (data, "nd", "ND/J");
  
  
  
  DrawImage(POT, PlaneInversion, data);
  DrawImageOfDEMJ_boundary(data, PlaneInversion, 255);
  SaveArray2PGMFile (data, "pot", "potential"); // name of the file is name.png 
  
  
  // test image
  DrawImage(DEM, PlaneInversion, data);
  CheckZPlaneOrientation(data);
  SaveArray2PGMFile (data, "defq", "boundary using DEM/J and first quadrant");
 
  
  // z window
  DrawImage(DEM, PlaneInversion, data);
  ShowWWindowOnZWindow(data);
  SaveArray2PGMFile (data, "wonz", "W Window On Z Window");
  
  
  
  
  // ------------------ inverterd plane = wplane = 1/z plane -------------------------------------
 
  PlaneInversion = 1; 
  
  DrawImage(DEM, PlaneInversion, data);
  SaveArray2PGMFile (data, "dei", "boundary using DEM/J inv");
  
  
  DrawImage(BD, PlaneInversion, data);
  SaveArray2PGMFile (data, "bdi", "BD/J inverted ");
  
 
  ComputeBoundaries(data, edge);
  SaveArray2PGMFile (edge, "bdbi", "boundaries of BD inv");
  
  
  DrawImage(BD, PlaneInversion, data); // data =  "bdi", "Binary Decomposition inv");
  ComputeBoundaries(data, edge); // edge =  "bdbi", "boundaries of BD inv");
  DrawImage(LSM, PlaneInversion, data); //(data = "lsi", "Level Set Method of Escape Time inv");
  ComputeBoundaries(data, edge2); //edge2 =  "lci", "LCM = boundaries of LSM = Level Curves Method inv");
  CopyBoundaries(edge, edge2); // lc+bdb
  SaveArray2PGMFile (edge2, "lcbdbi", "boundaries of LSM and BD inv");
  
  
  
  
  
  DrawImage(SAC, PlaneInversion, data);
  SaveArray2PGMFile (data, "sacdei", "SAC + DEM inverted");
  
  DrawImage(LSM, PlaneInversion, data);
  SaveArray2PGMFile (data, "lsi", "LSM inv");
 
 
  
  ComputeBoundaries(data, edge);
  SaveArray2PGMFile (edge, "lci", "boundaries of LSM inv");
    
  CopyBoundaries(edge, data);
  SaveArray2PGMFile (data, "lsci", "LSM + boundaries of LSM/J inv");
  
  
   
  DrawImage(NP, 1, data);
  SaveArray2PGMFile (data, "npi", "NP inverted");
   
  DrawImage(ND, 1, data);
  SaveArray2PGMFile (data, "ndi", "ND inverted");
 
  
  DrawImage(POT, PlaneInversion, data);
  DrawImageOfDEMJ_boundary(data, PlaneInversion, 255);
  SaveArray2PGMFile (data, "pot_i", "potential inverted"); // name of the file is name.png 
  
  DrawImage(DEM, PlaneInversion, data); // w window
  ShowZWindowOnWWindow(data);
  SaveArray2PGMFile (data, "zonw", "Z Window On W Window");
  
  
  
  //
  end();

  return 0;
}
