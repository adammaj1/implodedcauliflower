/*

gcc i.c -lm - Wall
./a.out


helper program for testing funcions



*/
#include <stdio.h>
#include <stdlib.h>		// malloc
#include <string.h>		// strcat
#include <math.h>		// M_PI; needs -lm also
#include <complex.h> 		// complex numbers : https://stackoverflow.com/questions/6418807/how-to-work-with-complex-numbers-in-c


// complex numbers of parametr plane 
double complex c = 0.35;		// parameter of function fc(z)=z^2 + c



// ------------------- ET -----------------------------------------
//static unsigned long int iterMax = 1000000;	//iHeight*100;
unsigned long int iterMax_LSM = 255;
const int iterMax_pot = 400; // potential 



double ER = 200.0;		// EscapeRadius for bailout test 
double EscapeRadius=1000000; // = ER big !!!!
double ER_LSM ; // see GiveER_LSM  // 27.764 =  manually find value such that level curves of escape time cross critical point and it's  preimages
double ER_DLD ; // see GiveER_LSM  // 27.764 =  manually find value such that level curves of escape time cross critical point and it's  preimages
double ER_NP = 100.0; 
double ER_pot = 100000.0;  // sqrt(1e24)


// -------------- DEM ---------------------------------------
double BoundaryWidth = 3.0; // % of image width  
double distanceMax; //distanceMax = BoundaryWidth*PixelWidth;

// --------------------  SAC/J  ------------------------------
double lnER; // ln(ER)
int i_skip = 2; // exclude (i_skip+1) elements from average
unsigned char s = 7; // stripe density






//  ------------- DLD  ----------------------
const int N = 20;		// fixed number : maximal number of iterations
double p = 0.180; //0.01444322;		//
double mi = 0.9;






// - -------------------- functions ------------------------------------------------------------







// find such ER for LSM/J that level curves croses critical point and it's preimages
double GiveER(int i_Max){

	complex double z= 0.0; // criical point
	int i;
	 ; // critical point escapes very fast here. Higher valus gives infinity
	for (i=0; i< i_Max; ++i ){
		z=z*z +c; 
	 
	 }
	 
	 return cabs(z);
	
	
}



int setup(){


	ER_LSM = GiveER(10); // find such ER for LSM/J that level curves croses critical point and it's preimages
  	ER_DLD = GiveER(15);
  	return 0;

}







// ***********************************************************************************************
//*************************************** SAC/J **************************************************
// *****************************************************************************************
// https://en.wikibooks.org/wiki/Fractals/Iterations_in_the_complex_plane/stripeAC
// SAC = Stripe Average Coloring

//

// the helper function
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
    if (de < distanceMax) A = FP_ZERO; //  boundary  !!!!! FP_ZERO is integer constant ( enum)
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
        reflection =  cdot(u, v) + h2;  // use the simplest model for the shading: Lambert, which consists in using the dot product of (x,y,1) with a constant vector indicating the direction of the light.
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









/*
 https://en.wikibooks.org/wiki/Fractals/mandelbrot-numerics
 https://code.mathr.co.uk/mandelbrot-numerics/blob/HEAD:/c/bin/m-describe.c
 
    double cre = 0;
    double cim = 0;
   
    double _Complex c = cre + I * cim;
    double _Complex z = 0;
    double _Complex dc = 0;
    printf("the input point was %+.18f + %+.18f i\n", creal(c), cimag(c));
    bool escaped = false;
    for (int i = 0; i < maxiters; ++i)
    {
      dc = 2 * z * dc + 1;
      z = z * z + c;
      if (!escaped && cabs2(z) > 1e10)
      {
        double mu = i + 1 - log2(log(cabs(z)));
        double de = 2 * cabs(z) * log(cabs(z)) / cabs(dc);
        printf("the point escaped with dwell %.5f\nand exterior distance estimate %.5g\n", mu, de);
        escaped = true;
      }
    }
    if (! escaped)
    {
      printf("the point didn't escape after %d iterations\n", maxiters);
    }
    
    
*/



int EscapeTimeAndDECheck(complex double z){

	int NotEscaped = 1; // bool flag
	
	int iMax = 200;
	double complex dc = 0;
	double ER = 1e5;
	
	for (int i = 0; i < iMax; ++i)
	{
		if (NotEscaped && cabs(z) > ER) // first test 
      			{
        			double mu = i + 1 - log2(log(cabs(z))); // http://linas.org/art-gallery/escape/escape.html  continous escape time = 
        			double de = 2 * cabs(z) * log(cabs(z)) / cabs(dc);
        			printf("the point escaped from circle with EscapeRadius = %.0f with after i = %d iterations , dwell mu =  %.5f and exterior distance estimate de = %.5f cabs(z) = %f\n", ER, i, mu, de, cabs(z));
        			NotEscaped = 0; //
      			}
      			
      			
      		// iteration after test		
      		dc = 2 * z * dc + 1;
      		z = z * z + c;
      		
    	}
    if (NotEscaped)
    	{printf("the point didn't escape from circle with EscapeRadius = %f after %d iterations, cabs(z) = %f\n",ER,  iMax, cabs(z)); }
	

	return 0;




}





int GiveEscapeTime(complex double z){


	int nMax = iterMax_LSM;
	double cabsz;
	int n;

  	for (n=0; n <= nMax; ++n){ //forward iteration
  	
		cabsz = cabs(z); // first test
    		if (cabsz > ER_LSM) break; // escaping
    		//if (cabsz< PixelWidth) break; // fails into finite attractor = interior, but not for disconnected Julia sets, then critical point and its preimages  !!!!
  		z = z*z +c ; /* forward iteration : complex quadratic polynomial */ 
  	}
  	
  	
  	if (n ==nMax) 
  		{ printf("the point didn't escape from circle with EscapeRadius = %f  after %d iterations, cabs(z) = %f\n", ER_LSM,  nMax, cabs(z));  }
  		else {printf("the point escaped from circle with EscapeRadius = %f  after %d iterations, cabs(z) = %f\n", ER_LSM,  n, cabs(z)); }
  		
	return n;

}


//


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










 
int PrintInfoAboutPoint(complex double z0){

	//unsigned int ix, iy;		// pixel coordinate
	// to do 
	
	double arg; //SAC/J
	unsigned char iColor;
	double reflection;//  = GiveReflection(z); 
	int et;
	double potential;
	
	
	printf ("input:  z = %.16f %+.16f*i\n", creal(z0), cimag(z0));
	printf ("input:  c = %.16f %+.16f*i\n", creal(c), cimag(c));
	
	// ET
	et = GiveEscapeTime(z0);
	printf("Escape Time = %d\n", et);
	// DE
	EscapeTimeAndDECheck(z0);
	
	
	
	// SAC/J
	arg = Give_Arg( z0, 2500); //   N in wiki
	iColor = ComputeColorOfSAC(z0);
	printf ("SAC/J : arg = %.16f ; iColor = %d  \n", arg, iColor);
	
	// NP
	reflection  = GiveReflection(z0);
	iColor = ComputeColorOfNP(z0);
	printf ("NP/J : reflection = %.16f ; iColor = %d  \n", reflection, iColor);
	
	// 
	potential = ComputePotential(z0);
	iColor = (unsigned char)(potential*255);
	printf ("potential = %.16f ; iColor = %d  \n", potential, iColor);

	return 0; 

}








int main(){

	
  
  	setup();
	PrintInfoAboutPoint(0.0);
	PrintInfoAboutPoint(-0.591607978309962*I);
	return 0;

}
