#include "rgig.h"
#include <cmath>
#include <algorithm>    //std::min, std::max

//include the std namespace (cmath and algorithm)
using namespace std;



//we need pi
const double pi = acos(-1.0);

//we also need a nan value
double gig::nan(){
  //does implementation have quite NaN
  if( numeric_limits<double>::has_quiet_NaN ){
    return numeric_limits<double>::quiet_NaN();
  }else{ //if not, create a NaN
    return 0.0/0.0;
  }
}//double gig::nan()

double gig::gig_propto(double x_in){
  return exp( (p_abs-1) * log(x_in)- (x_in + 1 / (x_in) )/two_d_beta);
}
	
double gig::sqrt_gig_propto(double x_in){
  return exp((p_abs-1)/2 *log(x_in) - (x_in + 1 / (x_in) )/(2*two_d_beta ) );
}

double gig::sqrt_gig_ratio(double x_in, double m_in){
	double ret = pow(x_in / m_in, (p_abs - 1)/2);
	ret *= exp( ( m_in + 1/m_in - x_in - 1/x_in  )/(2 * two_d_beta ) );
	return ret;
}	

double gig::scaleX(double x_in){
  if( p<0 ){ x_in = 1.0 / x_in; }
  x_in /= alpha;
  return  x_in;
}

double gig::algortihm_a0(){
  //a==0 -> invGamma(-p,b/2); p<0, b>0
  if( p>=0 ){ return  this->nan(); }
  gamma_distribution<double> gammaDist(-p, 1);
  return 1.0 / (gammaDist(rgen)*2.0/b);
}//double gig::algortihm_a0(double p, double b)

double gig::algortihm_b0(){
  //b==0 -> Gamma(p, a/2); p>0, a>0
  if( p<=0 ){ return  this->nan(); }
  gamma_distribution<double> gammaDist(p, 1);
  return gammaDist(rgen)*2.0/a;
}//double gig::algortihm_b0(double p, double a)


double gig::algortihm_region1(){
	two_d_beta = 2.0 / beta;
	
	// calc mode
  double mode = sqrt(pow(p_abs - 1,2) + pow(beta, 2) ) + (p_abs - 1);
	mode /= beta;
	
	
	// COMPUTING minimal bounding rectangle through
	//Setting up and solving Caradanos to get x_m,x_p
  double c1 =(- two_d_beta * (p_abs + 1)  - mode) / 3.0;
	double c2 =   two_d_beta * (p_abs - 1)  * mode - 1;
  double c3 = min(c2/3.0 - pow(c1, 2), 0.0); //c3 should be negative
  //c3 only ever used as sqrt(-c3)
  c3 = sqrt( -c3 );
  double c4 = c1 * (2.0 * pow(c1, 2) - c2) + mode;
  //computation of c5 -> ensure that -1 <= c5 <= 1 (numeric stability)
  double c5 = -c4/2.0 * pow( c3, -3);
  c5 = max( min(c5,1.0), -1.0);
  c5 = acos(c5) / 3.0;
  double c6 = 2.0 * c3;
  double x_m = c6 * cos(c5 + 4.0*pi/3.0) - c1;
  double x_p = c6 * cos(c5             ) - c1;
  //since c5 = acos(c5)/3.0 <= pi/3
  //then cos(c5) >= cos(c5+4*pi/3) -> x_m <= x_p;
  //numeric stability, ensure that x_p >= x_m
  x_p = max(x_p, x_m);
	
	// triangle to simulate
  double u_m_div_v = (x_m - mode) * sqrt_gig_ratio(x_m, mode);
  double u_p_div_v = (x_p - mode) * sqrt_gig_ratio(x_p, mode);
  double u_p_m = u_p_div_v - u_m_div_v;

	// generating
	bool new_value = true;
  double x;
	do{
    double U = U01_dist(rgen);
    double V = U01_dist(rgen);
		x  = (u_p_m*U + u_m_div_v) / V  + mode; 
		if(x < 0)
			continue;
		
		if(V <= sqrt_gig_ratio(x, mode)){
			new_value = false;
    }
	}while( new_value );
  
  return scaleX(x);
}//double gig::algortihm_region1()

double gig::algortihm_region2(){
	two_d_beta = 2.0 / beta;
	// calc mode
  double mode = beta;
	mode /= sqrt(pow(1 - p_abs, 2) + pow(beta, 2) ) + (1-p_abs);	
	//cout << "mode =" << mode <<"\n";
	// benerating constants and domains
	double x0 = beta / (1.0-p_abs);
  double xs = max(x0, two_d_beta);
	
  double c1 = gig_propto(mode);
  double c2 = 0;
  double A1 = c1 * x0;
  double A2 = 0;
  double x0_pow_p = 1;
  
	if( x0 < two_d_beta){
		c2 = exp(- beta);
		if(p_abs >0){
			x0_pow_p = pow(x0, p_abs); 
			A2 = c2 *( pow(two_d_beta, p_abs) - x0_pow_p)/ p_abs;
		}else{ // |p|=0
			A2 = c2 * log(two_d_beta / beta);
    }
	}else{
		c2 = 0;
		A2 = 0;
	}
  
  double c3 = pow(xs, p_abs - 1);
	double A3 = 2.0 * c3 * exp(-xs / two_d_beta) / beta;
	double A = A1 + A2 + A3;
	// generator
  double c4, U, x;
	do{
    U = U01_dist(rgen);
    double V = U01_dist(rgen) * A;
		
		if(V <= A1){ // region (0, x0)
			x = x0 * V / A1;
			c4 = c1;
		}else if(V <= A1 + A2){ // region (x0, two_d_beta)
			V = V - A1;
			if(p_abs > 0){
				x = pow(x0_pow_p + V * p_abs / c2, 1.0 / p_abs);
			}else{
				x = beta * exp(V * exp(beta));
      }
			c4 = c2 * pow(x, p_abs - 1);	
		}else{ // region (two_d_beta, \infinty)
			V = V - A1 - A2;
			x = - two_d_beta * log( exp( - xs / two_d_beta) - V  / (c3 *two_d_beta ));
			c4 = c3 * exp( - x / two_d_beta);
		}
	}while(U * c4 >= gig_propto(x));
  
  return scaleX(x);
}//double gig::algortihm_region2()

double gig::algortihm_region3(){
	two_d_beta = 2.0 / beta;
	// calc mode
  double mode = beta;
	mode /= sqrt(pow( 1 - p_abs ,2) + pow(beta, 2) ) + 1 - p_abs ;	
	// computing miniminal bounding rectangle
  double x_p = 1 + p_abs + sqrt( pow(1 + p_abs,2) + pow(beta, 2) );
	x_p /= beta;
  double v_p = sqrt_gig_propto(mode);
  double u_p = x_p * sqrt_gig_propto(x_p);

  double V, x;
	do{
    double U = U01_dist(rgen) * u_p;
    V = U01_dist(rgen) * v_p;
		x = U / V;
	}while(V >= sqrt_gig_propto(x));
  
  return scaleX(x);
}//double gig::algortihm_region3()


double gig::sample(const double p_in,const double a_in,const double b_in){
  //store the input parameters in the object
  p=p_in; a=a_in; b=b_in;

  double x_out=0;
  
  //invalid parameter space return NaN (or 0/0); cases:
  // 1) a<0
  // 2) b<0
  // 3) a==0 & b==0
  if( a<0 || b<0 || (a==0 && b==0) ){ return this->nan(); }
  
  //special boundary cases for a=0 or b=0
  if( a==0 ){
    //a==0 -> invGamma(-p,b/2); p<0, b>0
    x_out = algortihm_a0();

  }else if( b==0 ){
    //b==0 -> Gamma(p, a/2); p>0, a>0
    x_out = algortihm_b0();
    
  }else{
    alpha = sqrt(a / b);
    beta  = sqrt(a * b);
    //reformulating the model to:
    //\propto alpha^{-p} x^{p-1} e^{-beta/2 * ( (x/alpha)^-1 + (x/alpha)) }
	
    // utilize that x^-1 \propto GIG(-p,\alpha,\beta) 
    p_abs = fabs(p);
    if( beta>1 || p_abs>1){
      x_out = algortihm_region1();
    }else if( beta <= min(1.0/2.0, 2.0/3.0 * sqrt(1-p_abs)) ){
      x_out = algortihm_region2();
    }else{
      x_out = algortihm_region3();
    }
  }
  
  return x_out;
}
