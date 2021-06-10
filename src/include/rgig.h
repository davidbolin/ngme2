#ifndef __RGIG__H
#define __RGIG__H
// random generalized inverse Gaussian generator
// lambda, chi, psi using wiki notation
// *double where to put the r.v
// length of random variableser
// adding seed

#include <random>

class gig{
private:
  //	double mode;
  
  // parameters for model 
	double two_d_beta;
	double p, a, b; //original parameters
	double alpha, beta, p_abs; //parameter modification

  //A U(0,1)-distribution
  std::uniform_real_distribution<double> U01_dist;

  //A random-generator
  std::mt19937_64 rgen;
   
	/*
		Special case of a==0 -> invGamma(-p,b/2); p<0, b>0	
		@out return sample
	*/ 
	double algortihm_a0();

  /*
		Special case of b==0 -> Gamma(p, a/2); p>0, a>0
		@out return sample
	*/ 
	double algortihm_b0();

	/*
		Algorithm if the parameter region is
		sqrt(a * b)> 1 or abs(p) > 1
		Algorithm 3 from Hörmann Leydold (Created by Dagpunar)
		
		@out return sample
	*/ 
	double algortihm_region1();

	/*
		Algorithm if the parameter region is
		sqrt(a * b)< = 2/3 * \sqrt{1 - \lambda} or abs(p) a <= 1
		Algorithm 1 from Hörmann Leydold 
		
		@out return sample
	*/ 
	double algortihm_region2();

	/*
		Algorithm if the parameter region is
		1 >) sqrt(a * b) >  2/3 * \sqrt{1 - \lambda} or abs(p) a <= 1
		Algorithm 2 from Hörmann Leydold 
		
		@out return sample
	*/ 
	double algortihm_region3();

	/*
	 * Computes x^{pa - 1} exp( - \beta/2 (x + 1/x) )
	 */ 
	double gig_propto(double);
  
		/*
	 * Computes x^{(pa - 1)/2} exp( - \beta/4 (x + 1/x) )
	 */ 
	double sqrt_gig_propto(double);
  
	/* Computes the ratio between to square root of two gig distribtuion
	 * @ param x_in the denominator term 
	 * @ param m_in the  numerator
	 * Computes sqrt_gig_propto(x_in) / sqrt_gig_propto(m_in)
	 */ 
	double sqrt_gig_ratio(double, double);

	/* Inverts returned random number if p<0 and then (for all p) scales by 1/alpha 
	 * @param x_in random sample
	 * @out scaled sample
	 */ 
  double scaleX(double);
  
public:
	gig() : U01_dist(0.0, 1.0), rgen(){};
  explicit gig(uint_fast64_t seed) : U01_dist(0.0, 1.0), rgen(seed){};

	/*
		We need a nan value to return in cases of bad parameters
		@out returns a nan value
	*/ 
  double nan();
  
  //update seed of the random generator
  void seed(uint_fast64_t seed){ rgen.seed(seed); }

	/*
	 * sampling gig random variable
	 * @param p,a,b GIG-parameters
	 * @return the random sample
	 */
	double sample(const double,const double,const double);

}; //class gig
#endif /* GIG_PAR_H */
