#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <Rcpp.h>

// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace boost::numeric::odeint;

// const double sigma = 10.0, r = 28.0, b = 8.0 / 3.0;

typedef boost::array< double , 2> state_type;

class pz_ {

  double alpha, c, e, ml, mq;


  public:
    pz_( double alpha, double c, double e, double ml, double mq ) : alpha(alpha), c(c), e(e), ml(ml), mq(mq) { }

  void operator() ( const state_type &x , state_type &dxdt , const double /* t */ )
  {
    dxdt[0] = alpha * x[0] - c * x[0] * x[1];
    dxdt[1] = e * c * x[0] * x[1] - ml * x[1] - mq * x[1] * x[1];
  }
};

struct push_back_state_and_time
{
  std::vector< state_type >& m_states;
  std::vector< double >& m_times;

  push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
  : m_states( states ) , m_times( times ) { }

  void operator()( const state_type &x , double t )
  {
    m_states.push_back( x );
    m_times.push_back( t );
  }
};

using namespace Rcpp;

// // [[Rcpp::export]]
// NumericVector one_step_pz_(double P, double Z,  double t, double alpha, double c, double e, double ml, double mq ) {
//   state_type x = { P , Z }; // initial conditions
//   vector<state_type> x_vec;
//   vector<double> times;
//   pz_ pz_instance(alpha, c, e, ml, mq);
//   size_t steps = integrate( pz_instance , x , t , t + 1, 1.0 , push_back_state_and_time(x_vec, times));
//   NumericVector result = NumericVector::create(0, 0);
//   result(0) = x_vec[steps][0];
//   result(1) = x_vec[steps][1];
//   return result;
// }

// [[Rcpp::export]]
NumericMatrix one_step_pz_vector(NumericMatrix xparticles, NumericVector alphas, double t, NumericVector parameters){
  double c = parameters[0];
  double e = parameters[1];
  double  ml = parameters[2];
  double  mq = parameters[3];
  NumericMatrix result(2, xparticles.cols());
  for (int i = 0; i < xparticles.cols(); i++){
    double P = xparticles(0, i);
    double Z = xparticles(1, i);
    state_type x = { P , Z }; // initial conditions
    vector<state_type> x_vec;
    vector<double> times;
    pz_ pz_instance(alphas(i), c, e, ml, mq);
    size_t steps = integrate( pz_instance , x , t , t + 1, 1.0 , push_back_state_and_time(x_vec, times));
    result(0, i) = x_vec[steps][0];
    result(1, i) = x_vec[steps][1];
  }
  return result;
}

// [[Rcpp::export]]
NumericVector pz_generate_randomness_cpp(int nparticles, int datalength){
  RNGScope scope;
  NumericVector normal_draws = rnorm((2 + datalength) * nparticles, 0, 1);
  return normal_draws;
}

// [[Rcpp::export]]
NumericVector pz_perturb_randomness_cpp(const NumericVector & randomness, double rho){
  RNGScope scope;
  int l = randomness.size();
  NumericVector newrand(l);
  double v = sqrt(1.0 - rho*rho);
  NumericVector normal_draws = rnorm(l, 0, 1);
  for (int i = 0; i < l; i ++){
    newrand(i) =  rho * randomness(i) + v * normal_draws(i);
  }
  return newrand;
}
