// manual gradient for the calibration model

#define EIGEN_DONT_PARALLELIZE

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <typeinfo>
#include <cmath>

#include <boost/exception/all.hpp>

#include <stan/model/prob_grad.hpp>
#include <stan/prob/distributions.hpp>
#include <stan/math/matrix.hpp>
#include <stan/math.hpp>
#include <stan/io/dump.hpp>
#include <stan/io/reader.hpp>
#include <stan/io/writer.hpp>
#include <stan/io/csv_writer.hpp>

#include <unsupported/Eigen/KroneckerProduct>

//#include "timer.hpp"

namespace cal_pl_model_namespace {

using std::vector;
using std::string;
using std::stringstream;
using stan::model::prob_grad;
using stan::math::get_base1;
using stan::math::initialize;
using stan::math::stan_print;
using stan::math::lgamma;
using stan::io::dump;
using std::istream;
using namespace stan::math;
using namespace stan::prob;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

class cal_pl_model : public prob_grad {
private:
    int K;
    int N_cores;
    int N_cells;
    int N_pot;
    vector<vector<int> > y;
    vector<int> idx_cores;
    vector<vector_d> r;
    matrix_d d;
    matrix_d d_pot;
public:
    cal_pl_model(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad::prob_grad(0) {
        static const char* function__ = "cal_pl_model_namespace::cal_pl_model(%1%)";
        (void) function__; // dummy call to supress warning
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        context__.validate_dims("data initialization", "K", "int", context__.to_vec());
        K = int(0);
        vals_i__ = context__.vals_i("K");
        pos__ = 0;
        K = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N_cores", "int", context__.to_vec());
        N_cores = int(0);
        vals_i__ = context__.vals_i("N_cores");
        pos__ = 0;
        N_cores = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N_cells", "int", context__.to_vec());
        N_cells = int(0);
        vals_i__ = context__.vals_i("N_cells");
        pos__ = 0;
        N_cells = vals_i__[pos__++];
        context__.validate_dims("data initialization", "N_pot", "int", context__.to_vec());
        N_pot = int(0);
        vals_i__ = context__.vals_i("N_pot");
        pos__ = 0;
        N_pot = vals_i__[pos__++];
        context__.validate_dims("data initialization", "y", "int", context__.to_vec(N_cores,K));
        stan::math::validate_non_negative_index("y", "N_cores", N_cores);
        stan::math::validate_non_negative_index("y", "K", K);
        y = std::vector<std::vector<int> >(N_cores,std::vector<int>(K,int(0)));
        vals_i__ = context__.vals_i("y");
        pos__ = 0;
        size_t y_limit_1__ = K;
        for (size_t i_1__ = 0; i_1__ < y_limit_1__; ++i_1__) {
            size_t y_limit_0__ = N_cores;
            for (size_t i_0__ = 0; i_0__ < y_limit_0__; ++i_0__) {
                y[i_0__][i_1__] = vals_i__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "idx_cores", "int", context__.to_vec(N_cores));
        stan::math::validate_non_negative_index("idx_cores", "N_cores", N_cores);
        idx_cores = std::vector<int>(N_cores,int(0));
        vals_i__ = context__.vals_i("idx_cores");
        pos__ = 0;
        size_t idx_cores_limit_0__ = N_cores;
        for (size_t i_0__ = 0; i_0__ < idx_cores_limit_0__; ++i_0__) {
            idx_cores[i_0__] = vals_i__[pos__++];
        }
        stan::math::validate_non_negative_index("r", "N_cells", N_cells);
        stan::math::validate_non_negative_index("r", "K", K);
        r = std::vector<vector_d>(N_cells,vector_d(K));
        context__.validate_dims("data initialization", "r", "vector_d", context__.to_vec(N_cells,K));
        vals_r__ = context__.vals_r("r");
        pos__ = 0;
        size_t r_i_vec_lim__ = K;
        for (size_t i_vec__ = 0; i_vec__ < r_i_vec_lim__; ++i_vec__) {
            size_t r_limit_0__ = N_cells;
            for (size_t i_0__ = 0; i_0__ < r_limit_0__; ++i_0__) {
                r[i_0__][i_vec__] = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "d", "matrix_d", context__.to_vec(N_cells,N_cores));
        stan::math::validate_non_negative_index("d", "N_cells", N_cells);
        stan::math::validate_non_negative_index("d", "N_cores", N_cores);
        d = matrix_d(N_cells,N_cores);
        vals_r__ = context__.vals_r("d");
        pos__ = 0;
        size_t d_m_mat_lim__ = N_cells;
        size_t d_n_mat_lim__ = N_cores;
        for (size_t n_mat__ = 0; n_mat__ < d_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < d_m_mat_lim__; ++m_mat__) {
                d(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }
        context__.validate_dims("data initialization", "d_pot", "matrix_d", context__.to_vec(N_pot,2));
        stan::math::validate_non_negative_index("d_pot", "N_pot", N_pot);
        stan::math::validate_non_negative_index("d_pot", "2", 2);
        d_pot = matrix_d(N_pot,2);
        vals_r__ = context__.vals_r("d_pot");
        pos__ = 0;
        size_t d_pot_m_mat_lim__ = N_pot;
        size_t d_pot_n_mat_lim__ = 2;
        for (size_t n_mat__ = 0; n_mat__ < d_pot_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < d_pot_m_mat_lim__; ++m_mat__) {
                d_pot(m_mat__,n_mat__) = vals_r__[pos__++];
            }
        }

        // validate data
        try { 
            check_greater_or_equal(function__,K,0,"K");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of K: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,N_cores,0,"N_cores");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of N_cores: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,N_cells,0,"N_cells");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of N_cells: ") + std::string(e.what()));
        };
        try { 
            check_greater_or_equal(function__,N_pot,0,"N_pot");
        } catch (std::domain_error& e) { 
            throw std::domain_error(std::string("Invalid value of N_pot: ") + std::string(e.what()));
        };


        // validate transformed data

        // set parameter ranges
        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += K;
        ++num_params_r__;
        ++num_params_r__;
        ++num_params_r__;
    }

    ~cal_pl_model() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;


        if (!(context__.contains_r("phi")))
            throw std::runtime_error("variable phi missing");
        vals_r__ = context__.vals_r("phi");
        pos__ = 0U;
        context__.validate_dims("initialization", "phi", "vector_d", context__.to_vec(K));
        vector_d phi(K);
        for (int j1__ = 0U; j1__ < K; ++j1__)
            phi(j1__) = vals_r__[pos__++];
        try { writer__.vector_lub_unconstrain(0.01,300,phi); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable phi: ") + e.what()); }

        if (!(context__.contains_r("gamma")))
            throw std::runtime_error("variable gamma missing");
        vals_r__ = context__.vals_r("gamma");
        pos__ = 0U;
        context__.validate_dims("initialization", "gamma", "double", context__.to_vec());
        double gamma(0);
        gamma = vals_r__[pos__++];
        try { writer__.scalar_lub_unconstrain(0,1,gamma); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable gamma: ") + e.what()); }

        if (!(context__.contains_r("a")))
            throw std::runtime_error("variable a missing");
        vals_r__ = context__.vals_r("a");
        pos__ = 0U;
        context__.validate_dims("initialization", "a", "double", context__.to_vec());
        double a(0);
        a = vals_r__[pos__++];
        try { writer__.scalar_lub_unconstrain(0,500,a); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable a: ") + e.what()); }

        if (!(context__.contains_r("b")))
            throw std::runtime_error("variable b missing");
        vals_r__ = context__.vals_r("b");
        pos__ = 0U;
        context__.validate_dims("initialization", "b", "double", context__.to_vec());
        double b(0);
        b = vals_r__[pos__++];
        try { writer__.scalar_lub_unconstrain(2,100,b); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable b: ") + e.what()); }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


 inline double lub_transform(const double x, const double lb, const double ub,
				double &lja, double &ja, double &dj) const
    {
      double inv_logit_x;
      if (x > 0) {
        double exp_minus_x = exp(-x);
	double exp_minus_x_p1 = exp_minus_x + 1.0;
        inv_logit_x = 1.0 / (1.0 + exp_minus_x);
        lja = log(ub - lb) - x - 2 * log1p(exp_minus_x);
	ja  = (ub - lb) * exp_minus_x / (exp_minus_x_p1 * exp_minus_x_p1);
	dj  = -1.0 + 2.0 * exp_minus_x / (1 + exp_minus_x);
        if ((x < std::numeric_limits<double>::infinity())
            && (inv_logit_x == 1))
	  inv_logit_x = 1 - 1e-15;
      } else {
        double exp_x = exp(x);
	double exp_x_p1 = exp_x + 1.0;
        inv_logit_x = 1.0 - 1.0 / (1.0 + exp_x);
        lja = log(ub - lb) + x - 2 * log1p(exp_x);
	ja  = (ub - lb) * exp_x / (exp_x_p1 * exp_x_p1);
	dj  = -1.0 + 2.0 / (exp_x + 1.0);
        if ((x > -std::numeric_limits<double>::infinity())
            && (inv_logit_x== 0))
	  inv_logit_x = 1e-15;
      }
      return lb + (ub - lb) * inv_logit_x;
    }

    double normal_log_double(const vector_d& y, const double mu, const double sigma) const {
      double lp = 0.0;
      double inv_sigma = 1.0/sigma;
      double log_sigma = log(sigma);

      int size_y = y.size();

      for (int n = 0; n < size_y; n++) {
	const double y_minus_mu_over_sigma = (y[n] - mu) * inv_sigma;
	const double y_minus_mu_over_sigma_squared = y_minus_mu_over_sigma * y_minus_mu_over_sigma;
	lp -= 0.5 * y_minus_mu_over_sigma_squared;
      }

      return lp;
    }

    double normal_log_double(const double y, const double mu, const double sigma, const int sigma_fixed) const {
      double lp = 0.0;
      double inv_sigma = 1.0/sigma;
      double log_sigma = log(sigma);

      const double y_minus_mu_over_sigma = (y - mu) * inv_sigma;
      const double y_minus_mu_over_sigma_squared = y_minus_mu_over_sigma * y_minus_mu_over_sigma;

      lp -= 0.5 * y_minus_mu_over_sigma_squared;

      if (sigma_fixed != 1){
	lp -= log_sigma;
      }

      return lp;
    }

  ///////////////////////////////////////////////////////////////////////////////////////////////

    template <bool propto, bool jacobian>
    double log_prob_grad(Eigen::VectorXd& params_r,
                         Eigen::VectorXd& gradient,
                         std::ostream* pstream = 0) const {

      double lp = 0.0;

      vector<double> vec_params_r;
      vector<int>    vec_params_i;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      stan::io::reader<double> in(vec_params_r, vec_params_i);


        // // model parameters
        // stan::io::reader<T__> in__(params_r__,params_i__);
      vector_d phi(K), phi_lja(K), phi_ja(K), phi_dj(K);
      for (int i=0; i<K; ++i)
	phi[i] = lub_transform(in.scalar(), 0.01, 300, phi_lja[i], phi_ja[i], phi_dj[i]);
      
      if (jacobian)
	for (int i=0; i<K; ++i)
	  lp += phi_lja[i];

      double gamma,  gamma_lja, gamma_ja, gamma_dj;
      gamma = lub_transform(in.scalar(), 0, 1, gamma_lja, gamma_ja, gamma_dj);

      if (jacobian)
      	lp += gamma_lja;

      double a,  a_lja, a_ja, a_dj;
      a = lub_transform(in.scalar(), 0, 500, a_lja, a_ja, a_dj);

      if (jacobian)
      	lp += a_lja;

      double b,  b_lja, b_ja, b_dj;
      b = lub_transform(in.scalar(), 2, 100, b_lja, b_ja, b_dj);

      if (jacobian)
      	lp += b_lja;
      
      std::cout << "LP after constraints : " << lp << std::endl;

      //
      // compute log probability
      //

      matrix_d w;
      matrix_d w_pot;
      vector<vector_d> r_new(K);
      vector_d out_sum(K);
      double sum_w;
      double sum_w_pot;
      double max_r_new;
      int max_r_new_idx;
      
      vector_d N(N_cores);
      vector_d A(N_cores);;
      vector_d alpha(K);    

      N.fill(0);
      A.fill(0);

      sum_w_pot = 0;
      for (int v=0; v<N_pot; v++){
	sum_w_pot += d_pot(v,2) * (b-2) * (b-1) / ( 2 * pi() * a * a ) * pow( 1 + d_pot(v,1) / a, -b) ;
      }

      for (int i=0; i<N_cells; ++i)
  	for (int j; j<N_cores; ++j)
  	  w(i,j) = (b-2) * (b-1) / ( 2 * pi() * a * a ) * pow( 1 + d(i,j) / a, -b) ;
      
      out_sum.fill(0);
      sum_w = 0;

      for (int i; i<N_cores; ++i){
	for (int k; k<K; ++k){
	  // local piece
	  r_new[k][i] = gamma*r[k][idx_cores[i]];
       
	  for (int j; j<N_cells; ++j){ // change N_hood to N_cells
	    if (j != idx_cores[i]){
	      out_sum[k] += out_sum[k] + w(j,i)*r[k][j];
	      sum_w   += sum_w + w(j,i);
	    }  
	  }
  
      // for (int k; k<K; ++k){
	//local vs. non-local
	r_new[k][i] += out_sum[k] * (1-gamma) / sum_w_pot;
      // }
	}      

  //   // hacky!
  //   // find taxon with highest proportional value
  //   max_r_new <- 0;
  //   for (k in 1:K){
  //     if (r_new[i,k] > max_r_new){
  //       max_r_new     <- r_new[i,k];
  //       max_r_new_idx <- k;
  //      }
  //   }
    
	for (int k; k<K; ++k){
	  if (r_new[k][i] == 0){
	    //r_new[i,k] <- 0.0001;
	    //r_new[i,max_r_new_idx] <- r_new[i,max_r_new_idx] - 0.0001;
        
	    std::cout << "warning: zero proportion; core: " <<  i << "; taxon: " <<  k << " -> adjusting" << std::endl;
	  }
	}
      
	for (int k; k<K; ++k){
	  alpha[k] = phi[k] * r_new[k][i];
	  N[i] += y[i][k];     
      
	}

	A[i] = alpha.sum();
	lp += lgamma(N[i] + 1) + lgamma(A[i]) - lgamma(N[i] + A[i]);
	for (int k; k<K; ++k) 
	  lp += - lgamma(y[i][k] + 1) + lgamma(y[i][k] + alpha[k]) - lgamma(alpha[k]);

      }
    
	gradient.fill(0.0);

	for (int k=0; k<K; ++k){
	  for (int i=0; i<N_cores; i++){
	    gradient[k] = (-digamma(A[i] + N[i]) );
	  }
	}

	return lp;
    }


  template <bool propto, bool jacobian>
  double log_prob(vector<double>& params_r,
		  vector<int>& params_i,
		  std::ostream* pstream = 0) const {

    Eigen::VectorXd evec_params_r(params_r.size());
    Eigen::VectorXd evec_gradient(params_r.size());

    for (int i=0; i<params_r.size(); i++) evec_params_r[i] = params_r[i];
    double lp = log_prob_grad<propto, jacobian>(evec_params_r, evec_gradient, pstream);

    return lp;
  }

  template <bool propto__, bool jacobian__, typename T__>
  T__ log_prob(vector<T__>& params_r__,
	       vector<int>& params_i__,
	       std::ostream* pstream__ = 0) const {

    throw "log_prob called";

  }

  template <bool propto, bool jacobian, typename T>
  T log_prob(Eigen::Matrix<T,Eigen::Dynamic,1>& params_r,
	     std::ostream* pstream = 0) const {

    throw "log_prob called";

  }


    // template <bool propto__, bool jacobian__, typename T__>
    // T__ log_prob(vector<T__>& params_r__,
    //              vector<int>& params_i__,
    //              std::ostream* pstream__ = 0) const {



    //     T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    //     (void) DUMMY_VAR__;  // suppress unused var warning

    //     T__ lp__(0.0);
    //     stan::math::accumulator<T__> lp_accum__;

    //     // model parameters
    //     stan::io::reader<T__> in__(params_r__,params_i__);

    //     Eigen::Matrix<T__,Eigen::Dynamic,1>  phi;
    //     (void) phi;   // dummy to suppress unused var warning
    //     if (jacobian__)
    //         phi = in__.vector_lub_constrain(0.01,300,K,lp__);
    //     else
    //         phi = in__.vector_lub_constrain(0.01,300,K);

    //     T__ gamma;
    //     (void) gamma;   // dummy to suppress unused var warning
    //     if (jacobian__)
    //         gamma = in__.scalar_lub_constrain(0,1,lp__);
    //     else
    //         gamma = in__.scalar_lub_constrain(0,1);

    //     T__ a;
    //     (void) a;   // dummy to suppress unused var warning
    //     if (jacobian__)
    //         a = in__.scalar_lub_constrain(0,500,lp__);
    //     else
    //         a = in__.scalar_lub_constrain(0,500);

    //     T__ b;
    //     (void) b;   // dummy to suppress unused var warning
    //     if (jacobian__)
    //         b = in__.scalar_lub_constrain(2,100,lp__);
    //     else
    //         b = in__.scalar_lub_constrain(2,100);


    //     // transformed parameters

    //     // initialized transformed params to avoid seg fault on val access
        

    //     // validate transformed parameters

    //     const char* function__ = "validate transformed params %1%";
    //     (void) function__; // dummy to suppress unused var warning
    //     // model body
    //     {
    //         Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic>  w(N_cells,N_cores);
    //         (void) w;   // dummy to suppress unused var warning
    //         stan::math::fill(w,DUMMY_VAR__);
    //         Eigen::Matrix<T__,Eigen::Dynamic,1>  w_pot(N_pot);
    //         (void) w_pot;   // dummy to suppress unused var warning
    //         stan::math::fill(w_pot,DUMMY_VAR__);
    //         vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > r_new(N_cores, (Eigen::Matrix<T__,Eigen::Dynamic,1> (K)));
    //         stan::math::fill(r_new,DUMMY_VAR__);
    //         Eigen::Matrix<T__,Eigen::Dynamic,1>  out_sum(K);
    //         (void) out_sum;   // dummy to suppress unused var warning
    //         stan::math::fill(out_sum,DUMMY_VAR__);
    //         T__ sum_w;
    //         (void) sum_w;   // dummy to suppress unused var warning
    //         T__ sum_w_pot;
    //         (void) sum_w_pot;   // dummy to suppress unused var warning
    //         T__ max_r_new;
    //         (void) max_r_new;   // dummy to suppress unused var warning
    //         int max_r_new_idx(0);
    //         (void) max_r_new_idx;   // dummy to suppress unused var warning
    //         stan::math::initialize(w, DUMMY_VAR__);
    //         stan::math::initialize(w_pot, DUMMY_VAR__);
    //         stan::math::initialize(r_new, DUMMY_VAR__);
    //         stan::math::initialize(out_sum, DUMMY_VAR__);
    //         stan::math::initialize(sum_w, DUMMY_VAR__);
    //         stan::math::initialize(sum_w_pot, DUMMY_VAR__);
    //         stan::math::initialize(max_r_new, DUMMY_VAR__);
    //         lp_accum__.add(uniform_log<propto__>(phi, 0.01, 300));
    //         lp_accum__.add(uniform_log<propto__>(gamma, 0, 1));
    //         lp_accum__.add(uniform_log<propto__>(a, 0, 500));
    //         lp_accum__.add(uniform_log<propto__>(b, 2, 100));
    //         stan::math::assign(sum_w_pot, 0);
    //         for (int v = 1; v <= N_pot; ++v) {
    //             stan::math::assign(sum_w_pot, (sum_w_pot + ((((get_base1(d_pot,v,2,"d_pot",1) * (b - 2)) * (b - 1)) / (((2 * pi()) * a) * a)) * pow((1 + (get_base1(d_pot,v,1,"d_pot",1) / a)),-(b)))));
    //         }
    //         for (int i = 1; i <= N_cells; ++i) {
    //             for (int j = 1; j <= N_cores; ++j) {
    //                 stan::math::assign(get_base1_lhs(w,i,j,"w",1), ((((b - 2) * (b - 1)) / (((2 * pi()) * a) * a)) * pow((1 + (get_base1(d,i,j,"d",1) / a)),-(b))));
    //             }
    //         }
    //         for (int i = 1; i <= N_cores; ++i) {
    //             stan::math::assign(get_base1_lhs(r_new,i,"r_new",1), multiply(gamma,get_base1(r,get_base1(idx_cores,i,"idx_cores",1),"r",1)));
    //             for (int k = 1; k <= K; ++k) {
    //                 stan::math::assign(get_base1_lhs(out_sum,k,"out_sum",1), 0);
    //             }
    //             stan::math::assign(sum_w, 0);
    //             for (int j = 1; j <= N_cells; ++j) {
    //                 if (as_bool(logical_neq(j,get_base1(idx_cores,i,"idx_cores",1)))) {
    //                     stan::math::assign(out_sum, add(out_sum,multiply(get_base1(w,j,i,"w",1),get_base1(r,j,"r",1))));
    //                     stan::math::assign(sum_w, (sum_w + get_base1(w,j,i,"w",1)));
    //                 }
    //             }
    //             stan::math::assign(get_base1_lhs(r_new,i,"r_new",1), add(get_base1(r_new,i,"r_new",1),divide(multiply(out_sum,(1 - gamma)),sum_w_pot)));
    //             stan::math::assign(max_r_new, 0);
    //             for (int k = 1; k <= K; ++k) {
    //                 if (as_bool(logical_gt(get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2),max_r_new))) {
    //                     stan::math::assign(max_r_new, get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2));
    //                     stan::math::assign(max_r_new_idx, k);
    //                 }
    //             }
    //             for (int k = 1; k <= K; ++k) {
    //                 if (as_bool(logical_eq(get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2),0))) {
    //                     stan::math::assign(get_base1_lhs(get_base1_lhs(r_new,i,"r_new",1),k,"r_new",2), 0.0001);
    //                     stan::math::assign(get_base1_lhs(get_base1_lhs(r_new,i,"r_new",1),max_r_new_idx,"r_new",2), (get_base1(get_base1(r_new,i,"r_new",1),max_r_new_idx,"r_new",2) - 0.0001));
    //                     if (pstream__) {
    //                         stan_print(pstream__,"warning: zero proportion; core: ");
    //                         stan_print(pstream__,i);
    //                         stan_print(pstream__,"; taxon: ");
    //                         stan_print(pstream__,k);
    //                         stan_print(pstream__," -> adjusting");
    //                         *pstream__ << std::endl;
    //                     }
    //                 }
    //             }
    //             {
    //                 T__ N;
    //                 (void) N;   // dummy to suppress unused var warning
    //                 T__ A;
    //                 (void) A;   // dummy to suppress unused var warning
    //                 Eigen::Matrix<T__,Eigen::Dynamic,1>  alpha(K);
    //                 (void) alpha;   // dummy to suppress unused var warning
    //                 stan::math::fill(alpha,DUMMY_VAR__);
    //                 stan::math::initialize(N, DUMMY_VAR__);
    //                 stan::math::initialize(A, DUMMY_VAR__);
    //                 stan::math::initialize(alpha, DUMMY_VAR__);
    //                 stan::math::assign(alpha, elt_multiply(phi,get_base1(r_new,i,"r_new",1)));
    //                 stan::math::assign(A, sum(alpha));
    //                 stan::math::assign(N, sum(get_base1(y,i,"y",1)));
    //                 lp_accum__.add(((lgamma((N + 1)) + lgamma(A)) - lgamma((N + A))));
    //                 for (int k = 1; k <= K; ++k) {
    //                     lp_accum__.add(((-(lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + 1))) + lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + get_base1(alpha,k,"alpha",1)))) - lgamma(get_base1(alpha,k,"alpha",1))));
    //                 }
    //             }
    //         }
    //     }

    //     lp_accum__.add(lp__);
    //     return lp_accum__.sum();

    // } // log_prob()

    // template <bool propto, bool jacobian, typename T>
    // T log_prob(Eigen::Matrix<T,Eigen::Dynamic,1>& params_r,
    //            std::ostream* pstream = 0) const {
    //   std::vector<T> vec_params_r;
    //   vec_params_r.reserve(params_r.size());
    //   for (int i = 0; i < params_r.size(); ++i)
    //     vec_params_r.push_back(params_r(i));
    //   std::vector<int> vec_params_i;
    //   return log_prob<propto,jacobian,T>(vec_params_r, vec_params_i, pstream);
    // }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("phi");
        names__.push_back("gamma");
        names__.push_back("a");
        names__.push_back("b");
        names__.push_back("log_lik");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_cores);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "cal_pl_model_namespace::write_array(%1%)";
        (void) function__; // dummy call to supress warning
        // read-transform, write parameters
        vector_d phi = in__.vector_lub_constrain(0.01,300,K);
        double gamma = in__.scalar_lub_constrain(0,1);
        double a = in__.scalar_lub_constrain(0,500);
        double b = in__.scalar_lub_constrain(2,100);
        for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(phi[k_0__]);
        }
        vars__.push_back(gamma);
        vars__.push_back(a);
        vars__.push_back(b);

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__; // dummy call to supress warning
        stan::math::accumulator<double> lp_accum__;



        // validate transformed parameters

        // write transformed parameters

        if (!include_gqs__) return;
        // declare and define generated quantities
        vector_d log_lik(N_cores);
        (void) log_lik;   // dummy to suppress unused var warning

        {
            matrix_d w(N_cells,N_cores);
            (void) w;   // dummy to suppress unused var warning
            vector_d w_pot(N_pot);
            (void) w_pot;   // dummy to suppress unused var warning
            vector<vector_d> r_new(N_cores, (vector_d(K)));
            vector_d out_sum(K);
            (void) out_sum;   // dummy to suppress unused var warning
            double sum_w(0.0);
            (void) sum_w;   // dummy to suppress unused var warning
            double sum_w_pot(0.0);
            (void) sum_w_pot;   // dummy to suppress unused var warning
            double max_r_new(0.0);
            (void) max_r_new;   // dummy to suppress unused var warning
            int max_r_new_idx(0);
            (void) max_r_new_idx;   // dummy to suppress unused var warning
            double N(0.0);
            (void) N;   // dummy to suppress unused var warning
            double A(0.0);
            (void) A;   // dummy to suppress unused var warning
            vector_d alpha(K);
            (void) alpha;   // dummy to suppress unused var warning
            stan::math::initialize(w, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(w_pot, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(r_new, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(out_sum, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(sum_w, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(sum_w_pot, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(max_r_new, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(N, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(A, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(alpha, std::numeric_limits<double>::quiet_NaN());
            stan::math::assign(sum_w_pot, 0);
            for (int v = 1; v <= N_pot; ++v) {
                stan::math::assign(sum_w_pot, (sum_w_pot + ((((get_base1(d_pot,v,2,"d_pot",1) * (b - 2)) * (b - 1)) / (((2 * pi()) * a) * a)) * pow((1 + (get_base1(d_pot,v,1,"d_pot",1) / a)),-(b)))));
            }
            for (int i = 1; i <= N_cells; ++i) {
                for (int j = 1; j <= N_cores; ++j) {
                    stan::math::assign(get_base1_lhs(w,i,j,"w",1), ((((b - 2) * (b - 1)) / (((2 * pi()) * a) * a)) * pow((1 + (get_base1(d,i,j,"d",1) / a)),-(b))));
                }
            }
            for (int i = 1; i <= N_cores; ++i) {
                stan::math::assign(get_base1_lhs(r_new,i,"r_new",1), multiply(gamma,get_base1(r,get_base1(idx_cores,i,"idx_cores",1),"r",1)));
                for (int k = 1; k <= K; ++k) {
                    stan::math::assign(get_base1_lhs(out_sum,k,"out_sum",1), 0);
                }
                stan::math::assign(sum_w, 0);
                for (int j = 1; j <= N_cells; ++j) {
                    if (as_bool(logical_neq(j,get_base1(idx_cores,i,"idx_cores",1)))) {
                        stan::math::assign(out_sum, add(out_sum,multiply(get_base1(w,j,i,"w",1),get_base1(r,j,"r",1))));
                        stan::math::assign(sum_w, (sum_w + get_base1(w,j,i,"w",1)));
                    }
                }
                stan::math::assign(get_base1_lhs(r_new,i,"r_new",1), add(get_base1(r_new,i,"r_new",1),divide(multiply(out_sum,(1 - gamma)),sum_w_pot)));
                stan::math::assign(max_r_new, 0);
                for (int k = 1; k <= K; ++k) {
                    if (as_bool(logical_gt(get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2),max_r_new))) {
                        stan::math::assign(max_r_new, get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2));
                        stan::math::assign(max_r_new_idx, k);
                    }
                }
                for (int k = 1; k <= K; ++k) {
                    if (as_bool(logical_eq(get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2),0))) {
                        stan::math::assign(get_base1_lhs(get_base1_lhs(r_new,i,"r_new",1),k,"r_new",2), 0.0001);
                        stan::math::assign(get_base1_lhs(get_base1_lhs(r_new,i,"r_new",1),max_r_new_idx,"r_new",2), (get_base1(get_base1(r_new,i,"r_new",1),max_r_new_idx,"r_new",2) - 0.0001));
                        if (pstream__) {
                            stan_print(pstream__,"warning: zero proportion; core: ");
                            stan_print(pstream__,i);
                            stan_print(pstream__,"; taxon: ");
                            stan_print(pstream__,k);
                            stan_print(pstream__," -> adjusting");
                            *pstream__ << std::endl;
                        }
                    }
                }
                stan::math::assign(alpha, elt_multiply(phi,get_base1(r_new,i,"r_new",1)));
                stan::math::assign(A, sum(alpha));
                stan::math::assign(N, sum(get_base1(y,i,"y",1)));
                stan::math::assign(get_base1_lhs(log_lik,i,"log_lik",1), ((lgamma((N + 1)) + lgamma(A)) - lgamma((N + A))));
                for (int k = 1; k <= K; ++k) {
                    stan::math::assign(get_base1_lhs(log_lik,i,"log_lik",1), (((get_base1(log_lik,i,"log_lik",1) - lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + 1))) + lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + get_base1(alpha,k,"alpha",1)))) - lgamma(get_base1(alpha,k,"alpha",1))));
                }
            }
        }

        // validate generated quantities

        // write generated quantities
        for (int k_0__ = 0; k_0__ < N_cores; ++k_0__) {
            vars__.push_back(log_lik[k_0__]);
        }

    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }


    void write_csv_header(std::ostream& o__) const {
        stan::io::csv_writer writer__(o__);
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            writer__.comma();
            o__ << "phi" << '.' << k_0__;
        }
        writer__.comma();
        o__ << "gamma";
        writer__.comma();
        o__ << "a";
        writer__.comma();
        o__ << "b";
        for (int k_0__ = 1; k_0__ <= N_cores; ++k_0__) {
            writer__.comma();
            o__ << "log_lik" << '.' << k_0__;
        }
        writer__.newline();
    }

    template <typename RNG>
    void write_csv(RNG& base_rng__,
                   std::vector<double>& params_r__,
                   std::vector<int>& params_i__,
                   std::ostream& o__,
                   std::ostream* pstream__ = 0) const {
        stan::io::reader<double> in__(params_r__,params_i__);
        stan::io::csv_writer writer__(o__);
        static const char* function__ = "cal_pl_model_namespace::write_csv(%1%)";
        (void) function__; // dummy call to supress warning
        // read-transform, write parameters
        vector_d phi = in__.vector_lub_constrain(0.01,300,K);
        writer__.write(phi);
        double gamma = in__.scalar_lub_constrain(0,1);
        writer__.write(gamma);
        double a = in__.scalar_lub_constrain(0,500);
        writer__.write(a);
        double b = in__.scalar_lub_constrain(2,100);
        writer__.write(b);

        // declare, define and validate transformed parameters
        double lp__ = 0.0;
        (void) lp__; // dummy call to supress warning
        stan::math::accumulator<double> lp_accum__;




        // write transformed parameters

        // declare and define generated quantities
        vector_d log_lik(N_cores);
        (void) log_lik;   // dummy to suppress unused var warning

        {
            matrix_d w(N_cells,N_cores);
            (void) w;   // dummy to suppress unused var warning
            vector_d w_pot(N_pot);
            (void) w_pot;   // dummy to suppress unused var warning
            vector<vector_d> r_new(N_cores, (vector_d(K)));
            vector_d out_sum(K);
            (void) out_sum;   // dummy to suppress unused var warning
            double sum_w(0.0);
            (void) sum_w;   // dummy to suppress unused var warning
            double sum_w_pot(0.0);
            (void) sum_w_pot;   // dummy to suppress unused var warning
            double max_r_new(0.0);
            (void) max_r_new;   // dummy to suppress unused var warning
            int max_r_new_idx(0);
            (void) max_r_new_idx;   // dummy to suppress unused var warning
            double N(0.0);
            (void) N;   // dummy to suppress unused var warning
            double A(0.0);
            (void) A;   // dummy to suppress unused var warning
            vector_d alpha(K);
            (void) alpha;   // dummy to suppress unused var warning
            stan::math::initialize(w, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(w_pot, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(r_new, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(out_sum, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(sum_w, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(sum_w_pot, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(max_r_new, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(N, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(A, std::numeric_limits<double>::quiet_NaN());
            stan::math::initialize(alpha, std::numeric_limits<double>::quiet_NaN());
            stan::math::assign(sum_w_pot, 0);
            for (int v = 1; v <= N_pot; ++v) {
                stan::math::assign(sum_w_pot, (sum_w_pot + ((((get_base1(d_pot,v,2,"d_pot",1) * (b - 2)) * (b - 1)) / (((2 * pi()) * a) * a)) * pow((1 + (get_base1(d_pot,v,1,"d_pot",1) / a)),-(b)))));
            }
            for (int i = 1; i <= N_cells; ++i) {
                for (int j = 1; j <= N_cores; ++j) {
                    stan::math::assign(get_base1_lhs(w,i,j,"w",1), ((((b - 2) * (b - 1)) / (((2 * pi()) * a) * a)) * pow((1 + (get_base1(d,i,j,"d",1) / a)),-(b))));
                }
            }
            for (int i = 1; i <= N_cores; ++i) {
                stan::math::assign(get_base1_lhs(r_new,i,"r_new",1), multiply(gamma,get_base1(r,get_base1(idx_cores,i,"idx_cores",1),"r",1)));
                for (int k = 1; k <= K; ++k) {
                    stan::math::assign(get_base1_lhs(out_sum,k,"out_sum",1), 0);
                }
                stan::math::assign(sum_w, 0);
                for (int j = 1; j <= N_cells; ++j) {
                    if (as_bool(logical_neq(j,get_base1(idx_cores,i,"idx_cores",1)))) {
                        stan::math::assign(out_sum, add(out_sum,multiply(get_base1(w,j,i,"w",1),get_base1(r,j,"r",1))));
                        stan::math::assign(sum_w, (sum_w + get_base1(w,j,i,"w",1)));
                    }
                }
                stan::math::assign(get_base1_lhs(r_new,i,"r_new",1), add(get_base1(r_new,i,"r_new",1),divide(multiply(out_sum,(1 - gamma)),sum_w_pot)));
                stan::math::assign(max_r_new, 0);
                for (int k = 1; k <= K; ++k) {
                    if (as_bool(logical_gt(get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2),max_r_new))) {
                        stan::math::assign(max_r_new, get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2));
                        stan::math::assign(max_r_new_idx, k);
                    }
                }
                for (int k = 1; k <= K; ++k) {
                    if (as_bool(logical_eq(get_base1(get_base1(r_new,i,"r_new",1),k,"r_new",2),0))) {
                        stan::math::assign(get_base1_lhs(get_base1_lhs(r_new,i,"r_new",1),k,"r_new",2), 0.0001);
                        stan::math::assign(get_base1_lhs(get_base1_lhs(r_new,i,"r_new",1),max_r_new_idx,"r_new",2), (get_base1(get_base1(r_new,i,"r_new",1),max_r_new_idx,"r_new",2) - 0.0001));
                        if (pstream__) {
                            stan_print(pstream__,"warning: zero proportion; core: ");
                            stan_print(pstream__,i);
                            stan_print(pstream__,"; taxon: ");
                            stan_print(pstream__,k);
                            stan_print(pstream__," -> adjusting");
                            *pstream__ << std::endl;
                        }
                    }
                }
                stan::math::assign(alpha, elt_multiply(phi,get_base1(r_new,i,"r_new",1)));
                stan::math::assign(A, sum(alpha));
                stan::math::assign(N, sum(get_base1(y,i,"y",1)));
                stan::math::assign(get_base1_lhs(log_lik,i,"log_lik",1), ((lgamma((N + 1)) + lgamma(A)) - lgamma((N + A))));
                for (int k = 1; k <= K; ++k) {
                    stan::math::assign(get_base1_lhs(log_lik,i,"log_lik",1), (((get_base1(log_lik,i,"log_lik",1) - lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + 1))) + lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + get_base1(alpha,k,"alpha",1)))) - lgamma(get_base1(alpha,k,"alpha",1))));
                }
            }
        }

        // validate generated quantities

        // write generated quantities
        writer__.write(log_lik);

        writer__.newline();
    }

    template <typename RNG>
    void write_csv(RNG& base_rng,
                   Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                   std::ostream& o,
                   std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<int> params_i_vec;  // dummy
      write_csv(base_rng, params_r_vec, params_i_vec, o, pstream);
    }

    static std::string model_name() {
        return "cal_pl_model";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "phi" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "gamma";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "a";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "b";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
        for (int k_0__ = 1; k_0__ <= N_cores; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "phi" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "gamma";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "a";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "b";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
        for (int k_0__ = 1; k_0__ <= N_cores; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "log_lik" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
    }

}; // model

} // namespace

namespace stan {
  namespace model {

    template <bool propto, bool jacobian_adjust_transform>
    double log_prob_grad(const cal_pl_model_namespace::cal_pl_model& model,
                         Eigen::VectorXd& params_r,
                         Eigen::VectorXd& gradient,
                         std::ostream* msgs = 0) {
      double lp = model.template log_prob_grad<propto, jacobian_adjust_transform>(params_r, gradient, msgs);
      return lp;

    }

    template <bool propto, bool jacobian_adjust_transform>
    double log_prob_grad(const cal_pl_model_namespace::cal_pl_model& model,
                         std::vector<double>& params_r,
                         std::vector<int>& params_i,
                         std::vector<double>& gradient,
                         std::ostream* msgs = 0) {

      Eigen::VectorXd evec_params_r(params_r.size());
      Eigen::VectorXd evec_gradient(params_r.size());

      for (int i=0; i<params_r.size(); i++) evec_params_r[i] = params_r[i];
      double lp = model.template log_prob_grad<propto, jacobian_adjust_transform>(evec_params_r, evec_gradient, msgs);

      gradient.resize(params_r.size());
      for (int i=0; i<params_r.size(); i++) gradient[i] = evec_gradient[i];
      return lp;

    }

  }
}

#include <stan/gm/command.hpp>

int main(int argc, const char* argv[]) {
    try {
        return stan::gm::command<cal_pl_model_namespace::cal_pl_model>(argc,argv);
    } catch (std::exception& e) {
        std::cerr << std::endl << "Exception: " << e.what() << std::endl;
        std::cerr << "Diagnostic information: " << std::endl << boost::diagnostic_information(e) << std::endl;
        return -1;
    }
}

