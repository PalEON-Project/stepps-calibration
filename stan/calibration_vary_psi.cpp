// Code generated by Stan version 2.1

#include <stan/model/model_header.hpp>

namespace calibration_vary_psi_model_namespace {

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

class calibration_vary_psi_model : public prob_grad {
private:
    int K;
    int N_cores;
    int N_cells;
    vector<vector<int> > y;
    vector<int> idx_cores;
    vector<vector_d> r;
    matrix_d d2;
public:
    calibration_vary_psi_model(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad::prob_grad(0) {
        static const char* function__ = "calibration_vary_psi_model_namespace::calibration_vary_psi_model(%1%)";
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
        context__.validate_dims("data initialization", "d2", "matrix_d", context__.to_vec(N_cells,N_cores));
        stan::math::validate_non_negative_index("d2", "N_cells", N_cells);
        stan::math::validate_non_negative_index("d2", "N_cores", N_cores);
        d2 = matrix_d(N_cells,N_cores);
        vals_r__ = context__.vals_r("d2");
        pos__ = 0;
        size_t d2_m_mat_lim__ = N_cells;
        size_t d2_n_mat_lim__ = N_cores;
        for (size_t n_mat__ = 0; n_mat__ < d2_n_mat_lim__; ++n_mat__) {
            for (size_t m_mat__ = 0; m_mat__ < d2_m_mat_lim__; ++m_mat__) {
                d2(m_mat__,n_mat__) = vals_r__[pos__++];
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


        // validate transformed data

        // set parameter ranges
        num_params_r__ = 0U;
        param_ranges_i__.clear();
        num_params_r__ += K;
        num_params_r__ += K;
        ++num_params_r__;
    }

    ~calibration_vary_psi_model() { }


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

        if (!(context__.contains_r("psi")))
            throw std::runtime_error("variable psi missing");
        vals_r__ = context__.vals_r("psi");
        pos__ = 0U;
        context__.validate_dims("initialization", "psi", "vector_d", context__.to_vec(K));
        vector_d psi(K);
        for (int j1__ = 0U; j1__ < K; ++j1__)
            psi(j1__) = vals_r__[pos__++];
        try { writer__.vector_lub_unconstrain(0.10000000000000001,2,psi); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable psi: ") + e.what()); }

        if (!(context__.contains_r("gamma")))
            throw std::runtime_error("variable gamma missing");
        vals_r__ = context__.vals_r("gamma");
        pos__ = 0U;
        context__.validate_dims("initialization", "gamma", "double", context__.to_vec());
        double gamma(0);
        gamma = vals_r__[pos__++];
        try { writer__.scalar_lub_unconstrain(0,1,gamma); } catch (std::exception& e) {  throw std::runtime_error(std::string("Error transforming variable gamma: ") + e.what()); }
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


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        // model parameters
        stan::io::reader<T__> in__(params_r__,params_i__);

        Eigen::Matrix<T__,Eigen::Dynamic,1>  phi;
        (void) phi;   // dummy to suppress unused var warning
        if (jacobian__)
            phi = in__.vector_lub_constrain(0.01,300,K,lp__);
        else
            phi = in__.vector_lub_constrain(0.01,300,K);

        Eigen::Matrix<T__,Eigen::Dynamic,1>  psi;
        (void) psi;   // dummy to suppress unused var warning
        if (jacobian__)
            psi = in__.vector_lub_constrain(0.10000000000000001,2,K,lp__);
        else
            psi = in__.vector_lub_constrain(0.10000000000000001,2,K);

        T__ gamma;
        (void) gamma;   // dummy to suppress unused var warning
        if (jacobian__)
            gamma = in__.scalar_lub_constrain(0,1,lp__);
        else
            gamma = in__.scalar_lub_constrain(0,1);


        // transformed parameters

        // initialized transformed params to avoid seg fault on val access
        

        // validate transformed parameters

        const char* function__ = "validate transformed params %1%";
        (void) function__; // dummy to suppress unused var warning
        // model body
        {
            vector<Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> > w(K, (Eigen::Matrix<T__,Eigen::Dynamic,Eigen::Dynamic> (N_cells,N_cores)));
            stan::math::fill(w,DUMMY_VAR__);
            vector<Eigen::Matrix<T__,Eigen::Dynamic,1> > r_new(N_cores, (Eigen::Matrix<T__,Eigen::Dynamic,1> (K)));
            stan::math::fill(r_new,DUMMY_VAR__);
            Eigen::Matrix<T__,Eigen::Dynamic,1>  out_sum(K);
            (void) out_sum;   // dummy to suppress unused var warning
            stan::math::fill(out_sum,DUMMY_VAR__);
            T__ sum_w;
            (void) sum_w;   // dummy to suppress unused var warning
            T__ max_r_new;
            (void) max_r_new;   // dummy to suppress unused var warning
            int max_r_new_idx(0);
            (void) max_r_new_idx;   // dummy to suppress unused var warning
            stan::math::initialize(w, DUMMY_VAR__);
            stan::math::initialize(r_new, DUMMY_VAR__);
            stan::math::initialize(out_sum, DUMMY_VAR__);
            stan::math::initialize(sum_w, DUMMY_VAR__);
            stan::math::initialize(max_r_new, DUMMY_VAR__);
            lp_accum__.add(uniform_log<propto__>(phi, 0.01, 300));
            lp_accum__.add(uniform_log<propto__>(gamma, 0, 1));
            for (int k = 1; k <= K; ++k) {
                lp_accum__.add(uniform_log<propto__>(get_base1(psi,k,"psi",1), 0.10000000000000001, 2));
            }
            if (pstream__) {
                stan_print(pstream__,psi);
                *pstream__ << std::endl;
            }
            for (int k = 1; k <= K; ++k) {
                stan::math::assign(get_base1_lhs(w,k,"w",1), exp(divide(minus(d2),square(get_base1(psi,k,"psi",1)))));
            }
            for (int i = 1; i <= N_cores; ++i) {
                stan::math::assign(get_base1_lhs(r_new,i,"r_new",1), multiply(gamma,get_base1(r,get_base1(idx_cores,i,"idx_cores",1),"r",1)));
                for (int k = 1; k <= K; ++k) {
                    stan::math::assign(get_base1_lhs(out_sum,k,"out_sum",1), 0);
                }
                stan::math::assign(sum_w, 0);
                for (int k = 1; k <= K; ++k) {
                    for (int j = 1; j <= N_cells; ++j) {
                        if (as_bool(logical_neq(j,get_base1(idx_cores,i,"idx_cores",1)))) {
                            stan::math::assign(get_base1_lhs(out_sum,k,"out_sum",1), (get_base1(out_sum,k,"out_sum",1) + (get_base1(get_base1(w,k,"w",1),j,i,"w",2) * get_base1(get_base1(r,j,"r",1),k,"r",2))));
                        }
                    }
                }
                stan::math::assign(sum_w, sum(out_sum));
                stan::math::assign(get_base1_lhs(r_new,i,"r_new",1), add(get_base1(r_new,i,"r_new",1),divide(multiply(out_sum,(1 - gamma)),sum_w)));
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
                {
                    T__ N;
                    (void) N;   // dummy to suppress unused var warning
                    T__ A;
                    (void) A;   // dummy to suppress unused var warning
                    Eigen::Matrix<T__,Eigen::Dynamic,1>  alpha(K);
                    (void) alpha;   // dummy to suppress unused var warning
                    stan::math::fill(alpha,DUMMY_VAR__);
                    stan::math::initialize(N, DUMMY_VAR__);
                    stan::math::initialize(A, DUMMY_VAR__);
                    stan::math::initialize(alpha, DUMMY_VAR__);
                    stan::math::assign(alpha, elt_multiply(phi,get_base1(r_new,i,"r_new",1)));
                    stan::math::assign(A, sum(alpha));
                    stan::math::assign(N, sum(get_base1(y,i,"y",1)));
                    lp_accum__.add(((lgamma((N + 1)) + lgamma(A)) - lgamma((N + A))));
                    for (int k = 1; k <= K; ++k) {
                        lp_accum__.add(((-(lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + 1))) + lgamma((get_base1(get_base1(y,i,"y",1),k,"y",2) + get_base1(alpha,k,"alpha",1)))) - lgamma(get_base1(alpha,k,"alpha",1))));
                    }
                }
            }
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T>
    T log_prob(Eigen::Matrix<T,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("phi");
        names__.push_back("psi");
        names__.push_back("gamma");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(K);
        dimss__.push_back(dims__);
        dims__.resize(0);
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
        static const char* function__ = "calibration_vary_psi_model_namespace::write_array(%1%)";
        (void) function__; // dummy call to supress warning
        // read-transform, write parameters
        vector_d phi = in__.vector_lub_constrain(0.01,300,K);
        vector_d psi = in__.vector_lub_constrain(0.10000000000000001,2,K);
        double gamma = in__.scalar_lub_constrain(0,1);
        for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(phi[k_0__]);
        }
        for (int k_0__ = 0; k_0__ < K; ++k_0__) {
            vars__.push_back(psi[k_0__]);
        }
        vars__.push_back(gamma);

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__; // dummy call to supress warning
        stan::math::accumulator<double> lp_accum__;



        // validate transformed parameters

        // write transformed parameters

        if (!include_gqs__) return;
        // declare and define generated quantities


        // validate generated quantities

        // write generated quantities
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
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            writer__.comma();
            o__ << "psi" << '.' << k_0__;
        }
        writer__.comma();
        o__ << "gamma";
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
        static const char* function__ = "calibration_vary_psi_model_namespace::write_csv(%1%)";
        (void) function__; // dummy call to supress warning
        // read-transform, write parameters
        vector_d phi = in__.vector_lub_constrain(0.01,300,K);
        writer__.write(phi);
        vector_d psi = in__.vector_lub_constrain(0.10000000000000001,2,K);
        writer__.write(psi);
        double gamma = in__.scalar_lub_constrain(0,1);
        writer__.write(gamma);

        // declare, define and validate transformed parameters
        double lp__ = 0.0;
        (void) lp__; // dummy call to supress warning
        stan::math::accumulator<double> lp_accum__;




        // write transformed parameters

        // declare and define generated quantities


        // validate generated quantities

        // write generated quantities
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
        return "calibration_vary_psi_model";
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
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "psi" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "gamma";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
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
        for (int k_0__ = 1; k_0__ <= K; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "psi" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }
        param_name_stream__.str(std::string());
        param_name_stream__ << "gamma";
        param_names__.push_back(param_name_stream__.str());

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }

}; // model

} // namespace

int main(int argc, const char* argv[]) {
    try {
        return stan::gm::command<calibration_vary_psi_model_namespace::calibration_vary_psi_model>(argc,argv);
    } catch (std::exception& e) {
        std::cerr << std::endl << "Exception: " << e.what() << std::endl;
        std::cerr << "Diagnostic information: " << std::endl << boost::diagnostic_information(e) << std::endl;
        return -1;
    }
}

