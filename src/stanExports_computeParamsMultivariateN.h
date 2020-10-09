// Generated by rstantools.  Do not edit by hand.

/*
    BSure is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BSure is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with BSure.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_computeParamsMultivariateN_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_computeParamsMultivariateN");
    reader.add_event(22, 20, "end", "model_computeParamsMultivariateN");
    return reader;
}
#include <stan_meta_header.hpp>
class model_computeParamsMultivariateN
  : public stan::model::model_base_crtp<model_computeParamsMultivariateN> {
private:
        int M;
        std::vector<vector_d> y;
public:
    model_computeParamsMultivariateN(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_computeParamsMultivariateN(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_computeParamsMultivariateN_namespace::model_computeParamsMultivariateN";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 2;
            context__.validate_dims("data initialization", "M", "int", context__.to_vec());
            M = int(0);
            vals_i__ = context__.vals_i("M");
            pos__ = 0;
            M = vals_i__[pos__++];
            check_greater_or_equal(function__, "M", M, 0);
            current_statement_begin__ = 3;
            validate_non_negative_index("y", "2", 2);
            validate_non_negative_index("y", "M", M);
            context__.validate_dims("data initialization", "y", "vector_d", context__.to_vec(M,2));
            y = std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >(M, Eigen::Matrix<double, Eigen::Dynamic, 1>(2));
            vals_r__ = context__.vals_r("y");
            pos__ = 0;
            size_t y_j_1_max__ = 2;
            size_t y_k_0_max__ = M;
            for (size_t j_1__ = 0; j_1__ < y_j_1_max__; ++j_1__) {
                for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                    y[k_0__](j_1__) = vals_r__[pos__++];
                }
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 7;
            validate_non_negative_index("Omega", "2", 2);
            validate_non_negative_index("Omega", "2", 2);
            num_params_r__ += ((2 * (2 - 1)) / 2);
            current_statement_begin__ = 8;
            validate_non_negative_index("mu", "2", 2);
            num_params_r__ += 2;
            current_statement_begin__ = 9;
            validate_non_negative_index("sigma", "2", 2);
            num_params_r__ += 2;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_computeParamsMultivariateN() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 7;
        if (!(context__.contains_r("Omega")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable Omega missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("Omega");
        pos__ = 0U;
        validate_non_negative_index("Omega", "2", 2);
        validate_non_negative_index("Omega", "2", 2);
        context__.validate_dims("parameter initialization", "Omega", "matrix_d", context__.to_vec(2,2));
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Omega(2, 2);
        size_t Omega_j_2_max__ = 2;
        size_t Omega_j_1_max__ = 2;
        for (size_t j_2__ = 0; j_2__ < Omega_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < Omega_j_1_max__; ++j_1__) {
                Omega(j_1__, j_2__) = vals_r__[pos__++];
            }
        }
        try {
            writer__.corr_matrix_unconstrain(Omega);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable Omega: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 8;
        if (!(context__.contains_r("mu")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable mu missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("mu");
        pos__ = 0U;
        validate_non_negative_index("mu", "2", 2);
        context__.validate_dims("parameter initialization", "mu", "vector_d", context__.to_vec(2));
        Eigen::Matrix<double, Eigen::Dynamic, 1> mu(2);
        size_t mu_j_1_max__ = 2;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            mu(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(mu);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable mu: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 9;
        if (!(context__.contains_r("sigma")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable sigma missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("sigma");
        pos__ = 0U;
        validate_non_negative_index("sigma", "2", 2);
        context__.validate_dims("parameter initialization", "sigma", "vector_d", context__.to_vec(2));
        Eigen::Matrix<double, Eigen::Dynamic, 1> sigma(2);
        size_t sigma_j_1_max__ = 2;
        for (size_t j_1__ = 0; j_1__ < sigma_j_1_max__; ++j_1__) {
            sigma(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_lb_unconstrain(0, sigma);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable sigma: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 7;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Omega;
            (void) Omega;  // dummy to suppress unused var warning
            if (jacobian__)
                Omega = in__.corr_matrix_constrain(2, lp__);
            else
                Omega = in__.corr_matrix_constrain(2);
            current_statement_begin__ = 8;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu;
            (void) mu;  // dummy to suppress unused var warning
            if (jacobian__)
                mu = in__.vector_constrain(2, lp__);
            else
                mu = in__.vector_constrain(2);
            current_statement_begin__ = 9;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> sigma;
            (void) sigma;  // dummy to suppress unused var warning
            if (jacobian__)
                sigma = in__.vector_lb_constrain(0, 2, lp__);
            else
                sigma = in__.vector_lb_constrain(0, 2);
            // transformed parameters
            current_statement_begin__ = 12;
            validate_non_negative_index("Sigma", "2", 2);
            validate_non_negative_index("Sigma", "2", 2);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> Sigma(2, 2);
            stan::math::initialize(Sigma, DUMMY_VAR__);
            stan::math::fill(Sigma, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 13;
            stan::math::assign(Sigma, quad_form_diag(Omega, sigma));
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 12;
            size_t Sigma_j_1_max__ = 2;
            size_t Sigma_j_2_max__ = 2;
            for (size_t j_1__ = 0; j_1__ < Sigma_j_1_max__; ++j_1__) {
                for (size_t j_2__ = 0; j_2__ < Sigma_j_2_max__; ++j_2__) {
                    if (stan::math::is_uninitialized(Sigma(j_1__, j_2__))) {
                        std::stringstream msg__;
                        msg__ << "Undefined transformed parameter: Sigma" << "(" << j_1__ << ", " << j_2__ << ")";
                        stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable Sigma: ") + msg__.str()), current_statement_begin__, prog_reader__());
                    }
                }
            }
            stan::math::check_cov_matrix(function__, "Sigma", Sigma);
            // model body
            current_statement_begin__ = 16;
            lp_accum__.add(cauchy_log<propto__>(sigma, 0, 5));
            current_statement_begin__ = 17;
            lp_accum__.add(lkj_corr_log<propto__>(Omega, 1));
            current_statement_begin__ = 18;
            lp_accum__.add(normal_log<propto__>(mu, 0, 5));
            current_statement_begin__ = 19;
            lp_accum__.add(multi_normal_log<propto__>(y, mu, Sigma));
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("Omega");
        names__.push_back("mu");
        names__.push_back("sigma");
        names__.push_back("Sigma");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(2);
        dims__.push_back(2);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(2);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(2);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(2);
        dims__.push_back(2);
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
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_computeParamsMultivariateN_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Omega = in__.corr_matrix_constrain(2);
        size_t Omega_j_2_max__ = 2;
        size_t Omega_j_1_max__ = 2;
        for (size_t j_2__ = 0; j_2__ < Omega_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < Omega_j_1_max__; ++j_1__) {
                vars__.push_back(Omega(j_1__, j_2__));
            }
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> mu = in__.vector_constrain(2);
        size_t mu_j_1_max__ = 2;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            vars__.push_back(mu(j_1__));
        }
        Eigen::Matrix<double, Eigen::Dynamic, 1> sigma = in__.vector_lb_constrain(0, 2);
        size_t sigma_j_1_max__ = 2;
        for (size_t j_1__ = 0; j_1__ < sigma_j_1_max__; ++j_1__) {
            vars__.push_back(sigma(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 12;
            validate_non_negative_index("Sigma", "2", 2);
            validate_non_negative_index("Sigma", "2", 2);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Sigma(2, 2);
            stan::math::initialize(Sigma, DUMMY_VAR__);
            stan::math::fill(Sigma, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 13;
            stan::math::assign(Sigma, quad_form_diag(Omega, sigma));
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 12;
            stan::math::check_cov_matrix(function__, "Sigma", Sigma);
            // write transformed parameters
            if (include_tparams__) {
                size_t Sigma_j_2_max__ = 2;
                size_t Sigma_j_1_max__ = 2;
                for (size_t j_2__ = 0; j_2__ < Sigma_j_2_max__; ++j_2__) {
                    for (size_t j_1__ = 0; j_1__ < Sigma_j_1_max__; ++j_1__) {
                        vars__.push_back(Sigma(j_1__, j_2__));
                    }
                }
            }
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
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
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_computeParamsMultivariateN";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t Omega_j_2_max__ = 2;
        size_t Omega_j_1_max__ = 2;
        for (size_t j_2__ = 0; j_2__ < Omega_j_2_max__; ++j_2__) {
            for (size_t j_1__ = 0; j_1__ < Omega_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Omega" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        size_t mu_j_1_max__ = 2;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t sigma_j_1_max__ = 2;
        for (size_t j_1__ = 0; j_1__ < sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t Sigma_j_2_max__ = 2;
            size_t Sigma_j_1_max__ = 2;
            for (size_t j_2__ = 0; j_2__ < Sigma_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < Sigma_j_1_max__; ++j_1__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "Sigma" << '.' << j_1__ + 1 << '.' << j_2__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t Omega_j_1_max__ = ((2 * (2 - 1)) / 2);
        for (size_t j_1__ = 0; j_1__ < Omega_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "Omega" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t mu_j_1_max__ = 2;
        for (size_t j_1__ = 0; j_1__ < mu_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "mu" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t sigma_j_1_max__ = 2;
        for (size_t j_1__ = 0; j_1__ < sigma_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "sigma" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            size_t Sigma_j_1_max__ = (2 + ((2 * (2 - 1)) / 2));
            for (size_t j_1__ = 0; j_1__ < Sigma_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Sigma" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_computeParamsMultivariateN_namespace::model_computeParamsMultivariateN stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
