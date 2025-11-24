// This file contains code adapted from the SoftBart package by Antonio Linero.
// Source: https://CRAN.R-project.org/package=SoftBart

#ifndef GCSOFT_BART_H
#define GCSOFT_BART_H

#include <RcppArmadillo.h>
#include "GCfunctions.h"

struct Hypers;
struct Node;


// Extension for GcompBart:
// Initializes Opts object from an Rcpp::List.
// Includes new update flags for longitudinal and hierarchical components:
// - update_tvp, update_eta, update_phi, update_alpha_vec
// Also sets MCMC control parameters (num_burn, num_thin, num_save, num_print)
// and caching options.
struct Hypers {

  double alpha;
  double eta;
  arma::vec phi;
  arma::vec alpha_vec;
  double beta;
  double gamma;
  double sigma;
  double sigma_mu;
  double shape;
  double width;
  double tau_rate;
  double num_tree_prob;
  double temperature;
  int num_tree;
  int num_groups;
  int num_tgroups;
  arma::vec s;
  arma::vec logs;
  arma::uvec idx_t;
  arma::uvec idx_tv;
  double w;
  arma::vec logw;
  arma::vec v;
  arma::vec logv;
  arma::vec u;
  arma::vec logu;
  arma::vec pt;
  arma::vec logpt;
  arma::uvec group;
  arma::uvec tgroup;
  arma::uvec tgroup_size;

  arma::vec rho_propose;

  std::vector<std::vector<unsigned int> > group_to_vars;

  double sigma_hat;
  double sigma_mu_hat;
  double alpha_scale;
  double alpha_shape_1;
  double alpha_shape_2;
  arma::vec alpha_shape_1_vec;

  void UpdateSigma(const arma::vec& r, const arma::vec& weights);
  void UpdateSigmaMu(const arma::vec& means);
  void UpdateAlpha();
  void UpdateEta();
  void UpdatePhi();
  void UpdateAlpha_vec();
  void UpdateGamma(std::vector<Node*>& forest);
  void UpdateBeta(std::vector<Node*>& forest);
  void UpdateTauRate(const std::vector<Node*>& forest);

  // For updating tau
  // double loglik_tau(double tau,
  //                   const std::vector<Node*>& forest,
  //                   const arma::mat& X, const arma::vec& Y);
  // void update_tau(std::vector<Node*>& forest,
  //                 const arma::mat& X, const arma::vec& Y);

  int SampleVar() const;

  Hypers(Rcpp::List hypers);
  Hypers();

};

struct Node {

  bool is_leaf;
  bool is_root;
  Node* left;
  Node* right;
  Node* parent;

  // Branch parameters
  int var;
  double val;
  double lower;
  double upper;
  double tau;

  // Leaf parameters
  double mu;

  // Data for computing weights
  double current_weight;

  // Functions
  void Root(const Hypers& hypers);
  void GetLimits();
  void AddLeaves();
  void BirthLeaves(const Hypers& hypers);
  bool is_left();
  void GenTree(const Hypers& hypers);
  void GenBelow(const Hypers& hypers);
  void GetW(const arma::mat& X, int i);
  void DeleteLeaves();
  void UpdateMu(const arma::vec& Y, const arma::vec& weights,
                const arma::mat& X, const Hypers& hypers);
  void UpdateTau(const arma::vec& Y, const arma::vec& weights,
                 const arma::mat& X, const Hypers& hypers);
  void SetTau(double tau_new);
  double loglik_tau(double tau_new, const arma::mat& X, const arma::vec& Y,
                    const arma::vec& weights, const Hypers& hypers);

  Node();
  ~Node();

};


// Extension for GcompBart:
// Defines the Opts structure for MCMC and update options.
// Includes:
// - Iteration controls: num_burn, num_thin, num_save, num_print
// - Update flags for model components: sigma_mu, s, alpha, tvp, eta, phi, alpha_vec, beta, gamma, tau, tau_mean, num_tree, sigma
struct Opts {
  int num_burn;
  int num_thin;
  int num_save;
  int num_print;

  bool update_sigma_mu;
  bool update_s;
  bool update_alpha;
  bool update_tvp;
  bool update_eta;
  bool update_phi;
  bool update_alpha_vec;
  bool update_beta;
  bool update_gamma;
  bool update_tau;
  bool update_tau_mean;
  bool update_num_tree;
  bool update_sigma;
  bool cache_trees;

Opts() : update_sigma_mu(true), update_s(true), update_alpha(true),
         update_tvp(false), update_eta(false), update_phi(false), update_alpha_vec(false),
         update_beta(false), update_gamma(false), update_tau(true),
         update_tau_mean(false), update_num_tree(false), update_sigma(true), cache_trees(false) {

  num_burn = 1;
  num_thin = 1;
  num_save = 1;
  num_print = 100;

}

Opts(Rcpp::List opts_) {

  update_sigma_mu = opts_["update_sigma_mu"];
  update_s = opts_["update_s"];
  update_alpha = opts_["update_alpha"];
  update_tvp = opts_["update_tvp"];
  update_eta = opts_["update_eta"];
  update_phi = opts_["update_phi"];
  update_alpha_vec = opts_["update_alpha_vec"];
  update_beta = opts_["update_beta"];
  update_gamma = opts_["update_gamma"];
  update_tau = opts_["update_tau"];
  update_tau_mean = opts_["update_tau_mean"];
  update_num_tree = opts_["update_num_tree"];
  update_sigma = opts_["update_sigma"];
  num_burn = opts_["num_burn"];
  num_thin = opts_["num_thin"];
  num_save = opts_["num_save"];
  num_print = opts_["num_print"];
  cache_trees = opts_["cache_trees"];
}

};

class Forest {

 private:

  std::vector<Node*> trees;
  std::vector<std::vector<Node*>> saved_forests;
  Hypers hypers;
  Opts opts;

  arma::umat tree_counts;

 public:

  /* Forest(Rcpp::List hypers_); */
  Forest(Rcpp::List hypers_, Rcpp::List opts_);
  ~Forest();
  // arma::vec predict(const arma::mat& X);
  arma::mat do_gibbs(const arma::mat& X,
                     const arma::vec& Y,
                     const arma::mat& X_test, int num_iter);
  arma::vec get_s() {return hypers.s;}
  arma::mat do_gibbs_weighted(const arma::mat& X,
                              const arma::vec& Y, const arma::vec& weights,
                              const arma::mat& X_test, int num_iter);
  arma::uvec get_counts();
  arma::umat get_tree_counts();
  void set_sigma(double sigma);
  int num_gibbs;
  arma::vec predict_iteration(const arma::mat& X, std::size_t r_ter);
  double get_sigma();
  void set_s(const arma::vec& s_);
  arma::vec do_predict(const arma::mat& X);
  double get_sigma_mu();

};



// Extension for GcompBart:
// Initializes Opts object with all MCMC and update options.
// Includes new flags for longitudinal and hierarchical components:
// update_tvp, update_eta, update_phi, update_alpha_vec.
Opts InitOpts(int num_burn, int num_thin, int num_save, int num_print,
              bool update_sigma_mu, bool update_s, bool update_alpha,
              bool update_tvp, bool update_eta, bool update_phi, bool update_alpha_vec,
              bool update_beta, bool update_gamma, bool update_tau,
              bool update_tau_mean, bool update_num_tree, bool update_sigma);


// Extension for GcompBart:
// Initializes Hypers object with extended hyperparameters for longitudinal modeling.
// Adds support for eta, phi, alpha_vec, and time-varying group structures.
// Sets prior distributions and hierarchical parameters for Bayesian tree ensemble.
Hypers InitHypers(const arma::mat& X, double sigma_hat, double alpha, double eta,
                  const arma::vec& phi, const arma::vec& alpha_vec,
                  double beta, double gamma, double k, double width, double shape,
                  int num_tree, double alpha_scale, double alpha_shape_1,
                  double alpha_shape_2, double tau_rate, double num_tree_prob,
                  double temperature);

void GetSuffStats(Node* n, const arma::vec& y, const arma::vec& weights,
                  const arma::mat& X, const Hypers& hypers,
                  arma::vec& mu_hat_out, arma::mat& Omega_inv_out);

double LogLT(Node* n, const arma::vec& Y, const arma::vec& weights,
             const arma::mat& X, const Hypers& hypers);

double cauchy_jacobian(double tau, double sigma_hat);

double update_sigma(const arma::vec& r,
                    double sigma_hat, double sigma_old,
                    double temperature = 1.0);
double update_sigma(const arma::vec& r, const arma::vec& weights,
                    double sigma_hat, double sigma_old,
                    double temperature = 1.0);
arma::vec loglik_data(const arma::vec& Y, const arma::vec& weights, const arma::vec& Y_hat, const Hypers& hypers);
arma::vec predict(const std::vector<Node*>& forest,
                  const arma::mat& X,
                  const Hypers& hypers);

arma::vec predict(Node* node,
                  const arma::mat& X,
                  const Hypers& hypers);

bool is_left(Node* n);

double SplitProb(Node* node, const Hypers& hypers);
int depth(Node* node);
void leaves(Node* x, std::vector<Node*>& leafs);
std::vector<Node*> leaves(Node* x);
arma::vec get_means(std::vector<Node*>& forest);
void get_means(Node* node, std::vector<double>& means);
std::vector<Node*> init_forest(const arma::mat& X, const arma::vec& Y,
                               const Hypers& hypers);

Rcpp::List do_soft_bart(const arma::mat& X,
                        const arma::vec& Y,
                        const arma::vec& weights,
                        const arma::mat& X_test,
                        const Hypers& hypers,
                        const Opts& opts);

void IterateGibbsWithS(std::vector<Node*>& forest, arma::vec& Y_hat, const arma::vec& weights,
                       Hypers& hypers, const arma::mat& X, const arma::vec& Y,
                       const Opts& opts);
void IterateGibbsNoS(std::vector<Node*>& forest, arma::vec& Y_hat,
                     const arma::vec& weights,
                     Hypers& hypers, const arma::mat& X, const arma::vec& Y,
                     const Opts& opts);
void TreeBackfit(std::vector<Node*>& forest, arma::vec& Y_hat,
                 const arma::vec& weights,
                 Hypers& hypers, const arma::mat& X, const arma::vec& Y,
                 const Opts& opts);
double activation(double x, double c, double tau);
void birth_death(Node* tree, const arma::mat& X, const arma::vec& Y,
                 const arma::vec& weights, Hypers& hypers);
void node_birth(Node* tree, const arma::mat& X, const arma::vec& Y,
                const arma::vec& weights, Hypers& hypers);
void node_death(Node* tree, const arma::mat& X, const arma::vec& Y,
                const arma::vec& weights, Hypers& hypers);
Node* draw_prior(Node* tree, const arma::mat& X, const arma::vec& Y,
                 const arma::vec& weights, Hypers& hypers);
double growth_prior(int leaf_depth, const Hypers& hypers);
Node* birth_node(Node* tree, double* leaf_node_probability);
double probability_node_birth(Node* tree);
Node* death_node(Node* tree, double* p_not_grand);
std::vector<Node*> not_grand_branches(Node* tree);
void not_grand_branches(std::vector<Node*>& ngb, Node* node);
arma::uvec get_var_counts(std::vector<Node*>& forest, const Hypers& hypers);
void get_var_counts(arma::uvec& counts, Node* node, const Hypers& hypers);
arma::vec rdirichlet(const arma::vec& shape);
double alpha_to_rho(double alpha, double scale);
double rlgam(double shape);
double rho_to_alpha(double rho, double scale);
double logpdf_beta(double x, double a, double b);
double growth_prior(int node_depth, double gamma, double beta);
double forest_loglik(std::vector<Node*>& forest, double gamma, double beta);
double tree_loglik(Node* node, int node_depth, double gamma, double beta);
Node* rand(std::vector<Node*> ngb);
void UpdateS(std::vector<Node*>& forest, Hypers& hypers);
void UpdateTVP(std::vector<Node*>& forest, Hypers& hypers);
void copy_node(Node* nn, Node* n);
Node* copy_tree(Node* root, Hypers& hypers);
std::vector<Node*> copy_forest(std::vector<Node*> forest, Hypers& hypers);

// For tau
bool do_mh(double loglik_new, double loglik_old,
           double new_to_old, double old_to_new);
double logprior_tau(double tau, double tau_rate);
double tau_proposal(double tau);
double log_tau_trans(double tau_new);
arma::vec get_tau_vec(const std::vector<Node*>& forest);

// RJMCMC for trees
std::vector<Node*> TreeSwap(std::vector<Node*>& forest);
std::vector<Node*> TreeSwapLast(std::vector<Node*>& forest);
std::vector<Node*> AddTree(std::vector<Node*>& forest,
                           const Hypers& hypers, const Opts& opts);
std::vector<Node*> DeleteTree(std::vector<Node*>& forest);
void update_num_tree(std::vector<Node*>& forest, Hypers& hypers,
                     const Opts& opts,
                     const arma::vec& Y, const arma::vec& res,
                     const arma::mat& X);
double LogLF(const std::vector<Node*>& forest, const Hypers& hypers,
             const arma::vec& Y, const arma::mat& X);
double loglik_normal(const arma::vec& resid, const double& sigma);
void BirthTree(std::vector<Node*>& forest,
               Hypers& hypers,
               const Opts& opts,
               const arma::vec& Y,
               const arma::vec& res,
               const arma::mat& X);
void DeathTree(std::vector<Node*>& forest,
               Hypers& hypers,
               const arma::vec& Y,
               const arma::vec& res,
               const arma::mat& X);
double TPrior(const std::vector<Node*>& forest, const Hypers& hypers);
void RenormAddTree(std::vector<Node*>& forest,
                   std::vector<Node*>& new_forest,
                   Hypers& hypers);
void UnnormAddTree(std::vector<Node*>& forest,
                   std::vector<Node*>& new_forest,
                   Hypers& hypers);
void RenormDeleteTree(std::vector<Node*>& forest,
                      std::vector<Node*>& new_forest,
                      Hypers& hypers);
void UnnormDeleteTree(std::vector<Node*>& forest,
                      std::vector<Node*>& new_forest,
                      Hypers& hypers);

// Slice sampler

struct rho_loglik {
  double mean_log_s;
  double p;
  double alpha_scale;
  double alpha_shape_1;
  double alpha_shape_2;

  double operator() (double rho) {

    double alpha = rho_to_alpha(rho, alpha_scale);

    double loglik = alpha * mean_log_s
      + Rf_lgammafn(alpha)
      - p * Rf_lgammafn(alpha / p)
      + logpdf_beta(rho, alpha_shape_1, alpha_shape_2);

    return loglik;

  }
};

struct rho_loglik_tv {
  double logu_t;
  double p;
  double sum_lng_a_by_p_not_t;
  double sum_a_by_p_not_t;
  double sum_a_by_p_logs_not_t;
  double alpha_scale;
  double alpha_shape_1;
  double alpha_shape_2;

  double operator()(double rho) const {
    double alpha = rho_to_alpha(rho, alpha_scale);
    double alpha_over_p = alpha / p;

    // Update sums with current alpha contribution
    double sum_alpha_by_p = sum_a_by_p_not_t + alpha_over_p;
    double sum_alpha_by_p_logs = sum_a_by_p_logs_not_t + alpha_over_p * logu_t;
    double sum_lng_alpha_by_t = sum_lng_a_by_p_not_t + Rf_lgammafn(alpha_over_p);

    return sum_alpha_by_p_logs
    + Rf_lgammafn(sum_alpha_by_p)
      - sum_lng_alpha_by_t
      + logpdf_beta(rho, alpha_shape_1, alpha_shape_2);
  }
};

// struct rho_loglik_tv {
//   double logu_t;
//   double p;
//   double sum_lng_a_by_p_not_t;
//   double sum_a_by_p_not_t;
//   double sum_a_by_p_logs_not_t;
//   double alpha_scale;
//   double alpha_shape_1;
//   double alpha_shape_2;
//
//   double operator() (double rho) {
//
//     double alpha = rho_to_alpha(rho, alpha_scale);
//     double sum_lng_alpha_by_t = sum_lng_a_by_p_not_t + Rf_lgammafn(alpha / p);
//     double sum_alpha_by_p = sum_a_by_p_not_t + alpha / p;
//     double sum_alpha_by_p_logs = sum_a_by_p_logs_not_t + (alpha / p * logu_t);
//
//     double loglik = sum_alpha_by_p_logs
//     + Rf_lgammafn(sum_alpha_by_p)
//       - sum_lng_alpha_by_t
//       + logpdf_beta(rho, alpha_shape_1, alpha_shape_2);
//
//       return loglik;
//
//   }
// };
//
// Combines rho_loglik and rho_loglik_tv
// struct RhoLogLik {
//   double p;
//   double alpha_scale;
//   double alpha_shape_1;
//   double alpha_shape_2;
//
//   // Optional time-varying components
//   bool time_varying;
//   double logu_t;
//   double sum_lng_a_by_p_not_t;
//   double sum_a_by_p_not_t;
//   double sum_a_by_p_logs_not_t;
//   double mean_log_s;
//
//   double operator()(double rho) const {
//     double alpha = rho_to_alpha(rho, alpha_scale);
//
//     if (!time_varying) {
//       // Original formula
//       return alpha * mean_log_s
//       + Rf_lgammafn(alpha)
//       - p * Rf_lgammafn(alpha / p)
//       + logpdf_beta(rho, alpha_shape_1, alpha_shape_2);
//     } else {
//       // Extended formula
//       double sum_lng_alpha_by_t = sum_lng_a_by_p_not_t + Rf_lgammafn(alpha / p);
//       double sum_alpha_by_p = sum_a_by_p_not_t + alpha / p;
//       double sum_alpha_by_p_logs = sum_a_by_p_logs_not_t + (alpha / p * logu_t);
//
//       return sum_alpha_by_p_logs
//       + Rf_lgammafn(sum_alpha_by_p)
//         - sum_lng_alpha_by_t
//         + logpdf_beta(rho, alpha_shape_1, alpha_shape_2);
//     }
//   }
// };

double slice_sampler(double x0, rho_loglik& g, double w,
                     double lower, double upper) {


  /* Find the log density at the initial point, if not already known. */
  double gx0 = g(x0);

  /* Determine the slice level, in log terms. */

  double logy = gx0 - exp_rand();

  /* Find the initial interval to sample from */

  double u = w * unif_rand();
  double L = x0 - u;
  double R = x0 + (w-u);

  /* Expand the interval until its ends are outside the slice, or until the
     limit on steps is reached */

  do {

    if(L <= lower) break;
    if(g(L) <= logy) break;
    L -= w;

  } while(true);

  do {
    if(R >= upper) break;
    if(g(R) <= logy) break;
    R += w;
  } while(true);

  // Shrink interval to lower and upper bounds

  if(L < lower) L = lower;
  if(R > upper) R = upper;

  // Sample from the interval, shrinking it on each rejection

  double x1 = 0.0;

  do {

    x1 = (R - L) * unif_rand() + L;
    double gx1 = g(x1);

    if(gx1 >= logy) break;

    if(x1 > x0) {
      R = x1;
    }
    else {
      L = x1;
    }

  } while(true);

  return x1;

}

// PERTURB STUFF
void branches(Node* n, std::vector<Node*>& branch_vec);
std::vector<Node*> branches(Node* root);
double calc_cutpoint_likelihood(Node* node);
std::vector<double> get_perturb_limits(Node* branch);

void perturb_decision_rule(Node* tree,
                           const arma::mat& X,
                           const arma::vec& Y,
                           const arma::vec& weights,
                           Hypers& hypers);

#endif
