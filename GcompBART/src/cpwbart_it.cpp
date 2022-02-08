/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#include "tree.h"
#include "treefuns.h"

typedef std::vector<tree> vtree;

void getpred(int beg, int end, size_t p, size_t m, size_t np, xinfo& xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat);

RcppExport SEXP cpwbart_it(
    SEXP _itrees,		//treedraws list from fbart
    SEXP _ix,			//x matrix to predict at
    SEXP _iit
)
{
  //Rprintf("*****In main of C++ for bart prediction\n");
  
  //--------------------------------------------------
  //process trees
  Rcpp::List trees(_itrees);
  //get threadcount
  int it = Rcpp::as<int>(_iit);

  Rcpp::CharacterVector itrees(Rcpp::wrap(trees["trees"])); 
  std::string itv(itrees[0]);
  std::stringstream ttss(itv);
  
  size_t nd,m,p;
  //int it = 0;
  ttss >> nd >> m >> p;
  //cout << "number of bart draws: " << nd << endl;
  //cout << "number of trees in bart sum: " << m << endl;
  //cout << "number of x columns: " << p << endl;
  //--------------------------------------------------
  //process cutpoints (from trees)
  Rcpp::List  ixi(Rcpp::wrap(trees["cutpoints"]));
  size_t pp = ixi.size();
  if(p!=pp) cout << "WARNING: p from trees and p from x don't agree\n";
  xinfo xi;
  xi.resize(p);
  for(size_t i=0;i<p;i++) {
    Rcpp::NumericVector cutv(ixi[i]);
    xi[i].resize(cutv.size());
    std::copy(cutv.begin(),cutv.end(),xi[i].begin());
  }
  //--------------------------------------------------
  //process x
  Rcpp::NumericMatrix xpred(_ix);
  size_t np = xpred.ncol();
  //cout << "from x,np,p: " << xpred.nrow() << ", " << xpred.ncol() << endl;
  //read in trees
  std::vector<vtree> tmat(nd);
  for(size_t i=0;i<nd;i++) tmat[i].resize(m);
  for(size_t i=0;i<nd;i++) {
    for(size_t j=0;j<m;j++) ttss >> tmat[i][j];
  }
  //--------------------------------------------------
  //get predictions
  
  Rcpp::NumericMatrix yhat(nd,np);
  std::fill(yhat.begin(), yhat.end(), 0.0);
  double *px = &xpred(0,0);
  //--------------------------------------------------
  //get predictions
  //Rcpp::NumericVector yhat(np);
  //  Rcpp::NumericMatrix yhat(nd,np);
  //std::fill(yhat.begin(), yhat.end(), 0.0);
  //double *px = &xpred(0,0);

  getpred(it-1, it-1, p, m, np,  xi,  tmat, px,  yhat);
  //getpred(it, p, m, np,  xi,  tmat, px,  yhat);

  //--------------------------------------------------
  Rcpp::List ret;
  ret["yhat.test"] = yhat;
  return ret;
}

void getpred(int beg, int end, size_t p, size_t m, size_t np, xinfo& xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat)
{
  double *fptemp = new double[np];
  
  for(int i=beg;i<=end;i++) {
    for(size_t j=0;j<m;j++) {
      fit(tmat[i][j],xi,p,np,px,fptemp);
      for(size_t k=0;k<np;k++) yhat(i,k) += fptemp[k];
    }
  }
  
  delete [] fptemp;
}
