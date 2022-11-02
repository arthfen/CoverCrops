// Supplementary Material - From regional to parcel scale: a high-resolution map of cover crops across Europe combining satellite data with statistical surveys
// This file contains the C++ auxiliary codes for parameter estimation
// Author: Arthur Nicolaus Fendrich

#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::rowvec Arma_colSums(const arma::mat& x) {
  return arma::sum(x, 0);
}

// [[Rcpp::export]]
double fitCpp(arma::vec rho, const arma::colvec& Qty, const arma::mat& R, const arma::mat& RtR, const double& r, const double& N, const arma::vec& b0, const arma::mat& S1, const arma::mat& S2, const arma::mat& S3, const arma::mat& S4){
	rho = arma::exp(rho);
	if(min(rho) <= 1e-8) return arma::datum::inf;
	if(max(rho) == arma::datum::inf) return arma::datum::inf;

	arma::mat Sl = rho(1) * S1 + rho(2) * S2 + rho(3) * S3 + rho(4) * S4;
	double ldetSl;
	try {
		ldetSl = arma::log_det_sympd(Sl.rows(1, Sl.n_rows - 1).cols(1, Sl.n_cols - 1));
	} catch(...) {
		return arma::datum::inf;
	}

	arma::mat V = RtR/rho(0) + Sl;
	double ldetV;
	try {
		ldetV = arma::log_det_sympd(V);
	} catch(...) {
		return arma::datum::inf;
	}

	arma::mat g0 = R.t()*Qty/rho(0) - Sl*b0; // gradient
	arma::mat H0 = -V; // approximate Hessian
	arma::mat hat_beta;
	try {
		hat_beta = b0 - solve(H0, g0, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
	} catch(...) {
                return arma::datum::inf;
	}
	// this is test to approximate ||Y - A*g(Xb)||^2 + b'Sb

	arma::mat logLik = - Qty.t()*Qty/rho(0) - r/rho(0);
	logLik += - b0.t() * Sl * b0;
	logLik += - 2*b0.t()*Sl*(hat_beta - b0) + 2*Qty.t()*R/rho(0)*(hat_beta - b0);
	logLik += - (hat_beta - b0).t() * V * (hat_beta - b0);
	logLik += - N * std::log(rho(0)) + ldetSl - ldetV;
	
	return -1.0 * arma::as_scalar(logLik);
}

// [[Rcpp::export]]
arma::vec hbCpp(arma::vec rho, const arma::colvec& Qty, const arma::mat& R, const arma::mat& RtR, const double& r, const double& N, const arma::vec& b0, const arma::mat& S1, const arma::mat& S2, const arma::mat& S3, const arma::mat& S4, const double& w){
	rho = arma::exp(rho);

	arma::mat Sl = rho(1) * S1 + rho(2) * S2 + rho(3) * S3 + rho(4) * S4;
	arma::mat V = RtR/rho(0) + Sl;
	arma::mat g0 = R.t()*Qty/rho(0) - Sl*b0; // gradient
	arma::mat H0 = -V; // approximate Hessian
	arma::mat hat_beta = b0 - w * solve(H0, g0, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);

	return hat_beta;
}

// [[Rcpp::export]]
arma::vec getIndex(const arma::vec& v, const int b){
	arma::vec res(v.size());

	int j = 0;
	for (int i = 0; i != v.size(); i++) {
		if (v[i] == b){
			res[j] = i;
			j++;
		}
	}
	res = res.head(j);

	return res;
}

// [[Rcpp::export]]
arma::uvec getRowIndex_NA(const arma::mat& m){
	int n = m.n_rows;
	arma::uvec res = find_nonfinite(m);
	res -= floor(res/n)*n;

	return unique(res) + 1;
}


arma::vec g(arma::vec& x){
	arma::vec v = 1/(1 + arma::exp(x));
	return v;
}

arma::mat gmat(arma::mat& x){
	arma::mat v = 1/(1 + arma::exp(x));
	return v;
}

arma::vec dg(arma::vec& x){
	arma::vec v = - arma::exp(x)/arma::pow((arma::exp(x) + 1), 2);
	v.elem(find_nonfinite(v)).zeros();
	return v;
}

// [[Rcpp::export]]
Rcpp::List tpmm_agg_noblock(const arma::mat& Xs, const arma::mat& Xt, const arma::vec& groups, const arma::vec& blocks, const arma::vec& b, const arma::mat& qrQ){
	int ms = Xs.n_cols;
	int mt = Xt.n_cols;
	int m = 1 + ms * mt; // one for the intercept

	int n = groups.size();

	bool cstr;
	if(qrQ.n_rows == 1){
		cstr = false;
	} else {
		cstr = true;
		Rcout << "ERROR: tpmm_agg_noblock \n";
	}

	// declare the output vectors
	arma::vec gzd(n);
	arma::mat AGX;
	if(cstr){
		AGX.zeros(n, m-1);
	} else {
		AGX.zeros(n, m);
	}
	arma::vec MMt(n);

	// declare the temporary vectors
	arma::uvec idx;
	arma::mat Xs_tmp;
	arma::mat Xt_tmp;
	arma::vec zd_tmp;
	arma::vec gzd_tmp;
	arma::vec dgzd_tmp;
	int st;

	for(int i = 0; i != n; i++) {
		idx = arma::conv_to<arma::uvec>::from(getIndex(blocks, groups[i]));
		Xs_tmp = Xs.rows(idx);
		Xt_tmp = Xt.rows(idx);

		for(int j = 0; j != ms; j++){
			st = 1 + j * mt;

			if(j == 0){
				zd_tmp = (Xt_tmp * b.rows(st, st + mt - 1)) % Xs_tmp.col(j);
			} else {
				zd_tmp += (Xt_tmp * b.rows(st, st + mt - 1)) % Xs_tmp.col(j);
			}
		}
		zd_tmp += b[0];
		gzd_tmp = g(zd_tmp);
		dgzd_tmp = dg(zd_tmp);

		for(int j = 0; j != ms; j++){
			st = 1 + j * mt;

			AGX.row(i).cols(st, st + mt - 1) += (Xt_tmp.t() * (dgzd_tmp % Xs_tmp.col(j))).t();
		}
		gzd.row(i) = arma::sum(gzd_tmp);
		MMt.row(i) = idx.size();

		AGX.row(i).col(0) = arma::sum(dgzd_tmp);
	}

	return Rcpp::List::create(Rcpp::Named("gz.") = gzd, Rcpp::Named("AGX") = AGX, Rcpp::Named("MMt") = MMt);
}

// [[Rcpp::export]]
arma::vec tpmm_vec_q_noblock(const arma::mat& Xs, const arma::mat& Xt, const arma::mat& b, const arma::mat& qrQ, int type){
	int ms = Xs.n_cols;
	int mt = Xt.n_cols;
	int m = 1 + ms * mt; // one for the intercept

	int n = Xs.n_rows;
	int nb = b.n_cols;

	bool cstr;
	if(qrQ.n_rows == 1){
		cstr = false;
	} else {
		cstr = true;
		Rcout << "ERROR: tpmm_agg_noblock \n";
	}

	// declare the output vectors
	arma::mat gzd(n, nb, arma::fill::zeros);

	// declare the temporary vectors
	arma::mat X_tmp;
	int st;
	int en;

	X_tmp.ones(n, mt);
	for(int j = 0; j != ms; j++){
		X_tmp = Xs.col(j) % Xt.each_col();
		st = 1 + j * mt;
		gzd += X_tmp * b.rows(st, st + mt - 1);
	}
	gzd.each_row() += b.row(0); // intercept
	gzd = gmat(gzd);

	arma::vec out;
	if(type == 0){
		arma::vec P_005 = {0.50};
		out = arma::quantile(gzd, P_050, 1);
	} if(type == 1){
		out = arma::stddev(gzd, 0, 1);
	} else if(type == 2){
		arma::vec P_005 = {0.05};
		out = arma::quantile(gzd, P_005, 1);
	} else if(type == 3){
		arma::vec P_095 = {0.95};
		out = arma::quantile(gzd, P_095, 1);
	}	

	return out;
}


// [[Rcpp::export]]
Rcpp::List qrCpp(const arma::mat& m){
	arma::mat Q(m.n_rows, m.n_rows);
	arma::mat R(m.n_rows, m.n_cols);
	arma::qr(Q, R, m);

	return Rcpp::List::create(Rcpp::Named("Q") = Q, Rcpp::Named("R") = R);
}

