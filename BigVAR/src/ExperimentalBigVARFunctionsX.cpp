#include <numeric>      // std::iota
#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <algorithm>
#define NDEBUG 1

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// omp_set_num_threads(4);

// Soft thresholding
// [[Rcpp::export]]
double ST1a(double z, double gam) {
	if (z > 0 && gam < fabs(z)) return(z - gam);
	if (z < 0 && gam < fabs(z)) return(z + gam);
	if (gam >= fabs(z)) return(0);
	else return(0);
}

// Columnwise softthresholding

// [[Rcpp::export]]
colvec ST3a(colvec z, double gam) {
	int n = z.size();
	colvec z1(n);
	for (int i = 0; i < n; ++i) {
        double z11 = z(i);
        z1(i) = ST1a(z11, gam);
    }
	return(z1);
}

// Columnwise softthresholding pass by reference for Lasso-VAR
colvec ST3ar(const colvec& z, const double gam) {
	int n = z.size();
	colvec z1(n);
	for (int i = 0; i < n; ++i) {
        double z11 = z(i);
        z1(i) = ST1a(z11, gam);
    }
	return(z1);
}


uvec ind(int n2, int m) {
	std::vector<int> subs;
	for (int i = 0; i < n2; ++i) {
        subs.push_back(i);
    }
	subs.erase(subs.begin() + m);
	return(conv_to<uvec>::from(subs));
}

// Lasso Fista Function
mat FistaLV(const mat& Y, const mat& Z, mat& B, const rowvec gam, const double eps, double tk, int k, int p, bool sep_lambda = false){
	B = trans(B);
	colvec B1 = B.col(0);
	double j = 1;
	for (int i = 0; i < k; ++i) {
        B1 = B.col(i);
        colvec BOLD = B.col(i);
        colvec BOLDOLD = BOLD;
        double thresh = 10 * eps;
        j = 1;
        double tempgam;
        double maxiters = 1000;
        if (sep_lambda) {
            tempgam = gam(i);
        } else {
            tempgam = gam(0);
        }
        while ((thresh > eps) & (j < maxiters)) {
            colvec v = BOLD + ((j - 2) / (j + 1)) * (BOLD - BOLDOLD);
            B1 = ST3a(vectorise(v) + tk * vectorise((trans(Y.col(i)) - trans(v) * Z) * trans(Z)), tempgam * tk);
            thresh = max(abs(B1 - v));
            BOLDOLD = BOLD;
            BOLD = B1;
            j += 1;
        }
        B.col(i) = B1;
	}
	B = trans(B);
	return(B);
} 



// Lasso Fista Function
mat FistaLVEN(const mat& Y, const mat& Z, mat& B, const rowvec gam,double alpha, const double eps, double tk, int k,int p, bool sep_lambda = false) {
	B = trans(B);
	colvec B1 = B.col(0);
	double j = 1;
	mat I(Z.n_cols, Z.n_cols);
	I.eye();
	for (int i = 0; i < k; ++i) {
        B1 = B.col(i);
        colvec BOLD = B.col(i);
        colvec BOLDOLD = BOLD;
        double thresh = 10 * eps;
        j = 1;
        double tempgam;
        double maxiters = 1000;
        if (sep_lambda) {
            tempgam = gam(i);
        } else {
            tempgam = gam(0);
        }
        while ((thresh > eps) & (j < maxiters)) {
            colvec v = BOLD + ((j - 2) / (j + 1)) * (BOLD - BOLDOLD);
            //From Friedman et al
            B1 = ST3a(vectorise(v) + tk * vectorise((trans(Y.col(i)) - trans(v) * Z) * trans(Z)), tempgam * tk * alpha) / (1 + tempgam * tk * (1 - alpha));
            thresh = max(abs(B1 - v));
            BOLDOLD = BOLD;
            BOLD = B1;
            j += 1;
        }
        B.col(i) = B1;
    }
	B = trans(B);
	return(B);
} 


//Penalty Loop For FISTA

// [[Rcpp::export]]
cube gamloopFista(NumericVector beta_, const mat& Y, const mat& Z, const mat gammgrid, const double eps,const colvec& YMean2, const colvec& ZMean2, mat& B1, int k, int p,double tk, int k1,int s, bool sep_lambda=false){
	mat b2 = B1;
	mat B1F2 = B1;
	IntegerVector dims = beta_.attr("dim");
	cube bcube(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube bcube2(dims[0], dims[1] + 1, dims[2]);
	bcube2.fill(0);
	int nseries = dims[0];
	colvec nu = zeros<colvec>(dims[0]);
	//loop through candidate lambda values
	for (int i = 0; i < dims[2]; ++i) {
		rowvec gam = gammgrid.row(i);
		mat B1F2 = bcube.slice(i);
		B1 = FistaLV(Y, Z, B1F2, gam, eps, tk, nseries, p, sep_lambda);
		nu = YMean2 - B1 * ZMean2;
		bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}
    return(bcube2);
}


// [[Rcpp::export]]
cube gamloopFistaEN(NumericVector beta_, const mat& Y, const mat& Z, const mat gammgrid, vec& alpha, const double eps, const colvec& YMean2, const colvec& ZMean2, mat& B1, int k, int p, double tk, int k1, int s, bool sep_lambda = false) {
	mat b2 = B1;
	mat B1F2 = B1;
	IntegerVector dims = beta_.attr("dim");
	cube bcube(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube bcube2(dims[0], dims[1] + 1, dims[2]);
	bcube2.fill(0);
	colvec nu = zeros<colvec>(dims[0]);
	int nlambda = gammgrid.n_rows;
	int nalpha = alpha.n_elem;
	int nseries = dims[0];
	double a_temp;
	for (int i = 0; i < nlambda; ++i) {
		for (int j = 0; j < nalpha; ++j) {
			rowvec gam;
			if (sep_lambda) {
				gam = gammgrid.row(i);
				// can't have separate lambda and vector of alphas  
				a_temp = alpha(0);
			} else {
				rowvec gamm2(1);					   
				gamm2(0) = gammgrid(i, j);
				gam = gamm2;
				a_temp = alpha(j);
			}
			mat B1F2 = bcube.slice((i) * nalpha + j);
			B1 = FistaLVEN(Y, Z, B1F2, gam, a_temp, eps, tk, nseries, p, sep_lambda);
			nu = YMean2 - B1 * ZMean2;
			bcube2.slice((i) * nalpha + j) = mat(join_horiz(nu, B1));
		}
	}
    return(bcube2);
}


// Lag Group VARX-L functions

// Newton-Raphson Functions
inline double trust32(int k, const arma::mat& PP, double delta, double lambda, const arma::vec& EigVA, const arma::mat& EigVector) {
	double g = 0;
	for (int i = 0; i < k; ++i) {
        g += pow(arma::as_scalar(trans(EigVector.col(i)) * PP), 2) / pow(EigVA[i] * delta + lambda, 2);
    }
	return(g);
}

inline double fprime2(int k, const arma::mat& PP, double delta, double lambda, const arma::vec& EigVA, const arma::mat& EigVE) {
	double gg2 = 0;
	for (int i = 0; i < k; ++i) {
        gg2 += (pow(arma::as_scalar(trans(EigVE.col(i)) * PP), 2) * EigVA[i]) / pow(EigVA[i] * delta + lambda, 3);
    }
	double c1 = trust32(k, PP, delta, lambda, EigVA, EigVE);
	double res = -.5 * pow(c1, -1.5) * -2 * gg2;
	return(res);
}

//Newton raphson for trust region problem
inline double Newton2(int k, const arma::mat& P, double lambda, const arma::vec& EigVA, const arma::mat& EigVE) {
	double delta = 0;
	double threshold = 1;
	double phi = 0;
	double deltanew = delta;
	while (threshold > .0001) {
        phi = 1 - 1 / pow(trust32(k, P, delta, lambda, EigVA, EigVE), .5);
        deltanew += phi / fprime2(k, P, deltanew, lambda, EigVA, EigVE);
        threshold = fabs(delta - deltanew);
        delta = deltanew;
        phi = 0;
    }
	return(deltanew);
}

// Eigen-Decomposition
// [[Rcpp::export]]
List Eigencomp( mat& Z1, List groups, int n1, int k1) {
	List M2f(n1);
	List eigvalF(n1);
	List M3f(n1);
	List eigvecF(n1);
	int count = 0;
	for (int i = 0; i < n1; ++i) {
        NumericVector g1 = groups[i];
        count += g1.size();
        arma::uvec s4 = as<arma::uvec>(g1);
        arma::mat M1 = Z1.rows(s4);
        arma::mat M2 = M1 * trans(M1);
        M2f(i) = M2;
        arma::mat D(k1, k1);
        D.eye();
        arma::mat M3 = kron(M2, D);
        M3f(i) = M3;
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, M3);
        eigvalF(i) = eigval;
        eigvecF(i) = eigvec;
    }
    List results = List::create(Named("M3") = wrap(M3f), Named("eigval") = wrap(eigvalF), Named("eigvec") = wrap(eigvecF));
	return(results);
}

// Eigen Decomposition for own/other group VARX-L
// [[Rcpp::export]]
List EigencompOO( mat& ZZ1, List groups, int n1, int k) {
	List M2f(n1);
	List eigvecF(n1);
	List eigvalF(n1);
	int count = 0;
	for (int i = 0; i < n1; ++i) {
        arma::uvec s4 = as<arma::uvec>(groups[i]);
        count += s4.n_elem;
        arma::mat M1 = ZZ1.cols(s4);
        arma::mat M2 = trans(M1) * M1;
        M2f(i) = M2;
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, M2);
        eigvecF(i) = eigvec;
        eigvalF(i) = eigval;
    }
    List results = List::create(Named("M3") = wrap(M2f), Named("eigval") = wrap(eigvalF), Named("eigvec") = wrap(eigvecF));
	return(results);
}



// One group update: Lag VARX-L
List BlockUpdateGL(mat& beta, const mat& Z1, double lam, const mat& Y1, double eps, List groups, List fullgroups, List compgroups, int k, List M3f_, List eigvalF_, List eigvecF_, int k1){
	int n1 = groups.size();
	List active(n1);
	int n = beta.n_rows, m = beta.n_cols;
	arma::mat betaPrev = beta;
	int converge = 0;
	int count = 0;
 
	if (groups.size() == count) {
		beta.zeros(n, m);
		active = groups;
	} else {
		for (int i = 0; i < n1; ++i) {
            arma::uvec s45 = as<arma::uvec>(groups[i]);
            arma::uvec s45F = as<arma::uvec>(fullgroups[i]);

            //inactive groups are set to zero
            if (max(s45) == 0) {
                beta.cols(s45F) = arma::zeros(k1, s45F.n_elem);
                active(i) = 0;
            }
            if (max(s45) != 0) {
                arma::uvec scomp1 = as<arma::uvec>(compgroups[i]);
                mat r = beta.cols(scomp1) * Z1.rows(scomp1) - Y1;

                colvec p = vectorise((r) * trans(Z1.rows(s45)));
                double adjlam = sqrt(static_cast<double>(s45.n_elem)) * lam;

                // threshold to enforce sparsity
                if (arma::norm(p, "fro") <= adjlam) {
                        arma::mat astar = arma::zeros(k1, s45.n_elem);
                        active(i) = 0;
                } else {
                    arma::mat M3 = M3f_(i);
                    int k1a = M3.n_cols;
                    double deltfin = Newton2(k1a, p, adjlam, eigvalF_(i), eigvecF_(i));

                    M3.diag() += adjlam / deltfin;
                    arma::mat astar =- solve(M3, p);
                    astar.set_size(k1, s45.n_elem);
                    beta.cols(s45) = astar;
                    active(i) = s45;
                }
            }
        }
	}
	arma::mat thresh1 = arma::abs((betaPrev - beta) / (arma::ones(n, m) + arma::abs(betaPrev)));
	double thresh = arma::norm(thresh1,"inf");
	//convergence flag
	if (thresh < eps) {
		converge = 1;
	} else {
		converge = 0;
	}

	Rcpp::List results = Rcpp::List::create(Named("beta") = beta, Named("active") = wrap(active), Named("Converge") = wrap(converge));
	return(results);
}


//runs one group until convergence, for use with active set algorithm 
mat ThreshUpdate(mat& betaActive, const mat& Z1, double lam, const mat& Y1, double eps, List groups, List fullgroups, List compgroups, List M2f_, List eigvalF_, List eigvecF_, int k1) {
	int n = betaActive.n_rows, m = betaActive.n_cols;
	int n1 = groups.size();
	mat betaLast = betaActive;
	List active(n1);
	int count = 0;
	List betaActive2(3);
	for (int i = 0; i < n1; ++i) {
        NumericVector g1 = groups[i];
        count += max(g1);
    }
	if (count == 0) {
        betaActive.zeros(n, m);
        active = groups;
    } else {
		double threshold = 10 * eps;
		while (threshold > eps) {
            betaActive2 = BlockUpdateGL(betaActive, Z1, lam, Y1, eps, groups, fullgroups, compgroups, n, M2f_, eigvalF_, eigvecF_, k1);
            betaActive = as<mat>(betaActive2("beta"));
            arma::mat thresh1 = arma::abs((betaLast - betaActive) / (arma::ones(n, m) + arma::abs(betaLast)));
            threshold = arma::norm(thresh1, "inf");
            active = betaActive2("active");
            betaLast = betaActive;
        }
	}
	return(betaActive);
}


// [[Rcpp::export]]
List GamLoopGL2(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y1, const mat& Z1,List jj, List jjfull, List jjcomp, double eps,const colvec& YMean2, const colvec& ZMean2, int k, int pk, const List M2f_, const List eigvalF_, const List eigvecF_) {
	IntegerVector dims = beta_.attr("dim");

	int gran2 = dims[2];
	List activefinal(gran2);
	cube beta2(beta_.begin(), dims[0], dims[1], gran2, false);
	cube betafin(dims[0], dims[1] + 1, dims[2]);
	betafin.fill(0);
	List iterations(gran2);
	mat betaPrev = zeros<mat>(dims[0], dims[1]);
 
	// INDEX LISTS WITH PARENTHESES NOT BRACKETS
	// WHEN EXTRACTING FROM A LIST YOU NEED as<MAT>

	for (int i = 0; i < gran2; ++i) {
        double gam = gamm[i];
        betaPrev = beta2.slice(i);
        List Active = Activeset[i];
        int k2 = 0;
        int converge = 0;
        mat betaF = zeros(k, k);
        List betaFull(3);
        while (converge == 0) {
            betaPrev = ThreshUpdate(betaPrev, Z1, gam, Y1, eps, Active, jjfull, jjcomp, M2f_, eigvalF_, eigvecF_, k);
            betaFull = BlockUpdateGL(betaPrev, Z1, gam, Y1, eps, jjfull, jjfull, jjcomp, k, M2f_, eigvalF_, eigvecF_, k);
            betaF = as<mat>(betaFull("beta"));
            Active = betaFull("active");
            converge = betaFull("Converge");
            k2 += 1;
        }

        colvec nu = YMean2 - betaF * ZMean2;
        betafin.slice(i) = mat(join_horiz(nu, betaF));
        activefinal[i] = Active;
        iterations[i] = k2;
    }
	List Results = List::create(Named("beta") = betafin, Named("active") = wrap(activefinal), Named("iterations") = iterations);
	return(Results);
}


// *
// Group Lasso Own/Other Functions
// *

List BlockUpdate2(const mat& ZZ1, double lam, const mat& Y1, double eps, List groups, List fullgroups, List compgroups, int k, List M2f_, List eigvalF_, List eigvecF_, colvec& B, int k1) {
	int n1 = groups.size();
	List active(n1);
	colvec BPrev = B;
	int converge = 0;
	int count = 0;
 
	if (groups.size() == count) {
		B.zeros();
		active = groups;
	} else {
		for (int i = 0; i < n1; ++i) {
            NumericVector s1 = groups[i];
            IntegerVector s2 = fullgroups[i];
            NumericVector scomp = compgroups[i];
            arma::uvec s45 = as<arma::uvec>(s1);
            arma::uvec s45F = as<arma::uvec>(s2);
            if (max(s1) == 0) {
                B.elem(s45F) = arma::zeros(s2.size());
                active(i) = 0;
            }
            if (max(s1) != 0) {
                arma::uvec s4 = as<arma::uvec>(s1);
                arma::uvec scomp2 = as<arma::uvec>(scomp);
                arma::mat M2a = ZZ1.cols(scomp2);
                arma::colvec a1 = B.elem(scomp2);
                mat foo1 = (M2a * a1);
                arma::colvec r = M2a * a1 - arma::vectorise(Y1);
                arma::mat M1 = ZZ1.cols(s4);
                arma::mat M2 = M2f_(i);
                arma::vec eigval = eigvalF_(i);
                arma::mat eigvec = eigvecF_(i);
                arma::mat p = trans(M1) * r;

                double rho = sqrt(static_cast<double>(s1.size()));
                double adjlam = rho * lam;
                if (arma::norm(p, 2) <= adjlam) {
                    arma::colvec astar = arma::zeros(s1.size());
                    active(i) = 0;
                } else {
                    int k1a = M2.n_cols;
                    double deltfin = Newton2(k1a, p, adjlam, eigval, eigvec);
                    arma::mat D1(s1.size(), s1.size());
                    D1.eye();
                    //correct for rare occurrence where newton returns zero
                    if (deltfin == 0) {
                        deltfin += std::numeric_limits<double>::epsilon();
                    }
                    arma::mat astar = -solve(M2 + adjlam / deltfin * D1, p);
                    B.elem(s4) = astar;
                    active(i) = s4;
                }
            }
        }
	}
	double thresh = arma::norm(B - BPrev, "inf");
	if (thresh < eps) {
		converge = 1;
    } else {
		converge = 0;
	}
	Rcpp::List results = Rcpp::List::create(Named("beta") = B, Named("active") = active, Named("Converge") = converge);
	return(results);
}

colvec ThreshUpdateOO(const mat& ZZ, double lam, const mat& Y, double eps, List groups, List fullgroups, List compgroups, List M2f_, List eigvalF_, List eigvecF_, colvec& B, int n, int k1) {
    int kp = B.n_elem;
	int n1 = groups.size();
	colvec BPrev = B;
	List active(n1);
	List betaActive2(3);
	int count = 0;
	for (int i = 0; i < n1; ++i) {
        NumericVector g1 = groups[i];
        count += max(g1);
    }
	if (count == 0) {
        B.zeros(kp);
        active = groups;
    } else {
		double threshold = 10 * eps;
		while (threshold > eps) {
            betaActive2 = BlockUpdate2(ZZ, lam, Y, eps, groups, fullgroups, compgroups, n, M2f_, eigvalF_, eigvecF_, B, k1);
            B = Rcpp::as<arma::colvec>(betaActive2["beta"]);
            threshold = arma::norm(B - BPrev, "inf");
            active = betaActive2("active");
            BPrev = B;
        }
	}
	return(B);
}

// [[Rcpp::export]]
List GamLoopGLOO(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y, const mat& Z, List jj, List jjfull, List jjcomp, double eps, colvec& YMean2, colvec& ZMean2, int k, int pk, List M2f_, List eigvalF_, List eigvecF_, int k1) {
    IntegerVector dims = beta_.attr("dim");
	int gran2 = gamm.size();
	List activefinal(gran2);
	cube beta2(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube betafin(dims[0], dims[1] + 1, dims[2]);
	betafin.fill(0);
	List iterations(gran2);
	mat betaPrev = zeros<mat>(dims[0], dims[1]);

	arma::colvec B = arma::vectorise(betaPrev);
	NumericVector betaF2(dims[0] * dims[1]);
	for (int i = 0; i < gran2; ++i) {
        double gam = gamm[i];
        betaPrev = beta2.slice(i);
        B = arma::vectorise(betaPrev);
        List Active = Activeset[i];
        int k2 = 0;
        int converge = 0;
        List betaFull(3);
        while (converge == 0) {
            B = ThreshUpdateOO(Z, gam, Y, eps, Active, jjfull, jjcomp, M2f_, eigvalF_, eigvecF_, B, k1, k1);
            betaFull = BlockUpdate2(Z, gam, Y, eps, jjfull, jjfull, jjcomp, k, M2f_, eigvalF_, eigvecF_, B, k1);
            betaF2 = as<NumericVector>(betaFull("beta"));
            Active = betaFull("active");
            converge = betaFull("Converge");
            k2 += 1;
        }
        mat betaF(betaF2.begin(), dims[0], dims[1], false);
        colvec nu = YMean2 - betaF * ZMean2;
        betafin.slice(i) = mat(join_horiz(nu, betaF));
        activefinal[i] = Active;
        iterations[i] = k2;
    }
	List Results = List::create(Named("beta") = betafin, Named("active") = wrap(activefinal), Named("iterations") = iterations);
	return(Results);
}

// *
// Lag VARX-L
// *

mat sparseWLOO(const mat& M1a, const colvec& R1, const double ngroups, colvec& beta, const double t, const double alpha, const double lambda, const double eps, double rho) {
    double thresh = 10;
	colvec p = beta;
	colvec STS = beta;
	colvec thetaOLD = beta;
	double l = 1;
	colvec one = ones<vec>(beta.n_elem);
	while (thresh > eps) {
        p = trans(M1a) * (M1a * beta - R1) / ngroups;
        const colvec p1 = beta - t * vectorise(p);
        if (alpha > 0) {
            STS = ST3ar(p1, t * alpha * lambda);
        } else {
            STS = p1;
        }
        double denom2 = norm(STS, "fro") + sqrt(std::numeric_limits<double>::epsilon());
        double s3 = fmax(1 - (t * (1 - alpha) * lambda * rho) / denom2, 0);
        STS = s3 * STS;
        beta = thetaOLD + (l / (l + 3)) * (STS - thetaOLD);
        l = l + 1;
        //Relative thresholds seem to help with computation time
        thresh = max(abs(beta - STS) / (one + abs(STS)));
        thetaOLD = STS;
    }
	return(beta);
}

List blockUpdateSGLOO(colvec& beta, const mat& Z1, double lam, double alpha, const colvec& Y2, double eps, List groups_, const List fullgroups_, List compgroups_, const List M2f_, const NumericVector Eigs_, double k1, double m) {
	int n1 = groups_.size();
	List active(n1);
	colvec betaPrev = beta;
	int converge = 0;
	int count = 0;
	colvec one = ones<vec>(beta.n_elem);
	if (groups_.size() == count) {
        beta.zeros();
        active = groups_;
	} else {
		for (int i = 0; i < n1; ++i) {
            arma::uvec s4 = as<arma::uvec>(groups_[i]);
            arma::uvec s45F = as<arma::uvec>(fullgroups_[i]);
            uvec scomp2 = as<uvec>(compgroups_[i]);

            if (max(s4) == 0) {
                beta.elem(s45F) = arma::zeros(s45F.n_elem);
                active(i) = 0;
            } else {
                const arma::mat M2a = Z1.cols(scomp2);
                const arma::colvec a1 = beta.elem(scomp2);
                const arma::colvec beta2 = beta.elem(s4);
                const arma::colvec r = Y2 - M2a * a1;
                const mat M1 = Z1.cols(s4);
                const arma::mat p = -trans(M1) * (r - M1 * beta2);
                double rho = sqrt(static_cast<double>(s4.n_elem));
                colvec STS;
                if (alpha > 0) {
                    STS = ST3a(vectorise(p), lam * alpha);
                } else {
                    STS = vectorise(p);
                }

                double lamadj = lam * (1 - alpha) * rho;
                if (arma::norm(STS, "fro") <= lamadj) {
                    arma::colvec astar = arma::zeros(s4.n_elem);
                    active(i) = 0;
                    beta.elem(s4) = astar;
                } else {
                    colvec betaS = beta.elem(s4);
                    const double t = 1 / Eigs_(i);
                    double ngroups = (double) Y2.n_elem;
                    const  mat astar2 = sparseWLOO(M1, r, ngroups, betaS, t, alpha, lam, eps, rho);
                    beta.elem(s4) = astar2;
                    active(i) = s4;
                }
            }
        }
	}
	double thresh = max(abs(beta - betaPrev) / (one + abs(betaPrev)));
	if (thresh < eps) {
		converge = 1;
	} else {
		converge = 0;
	}

	Rcpp::List results = Rcpp::List::create(Named("beta") = wrap(beta), Named("active") = wrap(active), Named("Converge") = wrap(converge), Named("thresh") = wrap(thresh));
	return(results);
}


mat ThreshUpdateSGLOO(colvec& betaActive, const mat& Z, const double lam, const colvec& Y, const double eps, List groups_, const List fullgroups_, const List compgroups_, const List M2f_, const NumericVector eigs_, const double alpha, double k1, double m) {
    int n1 = groups_.size();
	colvec betaLast = betaActive;
	List active(n1);
	int count = 0;
	List betaActive2(3);
	for (int i = 0; i < n1; ++i) {
        NumericVector g1 = groups_[i];
        count += g1.size();
    }
	
	if (count == n1) {
        betaActive.zeros();
        active = groups_;
    } else {
		int converge = 0;
		double th = 10 * eps;
		int iters = 0;
		while ((converge == 0) && (th > eps)) {
            colvec betaOld = betaActive;
            betaActive2 = blockUpdateSGLOO(betaActive, Z, lam, alpha, Y, eps, groups_, fullgroups_, compgroups_, M2f_, eigs_, k1, m);
            betaActive = as<colvec>(betaActive2("beta"));
            converge = betaActive2("Converge");
            iters += 1;
            th = betaActive2("thresh");
            if (iters > 1000) {
                break;
            }
        }
	}
    return(betaActive);
}

//Loop through Lambda values, Sparse Own/Other VARX-L
// [[Rcpp::export]]
List GamLoopSGLOO(NumericVector beta_, const List Activeset_, const NumericVector gamm, const double alpha, const mat& Y, const mat& Z, List jj_, const List jjfull_, List jjcomp_, const double eps, const colvec& YMean2, const colvec& ZMean2, const int k1, const int pk, const List M2f_, const NumericVector eigs_, double m) {
	IntegerVector dims = beta_.attr("dim");
	int gran2 = dims[2];
	List activefinal(gran2);
	cube beta2(beta_.begin(), dims[0], dims[1], gran2, false);
	cube betafin(dims[0], dims[1] + 1, gran2);
	betafin.fill(0);
	List iterations(gran2);
	mat betaPrev = zeros<mat>(dims[0], dims[1]);
	NumericVector betaF2(dims[0] * dims[1]);
	const arma::colvec& Y2 = arma::vectorise(Y, 0);
	int converge;
	for (int i = 0; i < gran2; ++i) {
        double gam = gamm[i];
        betaPrev = beta2.slice(i);
        List Active = Activeset_[i];
        converge = 0;
        colvec B = vectorise(betaPrev);
        List betaFull(3);
        //Three components in the list
        int iters = 0;
        int maxiters = 100;
        while (converge == 0 && iters < maxiters) {
            B = ThreshUpdateSGLOO(B, Z, gam, Y2, eps, Active, jjfull_, jjcomp_, M2f_, eigs_, alpha, (double) k1, m);
            betaFull = blockUpdateSGLOO(B, Z, gam, alpha, Y2, eps, jjfull_, jjfull_, jjcomp_, M2f_, eigs_, (double) k1, m);
            betaF2 = as<NumericVector>(betaFull("beta"));
            Active = betaFull("active");
            converge = betaFull("Converge");
            iters += 1;
        }
        mat betaF(betaF2.begin(), k1, pk, false);
        colvec nu = YMean2 - betaF * ZMean2;
        betafin.slice(i) = mat(join_horiz(nu, betaF));
        activefinal[i] = Active;
        iterations[i] = iters;
    }
	List Results = List::create(Named("beta") = betafin, Named("active") = wrap(activefinal), Named("iterations") = iterations, Named("converge") = converge);
	return(Results);
}



// [[Rcpp::export]]
List GamLoopSGLOODP(NumericVector beta_, const List Activeset_, mat gamm, const colvec alpha, const mat& Y, const mat& Z, List jj_, const List jjfull_, List jjcomp_, const double eps, const colvec& YMean2, const colvec& ZMean2, const int k1, const int pk, const List M2f_, const NumericVector eigs_, double m) {
    int nlambda = gamm.n_rows;
	int nalpha = gamm.n_cols;
	IntegerVector dims = beta_.attr("dim");
	List activefinal(dims[2]);
	cube beta2(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube betafin(dims[0], dims[1] + 1, dims[2]);
	betafin.fill(0);
	List iterations(dims[2]);
	mat betaPrev = zeros<mat>(dims[0], dims[1]);
	NumericVector betaF2(dims[0] * dims[1]);
	const arma:: colvec& Y2 = arma::vectorise(Y, 0);
	int converge;

	for (int i = 0; i < nlambda; ++i) {
		for (int j = 0; j < nalpha; ++j) {
            double gam = gamm(i,j);
            double alpha1 = alpha(j);

            betaPrev = beta2.slice((i) * nalpha + j);
            List Active = Activeset_[(i) * nalpha + j];
            int k2 = 0;
            converge = 0;
            colvec B = vectorise(betaPrev);
            List betaFull(3);
            //Three components in the list
            int iters = 0;
            int maxiters = 1000;
            while ((converge == 0) && (iters < maxiters)) {
                B = ThreshUpdateSGLOO(B, Z, gam, Y2, eps, Active, jjfull_, jjcomp_, M2f_, eigs_, alpha1, (double) k1, m);
                betaFull = blockUpdateSGLOO(B, Z, gam, alpha1, Y2, eps, jjfull_, jjfull_, jjcomp_, M2f_, eigs_, (double) k1, m);
                betaF2 = as<NumericVector>(betaFull("beta"));
                Active = betaFull("active");
                converge = betaFull("Converge");
                iters += 1;
            }
            mat betaF(betaF2.begin(), dims[0], dims[1], false);
            colvec nu = YMean2 - betaF * ZMean2;
            betafin.slice((i) * nalpha + j) = mat(join_horiz(nu, betaF));
            activefinal[(i) * nalpha + j] = Active;
            iterations[(i) * nalpha + j] = k2;
        }
    }
	List Results=List::create(Named("beta") = betafin, Named("active") = wrap(activefinal), Named("iterations") = iterations, Named("converge") = converge);
	return(Results);
}



// *
// Componentwise HLAG
// *

//prox function 
rowvec proxcpp(colvec v2, int L, double lambda, int k, colvec w) {
	colvec r = v2;
	for (int q = (L - 1); q >= 0; --q) {
        std::vector<unsigned int> ivec(L * k - (q) * k);
        std::iota(ivec.begin(), ivec.end(), k * q);
        uvec res = conv_to<uvec>::from(ivec);
        if (norm(r(res) / (lambda * w(q)), "fro") < 1 + 1e-8) {
            r(res) = zeros(res.n_elem);
        } else {
            r(res) = (r(res) - lambda * w(q) * r(res) / norm(r(res), "fro"));
        }
    }
	return(trans(r));
}

//Fista algorithm
// [[Rcpp::export]]
mat Fistapar(const mat Y, const mat Z, const mat phi, const int L, const rowvec lambda, const double eps, const double tk, const int k, bool sep_lambda = false) {
	mat phiFIN = phi;
	colvec w(L);
	for (int r = 0; r < L; ++r) {
		w(r) = sqrt(static_cast<double>(k));
	}
	const colvec w2 = w;
	rowvec phiR = phi.row(0);
	rowvec v = phiR;
	rowvec phiOLD = phiR;
	rowvec phiOLDOLD = phiR;

	for (int i = 0; i < k; ++i) {
        phiR = phi.row(i);
        phiOLD = zeros(1, k * L);
        phiOLDOLD = phiOLD;
        double thresh = 10 * eps;
        double j = 1;
        double templambda;
        if (sep_lambda) {
            templambda = lambda(i);
        } else {
            templambda = lambda(0);
        }
        while (thresh > eps) {
            // Nesterov step
            v = phiOLD + ((j - 2) / (j + 1)) * (phiOLD - phiOLDOLD);
            phiR = proxcpp(vectorise(v) + tk * vectorise((trans(Y.col(i)) - v * Z) * trans(Z)), L, tk * templambda, k, w2);
            thresh = max(abs(phiR - v));

            phiOLDOLD = phiOLD;
            phiOLD = phiR;
            j += 1;
        }
        phiFIN.row(i) = phiR;
    }
	return(phiFIN);
}

// [[Rcpp::export]]
cube gamloopHLAG(NumericVector beta_, const mat& Y, const mat& Z, mat gammgrid, const double eps, const colvec& YMean2, const colvec& ZMean2, mat& B1, const int k, const int p, bool sep_lambda = false) {
	vec eigval;
	mat eigvec;
	const mat& Zt = Z * trans(Z);
	eig_sym(eigval, eigvec, Zt);

	double tk = 1 / max(eigval);
	IntegerVector dims = beta_.attr("dim");

	const int ngridpts = dims[2];
	cube bcube(beta_.begin(), dims[0], dims[1], ngridpts, false);
	cube bcube2(dims[0], dims[1] + 1, ngridpts);
	bcube2.fill(0);
	colvec nu = zeros<colvec>(dims[0]);
	for (int i = 0; i < ngridpts; ++i) {
		rowvec gamm = gammgrid.row(i);
		B1 = bcube.slice(i);
		B1 = Fistapar(Y, Z, B1, p, gamm, eps, tk, k);
		nu = YMean2 - B1 * ZMean2;
		bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}
    return(bcube2);
}

// *
// *
// *
// Own Other HLAG
// *
// *
// *


rowvec proxcppOO(colvec v2, int L, double lambda, List vsubs, int k, colvec w) {
	colvec r = v2;
	for (int i = (L - 1); i >= 0; --i) {
        uvec res = as<uvec>(vsubs(i));
        if (norm(r(res) / (lambda * w(i)), "fro") < 1 + 1e-8) {
            r(res) = zeros(res.n_elem);
        } else {
            r(res) = (r(res) - lambda * w(i) * r(res) / (norm(r(res), "fro")));
        }
    }
	return(trans(r));
}

mat FistaOO(const mat Y, const mat Z, mat phi, const int p, const int k, const rowvec lambda, List groups_, const double eps, const double tk, colvec w, bool sep_lambda) {
	double j = 1;
	mat phiFin = phi;
	rowvec phiR = phi.row(0);
	rowvec phiOLD = phiR;
	rowvec phiOLDOLD = phiOLD;
	rowvec v = phiOLD;
	uvec res1 = ind(p, 0);

	double thresh = 10;
	for (int i = 0; i < k; ++i) {
        j = 1;
        thresh = 10 * eps;
        phiR = phi.row(i);
        phiOLD = phiR;
        phiOLDOLD = phiOLD;
        double templambda;
        if (sep_lambda) {
            templambda = lambda(i);
        } else {
            templambda = lambda(0);
        }
        v = phiR;
        List vsubs = groups_[i];
        while (thresh > eps) {
            v = phiOLD + ((j - 2) / (j + 1)) * (phiOLD - phiOLDOLD);
            phiR = proxcppOO(vectorise(v) + tk * vectorise((trans(Y.col(i)) - v * Z) * trans(Z)), 2 * p, tk * templambda, vsubs, k, w);
            thresh = max(abs(phiR - v));
            phiOLDOLD = phiOLD;
            phiOLD = phiR;
            j += 1;
        }
        phiFin.row(i) = phiR;
    }
	return(phiFin);
}

// [[Rcpp::export]]
cube gamloopOO(NumericVector beta_, const mat Y, const mat Z, mat gammgrid, const double eps, const colvec YMean2, const colvec ZMean2, mat B1, const int k, const int p, colvec w, List groups_, bool sep_lambda = false) {
    mat B1F2 = B1;
	vec eigval;
	mat eigvec;
	const mat Zt = Z * trans(Z);
	eig_sym(eigval, eigvec, Zt);

	double tk = 1 / max(eigval);
    IntegerVector dims = beta_.attr("dim");
    const int ngridpts = dims[2];
	cube bcube(beta_.begin(), dims[0], dims[1], ngridpts, false);
	cube bcube2(dims[0], dims[1] + 1, ngridpts);
	bcube2.fill(0);
	colvec nu = zeros<colvec>(k);
	for (int i = 0; i < ngridpts; ++i) {
        rowvec gamm = gammgrid.row(i);
		B1F2 = bcube.slice(i);
		B1 = FistaOO(Y, Z, B1F2, p, k, gamm, groups_, eps, tk, w, sep_lambda);
		nu = YMean2 - B1 * ZMean2;
		bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}
    return(bcube2);
}

// *
// *
// // Elementwise HLAG
// *
// *

uvec vsubscppelem(int p, int pmax) {
	uvec vs(pmax - p + 1);
	for(int i = pmax; i >= p; --i) {
        vs(i - p) = i - 1;
    }
	return(vs);
}

uvec bbsubs(int j, int k, int p) {
	uvec bb(p);
	bb(0) = j;
	for (int i = 1; i < p; ++i) {
        bb(i) = j + k * (i);
    }
	return(bb);
}


rowvec proxcppelem(colvec v2, int L, double lambda, uvec res1, colvec w) {
	colvec r = v2;
	for (int i = (L - 1); i >= 0; --i) {
        uvec res = vsubscppelem(i + 1, L);
        if (norm(r(res) / (lambda * w(i)), "fro") < 1 + 1e-8) {
            r(res) = zeros(res.n_elem);
        } else {
            r(res) = r(res) - lambda * w(i) * r(res) / (norm(r(res), "fro"));
        }
    }
	return(trans(r));
}


rowvec prox2(colvec v, double lambda, int k, int p, uvec res1, colvec w) {
	rowvec v2(v.n_elem);
	rowvec v3(p);
	for (int i = 0; i < k; ++i) {
        uvec bb = bbsubs(i, k, p);
        colvec v1 = v(bb);
        v3 = proxcppelem(v1, p, lambda, res1, w);
        v2(bb) = v3;
    }
	return(v2);
}


// [[Rcpp::export]]
mat FistaElem(const mat& Y, const mat& Z, mat phi, const int p, const int k, rowvec lambda, const double eps, const double tk, bool sep_lambda = false) {
	double j = 1;
	mat phiFin = phi;
	rowvec phiR = phi.row(0);
	rowvec phiOLD = phiR;
	rowvec phiOLDOLD = phiOLD;
	rowvec v = phiOLD;
	uvec res1 = ind(p, 0);
	colvec w(p);
	w.ones();
    for (int i = 0; i < k; ++i) {
        j = 1;
        double thresh = 10 * eps;
        phiR = phi.row(i);
        phiOLD = phiR;
        phiOLDOLD = phiOLD;
        v = phiR;
        double templambda;
        if (sep_lambda) {
            templambda = lambda(i);
        } else {
            templambda = lambda(0);
        }
        while (thresh > eps) {
            v = phiOLD + ((j - 2) / (j + 1)) * (phiOLD - phiOLDOLD);
            phiR = prox2(vectorise(v) + tk * vectorise((trans(Y.col(i)) - v * Z) * trans(Z)), tk * templambda, k, p, res1, w);
            thresh = max(abs(phiR - v));
            phiOLDOLD = phiOLD;
            phiOLD = phiR;
            j += 1;
        }
        phiFin.row(i) = phiR;
    }
	return(phiFin);
}

// Lamba loop
// [[Rcpp::export]]
cube gamloopElem(NumericVector beta_, const mat& Y, const mat& Z, mat gammgrid, const double eps, const colvec YMean2, const colvec ZMean2, mat B1, const int k, const int p, bool sep_lambda = false) {
    mat B1F2 = B1;
	vec eigval;
	mat eigvec;
	const mat Zt = Z * trans(Z);
	eig_sym(eigval, eigvec, Zt);
	double tk = 1 / max(eigval);
    IntegerVector dims = beta_.attr("dim");
    const int ngridpts = dims[2];
	cube bcube(beta_.begin(), dims[0], dims[1], ngridpts, false);
	cube bcube2(dims[0], dims[1] + 1, ngridpts);
	bcube2.fill(0);
	colvec nu = zeros<colvec>(k);
	for (int i = 0; i < ngridpts; ++i) {
		rowvec gamm = gammgrid.row(i);
		B1F2 = bcube.slice(i);
		B1 = FistaElem(Y, Z, B1F2, p, k, gamm, eps, tk, sep_lambda);
		nu = YMean2 - B1 * ZMean2;
		bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}
    return(bcube2);
}

// *
// *
// // Power Method Algorithm for eigenvalue computation
// *
// *

// [[Rcpp::export]]
List powermethod(mat A, colvec x1) {
	double dd = 1.0;
	int nn = x1.n_elem;
	arma::mat x(x1.begin(), nn, 1, false);
	double eps = .001;
	arma::mat y = x;
	double theta = 0;
	while (dd > eps * fabs(theta)) {
		x = y / arma::norm(y, 2);
		y = A * x;
		theta = as_scalar(trans(x) * y);
		dd = arma::norm(y - theta * x, 2);
	}
	double lambda = theta;
	Rcpp::List results = Rcpp::List::create(Named("lambda") = as<double>(wrap(lambda)), Named("q1") = as<NumericVector>(wrap(x)));
	return(results);
}


// Useful 2-norm function
// [[Rcpp::export]]
double norm2(NumericVector x) {
	arma::vec xx = x;
	double g = arma::norm(xx, 2);
	return (as<double>(wrap(g)));
}



// QR Factorization for Relaxed VAR
mat QRF(const mat& K, mat R5, int i, int kp, int k1) {
	int RC = R5.n_cols;
	mat RA = zeros(kp + k1, RC + 1);
	for (int j = 0; j < RC; ++j) {
        for (int q = 0; q < kp; ++q) {
            RA(q, j) = R5(q, j);
        }
    }
	for (int j = kp; j < kp + k1; ++j) {
		if (j - kp == i) {
            RA(j, RC) = 1;
        }
    }
	mat K2 = K * RA;
	int q = K2.n_cols;
	double delta = (pow(static_cast<double>(q), 2) + q + 1) * sqrt(std::numeric_limits<double>::epsilon());
	colvec D1 = zeros(K2.n_cols);
	for (int i = 0; i < q; ++i) {
        D1(i) = norm(K2.col(i), "fro");
    }
	D1 = sqrt(delta) * D1;
	mat AA = diagmat(D1);
	mat K3 = mat(join_vert(K2, AA));
	mat Q1, R1;
	qr(Q1, R1, K3);
	mat R11 = R1.submat(0, 0, RC - 1, RC - 1);
	colvec R22 = vectorise(R1.submat(0, RC, RC - 1, RC));
	mat RLS = solve(R11, R22);
	return(RLS);
}

// Relaxed Least Squares 
// [[Rcpp::export]]
mat RelaxedLS(const mat K, mat B2) {
	mat B3 = B2.cols(1, B2.n_cols - 1);
	if (norm(B3, "inf") == 0) {
	    return(B2);
    } else {
		int kp = B3.n_cols;
		int k1 = B3.n_rows;
		colvec nu = B2.col(0);
		mat A = zeros(k1, kp);
		for (int i = 0; i < k1; ++i) {
            rowvec B3a = B3.row(i);
            double thresh = 1e-8;
            uvec R1a = find(abs(B3a) > thresh);
            if (R1a.n_elem < 2) {
                A.row(i) = B3a;
            } else {
                mat R5 = zeros(kp, R1a.n_elem);
                int jj = 0;
                std::vector<int> R1(R1a.begin(), R1a.end());
                for (int ii = 0; ii < kp; ++ii) {
                    if (std::find(R1.begin(), R1.end(), ii) != R1.end()) {
                        R5(ii, jj) = 1;
                        jj += 1;
                    }
                }
                mat RLS = QRF(K, R5, i, kp, k1);
                mat R6 = R5 * RLS;
                R6.reshape(1, kp);
                A.row(i) = R6;
            }
        }
		mat BR = mat(join_horiz(nu, A));
		return(BR);
	}
}
	

// Sparse Lag while loop
mat sparseWLX(const mat& M1a, const mat& R1, double ngroups, mat& beta, double t, double alpha, double lambda, double eps) {
	int n = M1a.n_rows;
	int n2 = beta.n_rows, k2 = beta.n_cols;
	double thresh = 10;
	mat p = zeros<mat>(n2, k2);
	mat STS = zeros<mat>(n, 1);
	mat thetaOLD = beta;
	mat thetaOLDOLD = beta;
	double l = 1;
	mat u = zeros<mat>(n2, k2);
	while (thresh > eps) {
        p = (beta * M1a - R1) * trans(M1a) / ngroups;
        if (alpha > 0) {
            STS = ST3a(vectorise(beta) - t * vectorise(p), t * alpha * lambda);
        } else {
            STS = vectorise(beta) - t * vectorise(p);
        }
        double denom2 = norm(STS, "fro");
        double s31 = fmax(1 - (t * (1 - alpha) * lambda) / denom2, 0);
        STS = s31 * STS;
        STS.set_size(n2, k2);
        beta = thetaOLD + (l / (l + 3)) * (STS - thetaOLD);
        l += 1;
        mat thresh1 = beta - STS;
        thresh = norm(thresh1, "inf");
        thetaOLD = STS;
    }
	return(beta);
}

// Very similar to group lasso case
List blockUpdateSGLX(mat& beta, const mat& Z1, double lam, double alpha, const mat& Y1, double eps, List groups, List fullgroups, List compgroups, int k1, List M2f_, NumericVector Eigs) {
	int n1 = groups.size();
	List active(n1);
	int n = beta.n_rows, m = beta.n_cols;
	arma::mat betaPrev = beta;
	int converge = 0;
 
	if (groups.size() == 0) {
        beta.zeros(n, m);
        active = groups;
    } else {
		for (int i = 0; i < n1; ++i) {
            arma::uvec s45=as<arma::uvec>(groups[i]);
            arma::uvec s45F=as<arma::uvec>(fullgroups[i]);

            if (max(s45) == 0) {
                beta.cols(s45F) = zeros(k1, s45F.n_elem);
                active(i) = 0;
            } else {
                arma::uvec scomp1 = as<arma::uvec>(compgroups[i]);
                arma::mat M2a = Z1.rows(scomp1);
                arma::mat a1 = beta.cols(scomp1);
                arma::mat beta2 = beta.cols(s45);
                arma::mat r = Y1 - a1 * M2a;
                arma::mat M1 = Z1.rows(s45);
                arma::mat M2 = M2f_(i);
                arma::mat p = (beta2 * M1 - r) * trans(M1);
                double rho = sqrt(static_cast<double>(s45.size()));
                colvec STS;
                if (alpha > 0) {
                    STS = ST3a(vectorise(p), alpha * lam);
                } else {
                    STS = vectorise(p);
                }
                double lamadj = lam * (1 - alpha) * rho;
                if (arma::norm(STS, "fro") <= lamadj) {
                    mat astar = zeros(k1, s45.n_elem);
                    active(i) = 0;
                } else {
                    mat betaS = beta.cols(s45);
                    double t = 1 / Eigs(i);
                    mat astar2 = sparseWLX(M1, r, k1, betaS, t, alpha, lam, eps);
                    beta.cols(s45) = astar2;
                    active(i) = s45;
                }
            }
        }
	}
	arma::mat thresh1 = arma::abs((betaPrev - beta) / (arma::ones(n, m) + arma::abs(betaPrev)));
	double thresh = arma::norm(thresh1, "inf");
	if (thresh < eps) {
		converge = 1;
	} else {
		converge = 0;
	}
	Rcpp::List results = Rcpp::List::create(Named("beta") = wrap(beta), Named("active") = wrap(active), Named("Converge") = wrap(converge));
	return(results);
}

mat ThreshUpdateSGLX(mat& betaActive, const mat& Z, double lam, const mat& Y, double eps, List groups, List fullgroups, List compgroups, List M2f, NumericVector eigs, double alpha, int k1) {
	int n = betaActive.n_rows, m = betaActive.n_cols;
	int n1 = groups.size();
	mat betaLast = betaActive;
	List active(n1);
	int count = 0;
	List betaActive2(3);
	for (int i = 0; i < n1; ++i) {
        NumericVector g1 = groups[i];
        count += g1.size();
    }
	if (groups.size() == count) {
        betaActive.zeros(n, m);
        active = groups;
    } else {
		double threshold = 10 * eps;
		while(threshold > eps) {
            betaActive2 = blockUpdateSGLX(betaActive, Z, lam, alpha, Y, eps, groups, fullgroups, compgroups, k1, M2f, eigs);
            betaActive = as<mat>(betaActive2("beta"));
            arma::mat thresh1 = arma::abs((betaLast - betaActive) / (arma::ones(n, m) + arma::abs(betaLast)));
            threshold = arma::norm(thresh1, "inf");
            active = betaActive2("active");
            betaLast = betaActive;
        }
	}
	return(betaActive);
}

// [[Rcpp::export]]
List GamLoopSGLX(NumericVector beta_, List Activeset, NumericVector gamm, double alpha, const mat& Y1, const mat& Z1, List jj, List jjfull, List jjcomp, double eps, colvec YMean2, colvec ZMean2, int k, int pk, List M2f_, NumericVector eigs, int k1) {
    IntegerVector dims = beta_.attr("dim");
	int gran2 = dims[2];
	List activefinal(gran2);
	cube beta2(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube betafin(dims[0], dims[1] + 1, gran2);
	betafin.fill(0);
	List iterations(gran2);
	mat betaPrev = zeros<mat>(dims[0], dims[1]);
	for (int i = 0; i < gran2; ++i) {
        double gam = gamm[i];
        betaPrev = beta2.slice(i);
        List Active = Activeset[i];
        int k2 = 0;
        int converge = 0;
        mat betaF = zeros(dims[0], dims[1]);
        List betaFull(3);
        while (converge == 0) {
            betaPrev = ThreshUpdateSGLX(betaPrev, Z1, gam, Y1, eps, jjfull, jjfull, jjcomp, M2f_, eigs, alpha, k1);
            betaFull = blockUpdateSGLX(betaPrev, Z1, gam, alpha, Y1, eps, jjfull, jjfull, jjcomp, k1, M2f_, eigs);
            betaF = as<mat>(betaFull("beta"));
            Active = betaFull("active");
            converge = betaFull("Converge");
            k2 += 1;
        }
        colvec nu = YMean2 - betaF * ZMean2;
        betafin.slice(i) = mat(join_horiz(nu, betaF));
        activefinal[i] = Active;
        iterations[i] = k2;
    }
	List Results = List::create(Named("beta") = betafin, Named("active") = wrap(activefinal), Named("iterations") = iterations);
	return(Results);
}



// Prox function for endogenous first VARX
// [[Rcpp::export]]
colvec proxvx2(colvec v2, int L, double lambda, int m, int k, int F1) {
	colvec r = v2;
	int start;
	if (F1 == 0) {
	    start = 1;
	} else {
	    start = 0;
	}
	for (int i = start; i >= 0; --i) {
        int s1;
        if (F1 == 0) {
            s1 = (k + m) - (i * k);
        } else {
            s1 = k;
        }
        std::vector<unsigned int> ivec(s1);
        if (F1 == 0) {
            std::iota(ivec.begin(), ivec.end(), i * k);
        } else {
            std::iota(ivec.begin(), ivec.end(), 0);
        }
        uvec res = conv_to<uvec>::from(ivec);
        if (norm(r(res) / (lambda), "fro") < 1 + 1e-8) {
            r(res) = zeros(res.n_elem);
        } else {
            r(res) = r(res) - lambda * r(res) / (norm(r(res), "fro"));
        }
    }
    return(r);
}


mat sparseWL(const mat& M1a, const mat& R1, double ngroups, mat& beta, double t, double alpha, double lambda, double eps) {
	int n = M1a.n_rows;
	int n2 = beta.n_rows;
	double thresh = 10;
	mat p = zeros<mat>(n, n);
	colvec STS(n);
	mat thetaOLD = beta;
	mat thetaOLDOLD = beta;
	double l = 1;
	mat u = zeros<mat>(n2, n2);
	while (thresh > eps) {
        p = (beta * M1a - R1) * trans(M1a) / ngroups;
        if (alpha > 0) {
            STS = ST3a(vectorise(beta) - t * vectorise(p), t * alpha * lambda);
        } else {
            STS = vectorise(beta) - t * vectorise(p);
        }
        // DON'T USE TWO NORM ON MATRICES IN ARMADILLO, IT IS MATRIX 2 NORM
        double denom2 = norm(STS, "fro");
        double s31 = fmax(1 - (t * (1 - alpha) * lambda) / denom2, 0);
        STS = s31 * STS;
        u = reshape(STS, n2, n2);
        // nesterov
        beta = thetaOLD + (l / (l + 3)) * (u - thetaOLD);
        l += 1;
        mat thresh1 = beta - u;
        thresh = norm(thresh1, "inf");
        thetaOLD = u;
    }
	return(beta);
}

// Very similar to group lasso case
List blockUpdateSGL(mat& beta, const mat& Z1, double lam, double alpha, const mat& Y1, double eps, List groups, List fullgroups, List compgroups, int k, List M1f_, List M2f_, NumericVector Eigs) {
    int n1 = groups.size();
	List active(n1);
	int n = beta.n_rows, m = beta.n_cols;

	arma::mat betaPrev = beta;
	int converge = 0;
	int count = 0;
 
	if (groups.size() == count) {
		beta.zeros(n, m);
		active = groups;
	} else {
		for (int i = 0; i < n1; ++i) {
            NumericVector s1 = groups[i];
            IntegerVector s2 = fullgroups[i];
            NumericVector scomp = compgroups[i];
            arma::uvec s45 = as<arma::uvec>(s1);
            arma::uvec s45F = as<arma::uvec>(s2);

            if (s1.size() == 1) {
                beta.cols(s45F) = arma::zeros(s2.size(), s2.size());
                active(i) = 0;
            }
            if (s1.size() != 1) {
                arma::colvec s3(s1.begin(), s1.size(), false);
                arma::colvec scomp1(scomp.begin(), scomp.size(), false);
                arma::uvec s4(s3.size());

                for (int j = 0; j < (int) s3.size(); ++j) {
                    s4(j) = s3(j);
                }

                arma::uvec scomp2(scomp1.size());
                for (int m = 0; m < (int) scomp1.size(); ++m) {
                    scomp2(m) = scomp1(m);
                }

                arma::mat M2a = Z1.rows(scomp2);
                arma::mat a1 = beta.cols(scomp2);
                arma::mat beta2 = beta.cols(s4);
                arma::mat r = Y1 - a1 * M2a;
                arma::mat M1 = M1f_(i);
                arma::mat M2 = M2f_(i);
                arma::mat p = (beta2 * M1 - r) * trans(M1);

                colvec STS;
                if (alpha > 0) {
                    STS = ST3a(vectorise(p), alpha * lam);
                } else {
                    STS = vectorise(p);
                }

                double lamadj = lam * (1 - alpha);
                if (arma::norm(STS, "fro") <= lamadj) {
                    arma::mat astar = arma::zeros(s3.size(), s3.size());
                    active(i) = 0;
                    beta.cols(s4) = astar;
                } else {
                    mat betaS = beta.cols(s4);
                    double t = 1 / Eigs(i);
                    mat astar2 = sparseWL(M1, r, k, betaS, t, alpha, lam, eps);
                    beta.cols(s4) = astar2;
                    active(i) = s4;
                }
            }
        }
	}
	arma::mat thresh1 = arma::abs((betaPrev - beta) / (arma::ones(n, m) + arma::abs(betaPrev)));
	double thresh = arma::norm(thresh1, "inf");
	if (thresh < eps) {
		converge = 1;
	} else {
		converge = 0;
	}
	Rcpp::List results = Rcpp::List::create(Named("beta") = wrap(beta), Named("active") = wrap(active), Named("Converge") = wrap(converge));
	return(results);
}

mat ThreshUpdateSGL(mat& betaActive, const mat& Z, double lam, const mat& Y, double eps, List groups, List fullgroups, List compgroups, List M1f, List M2f, NumericVector eigs, double alpha, int k) {
	int n = betaActive.n_rows, m = betaActive.n_cols;
	int n1 = groups.size();
	mat betaLast = betaActive;
	List active(n1);
	int count = 0;
	List betaActive2(3);
	int iters = 0;
	for (int i = 0; i < n1; ++i) {
        NumericVector g1 = groups[i];
        count += g1.size();
    }
	if (groups.size() == count) {
        betaActive.zeros(n, m);
        active = groups;
	} else {
		double threshold = 10 * eps;
		while (threshold > eps) {
            betaActive2 = blockUpdateSGL(betaActive, Z, lam, alpha, Y, eps, groups, fullgroups, compgroups, k, M1f, M2f, eigs);
            betaActive = as<mat>(betaActive2("beta"));
            arma::mat thresh1 = arma::abs((betaLast - betaActive) / (arma::ones(n, m) + arma::abs(betaLast)));
            threshold = arma::norm(thresh1, "inf");
            active = betaActive2("active");
            betaLast = betaActive;
            iters += 1;
        }
	}
	return(betaActive);
}

// [[Rcpp::export]]
List GamLoopSGL(NumericVector beta_, List Activeset, const NumericVector gamm, const double alpha, const mat& Y1, const mat& Z1,List jj, const List jjfull, const List jjcomp, const double eps, const colvec YMean2, const colvec ZMean2, const int k, const int pk, const List M1f_, const List M2f_, const NumericVector eigs_) {
	IntegerVector dims = beta_.attr("dim");
	int gran2 = dims[2];
	List activefinal(gran2);

	cube beta2(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube betafin(dims[0], dims[1] + 1, dims[2]);
	betafin.fill(0);
	List iterations(gran2);
	mat betaPrev = zeros<mat>(dims[0], dims[1]);

	for (int i = 0; i < gran2; ++i) {
        double gam = gamm[i];
        betaPrev = beta2.slice(i);
        List Active = Activeset[i];
        int k2 = 0;
        int converge = 0;
        mat betaF = zeros(dims[0], dims[1]);
        List betaFull(3);
        //Three components in the list
        while (converge == 0) {
            betaPrev = ThreshUpdateSGL(betaPrev, Z1, gam, Y1, eps, Active, jjfull, jjcomp, M1f_, M2f_, eigs_, alpha, k);
            betaFull = blockUpdateSGL(betaPrev, Z1, gam, alpha, Y1, eps, jj, jjfull, jjcomp, k, M1f_, M2f_, eigs_);
            betaF = as<mat>(betaFull("beta"));
            Active = betaFull("active");
            converge = betaFull("Converge");
            k2 += 1;
        }

        colvec nu = YMean2 - betaF * ZMean2;
        betafin.slice(i) = mat(join_horiz(nu, betaF));
        activefinal[i] = Active;
        iterations[i] = k2;
    }
	List Results = List::create(Named("beta") = betafin, Named("active") = wrap(activefinal), Named("iterations") = iterations);
	return(Results);
}


// [[Rcpp::export]]
List GamLoopSGLDP(NumericVector beta_, List Activeset, const mat gamm, const colvec alpha, const mat& Y1, const mat& Z1, List jj, const List jjfull, const List jjcomp, const double eps, const colvec YMean2, const colvec ZMean2, const int k, const int pk, const List M1f_, const List M2f_, const NumericVector eigs_) {
	int nlambda = gamm.n_rows;
	int nalpha = gamm.n_cols;
	List activefinal(nlambda * nalpha);
	IntegerVector dims = beta_.attr("dim");

	cube beta2(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube betafin(dims[0], dims[1] + 1, dims[2]);
	betafin.fill(0);
	List iterations(nlambda * nalpha);
	mat betaPrev = zeros<mat>(dims[0], dims[1]);

	for (int i = 0; i < nlambda; ++i) {
        for(int j = 0; j < nalpha; ++j) {
            double gam = gamm(i, j);
            double alpha1 = alpha(j);
            betaPrev = beta2.slice((i) * nalpha + j);
            List Active = Activeset[(i) * nalpha + j];
            int k2 = 0;
            int converge = 0;
            mat betaF = zeros(k, pk);
            List betaFull(3);
            //Three components in the list
            while (converge == 0) {
                betaPrev = ThreshUpdateSGL(betaPrev, Z1, gam, Y1, eps, Active, jjfull, jjcomp, M1f_, M2f_, eigs_, alpha1, k);
                betaFull = blockUpdateSGL(betaPrev, Z1, gam, alpha1, Y1, eps, jj, jjfull, jjcomp, k, M1f_, M2f_, eigs_);
                betaF = as<mat>(betaFull("beta"));
                Active = betaFull("active");
                converge = betaFull("Converge");
                k2 += 1;
            }
            colvec nu = YMean2 - betaF * ZMean2;
            betafin.slice((i) * nalpha + j) = mat(join_horiz(nu, betaF));
            activefinal[(i) * nalpha + j] = Active;
            iterations[(i) * nalpha + j] = k2;
        }
    }
	List Results = List::create(Named("beta") = betafin, Named("active") = wrap(activefinal), Named("iterations") = iterations);
	return(Results);
}




// [[Rcpp::export]]
List GamLoopSGLXDP(NumericVector beta_, List Activeset, mat gamm, colvec alpha, const mat& Y1, const mat& Z1,List jj, List jjfull, List jjcomp, double eps, colvec YMean2, colvec ZMean2, int k, int pk, List M2f_, NumericVector eigs, int k1) {
	int nlambda = gamm.n_rows;
	int nalpha = gamm.n_cols;
	List activefinal(nlambda * nalpha);
	IntegerVector dims = beta_.attr("dim");
	cube beta2(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube betafin(dims[0], dims[1] + 1, dims[2]);
	betafin.fill(0);
	List iterations(nlambda * nalpha);
	mat betaPrev = zeros<mat>(dims[0], dims[1]);
	for (int i = 0; i < nlambda; ++i) {
        for (int j = 0; j < nalpha; ++j) {
            double gam = gamm(i, j);
            double alpha1 = alpha(j);
            betaPrev = beta2.slice((i) * nalpha + j);
            List Active = Activeset[(i) * nalpha + j];
            int k2 = 0;
            int converge = 0;
            mat betaF = zeros(k1, pk);
            List betaFull(3);
            while (converge == 0) {
                betaPrev = ThreshUpdateSGLX(betaPrev, Z1, gam, Y1, eps, jjfull, jjfull, jjcomp, M2f_, eigs, alpha1, k1);
                betaFull = blockUpdateSGLX(betaPrev, Z1, gam, alpha1, Y1, eps, jjfull, jjfull, jjcomp, k1, M2f_, eigs);
                betaF = as<mat>(betaFull("beta"));
                Active = betaFull("active");
                converge = betaFull("Converge");
                k2 += 1;
            }
            colvec nu = YMean2 - betaF * ZMean2;
            betafin.slice((i) * nalpha + j) = mat(join_horiz(nu, betaF));
            activefinal[(i) * nalpha + j] = Active;
            iterations[(i) * nalpha + j] = k2;
        }
    }
	List Results = List::create(Named("beta") = betafin, Named("active") = wrap(activefinal), Named("iterations") = iterations);
	return(Results);
}


double MCP(double z, double l1, double l2, double gamma, double v) {
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma * l1 * (1 + l2)) return(s * (fabs(z) - l1) / (v * (1 + l2 - 1 / gamma)));
  else return(z / (v * (1 + l2)));
}

double SCAD(double z, double l1, double l2, double gamma, double v) {
  double s = 0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1 * (1 + l2) + l1)) return(s * (fabs(z) - l1) / (v * (1 + l2)));
  else if (fabs(z) <= gamma * l1 * (1 + l2)) return(s * (fabs(z) - gamma * l1 / (gamma - 1)) / (v * (1 - 1 / (gamma - 1) + l2)));
  else return(z / (v * (1 + l2)));
}



// We solve along the regularization path due to non-convexity concerns, need good starting vals  fit by series
// [[Rcpp::export]]
cube mcp_loop(mat Y, mat Z, cube B, const vec lambda, const double tol, double gamma, bool mcp = true) {
	int n1 = Z.n_rows;
	vec cps(Z.n_rows);
	double Z2 = Z.n_cols;
	for (int i = 0; i < n1; ++i) {
		cps(i) = as_scalar(Z.row(i) * trans(Z.row(i))) / Z2;
	}
	vec active = zeros(B.n_cols);
	double gloss;
	double T;
	double sdy;
	double pk = Z.n_rows;
	double k = B.n_rows;
	mat BOLD_FULL;
	int nlambdas = B.n_slices;
	for (int j = 0; j < nlambdas; ++j) {
		mat BPH = zeros(k, pk);
		if (j > 0) {
			BOLD_FULL = B.slice(j - 1);
		} else {
			BOLD_FULL = zeros(k, pk);
		}
		for (int i = 0; i < pk; ++i) {
			if (BOLD_FULL(i) != 0) {
				active(i) = 1;
			} else {
				active(i) = 0;
			}
		}
		mat BTemp_FULL = B.slice(j);
		int nrows = BTemp_FULL.n_rows;
		for (int jj = 0; jj < nrows; ++jj) {
			rowvec BTemp = BTemp_FULL.row(jj);
			rowvec BOLD = BOLD_FULL.row(jj);
			vec r = Y.col(jj);
			int len_r = r.n_elem;
			vec r2 = zeros(len_r);
			for (int ii = 0; ii < len_r; ++ii) {
				r2(ii) = pow(r(ii), 2);
			}
			gloss = sum(r2);
			T = Y.n_rows;
			sdy = sqrt(gloss / T);

			double thresh = 10 * tol;
			int counter = 0;

			while (counter < 100) {
                counter += 1;
				int counter2 = 0;
				while (counter2 < 100) {
                    thresh = 0;
                    for (int i = 0; i < pk; ++i) {
                        if (active(i) != 0) {
                            vec Z2 = trans(Z.row(i));
                            double G1 = as_scalar(trans(r) * Z2) / T + BOLD(i) * cps(i);
                            if (mcp) {
                                BTemp(i) = MCP(G1, lambda(j), lambda(j), gamma, cps[i]);
                            } else {
                                BTemp(i) = SCAD(G1, lambda(j), lambda(j), gamma, cps[i]);
                            }

                            double shift = BTemp(i) - BOLD(i);
                            if (shift != 0) {
                                r = r - trans(shift * Z.row(i));
                            }
                        }
                    }
                    vec thresh1 = abs(trans(BTemp) - trans(BOLD));
                    thresh = norm(thresh1, "inf");
                    BOLD = BTemp;

                    if (thresh < tol * sdy) {
                        break;
                    }
                    counter2 += 1;
                }

				int violations = 0;
				for (int i = 0; i < pk; ++i) {
					if (active(i) == 0) {
						vec Z2 = trans(Z.row(i));
						double G1 = as_scalar(trans(r) * Z2) / T;
						if (mcp) {
							BTemp(i) = MCP(G1, lambda(j), lambda(j), gamma, cps[i]);
						} else {
							BTemp(i) = SCAD(G1, lambda(j), lambda(j), gamma, cps[i]);
						}

						if (BTemp(i) != 0){
							active(i) = 1;
							r = r - trans(BTemp(i) * Z.row(i));
							violations += 1;
						}
					}
				}
				BOLD = BTemp;
				if (violations == 0) break;
			}
			BPH.row(jj) = BTemp;
		}
		B.slice(j) = BPH;
	}
	return(B);
}

// [[Rcpp::export]]
cube gamloopMCP(NumericVector beta_, const mat& Y, const mat& Z, vec lambda, const double eps, const colvec& YMean2, const colvec& ZMean2, double gamma, bool mcp) {

	IntegerVector dims = beta_.attr("dim");
	cube bcube(beta_.begin(), dims[0], dims[1], dims[2], false);
	cube bcube2(dims[0], dims[1] + 1, dims[2]);
	bcube2.fill(0);
	
	colvec nu = zeros<colvec>(dims[0]);

	bcube = mcp_loop(Y, Z, bcube, lambda, eps, gamma, mcp);
	for (int i = 0; i < dims[2]; ++i) {
		mat B1 = bcube.slice(i);
		nu = YMean2 - B1 * ZMean2;
		bcube2.slice(i) = mat(join_horiz(nu, B1));
	}

    return(bcube2);
}
