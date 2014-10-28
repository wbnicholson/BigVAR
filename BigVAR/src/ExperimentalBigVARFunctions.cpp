
#include <RcppArmadillo.h>
#include <math.h> 
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Soft thresholding
double ST1a(double z,double gam){
double sparse=0;
if(z>0 && gam<fabs(z)) return(z-gam);

if(z<0 && gam<fabs(z)) return(z+gam);
if(gam>=fabs(z)) return(sparse);
else return(0);
}

// Columnwise softthresholding
colvec ST3a(colvec z ,double gam)
{
int n=z.size();
colvec z1(n);
for( int i=0; i<n;++i)
{
double z11=z(i);
z1(i)=ST1a(z11,gam);
}
return(z1);
}
// Columnwise softthresholding pass by reference for Lasso-VAR
 colvec ST3ar(const colvec& z ,const double gam)
{
int n=z.size();
colvec z1(n);
for( int i=0; i<n;++i)
{
double z11=z(i);
z1(i)=ST1a(z11,gam);
}
return(z1);
}

// Indexing function for residuals
uvec ind(int n2,int m){
IntegerVector subs(n2);
for(int i =0 ; i<n2;++i)
{
subs(i)=i;
//	Rcpp::Rcout << "R" << subs[i] << std::endl;
 }
subs.erase(m);
return(as<uvec>(subs));
}

// Rowwise softthresholding
rowvec ST3bc(rowvec& z,double gam)
{
int n=z.size();
rowvec z1(n);
for( int i=0; i<n;++i)
{
double z11=z(i);
z1(i)=ST1a(z11,gam);
}
return(z1);
}

// Lasso Fista Function
// [[Rcpp::export]]
mat FistaLV(const mat& Y, const mat& Z, mat& B, const double gam, const double eps, double tk, int k,int p)
{
  B=trans(B);
  // rowvec B1=B.row(0);
  colvec B1=B.col(0);
 
  double j =1;

  for( int i =0; i<k; ++i)

    {
      B1=B.col(i);
      colvec BOLD=B.col(i);
      colvec BOLDOLD=BOLD;
	// v=B1;
     double thresh=10*eps;
      j=1;
      while(thresh>eps)
	{
	 colvec v=BOLD+((j-2)/(j+1))*(BOLD-BOLDOLD);

	 // If I pass by reference, I need to declare the varible first.  Otherwise it makes no sense to do so!
       	 B1=ST3a(vectorise(v)+tk*vectorise((trans(Y.col(i))-trans(v)*Z)*trans(Z)),gam*tk);
	  thresh=max(abs(B1-v));
	  BOLDOLD=BOLD;
	  BOLD=B1;
	  j+=1;
    }
      B.col(i)=B1;
    }

  B=trans(B);
  return(B);
} 

//Gam Loop For FISTA

// [[Rcpp::export]]
cube gamloopFista(NumericVector beta_, const mat& Y,const mat& Z,const  NumericVector gamm, const double eps,const colvec& YMean2, const colvec& ZMean2,mat& B1, int k, int p){

  //Data is read in from R in the form of NumericMatrix and NumericVector, but needs to be convered to the armadillo library

  //These are constant data vectors

vec eigval;
mat eigvec;
const mat Zt=Z*trans(Z);
eig_sym(eigval, eigvec, Zt);

double tk=1/max(eigval); 

//Penalty parameters that we iterate across
double gammgrid[gamm.size()];
// int gran2=gamm.size();

//Constant data matrices

//Coefficient matrices
 mat b2=B1;
 mat B1F2=B1;
 //These are typically defined inside of the "lassocore" function, but doing so causes problems with openMP, so they are defined out here and read in
 const int ngridpts=gamm.size();
 cube bcube(beta_.begin(),k,k*p,ngridpts,false);
 cube bcube2(k,k*p+1,ngridpts);
 bcube2.fill(0);
 
colvec nu=zeros<colvec>(k);
double gam =0;

for(int j=0;j<ngridpts;++j)
{
gammgrid[j]=gamm[j];
}
 int i;
  for (i=0; i<ngridpts;++i) {
            gam=gammgrid[i];

	mat B1F2=bcube.slice(i);
	B1 = FistaLV(Y,Z,B1F2,gam,eps,tk,k,p); 
	

  
	 nu = YMean2 - B1 *ZMean2;
         bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}

    return(bcube2);
}
/*
// Block Group Lasso Functions
*
*
*
*/

// Newton-Raphson Functions
double trust32(int k,const arma::mat PP,double delta, double lambda, arma::vec EigVA,arma::mat EigVector){
          // int k=M2.n_cols;
         double g=0;
         for(int i = 0; i < k; ++i)
         { 
    
	   g+=pow(arma::as_scalar(trans(EigVector.col(i))*PP),2)/pow(EigVA[i]*delta+lambda,2);

         }
	 return(g);
}

double fprime2(int k,const arma::mat PP,double delta, double lambda, arma::vec EigVA,arma::mat EigVE)
{
         // int k=M2.n_cols;
         double gg2=0;
         for(int i = 0; i < k; ++i)
         { 
    
	   gg2+=(pow(arma::as_scalar(trans(EigVE.col(i))*PP),2)*EigVA[i])/pow(EigVA[i]*delta+lambda,3);


         }
	 double c1=trust32(k, PP, delta, lambda, EigVA, EigVE);
	 double res=-.5*pow(c1,-1.5)*-2*gg2;
	 return(res);

}

double Newton2(int k,const  arma::mat P, double lambda, arma::vec EigVA,arma::mat EigVE){

  double delta=0;
  double threshold=1;
  double phi=0;
  double deltanew=delta;
  int iter=0;
  while(threshold>.0001)
    {
   phi=1-1/pow(trust32(k, P, delta, lambda, EigVA, EigVE),.5);
   deltanew+=phi/fprime2(k, P, deltanew, lambda, EigVA, EigVE);
  threshold = fabs(delta-deltanew);
  delta=deltanew;
  phi=0;
	  iter+=1;

    }
  

  return(deltanew);




}

// Eigen-Decomposition
 // [[Rcpp::export]]
List Eigencomp( mat Z1, List groups,int n1,int k) 
{
  List M2f(n1);
  List eigvalF(n1);
  List M3f(n1);
  List eigvecF(n1);
  int count=0;
 
  for(int i=0; i<n1; ++i)
    {
      NumericVector g1=groups[i];
      count+=g1.size();
      arma::uvec s4=as<arma::uvec>(g1);
      arma::mat M1=Z1.rows(s4);
      arma::mat M2=M1*trans(M1);
      M2f(i)=M2;
      arma::mat D(k,k);
      D.eye();
      arma::mat M3=kron(M2,D);
      M3f(i)=M3;
      arma::vec eigval;
      arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, M3);
      eigvalF(i)=eigval;
      eigvecF(i)=eigvec;

    }

 List results=List::create(Named("M3")=wrap(M3f),Named("eigval")=wrap(eigvalF),Named("eigvec")=wrap(eigvecF));
  return(results);

}

// Eigen Decomposition for own/other group lasso
 // [[Rcpp::export]]
List EigencompOO( mat ZZ1, List groups,int n1,int k) 
{
  List M2f(n1);
  List eigvecF(n1);
  List eigvalF(n1);
  int count=0;
  for(int i=0; i<n1; ++i)
    {
      NumericVector g1=groups[i];
      arma::uvec s4=as<arma::uvec>(g1);
      count+=g1.size();
      arma::mat M1=ZZ1.cols(s4);
      arma::mat M2=trans(M1)*M1;
      M2f(i)=M2;
	   
	  arma::vec eigval;
	  arma::mat eigvec;
	  arma::eig_sym(eigval, eigvec, M2);
	  eigvecF(i)=eigvec;
	  eigvalF(i)=eigval;

    }
 

 List results=List::create(Named("M3")=wrap(M2f),Named("eigval")=wrap(eigvalF),Named("eigvec")=wrap(eigvecF));
  return(results);

}



// Block Group Lasso Update
// [[Rcpp::export]]
List BlockUpdateGL(mat& beta,const mat& Z1, double lam, const mat& Y1,double eps, List groups, List fullgroups, List compgroups, int k,List M3f_,List eigvalF_, List eigvecF_){
  int n1=groups.size();
  List active(n1);
  int n=beta.n_rows, m=beta.n_cols;
  arma::mat betaPrev=beta;
  int converge=0;
  int count=0;
 
  if(groups.size()==count)
    {
      beta.zeros(n,m);
      active=groups;
    }
  else{
    for(int i=0; i<n1;++i)
 
      {
	// c++ knows nothing about data types in R for lists, need to be explicit!
      // Note: Cannot have null elements in a list, just set them equal to 0 and use length 1 as the criterion
	NumericVector s1=groups[i];
	IntegerVector s2=fullgroups[i];
	NumericVector scomp=compgroups[i];
	arma::uvec s45=as<arma::uvec>(s1);
	arma::uvec s45F=as<arma::uvec>(s2);


      

	
        if(s1.size()==1){
	  
	  beta.cols(s45F)=arma::zeros(s2.size(),s2.size());
	  active(i)=0;

	}
	if(s1.size()!=1){
	  
	  arma::colvec s3(s1.begin(),s1.size(),false);
	  arma::colvec scomp1(scomp.begin(),scomp.size(),false);

	  // Need index vectors to be uvecs, also need to make sure vecs start indexing at zero!
	  arma::uvec s4(s3.size());

	  for(int j=0; j<s3.size(); ++j)
	    {
	      s4(j)=s3(j);
	    }
	  
	  arma::uvec scomp2(scomp1.size());

	  for(int m=0; m<scomp1.size(); ++m)
	    {
	      scomp2(m)=scomp1(m);
	    }

          mat r = beta.cols(scomp2)*Z1.rows(scomp2)-Y1;
	  arma::mat M3=M3f_(i);
	  mat p=r*trans(Z1.rows(s4));


	  if(arma::norm(p,"fro")<=lam)
	    {
	      arma::mat astar=arma::zeros(s3.size(),s3.size());
	      active(i)=0;


	    }
	  else{
	    int k1=M3.n_cols;
	    double deltfin=  Newton2(k1,arma::vectorise(p),lam,eigvalF_(i),eigvecF_(i));
	
	    arma::mat D1(k1,k1);
	    D1.eye();

	arma::mat astar=-solve(M3+lam/deltfin*D1,arma::vectorise(p));
	arma::mat astar2(astar.begin(),s4.size(),s4.size(),false);
	beta.cols(s4)=astar2;
	
	 active(i)=s4;  

	  }

	
	


	}

      }



  }
  // arma::mat thresh1=arma::abs((betaPrev-beta)/(arma::ones(n,m)+arma::abs(betaPrev)));
  double thresh=arma::norm(betaPrev-beta,"inf");
  if(thresh<eps)
    {
      converge=1;

    }
  else{ 
converge=0;
}


  Rcpp::List results=Rcpp::List::create(Named("beta")=beta,Named("active")=wrap(active),Named("Converge")=wrap(converge));
  return(results);

}



// [[Rcpp::export]]
mat ThreshUpdate(mat& betaActive,const mat& Z1, double lam, const mat& Y1,double eps, List groups, List fullgroups, List compgroups,List M2f_, List eigvalF_, List eigvecF_)
  {

    int n=betaActive.n_rows, m=betaActive.n_cols;
  int n1=groups.size();
  mat betaLast=betaActive;
  List active(n1);
  int count=0;
  List betaActive2(3);
  for(int i=0; i<n1; ++i)
    {
      NumericVector g1=groups[i];
      count+=g1.size();

    }
 
  if(groups.size()==count)
    {

      betaActive.zeros(n,m);

      active=groups;
    }
  else{
   double threshold=10*eps;
    while(threshold>eps)
      {
        betaActive2=BlockUpdateGL(betaActive,Z1,lam,Y1,eps,groups,fullgroups,compgroups,n,M2f_,eigvalF_,eigvecF_);	 
        betaActive=as<mat>(betaActive2("beta"));
	// arma::mat thresh1=arma::abs((betaLast-betaActive)/(arma::ones(n,m)+arma::abs(betaLast)));
        threshold=arma::norm(betaLast-betaActive,"inf");
        active=betaActive2("active");
	betaLast=betaActive;

      }
  }
    return(betaActive);
  }

// [[Rcpp::export]]
List GamLoopGL(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y1, const mat& Z1,List jj, List jjfull, List jjcomp, double eps,const colvec& YMean2, const colvec&  ZMean2,int k,int pk,const List M2f_, const List eigvalF_, const List eigvecF_)
{
// arma::colvec YMean2=YMean;
// arma::colvec ZMean2=ZMean;
 int gran2=gamm.size();
 List activefinal(gran2);
 cube beta2(beta_.begin(),k,pk,gran2,false);
 cube betafin(k,pk+1,gran2);
 betafin.fill(0);
 List iterations(gran2);
 mat betaPrev=zeros<mat>(k,pk);
 
int n3=Z1.n_rows, m3=Z1.n_cols;
 // arma::mat Z1(Z.begin(),n3,m3,false);
int n4=Y1.n_rows, m4=Y1.n_rows;
// arma::mat Y1(Y.begin(),n4,m4,false);
 
 //INDEX LISTS WITH PARENTHESES NOT BRACKETS 
 // WHEN EXTRACTING FROM A LIST YOU NEED as<MAT>

 for(int i=0; i<gran2;++i)
    {
        double gam=gamm[i];
	betaPrev=beta2.slice(i);
	List Active = Activeset[i];
	int k2=0;
	int converge=0;
        mat betaF=zeros(k,k);
	List betaFull(3);
	int thresh=0;
	while(converge==0)
	  {
	    
       	    betaPrev = ThreshUpdate(betaPrev, Z1, gam, Y1, eps, Active, jjfull, jjcomp,M2f_,eigvalF_,eigvecF_);
	 
	    betaFull=BlockUpdateGL(betaPrev,Z1,gam,Y1,eps,jjfull,jjfull,jjcomp,k,M2f_,eigvalF_,eigvecF_);
	     betaF=as<mat>(betaFull("beta"));
	     Active=betaFull("active");
	     converge =betaFull("Converge");
             k2+=1;

	  }
  colvec nu= YMean2 - betaF *ZMean2;
  betafin.slice(i)=mat(join_horiz(nu, betaF));
    activefinal[i]=Active;
  iterations[i]=k2; 
    }
 List Results=List::create(Named("beta")=betafin,Named("active")=wrap(activefinal),Named("iterations")=iterations);
   return(Results);
    }
// *
// *
// *
// *
// *
// *
// *
// *
// Group Lasso Own/Other Functions
// *
// *
// *
// *
// *

// [[Rcpp::export]]
List BlockUpdate2(const mat& ZZ1, double lam,const mat& Y1,double eps, List groups, List fullgroups, List compgroups, int k, List M2f_, List eigvalF_, List eigvecF_,colvec& B){
  int n1=groups.size();
  List active(n1);
  // int n=beta.n_rows, m=beta.n_cols;
  // arma::colvec B=arma::vectorise(beta);
  // int n3=ZZ.nrow(), m3=ZZ.ncol();
  // arma::mat ZZ1(ZZ.begin(),n3,m3,false);
  // int n4=Y.nrow(), m4=Y.ncol();
  // arma::mat Y1(Y.begin(),n4,m4,false);
  colvec  BPrev=B;
  int converge=0;
  int count=0;
 
  if(groups.size()==count)
    {
      // beta.zeros(n,m);
      B.zeros();
      active=groups;
    }
  else{
    for(int i=0; i<n1;++i)
 
      {
	NumericVector s1=groups[i];
	IntegerVector s2=fullgroups[i];
	NumericVector scomp=compgroups[i];
	arma::uvec s45=as<arma::uvec>(s1);
	arma::uvec s45F=as<arma::uvec>(s2);


      

	
        if(s1.size()==1){
	  
	  B.elem(s45F)=arma::zeros(s2.size());
	  active(i)=0;

	}
	if(s1.size()!=1){
	  
	  arma::colvec s3(s1.begin(),s1.size(),false);
	  // Need index vectors to be uvecs, also need to make sure vecs start indexing at zero!
	  arma::uvec s4=as<arma::uvec>(s1);

	  
	  arma::uvec scomp2=as<arma::uvec>(scomp);

	  arma::mat M2a= ZZ1.cols(scomp2);
	  arma::colvec a1= B.elem(scomp2);

	  arma::colvec r=M2a*a1-arma::vectorise(Y1);

	  arma::mat M1=ZZ1.cols(s4);
	  arma::mat M2=M2f_(i);
	  arma::vec eigval=eigvalF_(i);
	  arma::mat eigvec=eigvecF_(i);
	  arma::mat p=trans(M1)*r;

	  double rho=sqrt(s3.size());
	  double adjlam=rho*lam;

	  if(arma::norm(p,2)<=adjlam)
	    {
	      arma::colvec astar=arma::zeros(s3.size());
	      active(i)=0;


	    }
	  else{
        int k1=M2.n_cols;
    
	double deltfin=  Newton2(k1,p,adjlam,eigval,eigvec);
	
	arma::mat D1(s3.size(),s3.size());
	D1.eye();
	arma::mat astar=-solve(M2+adjlam/deltfin*D1,p);

	 B.elem(s4)=astar;
	
	 active(i)=s4;  

	  }

	
	


	}

      }



  }
  // arma::mat betafin(B.begin(),n,m,false);
  // arma::mat thresh1=arma::abs((beta-betafin)/(arma::ones(n,m)+arma::abs(beta)));
  // double thresh=arma::norm(thresh1,"inf");
  double thresh=arma::norm(B-BPrev,"inf");
  
  if(thresh<eps)
    {
      converge=1;

    }
  else{ 
converge=0;
}


  Rcpp::List results=Rcpp::List::create(Named("beta")=B,Named("active")=active,Named("Converge")=converge);
  return(results);

}
// [[Rcpp::export]]
colvec ThreshUpdateOO(const mat& ZZ, double lam,const mat& Y,double eps, List groups, List fullgroups, List compgroups,List M2f_, List eigvalF_, List eigvecF_,colvec& B, int n)
  {

  int kp=B.n_elem;
  int n1=groups.size();
  colvec BPrev=B;
  List active(n1);
  int count=0;
 List betaActive2(3);

  for(int i=0; i<n1; ++i)
    {
      NumericVector g1=groups[i];
      count+=g1.size();

    }
 
  if(groups.size()==count)
    {

      B.zeros(kp);

      active=groups;
    }
  else{
   double threshold=10*eps;
    while(threshold>eps)
      {
	betaActive2=BlockUpdate2(ZZ,lam,Y,eps,groups,fullgroups,compgroups,n,M2f_,eigvalF_,eigvecF_,B);
        B=Rcpp::as<arma::colvec>(betaActive2["beta"]);
	// B=vectorise(betaActive);			  
	// arma::mat thresh1=arma::abs((betaLast-betaActive)/(arma::ones(n,m)+arma::abs(betaLast)));
	
        threshold=arma::norm(B-BPrev,"inf");
        active=betaActive2("active");
	BPrev=B;

      }
  }
  return(B);
  }

// [[Rcpp::export]]
List GamLoopGLOO(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y, const mat& Z,List jj, List jjfull, List jjcomp, double eps, colvec& YMean2, colvec& ZMean2,int k,int pk,List M2f_, List eigvalF_, List eigvecF_)
{
// arma::colvec YMean2=YMean;
// arma::colvec ZMean2=ZMean;
 int gran2=gamm.size();
 List activefinal(gran2);
 cube beta2(beta_.begin(),k,pk,gran2,false);
 cube betafin(k,pk+1,gran2);
 betafin.fill(0);
 List iterations(gran2);
 mat betaPrev=zeros<mat>(k,pk);
 arma::colvec B=arma::vectorise(betaPrev);
 NumericVector betaF2(k*pk);

 for(int i=0; i<gran2;++i)
    {
        double gam=gamm[i];
	betaPrev=beta2.slice(i);
	B=arma::vectorise(betaPrev);
	List Active = Activeset[i];
	int k2=0;
	int converge=0;
        // mat betaF=zeros(k,k);
	List betaFull(3);
	int thresh=0;
	while(converge==0)
	  {
	    
       	    B = ThreshUpdateOO(Z, gam, Y, eps, Active, jjfull, jjcomp, M2f_,eigvalF_,eigvecF_,B,k);
	 
	     betaFull=BlockUpdate2(Z,gam,Y,eps,jjfull,jjfull,jjcomp,k,M2f_,eigvalF_,eigvecF_,B);
	     betaF2=as<NumericVector>(betaFull("beta"));
	     Active=betaFull("active");
	     converge =betaFull("Converge");
             k2+=1;

	  }
  // arma::mat betaF(betaF2,k,pk,false);	
	mat betaF(betaF2.begin(),k,pk,false);
  colvec nu= YMean2 - betaF *ZMean2;
  betafin.slice(i)=mat(join_horiz(nu, betaF));
    activefinal[i]=Active;
  iterations[i]=k2; 
    }
 List Results=List::create(Named("beta")=betafin,Named("active")=wrap(activefinal),Named("iterations")=iterations);
   return(Results);
    }

// *
// *
// *
// Block Sparse Group Lasso Functions
// *
// *
// *

// Inner while loop from simon et al (2013)
mat createmat(colvec U, int n)
{

  mat u(U.begin(),n,n,false);
  return(u) ;

}

mat sparseWL(const mat& M1a,const  mat& R1, double ngroups, mat& beta,  double t,  double alpha,  double lambda, double eps)
{
  int n=M1a.n_rows,k=M1a.n_cols;
  
int n2=beta.n_rows,k2=beta.n_cols;

double thresh=10;

 mat p=zeros<mat>(n,n);

 colvec STS(n);
 mat thetaOLD=beta;
 mat thetaOLDOLD=beta;
 double l=1;
 mat u=zeros<mat>(n2,n2);
// while(thresh>eps)
//   {
//     beta=thetaOLD+((l-2)/(l+1))*(thetaOLDOLD-thetaOLD);
//     p = (beta*M1a-R1)*trans(M1a)/ngroups;   
//     STS=ST3a(vectorise(beta)-t*vectorise(p),t*alpha*lambda);
//     // DON'T USE TWO NORM IN ARMADILLO, IT IS MATRIX 2 NORM
//     double denom2= norm(STS,"fro");
//     double s31=fmax(1-(t*(1-alpha)*lambda)/denom2,0);
//     STS=s31*STS;
 
//     // mat u(STS.begin(),n2,n2,false);
//     u=createmat(STS,n2);
   
    
//     // beta=thetaOLD+(l/(l+3))*(u-thetaOLD);
 
//     l+=1;
//     mat thresh1=abs(beta-u)/(ones(n2,k2)+abs(u));
//     thresh=norm(thresh1,"inf");
//     thetaOLDOLD=thetaOLD;
//     thetaOLD=u;
//   }

while(thresh>eps)
  {
    p = (beta*M1a-R1)*trans(M1a)/ngroups;
   
    STS=ST3a(vectorise(beta)-t*vectorise(p),t*alpha*lambda);
    // DON'T USE TWO NORM IN ARMADILLO, IT IS MATRIX 2 NORM
    double denom2= norm(STS,"fro");
    double s31=fmax(1-(t*(1-alpha)*lambda)/denom2,0);
    STS=s31*STS;
 
    u= createmat(STS,n2);
   
    
    beta=thetaOLD+(l/(l+3))*(u-thetaOLD);
    l+=1;
    mat thresh1=beta-u;
     thresh=norm(thresh1,"inf");
 
    thetaOLD=u;
  }
  
return(beta);
  }

// Very similar to group lasso case
List blockUpdateSGL(mat& beta,const mat& Z1, double lam, double alpha,const mat& Y1, double eps, List groups, List fullgroups, List compgroups, int k, List M1f, List M2f,NumericVector Eigs)
{

 int n1=groups.size();
  List active(n1);
  int n=beta.n_rows, m=beta.n_cols;

  arma::mat betaPrev=beta;
  int converge=0;
  int count=0;
 
  if(groups.size()==count)
    {
      beta.zeros(n,m);
      active=groups;
    }
  else{
    for(int i=0; i<n1;++i)
 
      {


	NumericVector s1=groups[i];
	IntegerVector s2=fullgroups[i];
	NumericVector scomp=compgroups[i];
	arma::uvec s45=as<arma::uvec>(s1);
	arma::uvec s45F=as<arma::uvec>(s2);


      

	
        if(s1.size()==1){
	  
	  beta.cols(s45F)=arma::zeros(s2.size(),s2.size());
	  active(i)=0;

	}
	if(s1.size()!=1){
	  
	  arma::colvec s3(s1.begin(),s1.size(),false);
	  arma::colvec scomp1(scomp.begin(),scomp.size(),false);

	  arma::uvec s4(s3.size());

	  for(int j=0; j<s3.size(); ++j)
	    {
	      s4(j)=s3(j);
	    }
	  
	  arma::uvec scomp2(scomp1.size());

	  for(int m=0; m<scomp1.size(); ++m)
	    {
	      scomp2(m)=scomp1(m);
	    }
	
	  arma::mat M2a= Z1.rows(scomp2);
	  arma::mat a1= beta.cols(scomp2);
	  arma::mat beta2=beta.cols(s4);

	  arma::mat r=Y1-a1*M2a;

	  arma::mat M1=M1f(i);
	  arma::mat M2=M2f(i);
	
	  arma::mat p=(beta2*M1-r)*trans(M1);
	  colvec STS = ST3a(vectorise(p),alpha*lam);
	    

	  double lamadj=lam*(1-alpha);
	  if(arma::norm(STS,"fro")<=lamadj)
	    {
	      arma::mat astar=arma::zeros(s3.size(),s3.size());
	      active(i)=0;
	      beta.cols(s4)=astar;


	    }
	  else{
	    mat betaS=beta.cols(s4);
	    double t=1/Eigs(i);
	    mat astar2= sparseWL(M1, r, k, betaS, t,  alpha, lam, eps);
	    beta.cols(s4)=astar2;
	    active(i)=s4;  

	  }

	
	


	}

      }



  }
  arma::mat thresh1=arma::abs((betaPrev-beta)/(arma::ones(n,m)+arma::abs(betaPrev)));
  double thresh=arma::norm(thresh1,"inf");
  if(thresh<eps)
    {
      converge=1;

    }
  else{ 
converge=0;
}


  Rcpp::List results=Rcpp::List::create(Named("beta")=wrap(beta),Named("active")=wrap(active),Named("Converge")=wrap(converge));
  return(results);

}

mat ThreshUpdateSGL(mat& betaActive,const mat& Z, double lam,const mat& Y,double eps, List groups, List fullgroups, List compgroups,List M1f,List M2f, NumericVector eigs, double alpha, int k)
  {

    int n=betaActive.n_rows, m=betaActive.n_cols;
  int n1=groups.size();
  mat betaLast=betaActive;
  List active(n1);
  int count=0;
  List betaActive2(3);
  for(int i=0; i<n1; ++i)
    {
      NumericVector g1=groups[i];
      count+=g1.size();

    }
 
  if(groups.size()==count)
    {

      betaActive.zeros(n,m);

      active=groups;
    }
  else{
   double threshold=10*eps;
    while(threshold>eps)
      {
        betaActive2=blockUpdateSGL(betaActive,Z,lam,alpha,Y,eps,groups,fullgroups,compgroups,k, M1f,M2f,eigs);	 
        betaActive=as<mat>(betaActive2("beta"));
	arma::mat thresh1=arma::abs((betaLast-betaActive)/(arma::ones(n,m)+arma::abs(betaLast)));
        threshold=arma::norm(thresh1,"inf");
        active=betaActive2("active");
	betaLast=betaActive;

      }
  }
    return(betaActive);
  }

// [[Rcpp::export]]
List GamLoopSGL(NumericVector beta_, List Activeset,const NumericVector gamm,const double alpha, const mat& Y1, const mat& Z1,List jj,const List jjfull,const List jjcomp,const double eps,const colvec YMean2,const colvec ZMean2,const int k,const int pk,const List M1f_,const List M2f_, const NumericVector eigs_)
{
// arma::colvec YMean2=YMean;
// arma::colvec ZMean2=ZMean;
 int gran2=gamm.size();
 List activefinal(gran2);
  // int n3=Z1.nrow(), m3=Z1.ncol();
  // arma::mat Z(Z1.begin(),n3,m3,false);
  // int n4=Y1.nrow(), m4=Y1.ncol();
  // arma::mat Y(Y1.begin(),n4,m4,false);

 cube beta2(beta_.begin(),k,pk,gran2,false);
 cube betafin(k,pk+1,gran2);
 betafin.fill(0);
 List iterations(gran2);
 mat betaPrev=zeros<mat>(k,pk);

 //INDEX LISTS WITH PARENTHESES NOT BRACKETS :O 
 // WHEN EXTRACTING FROM A LIST YOU NEED as<MAT>
 for(int i=0; i<gran2;++i)
    {
        double gam=gamm[i];
	betaPrev=beta2.slice(i);
	List Active = Activeset[i];
	int k2=0;
	int converge=0;
        mat betaF=zeros(k,pk);
	List betaFull(3);
	//Three components in the list
	int thresh=0;
	while(converge==0)
	  {
	    
       	    betaPrev = ThreshUpdateSGL(betaPrev, Z1, gam, Y1, eps, Active, jjfull, jjcomp, M1f_, M2f_, eigs_, alpha,k);
	 
	    betaFull=blockUpdateSGL(betaPrev,Z1,gam,alpha,Y1,eps,jj,jjfull,jjcomp,k,M1f_,M2f_,eigs_);
	     betaF=as<mat>(betaFull("beta"));

	     Active=betaFull("active");
	     converge =betaFull("Converge");
             k2+=1;

	  }

  colvec nu= YMean2 - betaF *ZMean2;
  betafin.slice(i)=mat(join_horiz(nu, betaF));
    activefinal[i]=Active;
  iterations[i]=k2; 
    }
 List Results=List::create(Named("beta")=betafin,Named("active")=wrap(activefinal),Named("iterations")=iterations);
   return(Results);
    }

// *
// *
// *
// Own-Other Sparse Group Lasso
// *
// *
// *

// [[Rcpp::export]]
mat sparseWLOO(const mat& M1a, const colvec& R1, const int ngroups, colvec& beta, const double t, const double alpha, const double lambda,const double eps,double rho)
{
int n=M1a.n_rows,k=M1a.n_cols;
int n1=R1.n_rows,k1=R1.n_cols;

double thresh=10;

colvec p=beta;

 colvec STS=beta;

 colvec thetaOLD=beta;
 double l=1;
 colvec one=ones<vec>(beta.n_elem);
while(thresh>eps)
  {
    p = trans(M1a)*(M1a*beta-R1)/ngroups;
   
    const colvec& p1=beta-t*vectorise(p);
    STS=ST3ar(p1,rho*t*alpha*lambda);

    double denom2= norm(STS,2);

    double s3=fmax(1-(t*(1-alpha)*lambda*rho)/denom2,0);
    STS=s3*STS;
 
    beta=thetaOLD+(l/(l+3))*(STS-thetaOLD);
 
    l=l+1;
    thresh=max(abs(beta-STS)/(one+abs(STS)));
	// thresh=
    thetaOLD=STS;
  }
 
return(beta);
  }

List blockUpdateSGLOO( colvec& beta,const mat& Z1, double lam, double alpha,const colvec& Y2, double eps, List groups_, const List fullgroups_, List compgroups_, int k,const List M1f_,const List M2f_,const NumericVector Eigs_)
{
 
  int n1=groups_.size();
  List active(n1);
  // int n=beta2.n_rows, m=beta2.n_cols;

  // colvec beta=vectorise(beta2);

  colvec betaPrev=beta;

  int converge=0;
  int count=0;
  colvec one=ones<vec>(beta.n_elem);
  // arma:: colvec Y2=arma::vectorise(Y1,0);

  if(groups_.size()==count)
    {
      // beta.zeros(n,m);
      beta.zeros();
      active=groups_;
    }
  else{
    for(int i=0; i<n1;++i)
 
      {


	// NumericVector s1=groups_[i];
	// IntegerVector s2=fullgroups_[i];
	// NumericVector scomp=compgroups_[i];
	// arma::uvec s45=as<arma::uvec>(s1);
	// arma::uvec s45F=as<arma::uvec>(s2);


	// NumericVector s1=groups_[i];
	// IntegerVector s2=fullgroups_[i];
	// NumericVector scomp=compgroups_[i];
	arma::uvec s4=as<arma::uvec>(groups_[i]);
	arma::uvec s45F=as<arma::uvec>(fullgroups_[i]);
	uvec scomp2=as<uvec>(compgroups_[i]);
      

	
        // if(s1.size()==1){
	    if(s4.n_elem==1){
	  
	  // beta.elem(s45F)=arma::zeros(s2.size());
			  beta.elem(s45F)=arma::zeros(s45F.n_elem);
	
	  active(i)=0;

	}
	if(s4.n_elem!=1){
	  
	  // arma::colvec s3(s1.begin(),s1.size(),false);
	  // arma::colvec scomp1(scomp.begin(),scomp.size(),false);

	  // arma::uvec s4(s3.size());

	  // for(int j=0; j<s3.size(); ++j)
	  //   {
	  //     s4(j)=s3(j);
	  //   }
	  
	  // arma::uvec scomp2(scomp1.size());

	  // for(int m=0; m<scomp1.size(); ++m)
	  //   {
	  //     scomp2(m)=scomp1(m);
	  //   }
	
	  const arma::mat& M2a= Z1.cols(scomp2);
	  const arma::colvec& a1= beta.elem(scomp2);
	  const arma::colvec& beta2=beta.elem(s4);
	  
	  const arma::colvec& r=Y2-M2a*a1;
	  // colvec r=Y2-Z1.cols(scomp2)*beta.elem(scomp2);
	  
	  const arma::mat& M1=M1f_(i);
	  const arma::mat& M2=M2f_(i);

 	  const arma::mat& p=-trans(M1)*(r-M1*beta2);
	  // arma::mat p=-trans(as<mat>(M1f_(i)))*(r-M1*beta2);


	  double rho=sqrt(s4.n_elem);
	  const colvec& STS = ST3a(vectorise(p),rho*lam*alpha);

	  double lamadj=lam*(1-alpha)*rho;
	  if(arma::norm(STS,"fro")<=lamadj)
	    {
	      arma::colvec astar=arma::zeros(s4.n_elem);	      
	      active(i)=0;
	      beta.elem(s4)=astar;


	    }
	  else{

	    colvec betaS=beta.elem(s4);
	   const double t=1/Eigs_(i);

	   const  mat astar2= sparseWLOO(M1, r,2* k, betaS, t,  alpha, lam, eps,rho);
	    
	    beta.elem(s4)=astar2;
	    active(i)=s4;  

	  }

	
	


	}

      }



  }
  double thresh=max(abs(beta-betaPrev)/(one+abs(betaPrev)));
  if(thresh<eps)
    {
      converge=1;

    }
  else{ 
converge=0;
}
  // mat beta3(beta.begin(),n,m,false);

  Rcpp::List results=Rcpp::List::create(Named("beta")=wrap(beta),Named("active")=wrap(active),Named("Converge")=wrap(converge));
  return(results);

}


mat ThreshUpdateSGLOO(colvec& betaActive,const mat& Z,const double lam,const colvec& Y,const double eps, List groups_, const List fullgroups_,const List compgroups_,const List M1f_,const List M2f_,const NumericVector eigs_,const double alpha,const int k)
  {
    
  int kp = betaActive.n_elem;
    // int n=betaActive.n_rows, m=betaActive.n_cols;
  int n1=groups_.size();
  colvec betaLast=betaActive;
  List active(n1);
  int count=0;
  List betaActive2(3);
  for(int i=0; i<n1; ++i)
    {
      NumericVector g1=groups_[i];
      count+=g1.size();

    }
 
  if(groups_.size()==count)
    {

      betaActive.zeros();

      active=groups_;
    }
  else{
	int converge=0;
   // double threshold=10*eps;
      while(converge==0)	
      {
        betaActive2=blockUpdateSGLOO(betaActive,Z,lam,alpha,Y,eps,groups_,fullgroups_,compgroups_,k, M1f_,M2f_,eigs_);	 
        betaActive=as<colvec>(betaActive2("beta"));
        converge =betaActive2("Converge");

      }
  }
    return(betaActive);
  }

// [[Rcpp::export]]
List GamLoopSGLOO(NumericVector beta_,const List Activeset_,const NumericVector gamm,const double alpha,const mat& Y,const mat& Z,List jj_,const List jjfull_, List jjcomp_,const double eps,const colvec& YMean2,const colvec& ZMean2,const int k,const int pk,const List M1f_,const List M2f_,const NumericVector eigs_)
{
 int gran2=gamm.size();
 List activefinal(gran2);

 cube beta2(beta_.begin(),k,pk,gran2,false);
 cube betafin(k,pk+1,gran2);
 betafin.fill(0);
 List iterations(gran2);
 mat betaPrev=zeros<mat>(k,pk);
 NumericVector betaF2(k*pk);
const arma:: colvec& Y2=arma::vectorise(Y,0);

 for(int i=0; i<gran2;++i)
    {
        double gam=gamm[i];
	betaPrev=beta2.slice(i);
	List Active = Activeset_[i];
	int k2=0;
	int converge=0;
        // mat betaF=zeros(k,pk);
	colvec B=vectorise(betaPrev);
	List betaFull(3);
	//Three components in the list
	while(converge==0)
	  {
	    
       	    B = ThreshUpdateSGLOO(B, Z, gam, Y2, eps, Active, jjfull_, jjcomp_, M1f_, M2f_, eigs_, alpha,k);
	 
	    betaFull=blockUpdateSGLOO(B,Z,gam,alpha,Y2,eps,jjfull_,jjfull_,jjcomp_,k,M1f_,M2f_,eigs_);
	    betaF2=as<NumericVector>(betaFull("beta"));

	     Active=betaFull("active");
	     converge =betaFull("Converge");
             k2+=1;

	  }
  mat betaF(betaF2.begin(),k,pk,false);
  colvec nu= YMean2 - betaF *ZMean2;
  betafin.slice(i)=mat(join_horiz(nu, betaF));
    activefinal[i]=Active;
  iterations[i]=k2; 
    }
 List Results=List::create(Named("beta")=betafin,Named("active")=wrap(activefinal),Named("iterations")=iterations);
   return(Results);
    }


// *
// *
// *
// *
// // Power Rule Algorithm for eigenvalue generation
// *
// *
// *
// *



// [[Rcpp::export]]

List powermethod(NumericMatrix A1, NumericVector x1) {
   double dd = 1.0;
   int nm=A1.nrow(), km=A1.ncol();
arma::mat A(A1.begin(),nm,km,false);
int nn=x1.length();
arma::mat x(x1.begin(),nn,1,false);
   double eps = .001;
   arma::mat y=x;
double alpha=0;
double theta=0;
   while (dd > eps*fabs(theta)) {
    x=y/arma::norm(y,"fro");
    y = A*x;
    theta=as_scalar(trans(x)*y);
     dd=arma::norm(y-theta*x,2);


   }
double lambda=theta;
Rcpp::List results = Rcpp::List::create(Named("lambda")=as<double>(wrap(lambda)),Named("q1")=as<NumericVector>(wrap(x)));
return( results);
}


// Useful 2-norm function

// [[Rcpp::export]]
double norm2(NumericVector x){
       arma::vec xx = x;
       double g=arma::norm(xx,2);
       return (as<double>(wrap(g)));
}

