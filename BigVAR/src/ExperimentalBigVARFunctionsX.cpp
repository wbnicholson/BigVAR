
#include <RcppArmadillo.h>
#include <math.h> 
#include <vector>
#include <limits>
#include <algorithm>
#include <string>
#include <numeric>      // std::iota

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]





using namespace Rcpp;
using namespace arma;

// Soft thresholding
// [[Rcpp::export]]
double ST1a(double z,double gam){
//double sparse=0;
if(z>0 && gam<fabs(z)) return(z-gam);

if(z<0 && gam<fabs(z)) return(z+gam);
if(gam>=fabs(z)) return(0);
else return(0);
}

// Columnwise softthresholding
// [[Rcpp::export]]
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

uvec ind(int n2,int m){
  std::vector<int> subs;
for(int i =0 ; i<n2;++i)
{
  subs.push_back(i);
//	Rcpp::Rcout << "R" << subs[i] << std::endl;
 }
subs.erase(subs.begin()+m);
 // uvec subs2=conv_to<uvec>::from(subs);
return(conv_to<uvec>::from(subs));
}


// Rowwise softthresholding
// [[Rcpp::export]]
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
  colvec B1=B.col(0);
 
  double j =1;

  for( int i =0; i<k; ++i)

    {
      B1=B.col(i);
      colvec BOLD=B.col(i);
      colvec BOLDOLD=BOLD;
     double thresh=10*eps;
      j=1;
      while(thresh>eps)
	{
	 colvec v=BOLD+((j-2)/(j+1))*(BOLD-BOLDOLD);

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


//Penalty Loop For FISTA

// [[Rcpp::export]]
cube gamloopFista(NumericVector beta_, const mat& Y,const mat& Z,const  colvec gammgrid, const double eps,const colvec& YMean2, const colvec& ZMean2,mat& B1, int k, int p,double tk, int k1,int s){

 mat b2=B1;
 mat B1F2=B1;
 const int ngridpts=gammgrid.size();
 cube bcube(beta_.begin(),k1,k1*p+(k-k1)*s,ngridpts,false);
 cube bcube2(k1,k1*p+(k-k1)*s+1,ngridpts);
 bcube2.fill(0);
 
colvec nu=zeros<colvec>(k1);
double gam =0;

 int i;
  for (i=0; i<ngridpts;++i) {
            gam=gammgrid[i];

	mat B1F2=bcube.slice(i);
	B1 = FistaLV(Y,Z,B1F2,gam,eps,tk,k1,p); 
	

  
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
inline double trust32(int k,const arma::mat& PP,double delta, double lambda,const arma::vec& EigVA,const arma::mat& EigVector){
         double g=0;
         for(int i = 0; i < k; ++i)
         { 
    
	   g+=pow(arma::as_scalar(trans(EigVector.col(i))*PP),2)/pow(EigVA[i]*delta+lambda,2);

         }
	 return(g);
}

inline double fprime2(int k,const arma::mat& PP,double delta, double lambda, const arma::vec& EigVA,const arma::mat& EigVE)
{
         double gg2=0;
         for(int i = 0; i < k; ++i)
         { 
    
	   gg2+=(pow(arma::as_scalar(trans(EigVE.col(i))*PP),2)*EigVA[i])/pow(EigVA[i]*delta+lambda,3);


         }
	 double c1=trust32(k, PP, delta, lambda, EigVA, EigVE);
	 double res=-.5*pow(c1,-1.5)*-2*gg2;
	 return(res);

}

inline double Newton2(int k,const  arma::mat& P, double lambda, const arma::vec& EigVA,const arma::mat& EigVE){
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
List Eigencomp( mat& Z1, List groups,int n1,int k1) 
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
      arma::mat D(k1,k1);
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
List EigencompOO( mat& ZZ1, List groups,int n1,int k) 
{
  List M2f(n1);
  List eigvecF(n1);
  List eigvalF(n1);
  int count=0;
  for(int i=0; i<n1; ++i)
    {
//      NumericVector g1=groups[i];
      arma::uvec s4=as<arma::uvec>(groups[i]);
      count+=s4.n_elem;
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
List BlockUpdateGL(mat& beta,const mat& Z1, double lam, const mat& Y1,double eps, List groups, List fullgroups, List compgroups, int k,List M3f_,List eigvalF_, List eigvecF_,int k1){
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
	arma::uvec s45=as<arma::uvec>(groups[i]);
	arma::uvec s45F=as<arma::uvec>(fullgroups[i]);


      

	
	if(max(s45)==0){
	  
	  beta.cols(s45F)=arma::zeros(k1,s45F.n_elem);
	 
	  active(i)=0;

	}
	if(max(s45)!=0){
	  
	   arma::uvec scomp1=as<arma::uvec>(compgroups[i]);


          mat r =beta.cols(scomp1)*Z1.rows(scomp1)-Y1 ;

	  colvec p=vectorise((r)*trans(Z1.rows(s45)));
	  double adjlam=sqrt(s45.n_elem)*lam;

	  if(arma::norm(p,"fro")<=adjlam)
	    {
	      arma::mat astar=arma::zeros(k1,s45.n_elem);
	      active(i)=0;


	    }
	  else{
	    arma::mat M3=M3f_(i);
	    int k1a=M3.n_cols;
	    double deltfin=  Newton2(k1a,p,adjlam,eigvalF_(i),eigvecF_(i));
	
		M3.diag()+=adjlam/deltfin;
	arma::mat astar=-solve(M3,p);
	astar.set_size(k1,s45.n_elem);
	beta.cols(s45)=astar;
	
	 active(i)=s45;  

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


  Rcpp::List results=Rcpp::List::create(Named("beta")=beta,Named("active")=wrap(active),Named("Converge")=wrap(converge));
  return(results);

}



// [[Rcpp::export]]
mat ThreshUpdate(mat& betaActive,const mat& Z1, double lam, const mat& Y1,double eps, List groups, List fullgroups, List compgroups,List M2f_, List eigvalF_, List eigvecF_,int k1)
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
      count+=max(g1);

    }
 
  if(count==0)
    {

      betaActive.zeros(n,m);

      active=groups;
    }
  else{
   double threshold=10*eps;
    while(threshold>eps)
      {
		  betaActive2=BlockUpdateGL(betaActive,Z1,lam,Y1,eps,groups,fullgroups,compgroups,n,M2f_,eigvalF_,eigvecF_,k1);	 
        betaActive=as<mat>(betaActive2("beta"));
	arma::mat thresh1=arma::abs((betaLast-betaActive)/(arma::ones(n,m)+arma::abs(betaLast)));
        threshold=arma::norm(thresh1,"inf");
        active=betaActive2("active");
	betaLast=betaActive;

      }
  }
    return(betaActive);
  }

List GamLoopGL(cube& betafin, List Activeset, colvec gamm, const mat& Y1, const mat& Z1,List jj, List jjfull, List jjcomp, double eps,const colvec& YMean2, const colvec&  ZMean2,int k,int pk,const List M2f_, const List eigvalF_, const List eigvecF_,int k1,bool MN)
{

    mat C(k1,pk);
    C.zeros();

  if(MN){
    // mat C(k1,k1*p+(k-k1)*s);
    // C.zeros();
    C.diag(1);
    // Y=trans(Y);
    // Y=Y-C*Z;
    // Y=trans(Y);
  }


 int gran2=gamm.size();
 List activefinal(gran2);
 cube beta2=betafin.subcube(0,1,0,k1-1,pk,gran2-1);

 List iterations(gran2);
 mat betaPrev=zeros<mat>(k1,pk);
 
//int n3=Z1.n_rows;// m3=Z1.n_cols;
//int n4=Y1.n_rows;// m4=Y1.n_rows;
 
 //INDEX LISTS WITH PARENTHESES NOT BRACKETS
 // WHEN EXTRACTING FROM A LIST YOU NEED as<MAT>

 for(int i=0; i<gran2;++i)
    {
        double gam=gamm[i];
	betaPrev=beta2.slice(i);
	List Active = Activeset[i];
	int k2=0;
	int converge=0;
        mat betaF=zeros(k1,k);
	List betaFull(3);
	while(converge==0)
	  {
	    
		  betaPrev = ThreshUpdate(betaPrev, Z1, gam, Y1, eps, Active, jjfull, jjcomp,M2f_,eigvalF_,eigvecF_,k1);
	 
		  betaFull=BlockUpdateGL(betaPrev,Z1,gam,Y1,eps,jjfull,jjfull,jjcomp,k,M2f_,eigvalF_,eigvecF_,k1);
	     betaF=as<mat>(betaFull("beta"));
	     Active=betaFull("active");
	     converge =betaFull("Converge");
             k2+=1;

	  }

	if(MN)
	  {
	    for(int i=0;i<gran2;++i)
	      {
		betaF=betaF+C;
	      }
	  }


  colvec nu= YMean2 - betaF *ZMean2;
  betafin.slice(i)=mat(join_horiz(nu, betaF));
    activefinal[i]=Active;
  iterations[i]=k2; 
    }
 List Results=List::create(Named("beta")=wrap(betafin),Named("active")=wrap(activefinal),Named("iterations")=iterations);
	 return(Results);
    }

// [[Rcpp::export]]
List GamLoopGL2(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y1, const mat& Z1,List jj, List jjfull, List jjcomp, double eps,const colvec& YMean2, const colvec&  ZMean2,int k,int pk,const List M2f_, const List eigvalF_, const List eigvecF_)
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
 
//int n3=Z1.n_rows, m3=Z1.n_cols;
 // arma::mat Z1(Z.begin(),n3,m3,false);
//int n4=Y1.n_rows, m4=Y1.n_rows;
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
//	int thresh=0;
	while(converge==0)
	  {
	    
       	    betaPrev = ThreshUpdate(betaPrev, Z1, gam, Y1, eps, Active, jjfull, jjcomp,M2f_,eigvalF_,eigvecF_,k);
	 
	    betaFull=BlockUpdateGL(betaPrev,Z1,gam,Y1,eps,jjfull,jjfull,jjcomp,k,M2f_,eigvalF_,eigvecF_,k);
	     betaF=as<mat>(betaFull("beta"));
	     Active=betaFull("active");
	     converge =betaFull("Converge");
             k2+=1;

	  }
	// if(MN)
	//   {
	//     cube beta2f=as<cube>(beta3("beta"));
	//     for(int i=0;i<gran2;++i)
	//       {
	// 	beta2f.slice(i)=beta2f.slice(i)+C;
	//       }
	//     beta3("beta")=beta2f;
	//   }


  colvec nu= YMean2 - betaF *ZMean2;
  betafin.slice(i)=mat(join_horiz(nu, betaF));
    activefinal[i]=Active;
  iterations[i]=k2; 
    }



 List Results=List::create(Named("beta")=betafin,Named("active")=wrap(activefinal),Named("iterations")=iterations);
   return(Results);
    }

// [[Rcpp::export]]
List GroupLassoVAR(NumericVector beta_, mat Y, mat Z, const colvec& gamm, double eps, int k, int p,List INIActive_, List jj,List Compgroups, int k1,int s,bool MN)
{
    mat C(k1,k1*p+(k-k1)*s);
    C.zeros();

  if(MN){
    // mat C(k1,k1*p+(k-k1)*s);
    // C.zeros();
    C.diag(1);
    Y=trans(Y);
    Y=Y-C*Z;
    Y=trans(Y);
  }
  
	Y=trans(Y);
	colvec YMean=arma::mean(Y,1);

    colvec ZMean=mean(Z,1);
	colvec zones=ones(Y.n_cols,1);
	colvec yones=ones(Z.n_cols,1);
	List beta3;
	// if(LG==false){
	Y=Y-YMean*trans(yones);
	Z=Z-ZMean*trans(zones);
	List jjfull = jj;
	List Eigsys_=Eigencomp(Z,jj,jj.size(),k1);
	List M2_= Eigsys_("M3");
	List eigvals_= Eigsys_("eigval");
	List eigvecs_= Eigsys_("eigvec");
	int gran2=gamm.n_elem;
	cube beta2(beta_.begin(),k1,p*k1+s*(k-k1)+1,gran2,false);
	beta3=	GamLoopGL(beta2,INIActive_,gamm,Y,Z,jj,jjfull,Compgroups,eps,YMean,ZMean,k,p*k1+s*(k-k1),M2_,eigvals_,eigvecs_,k1,MN);
	// }
	// else{
	// mat Y1=Y-YMean*trans(yones);
	// mat Z1=Z-ZMean*trans(zones);
	// // List jj = Groups;
	// List jjfull = jj;
	// List Eigsys_=Eigencomp(Z,jj,jj.size(),k1);
	// // Rcout<<"Foo"<<std::endl
	// List M2_= Eigsys_("M3");
	// List eigvals_= Eigsys_("eigval");
	// List eigvecs_= Eigsys_("eigvec");
	// int gran2=gamm.n_elem;
    // cube beta2(beta_.begin(),k1,p*k1+s*(k-k1)+1,gran2,false);
	// beta3=	GamLoopGL(beta2,INIActive_,gamm,Y1,Z1,jj,jjfull,Compgroups,eps,YMean,ZMean,k,p*k,M2_,eigvals_,eigvecs_,k1);
	// }

	// if(MN)
	//   {
	//     cube beta2f=as<cube>(beta3("beta"));
	//     for(int i=0;i<gran2;++i)
	//       {
	// 	beta2f.slice(i)=beta2f.slice(i)+C;
	//       }
	//     beta3("beta")=beta2f;
	//   }

	return(beta3);

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
List BlockUpdate2(const mat& ZZ1, double lam,const mat& Y1,double eps, List groups, List fullgroups, List compgroups, int k, List M2f_, List eigvalF_, List eigvecF_,colvec& B,int k1){
  int n1=groups.size();
  List active(n1);
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


      

	
        if(max(s1)==0){
	  
	  B.elem(s45F)=arma::zeros(s2.size());
	  active(i)=0;

	}
	if(max(s1)!=0){
	  
//	  arma::colvec s3(s1.begin(),s1.size(),false);
	  // Need index vectors to be uvecs, also need to make sure vecs start indexing at zero!
	  arma::uvec s4=as<arma::uvec>(s1);

	  
	  arma::uvec scomp2=as<arma::uvec>(scomp);

	  arma::mat M2a= ZZ1.cols(scomp2);
	  arma::colvec a1= B.elem(scomp2);
	  mat foo1=(M2a*a1);
	   
	  arma::colvec r=M2a*a1-arma::vectorise(Y1);

	  arma::mat M1=ZZ1.cols(s4);
	  arma::mat M2=M2f_(i);
	  arma::vec eigval=eigvalF_(i);
	  arma::mat eigvec=eigvecF_(i);
	  arma::mat p=trans(M1)*r;

	  double rho=sqrt(s1.size());
	  double adjlam=rho*lam;

	  if(arma::norm(p,2)<=adjlam)
	    {
	      arma::colvec astar=arma::zeros(s1.size());
	      active(i)=0;


	    }
	  else{
        int k1a=M2.n_cols;
    
	double deltfin=  Newton2(k1a,p,adjlam,eigval,eigvec);
	
	arma::mat D1(s1.size(),s1.size());
	D1.eye();
	arma::mat astar=-solve(M2+adjlam/deltfin*D1,p);

	 B.elem(s4)=astar;
	
	 active(i)=s4;  

	  }

	
	


	}

      }



  }
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
colvec ThreshUpdateOO(const mat& ZZ, double lam,const mat& Y,double eps, List groups, List fullgroups, List compgroups,List M2f_, List eigvalF_, List eigvecF_,colvec& B,int n, int k1)
  {

  int kp=B.n_elem;
  int n1=groups.size();
  colvec BPrev=B;
  List active(n1);
//  int count=0;
 List betaActive2(3);

  // for(int i=0; i<n1; ++i)
  //   {
  //     NumericVector g1=groups[i];
  //     count+=g1.size();

  //   }
 
  if(max(groups)==0)
    {

      B.zeros(kp);

      active=groups;
    }
  else{
   double threshold=10*eps;
    while(threshold>eps)
      {
		  betaActive2=BlockUpdate2(ZZ,lam,Y,eps,groups,fullgroups,compgroups,n,M2f_,eigvalF_,eigvecF_,B,k1);
        B=Rcpp::as<arma::colvec>(betaActive2["beta"]);
	
        threshold=arma::norm(B-BPrev,"inf");
        active=betaActive2("active");
	BPrev=B;

      }
  }
  return(B);
  }

// [[Rcpp::export]]
List GamLoopGLOO(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y, const mat& Z,List jj, List jjfull, List jjcomp, double eps, colvec& YMean2, colvec& ZMean2,int k,int pk,List M2f_, List eigvalF_, List eigvecF_,int k1)
{
 int gran2=gamm.size();
 List activefinal(gran2);
 cube beta2(beta_.begin(),k1,pk,gran2,false);
 cube betafin(k1,pk+1,gran2);
 betafin.fill(0);
 List iterations(gran2);
 mat betaPrev=zeros<mat>(k1,pk);

 arma::colvec B=arma::vectorise(betaPrev);
 // B.print();
 NumericVector betaF2(k1*pk);

 for(int i=0; i<gran2;++i)
    {
        double gam=gamm[i];
	betaPrev=beta2.slice(i);
	B=arma::vectorise(betaPrev);
	List Active = Activeset[i];
	int k2=0;
	int converge=0;
	List betaFull(3);
	while(converge==0)
	  {
	    
		 B = ThreshUpdateOO(Z, gam, Y, eps, Active, jjfull, jjcomp, M2f_,eigvalF_,eigvecF_,B,k1,k1);
		 betaFull=BlockUpdate2(Z,gam,Y,eps,jjfull,jjfull,jjcomp,k,M2f_,eigvalF_,eigvecF_,B,k1);
	     betaF2=as<NumericVector>(betaFull("beta"));
	     Active=betaFull("active");
	     converge =betaFull("Converge");
         k2+=1;

	  }
	// Rcout<<betaF2<<std::endl;
	mat betaF(betaF2.begin(),k1,pk,false);
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


// *
// *
// *
// Own-Other Sparse Group Lasso
// *
// *
// *

mat sparseWLOO(const mat& M1a, const colvec& R1, const int ngroups, colvec& beta, const double t, const double alpha, const double lambda,const double eps,double rho)
{
//int n=M1a.n_rows,k=M1a.n_cols;
//int n1=R1.n_rows,k1=R1.n_cols;

 double thresh=10;

 colvec p=beta;

 colvec STS=beta;

 colvec thetaOLD=beta;
 double l=1;
 colvec one=ones<vec>(beta.n_elem);
while(thresh>eps)
  {
    p = trans(M1a)*(M1a*beta-R1)/ngroups;
   
    const colvec p1=beta-t*vectorise(p);
    STS=ST3ar(p1,rho*t*alpha*lambda);

    double denom2= norm(STS,2);

    double s3=fmax(1-(t*(1-alpha)*lambda*rho)/denom2,0);
    STS=s3*STS;
 
    beta=thetaOLD+(l/(l+3))*(STS-thetaOLD);
 
    l=l+1;
	//Relative thresholds seem to help with computation time
    thresh=max(abs(beta-STS)/(one+abs(STS)));
	// Rcout<<thresh<<std::endl;
	// thresh=max(abs(beta-STS));
	// thresh=
    thetaOLD=STS;
  }
 
return(beta);
  }

List blockUpdateSGLOO( colvec& beta,const mat& Z1, double lam, double alpha,const colvec& Y2, double eps, List groups_, const List fullgroups_, List compgroups_,const List M2f_,const NumericVector Eigs_,double k1)
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
	if(max(s4)==0){
	  
	  // beta.elem(s45F)=arma::zeros(s2.size());
			  beta.elem(s45F)=arma::zeros(s45F.n_elem);
	
	  active(i)=0;

	}
	if(max(s4)!=0){

	
	  const arma::mat& M2a= Z1.cols(scomp2);
	  const arma::colvec& a1= beta.elem(scomp2);
	  const arma::colvec& beta2=beta.elem(s4);
	  
	  const arma::colvec& r=Y2-M2a*a1;
	  // colvec r=Y2-Z1.cols(scomp2)*beta.elem(scomp2);
	  
	  // const arma::mat& M1=M1f_(i);
	  const mat& M1=Z1.cols(s4);
	//  const arma::mat& M2=M2f_(i);

 	  const arma::mat& p=-trans(M1)*(r-M1*beta2);
	  // arma::mat p=-trans(as<mat>(M1f_(i)))*(r-M1*beta2);


	  double rho=sqrt(s4.n_elem);
	  const colvec& STS = ST3a(vectorise(p),lam*alpha);

	  double lamadj=lam*(1-alpha)*rho;
	  if(arma::norm(STS,"fro")<=lamadj)
	    {
	      arma::colvec astar=arma::zeros(s4.n_elem);	      
	      active(i)=0;
	      beta.elem(s4)=astar;


	    }
	  else{

	    colvec betaS=beta.elem(s4);
	    // Rcout<<betaS.n_elem<<std::endl;
		// Rcout<<"Foo"<<std::endl;
	   const double t=1/Eigs_(i);

	   const  mat astar2= sparseWLOO(M1, r,k1, betaS, t,  alpha, lam, eps,rho);

	    beta.elem(s4)=astar2;
	    active(i)=s4;  

	  }

	
	


	}

      }



  }
  double thresh=max(abs(beta-betaPrev)/(one+abs(betaPrev)));
  // Rcout<< thresh<<std::endl;
  // double thresh=max(abs(beta-betaPrev));
  if(thresh<eps)
    {
		// Rcout<<"CONVERGE"<<std::endl;
      converge=1;

    }
  else{ 
converge=0;
}
  // mat beta3(beta.begin(),n,m,false);

  Rcpp::List results=Rcpp::List::create(Named("beta")=wrap(beta),Named("active")=wrap(active),Named("Converge")=wrap(converge));
  return(results);

}


mat ThreshUpdateSGLOO(colvec& betaActive,const mat& Z,const double lam,const colvec& Y,const double eps, List groups_, const List fullgroups_,const List compgroups_,const List M2f_,const NumericVector eigs_,const double alpha,double k1)
  {
    
 // int kp = betaActive.n_elem;
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
		  betaActive2=blockUpdateSGLOO(betaActive,Z,lam,alpha,Y,eps,groups_,fullgroups_,compgroups_,M2f_,eigs_,k1);	 
        betaActive=as<colvec>(betaActive2("beta"));
        converge =betaActive2("Converge");

      }
  }
    return(betaActive);
  }

// [[Rcpp::export]]
List GamLoopSGLOO(NumericVector beta_,const List Activeset_,const NumericVector gamm,const double alpha,const mat& Y,const mat& Z,List jj_,const List jjfull_, List jjcomp_,const double eps,const colvec& YMean2,const colvec& ZMean2,const int k1,const int pk,const List M2f_,const NumericVector eigs_)
{

 int gran2=gamm.size();
 List activefinal(gran2);
 cube beta2(beta_.begin(),k1,pk,gran2,false);
 cube betafin(k1,pk+1,gran2);
 betafin.fill(0);
 List iterations(gran2);
 mat betaPrev=zeros<mat>(k1,pk);
 NumericVector betaF2(k1*pk);
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
 // Rcout<<"TEST@"<<std::endl;
	    
		  B = ThreshUpdateSGLOO(B, Z, gam, Y2, eps, Active, jjfull_, jjcomp_, M2f_, eigs_, alpha,(double) k1);
 // Rcout<<"TEST"<<std::endl;
	 
		  betaFull=blockUpdateSGLOO(B,Z,gam,alpha,Y2,eps,jjfull_,jjfull_,jjcomp_,M2f_,eigs_,(double) k1);
	    betaF2=as<NumericVector>(betaFull("beta"));

	     Active=betaFull("active");
	     converge =betaFull("Converge");
             k2+=1;

	  }
  mat betaF(betaF2.begin(),k1,pk,false);
  // betaF.print();
  colvec nu= YMean2 - betaF *ZMean2;
  betafin.slice(i)=mat(join_horiz(nu, betaF));
    activefinal[i]=Active;
  iterations[i]=k2; 
    }
 // betafin.print();
 List Results=List::create(Named("beta")=betafin,Named("active")=wrap(activefinal),Named("iterations")=iterations);
   return(Results);
    }

// *
// *
// *
// *
// Block HVAR
// *
// *
// *
// *


uvec vs2(int p,int k,int j)
{
  uvec vs(p*k-(j-1)*k);
  for(int i=(k)*(j-1);i<(p*k);++i)
    {
      vs(i-(k)*(j-1))=i;
    }
	return(vs);
}

colvec zeropad(colvec& r, int L,int k)
{
  colvec U=zeros(k*L);
  int n=r.n_elem;
  // for(int i=(L-1);i>-1;--i)
  //   {
      for(int j=0;j<n;++j)
	{
	  U(k*L-j-1)=r(n-j-1);
	}
    //   if(i==n)
    // 	{break;}
    // }

  return(U);
}


// [[Rcpp::export]]
rowvec proxcpp(colvec v2,int L,double lambda,int k,colvec w)
{

 colvec r=v2;

      for(int q=(L-1); q>=0;--q)

	{

	  
	    
       	      uvec res=vs2(L,k,q+1);
	    
	    
			  if(norm(r(res)/(lambda*w(q)),"fro")<1+pow(10,-8))
	    {
	  		  r(res)=zeros(res.n_elem);
	
	      }
	  else{
		  r(res)=r(res)-lambda*w(q)*r(res)/norm(r(res),"fro");
	  
	  }

	
	
	    }
      return(trans(r));


}

// [[Rcpp::export]]
mat Fistapar(const mat Y,const mat Z,const mat phi, const int L,const double lambda,const List vsubs_,const double eps,const double tk,const int k)
{
  mat phiFIN=phi;



  colvec w(L);
  for(int r=0;r<L;++r)

    {
      w(r)=sqrt(k);

    }
  const colvec w2=w;



 rowvec phiR=phi.row(0);
 rowvec v=phiR;
 rowvec phiOLD=phiR;
 rowvec phiOLDOLD=phiR;
 int i;
  for( i=0;i<k;++i)
    {
  phiR=phi.row(i);
  phiOLD=zeros(1,k*L);
  phiOLDOLD=phiOLD;
   double thresh=10*eps;
   double   j=1;
  while(thresh>eps)
    {

      v=phiOLD+((j-2)/(j+1))*(phiOLD-phiOLDOLD);
      phiR=proxcpp(vectorise(v)+tk*vectorise((trans(Y.col(i))-v*Z)*trans(Z)),L,tk*lambda,k,w2);
      thresh=max(abs(phiR-v));

      phiOLDOLD=phiOLD;
      phiOLD=phiR;
      j+=1;
  
 }
   phiFIN.row(i)=phiR;

    } 

  return(phiFIN);
}

// [[Rcpp::export]]
cube gamloopHVAR(NumericVector beta_, const mat& Y,const mat& Z, colvec gammgrid, const double eps,const colvec& YMean2, const colvec& ZMean2,mat& B1, const int k, const int p,List vsubs_){


vec eigval;
mat eigvec;
 const mat& Zt=Z*trans(Z);
eig_sym(eigval, eigvec, Zt);

 double tk=1/max(eigval); 

 const int ngridpts=gammgrid.n_elem;
 cube bcube(beta_.begin(),k,k*p,ngridpts,false);
 cube bcube2(k,k*p+1,ngridpts);
 bcube2.fill(0);
colvec nu=zeros<colvec>(k);

 int i;

  for (i=0; i<ngridpts;++i) {
            

	 B1=bcube.slice(i);
	 B1 = Fistapar(Y,Z,B1,p,gammgrid[i],vsubs_,eps,tk,k); 
	  
	 nu = YMean2 - B1 *ZMean2;
         bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}

    return(bcube2);
}

// *
// *
// *
// *
// *
// *
// Own Other HVAR
// *
// *
// *
// *
// *
// *
// *
// *

colvec zeropadOO(colvec& r, int& L,int& k,uvec& vs)
{
  colvec U=zeros(k*L/2);
  int n=r.n_elem;
  int n2=vs.n_elem;

      for(int j=0;j<n2;++j)
	{
	  U(vs(n2-j-1))=r(n-j-1);
	}
    

  return(U);
}

// [[Rcpp::export]]
rowvec proxcppOO(colvec v2,int L,double lambda,List vsubs,int k,colvec w)
{

	colvec r=v2;
      for(int i=(L-1); i>=0;--i)

	{

	  uvec res=as<uvec>(vsubs(i));

	  if(norm(r(res)/(lambda*w(i)),"fro")<1+pow(10,-8))
	    {
			r(res)=zeros(res.n_elem);
	    }
	  else{
		  r(res)=r(res)-lambda*w(i)*r(res)/(norm(r(res),"fro"));
	  }

	}
  

      return(trans(r));
}
// [[Rcpp::export]]
mat FistaOO(const mat Y, const mat Z, mat phi, const int p, const int k, double lambda, List groups_, const double eps, const double tk,colvec w)
{

  double j=1;
  mat phiFin=phi;
  rowvec phiR=phi.row(0); 
  rowvec phiOLD=phiR;
  rowvec phiOLDOLD=phiOLD;
  rowvec v=phiOLD;
  //changed 5-13 used to be res1=ind(p,1)
  uvec res1=ind(p,0);

  double thresh=10;
  for(int i=0;i<k;++i)
    {
   j=1;
   thresh=10*eps;
   phiR=phi.row(i);
   phiOLD=phiR;
   phiOLDOLD=phiOLD;
   v=phiR;
   List vsubs=groups_[i];
  while(thresh>eps)
    {
       v=phiOLD+((j-2)/(j+1))*(phiOLD-phiOLDOLD);

       phiR=proxcppOO(vectorise(v)+tk*vectorise((trans(Y.col(i))-v*Z)*trans(Z)),2*p,tk*lambda,vsubs,k,w);
       
      thresh=norm(phiR-v,"inf");
      phiOLDOLD=phiOLD;
      phiOLD=phiR;
      j+=1;
    }
  phiFin.row(i)=phiR;
    }
  return(phiFin);

}

// [[Rcpp::export]]
cube gamloopOO(NumericVector beta_, const mat Y,const mat Z, colvec gammgrid, const double eps,const colvec YMean2, const colvec ZMean2,mat B1, const int k, const int p,colvec w, List groups_){


 mat B1F2=B1;

vec eigval;
mat eigvec;
 const mat Zt=Z*trans(Z);
eig_sym(eigval, eigvec, Zt);

 double tk=1/max(eigval); 


 const int ngridpts=gammgrid.n_elem;
 cube bcube(beta_.begin(),k,k*p,ngridpts,false);
 cube bcube2(k,k*p+1,ngridpts);
 bcube2.fill(0);
colvec nu=zeros<colvec>(k);

 int i;

  for (i=0; i<ngridpts;++i) {
            

	 B1F2=bcube.slice(i);
	 B1 = FistaOO(Y,Z,B1F2,p,k,gammgrid[i],groups_,eps,tk,w); 
	  
	 nu = YMean2 - B1 *ZMean2;
         bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}

    return(bcube2);
}




	
	








// *
// *
// *
// *
// // Elementwise HVAR
// *
// *
// *
// *

uvec vsubscppelem(int p,int pmax)
{
  uvec vs(pmax-p+1);
  for(int i=pmax;i>=p;--i)
    {
      vs(i-p)=i-1;
    }
  return(vs);
}
uvec bbsubs(int j,int k,int p)
{
  uvec bb(p);
  bb(0)=j;
  for(int i=1;i<p;++i)
    {
      bb(i)=j+k*(i);

    }
  return(bb);

}


rowvec proxcppelem(colvec v2,int L,double lambda,uvec res1,colvec w)
{
	colvec r =v2;
      for(int i=(L-1); i>=0;--i)

	{
           
	  
	  uvec res=vsubscppelem(i+1,L);
	  


	  if(norm(r(res)/(lambda*w(i)),"fro")<1+pow(10,-8))
	    {
			r(res)=zeros(res.n_elem);
	    }
	  else{
		  r(res)=r(res)-lambda*w(i)*r(res)/(norm(r(res),"fro"));
	  }

	}
  

      return(trans(r));


}


rowvec prox2(colvec v,double lambda, int k,int p,uvec res1,colvec w)
{
  rowvec v2(v.n_elem);
  rowvec v3(p);
  for(int i=0;i<k;++i)
    {
      uvec bb=bbsubs(i,k,p);
      colvec v1=v(bb);
      v3=proxcppelem(v1,p,lambda,res1,w);
      v2(bb)=v3;
   }
  return(v2);
}


// [[Rcpp::export]]
mat FistaElem(const mat& Y,const mat& Z, mat phi, const int p,const int k,double lambda, const double eps,const double tk)
{
  double j=1;
  mat phiFin=phi;
  rowvec phiR=phi.row(0); 
  rowvec phiOLD=phiR;
  rowvec phiOLDOLD=phiOLD;
  rowvec v=phiOLD;
  //changed to allow for maxlag of 1 from res1=ind(p,1
  uvec res1=ind(p,0);
  colvec w(p);
  w.ones();


  for(int i=0;i<k;++i)
    {
   j=1;
   double thresh=10*eps;
   phiR=phi.row(i);
   phiOLD=phiR;
   phiOLDOLD=phiOLD;
   v=phiR;
  while(thresh>eps)
    {
       v=phiOLD+((j-2)/(j+1))*(phiOLD-phiOLDOLD);

       phiR=prox2(vectorise(v)+tk*vectorise((trans(Y.col(i))-v*Z)*trans(Z)),tk*lambda,k,p,res1,w);
       
      thresh=norm(phiR-v,"fro");
      phiOLDOLD=phiOLD;
      phiOLD=phiR;
      j+=1;
    }
  phiFin.row(i)=phiR;
    }
  return(phiFin);





}

// Lamba loop
// [[Rcpp::export]]
cube gamloopElem(NumericVector beta_, const mat& Y,const mat& Z, colvec gammgrid, const double eps,const colvec YMean2, const colvec ZMean2,mat B1, const int k, const int p){

 mat B1F2=B1;

vec eigval;
mat eigvec;
 const mat Zt=Z*trans(Z);
eig_sym(eigval, eigvec, Zt);

 double tk=1/max(eigval); 


 const int ngridpts=gammgrid.n_elem;
 cube bcube(beta_.begin(),k,k*p,ngridpts,false);
 cube bcube2(k,k*p+1,ngridpts);
 bcube2.fill(0);
colvec nu=zeros<colvec>(k);

 int i;

  for (i=0; i<ngridpts;++i) {
            

	 B1F2=bcube.slice(i);
	 B1 = FistaElem(Y,Z,B1F2,p,k,gammgrid[i],eps,tk); 
	  
	 nu = YMean2 - B1 *ZMean2;
         bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}

    return(bcube2);
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

List powermethod(mat A, colvec x1) {
   double dd = 1.0;
int nn=x1.n_elem;
arma::mat x(x1.begin(),nn,1,false);
   double eps = .001;
   arma::mat y=x;
//double alpha=0;
double theta=0;
   while (dd > eps*fabs(theta)) {
    x=y/arma::norm(y,2);
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



cube FistaparDL(const mat Y,const mat Z,cube phi, const int L,const cube lambda,const List vsubs_,const double eps,const double tk,const int k, int ngridpts)
{
	mat phiFIN=zeros(k,k*L);






  colvec w(L);
  for(int r=0;r<L;++r)

    {
      w(r)=sqrt(k);

    }
  const colvec w2=w;

mat phiC=phi.slice(0); 
 rowvec phiR=phiC.row(0);
 rowvec v=phiR;
 rowvec phiOLD=phiR;
 rowvec phiOLDOLD=phiR;
 int i;
		  
for( i=0;i<k;++i)
    {		  
  colvec gammgrid= lambda.slice(i);
    for(int q=0;q<ngridpts;++q)
	  {
 phiC=phi.slice(q);
  phiR=phiC.row(i);
 phiOLD=zeros(1,k*L);
  phiOLDOLD=phiOLD;
 v=phiR;
   double thresh=10*eps;
   double   j=1;
   double lam=gammgrid[q];
  while(thresh>eps)
    {

      v=phiOLD+((j-2)/(j+1))*(phiOLD-phiOLDOLD);
      phiR=proxcpp(vectorise(v)+tk*vectorise((trans(Y.col(i))-v*Z)*trans(Z)),L,tk*lam,k,w2);
      thresh=max(abs(phiR-v));
      phiOLDOLD=phiOLD;
      phiOLD=phiR;
      j+=1;
  
 }
  phi.slice(q).row(i)=phiR;
    }
	}
  return(phi);
}

					 
// [[Rcpp::export]]
cube gamloopHVARDL(NumericVector beta_, const mat& Y,const mat& Z, NumericVector gammgrid_, const double eps,const colvec& YMean2, const colvec& ZMean2, const int k, const int p,List vsubs_,int gran2){

// omp_set_num_threads(4);

//  setenv("OMP_STACKSIZE","64M",1);

// mat B1F2=B1;

vec eigval;
mat eigvec;
 const mat& Zt=Z*trans(Z);
eig_sym(eigval, eigvec, Zt);

 double tk=1/max(eigval); 

// omp_set_num_threads(6);

 const int ngridpts=gran2;
 cube bcube(beta_.begin(),k,k*p,ngridpts,false);
 cube bcube2(k,k*p+1,ngridpts);
 bcube2.fill(0);
 cube lambdas(gammgrid_.begin(),gran2,1,k);
colvec nu=zeros<colvec>(k);

 int i;

 bcube = FistaparDL(Y,Z,bcube,p,lambdas,vsubs_,eps,tk,k,ngridpts); 

  for (i=0; i<ngridpts;++i) {
            

mat	 B1=bcube.slice(i);
 // B1.print();
	 nu = YMean2 - B1 *ZMean2;
         bcube2.slice(i) = mat(join_horiz(nu, B1)); 
	}

    return(bcube2);
}

// [[Rcpp::export]]
NumericVector HLassoVAR(NumericVector betaL_, mat& Y, mat& Z, const NumericVector gamm_, double eps, int k, int p,List groups,int gran2)
{
	Y=trans(Y);
	colvec YMean=mean(Y,1);

    colvec ZMean=mean(Z,1);
	colvec zones=ones(Y.n_cols,1);
	colvec yones=ones(Z.n_cols,1);	
	Y=Y-YMean*trans(yones);
	Z=Z-ZMean*trans(zones);
	Y=trans(Y);
	
	cube betafin = gamloopHVARDL(betaL_,Y,Z,gamm_,eps,YMean,ZMean,k,p,groups,gran2);
    
	return(as<NumericVector>(wrap(betafin)));
        }



// [[Rcpp::export]]
cube HVARCVAL(NumericVector beta_,const List Zfull_, const NumericVector gamm_, double eps,int T0, int T1, int p, int k,const mat& Z2,int gran2,List groups)
{
	mat Z1=Zfull_("Z");
	mat Y1=Zfull_("Y");
	cube MSFE(T1-T0+1,gran2,k);
	MSFE.fill(0.0);
	// MSFE.print();
	for(int i=T0;i<T1+1;++i)
		{
			mat Z1a=Z1.cols(0,i-p-2);
			mat Y1a=Y1.rows(0,i-p-2);
			// Rcout<< Z1a.col(i-p-2) <<std::endl;
            // Rcout<< Y1a <<std::endl;

			beta_=HLassoVAR(beta_,Y1a,Z1a,gamm_,eps,k,p,groups,gran2);
			colvec eZ=Z2.col(i-1);
				    // Rcout<< eZ <<std::endl;
			cube beta1(beta_.begin(),k,k*p+1,gran2,false);
					for(int q=0; q<k; ++q)
				{
			for(int j=0;j<gran2;++j)
				{
					mat B1=beta1.slice(j);
					
					rowvec B1a=B1.row(q);
					// B1a.print();
					double loss=as_scalar(Y1(i-1,q)-B1a*eZ);
				    MSFE.slice(q)(i-T0,j)=pow(loss,2);
				    // Rcout<< loss <<std::endl;
					 
					// MSFE(i-T0,j)=pow(loss,2);

				}
			
				}

				 Rcout<< i <<std::endl;

        }
	
	return(MSFE);
}

// [[Rcpp::export]]
mat HVAREval(NumericVector beta_,const List Zfull_, const NumericVector gamm_, double eps,int T1, int T2, int p, int k,const mat& Z2, int gran2, List groups)
{
	mat Z1=Zfull_("Z");
	mat Y1=Zfull_("Y");
	// int ngridpts=alpha.n_elem*gamm.n_elem;
	mat MSFE(T2-T1+1,gran2);
	for(int i=T1;i<T2+1;++i)
		{
			mat Z1a=Z1.cols(0,i-p-2);
			mat Y1a=Y1.rows(0,i-p-2);
			// Rcout<< Z1a.col(i-p-2) <<std::endl;
            // Rcout<< Y1a <<std::endl;

			beta_=HLassoVAR(beta_,Y1a,Z1a,gamm_,eps,k,p,groups,gran2);
			colvec eZ=Z2.col(i-1);
				    // Rcout<< eZ <<std::endl;
			cube beta1(beta_.begin(),k,k*p+1,gran2,false);
				// 	for(int q=0; q<k; ++q)
				// {
			for(int j=0;j<gran2;++j)
				{
					mat B1=beta1.slice(j);
					
					// rowvec B1a=B1.row(q);
					// B1a.print();
					double loss=norm(trans(Y1.row(i-1))-B1*eZ,"fro");
				    MSFE(i-T1,j)=pow(loss,2);
				    // Rcout<< loss <<std::endl;
					 
					// MSFE(i-T0,j)=pow(loss,2);

				}
			
				

				 Rcout<< i <<std::endl;
				 // Rcout<< i <<std::endl;

        }
	
	return(MSFE);
}

mat QRF(const mat& K, mat R5, int i, int kp, int k1, int p)
{
	int RC=R5.n_cols;

	mat RA=zeros(kp+k1,RC+1);
			for(int j=0; j<RC;++j)
				{
		    for(int q=0; q<kp;++q)
		        {

					RA(q,j)=R5(q,j);
				}
			    }

			for(int j = kp; j<kp+k1;++j)
				{
					if(j-kp==i)
						{
							RA(j,RC)=1;
						}


				}
			mat K2=K*RA;
			int q=K2.n_cols;
			double delta=(pow(q,2)+q+1)*sqrt(std::numeric_limits<double>::epsilon());
			colvec D1=zeros(K2.n_cols);
			for(int i=0;i<q;++i)
			  {
			    D1(i)=norm(K2.col(i),"fro");
			  }

			D1=sqrt(delta)*D1;
			mat AA=diagmat(D1);
			mat K3=mat(join_vert(K2,AA));
			mat Q1, R1;

			qr(Q1,R1,K3);
			mat R11=R1.submat(0,0,RC-1,RC-1);
			colvec R22=vectorise(R1.submat(0,RC,RC-1,RC));
			mat RLS=solve(R11,R22);

			return(RLS);
								 
					



}

// template <typename T> 
// const bool Contains( std::vector<T>& Vec, const T& Element ) 
// {
//     if (std::find(Vec.begin(), Vec.end(), Element) != Vec.end())
//         return true;

//     return false;
// }
 
// [[Rcpp::export]]
mat RelaxedLS(const mat K,  mat B2, int k, int p,int k1, int s)
{
	mat B3=B2.cols(1,B2.n_cols-1);
	// mat B3=B2;
	if(norm(B3,"inf")==0){return(B2);}
	else{

		int kp=B3.n_cols;
		colvec nu=B2.col(0);
		mat A=zeros(k1,kp);
	
		for(int i=0;i<k1;++i)
			{
				rowvec B3a=B3.row(i);
				unsigned int thresh=pow(10,-8);

				uvec R1a=find(abs(B3a)>thresh);
				if(R1a.n_elem<2){A.row(i)=B3a;}
				else{
				mat R5=zeros(kp,R1a.n_elem);
					
				int jj=0;
				std::vector<int> R1(R1a.begin(),R1a.end());
				
				for(int ii=0; ii<kp;++ii)
					{
						// IntIterator q = ;
						// std::find(R1.begin(),R1.end(),ii)!=R1.end()
						// Rcout<< Contains(R1,ii) <<std::endl;
		
						if(std::find(R1.begin(),R1.end(),ii)!=R1.end()){
						  
							R5(ii,jj)=1;
						    jj+=1;
						}
						
					}
				mat RLS=QRF(K,R5,i,kp,k1,p);
				mat R6=R5*RLS	;
				R6.reshape(1,kp);
			    A.row(i)=R6;

				}
			}
		mat BR=mat(join_horiz(nu, A));
		return(BR);		
				



	}
}
	

// [[Rcpp::export]]
List ARFit(const mat& K2, int k, int p)
{
                        int RC=k*p+1;   
                        int T=K2.n_rows;
			int q=K2.n_cols;
			double delta=(pow(q,2)+q+1)*sqrt(std::numeric_limits<double>::epsilon());
			colvec D1=zeros(K2.n_cols);
			for(int i=0;i<q;++i)
			  {
			    D1(i)=norm(K2.col(i),"fro");
			  }

			D1=sqrt(delta)*D1;
			// D1.print();
			mat AA=diagmat(D1);
			mat K3=mat(join_vert(K2,AA));
			mat Q1, R1;

			qr(Q1,R1,K3);
			mat R11=R1.submat(0,0,RC-1,RC-1);
			mat R22=R1.submat(RC,RC,RC+k-1,RC+k-1);
			// colvec R12=vectorise(R1.submat(0,RC,RC-1,RC));
			mat R12=R1.submat(0,RC,RC-1,RC+k-1);
			
			mat RLS=solve(R11,R12);
			mat Sigma=trans(R22)*R22/(T);
			List Results=List::create(Named("Rhat")=trans(RLS),Named("SigmaU")=Sigma);
			return(Results);
								 
					
			






}

// [[Rcpp::export]]
mat Zmat1(mat& Y, int p, int k)
{
	mat Y2=fliplr(Y);

	int T=Y.n_rows;
	colvec Y1=vectorise(trans(Y2.submat(0,0,p-1,k-1)));
	std::vector<double> Y1a(Y1.begin(),Y1.end());
	Y1a.push_back(1);
	std::reverse(Y1a.begin(),Y1a.end());

	mat Z;
	Z.zeros(k*p+1,T-p);
	colvec Y1c = conv_to<colvec>::from(Y1a);
	Z.col(0)=Y1c;
	for(int i=1; i<T-p;++i)
		{
	colvec Y1=vectorise(trans(Y2.submat(i,0,p+i-1,k-1)));
	std::vector<double> Y1b(Y1.begin(),Y1.end());
	Y1b.push_back(1);
	std::reverse(Y1b.begin(),Y1b.end());
    Y1c = conv_to<colvec>::from(Y1b);
	Z.col(i)=Y1c;
			


		}
	
	return(Z);

}

// [[Rcpp::export]]
List AIC(mat& Y, double k, int pmax)
{
         std::vector<double> AIC;
	// colvec AIC;
	// AIC.zeros(pmax);
	int T=Y.n_rows;
	mat ZF=Zmat1(Y,pmax,k);
	
	for(int i=1;i<=pmax;++i)
		{
			if((int) k*i>(T))
				{break;}
			mat sigmaU;
			mat Z=Zmat1(Y,i,k);


			mat Y2=Y.rows(i,T-1);

		
			mat K=join_horiz(trans(Z),Y2);
		    List BRes=ARFit(K,k,i);
			sigmaU = as<mat>(BRes("SigmaU"));
			double T2=Y2.n_rows;
			double p=(double) i;

			double AICi=log(det(sigmaU))+2*pow(k,2)*p/T2;

			AIC.push_back(AICi);
			


		}
	int min= std::distance(AIC.begin(), std::min_element(AIC.begin(), AIC.end()));
	
	min+=1; //to match R indexing convention

	mat ZFF=Zmat1(Y,min,k);


	mat Y2F=Y.rows(min,T-1);

	mat K2=join_horiz(trans(ZFF),Y2F);

	mat BResF=ARFit(K2,k,min)("Rhat");
	
	List Results=List::create(Named("B")=BResF,Named("p")=min);	
	// return(min);
	return(Results);

}


// [[Rcpp::export]]
List BIC(mat& Y, double k, int pmax)
{
         std::vector<double> BIC;
	int T=Y.n_rows;
	
	for(int i=1;i<=pmax;++i)
		{
			// if((int)k*i>T-10)
		    mat Z=Zmat1(Y,i,k);
				  // Rcout<<cond(Z*trans(Z))<<std::endl;

		    if(cond(Z*trans(Z))>pow(10,6))
				{
break;}



			mat Y2=Y.rows(i,T-1);

			mat K=join_horiz(trans(Z),Y2);
		    List BRes=ARFit(K,k,i);
			mat sigmaU = as<mat>(BRes("SigmaU"));
			double T2=Y2.n_rows;
			double p=(double) i;
			double BICi=log(det(sigmaU))+log(T2)*pow(k,2)*p/T2;
			// Rcout <<BICi<<std::endl;

			BIC.push_back(BICi);
			


		}
	int min= std::distance(BIC.begin(), std::min_element(BIC.begin(), BIC.end()));

		
	min+=1; //to match R indexing convention

	mat ZFF=Zmat1(Y,min,k);

	mat Y2F=Y.rows(min,T-1);

	mat K2=join_horiz(trans(ZFF),Y2F);

	mat BResF=ARFit(K2,k,min)("Rhat");
	
	List Results=List::create(Named("B")=BResF,Named("p")=min);	

	return(Results);

}



// [[Rcpp::export]]
List EvalIC(mat& Y, int T1, int k, int pmax,std::string IC)
{
	       int T2=Y.n_rows;

	colvec MSFE=zeros(T2-T1);
	colvec order=zeros(T2-T1);

	for(int i=T1;i<T2;++i)
		{
			mat Y1a=Y.rows(0,i-1);

			List popt;
            if(IC=="AIC"){
				 popt=AIC(Y1a,k,pmax);}
			            if(IC=="BIC"){
				 popt=BIC(Y1a,k,pmax);
				 		}

			int p=popt("p");
			mat B1=as<mat>(popt("B"));
			mat evalY=Y.rows(i-p,i);

			colvec eZ=Zmat1(evalY,p,k);
			order(i-T1)=p;
			colvec loss=trans(Y.row(i))-B1*eZ;
					 
			MSFE(i-T1)=pow(norm(loss,"fro"),2);

				}
			
        
	List Results=List::create(Named("MS")=MSFE,Named("p")=order);	
	return(Results);
}



// [[Rcpp::export]]
List QRCons(mat& Z, int T, int k,mat B2,int R2, int p)
{
	int kp=k*p+1;
	int K1OLD=0;
	int K2OLD=0;
	mat Q1AA=zeros(k*T,k*T);
	mat Q1AR=zeros(k*T,R2);
	for(int i=0; i<k;++i)
		{
  			rowvec B3a=B2.row(i);
			unsigned int thresh=pow(10,-8);
				uvec R1a=find(abs(B3a)>thresh);
				mat R5=zeros(kp,R1a.n_elem);
				int jj=0;
				std::vector<int> R1(R1a.begin(),R1a.end());
				for(int ii=0; ii<kp;++ii)
					{
		
						if(std::find(R1.begin(),R1.end(),ii)!=R1.end()){
						  
							R5(ii,jj)=1;
						    jj+=1;
						}
						
					}
	    				mat Q, R;
						mat Z2=Z*R5;
					   
						qr(Q,R,Z2);
						int K1=R5.n_cols;
						R=R.rows(0,Z2.n_cols-1);
						Q1AA.submat(i*T,K1OLD,(i+1)*T-1,K1+K1OLD-1)=Q.cols(0,K1-1);
						Q1AR.submat(K1OLD,K1OLD,K1+K1OLD-1,K1+K1OLD-1)=R;
						K1OLD=K1+K1OLD;
						mat Q1C=Q.cols(K1,Q.n_cols-1);
						int K2=Q1C.n_cols;
						Q1AA.submat(i*T,K2OLD+R2,(i+1)*T-1,K2+K2OLD+R2-1)=Q1C;
						K2OLD=K2+K2OLD;
						
		}
	List Q1F=List::create(Named("Q1A")=Q1AA,Named("Q1R")=Q1AR);
						
return(Q1F);
		}

mat sparseWLX(const mat& M1a,const  mat& R1, double ngroups, mat& beta,  double t,  double alpha,  double lambda, double eps)
{
  int n=M1a.n_rows;//k=M1a.n_cols;
  
int n2=beta.n_rows,k2=beta.n_cols;

double thresh=10;

 mat p=zeros<mat>(n2,k2);

 mat STS=zeros<mat>(n,1);
 mat thetaOLD=beta;
 mat thetaOLDOLD=beta;
 double l=1;
 mat u=zeros<mat>(n2,k2);

while(thresh>eps)
  {
    p = (beta*M1a-R1)*trans(M1a)/ngroups;
   
    STS=ST3a(vectorise(beta)-t*vectorise(p),t*alpha*lambda);
    double denom2= norm(STS,"fro");
    double s31=fmax(1-(t*(1-alpha)*lambda)/denom2,0);
    STS=s31*STS;
    STS.set_size(n2,k2);
    
    beta=thetaOLD+(l/(l+3))*(STS-thetaOLD);
    l+=1;
    if(l>1000){
      Rcout<< "Warning: Sparse GL Did Not Converge" <<std::endl;
      break;
    }
    mat thresh1=beta-STS;
     thresh=norm(thresh1,"inf");
 
    thetaOLD=STS;
  }
  
return(beta);
  }

// Very similar to group lasso case
List blockUpdateSGLX(mat& beta,const mat& Z1, double lam, double alpha,const mat& Y1, double eps, List groups, List fullgroups, List compgroups, int k1, List M2f_,NumericVector Eigs)
{

  int iters=0;
 int n1=groups.size();
  List active(n1);
  int n=beta.n_rows, m=beta.n_cols;

  arma::mat betaPrev=beta;
  int converge=0;
 
  if(groups.size()==0)
    {
      beta.zeros(n,m);
      active=groups;
    }
  else{
    for(int i=0; i<n1;++i)
 
      {


	arma::uvec s45=as<arma::uvec>(groups[i]);
	arma::uvec s45F=as<arma::uvec>(fullgroups[i]);


 // Rcpp::Rcout << "beta rows" <<n <<  std::endl;
      

	
 if(max(s45)==0){
   beta.cols(s45F)=zeros(k1,s45F.n_elem);
	  active(i)=0;

 }	
  else{
	  
	   arma::uvec scomp1=as<arma::uvec>(compgroups[i]);


	
	  arma::mat M2a= Z1.rows(scomp1);
	  arma::mat a1= beta.cols(scomp1);
	  arma::mat beta2=beta.cols(s45);

	  arma::mat r=Y1-a1*M2a;

	  arma::mat M1=Z1.rows(s45);
	  arma::mat M2=M2f_(i);
	  arma::mat p=(beta2*M1-r)*trans(M1);


	  double rho=sqrt(s45.size());
	  colvec STS = ST3a(vectorise(p),alpha*lam);
	    

	  double lamadj=lam*(1-alpha)*rho;
	  if(arma::norm(STS,"fro")<=lamadj)
	    {

	      mat astar = zeros(k1,s45.n_elem);
	      active(i)=0;


	    }
	  else{
	    mat betaS=beta.cols(s45);
	    double t=1/Eigs(i);
	    mat astar2= sparseWLX(M1, r,k1, betaS, t,  alpha, lam, eps);
	    beta.cols(s45)=astar2;
	    active(i)=s45;  
	    
	    
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

mat ThreshUpdateSGLX(mat& betaActive,const mat& Z, double lam,const mat& Y,double eps, List groups, List fullgroups, List compgroups,List M2f, NumericVector eigs, double alpha, int k1)
  {
    int iter=0;
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
        betaActive2=blockUpdateSGLX(betaActive,Z,lam,alpha,Y,eps,groups,fullgroups,compgroups,k1,M2f,eigs);	    
        betaActive=as<mat>(betaActive2("beta"));
	arma::mat thresh1=arma::abs((betaLast-betaActive)/(arma::ones(n,m)+arma::abs(betaLast)));
        threshold=arma::norm(thresh1,"inf");
        active=betaActive2("active");
	betaLast=betaActive;
	if(iter>500){
	  Rcout<<"Warning: SGL ThreshUpdate Did not Converge"<<std::endl;
	  break;
	    }
      }
  }
    return(betaActive);
  }

// [[Rcpp::export]]
List GamLoopSGLX(NumericVector beta_, List Activeset, NumericVector gamm,double alpha, const mat& Y1, const mat& Z1,List jj, List jjfull, List jjcomp, double eps, colvec YMean2, colvec ZMean2,int k,int pk, List M2f_, NumericVector eigs,int k1)
{
 int gran2=gamm.size();
 List activefinal(gran2);

 cube beta2(beta_.begin(),k1,pk,gran2,false);
 cube betafin(k1,pk+1,gran2);
 betafin.fill(0);
 List iterations(gran2);
 mat betaPrev=zeros<mat>(k1,pk);

 //INDEX LISTS WITH PARENTHESES NOT BRACKETS :O 
 // WHEN EXTRACTING FROM A LIST YOU NEED as<MAT>
 for(int i=0; i<gran2;++i)
    {
        double gam=gamm[i];
	betaPrev=beta2.slice(i);
	List Active = Activeset[i];
	int k2=0;
	int converge=0;
        mat betaF=zeros(k1,pk);
	List betaFull(3);
	//Three components in the list
//	int thresh=0;
	while(converge==0)
	  {
	    
       	    betaPrev = ThreshUpdateSGLX(betaPrev, Z1, gam, Y1, eps, jjfull, jjfull, jjcomp, M2f_, eigs, alpha,k1);
	 
	    betaFull=blockUpdateSGLX(betaPrev,Z1,gam,alpha,Y1,eps,jjfull,jjfull,jjcomp,k1,M2f_,eigs);
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




// [[Rcpp::export]]
colvec proxvx2(colvec v2,int L,double lambda,int m,int k,int F1)
{
	colvec r =v2;
	int start;
	if(F1==0){start=1;}
	else{start=0;}
      for(int i=start; i>=0;--i)

	{
           
	  
	  // uvec res=vsubscppelem(i+1,L);
		int s1;
	//	int t1=0;
		if(F1==0){
			 s1 = (k+m)-(i*k);
	//		 t1 = i*k;
		}
		else{ s1=k;
	//		t1=k;
		}
	   std::vector<unsigned int> ivec(s1);
	   if(F1==0){ 
		   std::iota(ivec.begin(), ivec.end(), i*k);}
	   else{std::iota(ivec.begin(), ivec.end(), 0);}
		uvec res = conv_to<uvec>::from(ivec);	
		// res.print();

	  if(norm(r(res)/(lambda),"fro")<1+pow(10,-8))
	    {
			r(res)=zeros(res.n_elem);
	    }
	  else{
		  r(res)=r(res)-lambda*r(res)/(norm(r(res),"fro"));
	  }

	}
  

      return(r);


}

//These functions need to be updated/incorporated with VARX case

// Inner while loop from simon et al (2013)


mat sparseWL(const mat& M1a,const  mat& R1, double ngroups, mat& beta,  double t,  double alpha,  double lambda, double eps)
{
  int n=M1a.n_rows;//,k=M1a.n_cols;
  
int n2=beta.n_rows;//k2=beta.n_cols;

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
	//int thresh=0;
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



