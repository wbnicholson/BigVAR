#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <math.h> 
#include <vector>
#include <limits>
#include <algorithm>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppEigen)]]

// using namespace std;
using namespace Rcpp;
using namespace arma;

using namespace Eigen;


// [[Rcpp::export]]
MatrixXd backsolve_rcpp(const MatrixXd& R2, const MatrixXd R)
{

	int k = R.rows();
	MatrixXd RI(k,1);

	RI=R.triangularView<Upper>().solve(R2);
	
	return(RI);
}
						

// [[Rcpp::export]]
MatrixXd RestMat(mat& BHAT, int i, int kp,double zerothresh)
			{
				rowvec B3a=BHAT.row(i);
				uvec R1a=find(abs(B3a)>zerothresh);
				// R1a.print();
				// if(R1a.n_elem<2){A.row(i)=B3a;}
				// else{
				MatrixXd R5(kp,R1a.n_elem);
				R5.setZero();
				int jj=0;
				// R1a.print();
			    // std::vector<int> R1=conv_to<std::vector<int> >::from(R1a);
				std::vector<int> R1(R1a.begin(),R1a.end());
				// printf("%d\n",R1[0]);
				// int K1OLD=0;
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

				return(R5);

			}


// [[Rcpp::export]]
MatrixXd RestMatFull(mat& BHAT, int kpp,double zerothresh)
			{
				colvec B3a=vectorise(BHAT);
				uvec R1a=find(abs(B3a)>zerothresh);
				// R1a.print();
				// if(R1a.n_elem<2){A.row(i)=B3a;}
				// else{
				MatrixXd R5(kpp,R1a.n_elem);
				R5.setZero();
				int jj=0;
				// R1a.print();
			    // std::vector<int> R1=conv_to<std::vector<int> >::from(R1a);
				std::vector<int> R1(R1a.begin(),R1a.end());
				// printf("%d\n",R1[0]);
				// int K1OLD=0;
				for(int ii=0; ii<kpp;++ii)
					{
						// IntIterator q = ;
						// std::find(R1.begin(),R1.end(),ii)!=R1.end()
						// Rcout<< Contains(R1,ii) <<std::endl;
		
						if(std::find(R1.begin(),R1.end(),ii)!=R1.end()){
						  
							R5(ii,jj)=1;
						    jj+=1;
						}
						
					}

				return(R5);

			}



// [[Rcpp::export]]
List QRConsM2(NumericMatrix Z1, int T1, int k, mat B2, int R2, int p, double zerothresh)
{
	// Rcout<<"starting"<<std::endl;
	const Map<MatrixXd>  Z(as<Map<MatrixXd> >(Z1));

	typedef Triplet<double> T;
	std::vector<T> tripletList;
	std::vector<T> tripletList2;
	List QQ(k);
	int kp=k*p+1;
	int K1OLD=0;
	int K2OLD=0;
	SparseMatrix<double> Q1AA(k*T1,k*T1);
	Q1AA.reserve(VectorXd::Constant(k*T1,T1));
	SparseMatrix<double> Q1AR(k*T1,R2);

	for(int i=0;i<k;++i)
		{
			MatrixXd R5=RestMat(B2, i, kp, zerothresh);
			MatrixXd Z2=Z*R5;
			HouseholderQR<MatrixXd> QR1(Z2);
			MatrixXd R = QR1.matrixQR().triangularView<Upper>();
			MatrixXd Q=QR1.householderQ();
			int K1=R5.cols();
			QQ(i)=Q;
			MatrixXd QSS=Q.leftCols(K1);
			int q=0;
			int r=0;


			for(int jj=K1OLD;jj<K1+K1OLD;++jj)
				{
					r=0;
					for(int ii=i*T1;ii<(i+1)*T1;++ii)
						{
							tripletList.push_back(T(jj,ii,QSS(r,q)));
							r+=1;
							
						}
					q+=1;

				}
			MatrixXd Q1C=Q.rightCols(Q.cols()-K1);
			int K2=Q1C.cols();
			q=0;
			r=0;
				for(int jj=K2OLD+R2;jj<K2+K2OLD+R2;++jj)
				{
					r=0;
					for(int ii=i*T1;ii<(i+1)*T1;++ii)
						{
							tripletList.push_back(T(jj,ii,Q1C(r,q)));
							r+=1;
							
						}
					q+=1;

				}
				q=0;
				r=0;
				for(int jj=K1OLD;jj<K1+K1OLD;++jj)
				{
					r=0;
					for(int ii=K1OLD;ii<K1OLD+K1;++ii)
						{
							tripletList2.push_back(T(ii,jj,R(r,q)));
							r+=1;
							
						}
					q+=1;

				}
			K1OLD+=K1;
			K2OLD+=K2;
		}
	Q1AA.setFromTriplets(tripletList.begin(),tripletList.end());
	Q1AR.setFromTriplets(tripletList2.begin(),tripletList2.end());
	Q1AA.makeCompressed();
	List Results=List::create(Named("QQ")=QQ,Named("Q1A")=Q1AA,Named("Q1R")=Q1AR);
	return(Results);
}


// // [[Rcpp::export]]
// int RMatCol1(mat B2, int i,double zerothresh)
// 			 {
// 				 if(i<0){return(0);}
// 				 else{
// 			   int k1=0;
// 			   for(int j=0;j<(i+1);++j)
// 				   {
//   			   rowvec B3a=B2.row(j);					   
// 			   uvec R1a=find(abs(B3a)>pow(10,-8));
// 			   k1+=R1a.n_elem;
// 			   // Rcout<<R1a.n_elem;
// 				   }
// 			   return(k1);

// 			 }
// 			 }
// // [[Rcpp::export]]
// int RMatCol1Comp(mat B2, int i, int T,double zerothresh)
// 			 {
// 				 if(i<0){return(0);}
// 				 else{
// 			   int k1=0;
// 			   for(int j=0;j<i+1;++j)
// 				   {
//   			   rowvec B3a=B2.row(j);					   
// 			   uvec R1a=find(abs(B3a)>zerothresh);
// 			   k1+=T-R1a.n_elem;
// 				   }
// 			   return(k1);
// 				 }
// 			 }
// // // [[Rcpp::export]]
// // mat RestMat(mat BHAT, int i, int kp,double zerothresh)
// // 			{
// // 				rowvec B3a=BHAT.row(i);
// // 				uvec R1a=find(abs(B3a)>zerothresh);
// // 				// R1a.print();
// // 				// if(R1a.n_elem<2){A.row(i)=B3a;}
// // 				// else{
// // 				mat R5=zeros(kp,R1a.n_elem);
// // 				int jj=0;
// // 				// R1a.print();
// // 			    // std::vector<int> R1=conv_to<std::vector<int> >::from(R1a);
// // 				std::vector<int> R1(R1a.begin(),R1a.end());
// // 				// printf("%d\n",R1[0]);
// // 				// int K1OLD=0;
// // 				for(int ii=0; ii<kp;++ii)
// // 					{
// // 						// IntIterator q = ;
// // 						// std::find(R1.begin(),R1.end(),ii)!=R1.end()
// // 						// Rcout<< Contains(R1,ii) <<std::endl;
		
// // 						if(std::find(R1.begin(),R1.end(),ii)!=R1.end()){
						  
// // 							R5(ii,jj)=1;
// // 						    jj+=1;
// // 						}
						
// // 					}

// // 				return(R5);

// // 			}



// [[Rcpp::export]]
int RMatCol1(mat& B2, int i,double zerothresh)
			 {
				 if(i<0){return(0);}
				 else{
			   int k1=0;
			   for(int j=0;j<(i+1);++j)
				   {
  			   const rowvec& B3a=B2.row(j);					   
			   const uvec& R1a=find(abs(B3a)>zerothresh);
			   k1+=R1a.n_elem;
			   // Rcout<<R1a.n_elem;
				   }
			   return(k1);

			 }
			 }
// [[Rcpp::export]]
int RMatCol1Comp(mat& B2, int i, int T1,double zerothresh)
			 {
				 if(i<0){return(0);}
				 else{
			   int k1=0;
			   for(int j=0;j<i+1;++j)
				   {
  			 const  rowvec& B3a=B2.row(j);					   
			 const  uvec& R1a=find(abs(B3a)>zerothresh);
			   k1+=T1-R1a.n_elem;
				   }
			   return(k1);
				 }
			 }


// [[Rcpp::export]]
List KronMatcppEigen(int k, int T1, mat Z,int RR, mat BHAT, mat C,int p, int R2i,List QQ_,double zerothresh,NumericMatrix Ystar1, NumericMatrix Q1AR12)
{


	// Rcout<<"Test1"<<endl;

	const Map<MatrixXd>  Ystar(as<Map<MatrixXd> >(Ystar1));

	const Map<MatrixXd>  Q1AR(as<Map<MatrixXd> >(Q1AR12));

	int kp=k*p+1;

	// Rcout<<"Test1"<<std::endl;
	typedef Triplet<double> T;
	std::vector<T> tripletList;
	std::vector<T> tripletList2;
	std::vector<T> tripletList3;
	std::vector<T> tripletList4;
	// Rcout<<(RR*(k*T1-RR))/2<<std::endl;
	tripletList3.reserve((RR*(k*T1-RR))/2);
	tripletList4.reserve((RR*(k*T1-RR))/2);
	int nreserve=pow(k*T1-RR,2);
	tripletList.reserve(nreserve/2);

	
	// Rcout<<"Test2"<<std::endl;
	
	// mat W22=zeros(k*T-RR,k*T-RR);
    int RBOLD=0;
	for(int i=0;i<k;++i)
		{
			rowvec B3a=BHAT.row(i);
 		    uvec R1a=find(abs(B3a)>zerothresh);
			int RB=T1-R1a.n_elem;
			int q=0;
			int r=0;
			for(int jj=RBOLD;jj<RBOLD+RB;++jj)
				{
					r=0;
					for(int ii=RBOLD;ii<RBOLD+RB;++ii)
						{
							if(r==q){
							tripletList.push_back(T(ii,jj,C(i,i)));
							}
						    r+=1;
						}
					q+=1;
				}			
			RBOLD+=RB;
		}
    RBOLD=0;
	for(int i=0;i<k;++i)
		{

			// Rcout<<"Test3"<<std::endl;
		   rowvec B3a=BHAT.row(i);

		   uvec R1a=find(abs(B3a)>zerothresh);

		    int RB=R1a.n_elem;
			int q=0;
			int r=0;
			for(int jj=RBOLD;jj<RBOLD+RB;++jj)
				{
					r=0;
					for(int ii=RBOLD;ii<RBOLD+RB;++ii)
						{
							if(r==q){
						tripletList2.push_back(T(ii,jj,C(i,i)));
							}
						    r+=1;
						}
					q+=1;
				}			
			
			// W11.submat(RBOLD,RBOLD,RB+RBOLD-1,RB+RBOLD-1)=C(i,i)*D1a.eye();
			RBOLD+=RB;
		}
	// Rcout<<"Test4"<<std::endl;

	MatrixXd Q1B= QQ_(0);
	MatrixXd Q2B=QQ_(0);
	MatrixXd QTild1=Q1B;
	MatrixXd QTild2=Q1B;
	MatrixXd QHat1=Q1B;
	MatrixXd QHat2=Q1B;
	MatrixXd WP=Q1B;

	MatrixXd R1(k*k*p,k*k*p);
	MatrixXd R2(k*k*p,k*k*p);
	MatrixXd Atild=Q1B;
	MatrixXd Atild2=Q1B;
	MatrixXd Ahat=Q1B;

	for(int j=0;j<(k-1);++j)
		{
			int ROLD=RMatCol1(BHAT,j,zerothresh);
		   			
		}
		// mat W21=zeros(RR,k*T-RR);
		// mat W12=zeros(RR,k*T-RR);
		for(int j=0;j<(k-1);++j)
			{
				int ROLD=RMatCol1(BHAT,j,zerothresh);

				for(int i=j+1;i<k;++i)
					{


						R1= RestMat(BHAT,i,kp,zerothresh);
						R2= RestMat(BHAT,j,kp,zerothresh);
						Q1B= QQ_(i);
						Q2B=QQ_(j);
						int K1=R1.cols();
						int K2=R2.cols();
						QTild1=Q1B.leftCols(K1);
						QTild2=Q2B.leftCols(K2);
						QHat1=Q1B.rightCols(Q1B.cols()-K1);
						QHat2=Q2B.rightCols(Q2B.cols()-K2);
						WP=QTild2.transpose()*QTild1;
						WP*=C(j,i);
						// R1.resize(0,0);
						// R2.resize(0,0);
						// Q1B.resize(0,0);
						// Q2B.resize(0,0);
						
						// Rcout<<"Max Qhat1"<<WP.lpNorm<Infinity>()<<std::endl;
						// Rcout<<"Max Qhat2"<<QHat2.lpNorm<Infinity>()<<endl;
						// Rcout<<"Max QTild1"<<QTild1.lpNorm<Infinity>()<<endl;
						// Rcout<<"Max QTild2"<<QTild2.lpNorm<Infinity>()<<endl;
						
						int Rrows1=RMatCol1(BHAT,j-1,zerothresh);
						int Rcols1=RMatCol1(BHAT,i-1,zerothresh);
						int Rrows2=RMatCol1(BHAT,j-1,zerothresh)+K2;
						int Rcols2=RMatCol1(BHAT,i-1,zerothresh)+K1;

						int q=0;
						int r=0;

						// Rcout<<WP<<std::endl;

						for(int jj=Rcols1;jj<Rcols2;++jj)
							{
								r=0;
								for(int ii=Rrows1;ii<Rrows2;++ii)
									{
										tripletList2.push_back(T(ii,jj,WP(r,q)));
								r+=1;							   
							        }
							q+=1;

							}


						
						// Rcout<<"j"<<j<<std::endl;
						// Rcout<<RMatCol1(BHAT,j-1)<<std::endl;
						// W11.submat(RMatCol1(BHAT,j-1,zerothresh),RMatCol1(BHAT,i-1,zerothresh),RMatCol1(BHAT,j-1,zerothresh)+K2-1,RMatCol1(BHAT,i-1,zerothresh)+K1-1)=WP;
						int K3=QHat1.cols();
						Atild2= QTild2.transpose()*QHat1;
						Atild2*=C(j,i);
						// Rcout<<"Max Atild2"<<Atild2.lpNorm<Infinity>()<<endl;
						int QF=Atild2.rows();
						 Rrows1=RMatCol1(BHAT,j-1,zerothresh);
						 Rcols1=RMatCol1Comp(BHAT,i-1,T1,zerothresh);
						 Rrows2=RMatCol1(BHAT,j-1,zerothresh)+QF;
						 Rcols2=RMatCol1Comp(BHAT,i-1,T1,zerothresh)+K3;

						 q=0;
						 r=0;
						 for(int jj=Rcols1;jj<Rcols2;++jj)
						 	{
						 		r=0;
						 		for(int ii=Rrows1;ii<Rrows2;++ii)
						 			{
						 				tripletList4.push_back(T(ii,jj,Atild2(r,q)));
						 		r+=1;							   
						 	        }
						 	q+=1;

						 	}
						
						 
						 
						
						// int QF=Atild2.n_rows;
						// W12.submat(RMatCol1(BHAT,j-1,zerothresh),RMatCol1Comp(BHAT,i-1,T,zerothresh),RMatCol1(BHAT,j-1,zerothresh)+QF-1,RMatCol1Comp(BHAT,i-1,T,zerothresh)+K3-1)=Atild2;
						Atild=QTild1.transpose()*QHat2;
						Atild*=C(j,i);
						// MatrixXd Atild2=C(j,i)*QTild2.transpose()*QHat1;
						// Rcout<<"Max Atild"<<Atild.lpNorm<Infinity>()<<endl;

						// Rcout<<Atild2<<std::endl;
						K3=Atild.rows();
					    QF=Atild.cols();


						 Rrows1=ROLD;
						 Rcols1=RMatCol1Comp(BHAT,j-1,T1,zerothresh);
						 Rrows2=ROLD+K3;
						 Rcols2=RMatCol1Comp(BHAT,j-1,T1,zerothresh)+QF;

						 q=0;
						 r=0;
						 for(int jj=Rcols1;jj<Rcols2;++jj)
							{
								r=0;
								for(int ii=Rrows1;ii<Rrows2;++ii)
									{
										tripletList3.push_back(T(ii,jj,Atild(r,q)));
										// tripletList4.push_back(T(ii,jj,Atild(r,q)));
								r+=1;							   
							        }
							q+=1;

							}
											   
						ROLD+=K3;
						Ahat=QHat2.transpose()*QHat1;
						Ahat*=C(j,i);
						// Rcout<<"Max Ahat "<<Atild.lpNorm<Infinity>()<<endl;

						int A1 = Ahat.rows();
						int A2 = Ahat.cols();

						
						Rrows1=RMatCol1Comp(BHAT,j-1,T1,zerothresh);
						Rcols1=RMatCol1Comp(BHAT,i-1,T1,zerothresh);
						Rrows2=RMatCol1Comp(BHAT,j-1,T1,zerothresh)+A1;
						Rcols2=RMatCol1Comp(BHAT,i-1,T1,zerothresh)+A2;

						 q=0;
						 r=0;
						 for(int jj=Rcols1;jj<Rcols2;++jj)
							{
								r=0;
								for(int ii=Rrows1;ii<Rrows2;++ii)
									{
										tripletList.push_back(T(ii,jj,Ahat(r,q)));
								r+=1;							   
							        }
							q+=1;

							}
						 // QTild1.resize(0,0);
						 // QTild2.resize(0,0);
						 // QHat1.resize(0,0);
						 // QHat2.resize(0,0);
						 // WP.resize(0,0);
						 // Ahat.resize(0,0);
						 // Atild.resize(0,0);
						 // Atild2.resize(0,0);

						 
					}

				// Rcout<<j<<std::endl;				
			}

		// Rcout<<"Finished Indexing"<<std::endl;
						//  QTild1.resize(0,0);
						//  QTild2.resize(0,0);
						//  QHat1.resize(0,0);
						//  QHat2.resize(0,0);
						//  WP.resize(0,0);
						//  Ahat.resize(0,0);

						// R1.resize(0,0);
						// R2.resize(0,0);
						// Q1B.resize(0,0);
						// Q2B.resize(0,0);
		 // Rcout<<"Removed matrices from memory"<<std::endl;

		 // 				for(int i=0;i<k;++i){
		 // 					QQ_(i)=0;
		 // 				}
		 // Rcout<<"Removed Qlist from memory"<<std::endl;
	
						
		 // (gdb) b foo
	SparseMatrix<double> W11(RR,RR);
	SparseMatrix<double> W21(RR,k*T1-RR); 
	SparseMatrix<double> W12(RR,k*T1-RR);
	SparseMatrix<double> W22(k*T1-RR,k*T1-RR);

		 
		W11.setFromTriplets(tripletList2.begin(),tripletList2.end());
		W22.setFromTriplets(tripletList.begin(),tripletList.end());
		W12.setFromTriplets(tripletList3.begin(),tripletList3.end());
		W21.setFromTriplets(tripletList4.begin(),tripletList4.end());
		// tripletList.clear();
		// tripletList2.clear();
		// tripletList3.clear();
		// tripletList4.clear();
						
		// Rcout<<"W Constructed As Sparse Matrix"<<std::endl;

	    MatrixXd w1(W11.cols(),W11.rows()+W12.cols());
				// Rcout<<"combine W1"<<std::endl;

		w1.leftCols(W11.rows())=W11.transpose();
		// Rcout<<"combine W1"<<std::endl;
		w1.rightCols(w1.cols()-W11.rows())=W12;
		// Rcout<<"combine W1a"<<std::endl;
		MatrixXd w2(W21.cols(),W21.rows()+W22.rows());
		// MatrixXd w2(W21.cols(),W21.rows());
		w2.setZero();
		// Rcout<<"combine W1b"<<std::endl;
		// Rcout<<W21.rows()<<std::endl;
		// Rcout<<W21.cols()<<std::endl;
		// Rcout<<w2.cols()<<std::endl;
		bool foo= w2.rows()==W21.cols();
		// Rcout<<foo<<std::endl;

		MatrixXd W2a=W21.transpose();
		w2.leftCols(W2a.cols())=W2a;

		// w2.leftCols(W21.cols())=W21.transpose();
// Rcout<<"combine W1c"<<std::endl;
		w2.rightCols(w2.cols()-W21.rows())=W22.transpose();
		// Rcout<<"combine W1d"<<std::endl;
		MatrixXd W(w1.rows()+w2.rows(),w1.cols());
		// W.setZero();
		// Rcout<< w1.rows()<<std::endl;
		// Rcout<<w1.cols()<<std::endl;
		// Rcout<<w2.rows()<<std::endl;
		// THIS IS THE PROBLEM
		W<< w1,
			w2;
		// MatrixXd WTEST=W;
				// Rcout<<"Constructed W as dense matrix"<<std::endl;
				// w1.resize(0,0);
				// w2.resize(0,0);
				
		W = W.colwise().reverse().eval();
		// W.colwise().reverse().eval();
		W.transposeInPlace();
		HouseholderQR<MatrixXd> qr(W);
		// MatrixXd Q1 = qr.householderQ();
		// W.resize(0,0);
		MatrixXd R1a = qr.matrixQR().triangularView<Upper>();
		// R1a.transposeInPlace();
		// Rcout<< "QR Decomposition" <<std::endl;
		// R1a.colwise().reverse().eval();
		// R1a.rowwise().reverse().eval();
		R1a.transposeInPlace();
		R1a=R1a.colwise().reverse().eval();
		R1a=R1a.rowwise().reverse().eval();
		
	// MatrixXd R1ab2=R1a.transpose();
	// MatrixXd	R1ab=R1ab2.colwise().reverse().eval();
	// MatrixXd    R1ac=R1ab.rowwise().reverse().eval();

		// Q1.transposeInPlace();
		// Q1.rowwise().reverse().eval();

		MatrixXd YStarF=Ystar.bottomRows(Ystar.rows()-R2i);

		MatrixXd W22a = R1a.block(R2i,R2i,R1a.rows()-R2i,R1a.cols()-R2i);

		MatrixXd W12a = R1a.rightCols(R1a.cols()-R2i);

		// Rcout<<YStarF.rows()<<std::endl;
	    // Rcout<<W22a.rows()<<std::endl;

		//Issue here
		
		MatrixXd vjstar= W12a*backsolve_rcpp(YStarF,W22a);

		MatrixXd diffmat = Ystar-vjstar;
		MatrixXd BFGLS(RR,1);

		MatrixXd diffmat2= diffmat.topRows(RR);

		MatrixXd Q1AR2=Q1AR.block(0,0,RR,RR);

		// Rcout<<"TEST2"<<std::endl;

		
		BFGLS=Q1AR2.triangularView<Upper>().solve(diffmat2);
		
		
		
		// List Results= List::create(Named("W11")=W11,Named("W12")=W12,Named("W21")=W21,Named("W22")=W22,Named("W2")=W,Named("BFGLS")=BFGLS,Named("W22a")=W22a,Named("W12a")=W12a,Named("diffmat")=vjstar);
		List Results =List::create(Named("BFGLS")=BFGLS,Named("vstar")=vjstar,Named("Ystar")=Ystar,Named("W")=W);
		
		return(Results);
		// return(W);






}


// [[Rcpp::export]]
double onorm(mat S)
{
	return(norm(S,2));

}
