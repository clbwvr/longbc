#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
using namespace Rcpp;
using namespace RcppEigen;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double sqrt_double( double x ){ return ::sqrt( x ); }

// [[Rcpp::export]]
NumericVector cppconcat(NumericVector x, NumericVector y) {
  int nx=x.size();
  int ny=y.size();
  // allocate the matrix we will return
  NumericVector z(nx+ny);
  z[Rcpp::Range(0, nx)] = x;
  z[Rcpp::Range(nx, nx+ny)] = y;
  return(z);
}

// [[Rcpp::export]]
NumericVector cppunlist(const List& list){
  int n1 = list.length();
  List l0 = list[1];
  int n2 = l0.length();
  NumericVector ll0 = l0[2];
  int n3 = ll0.length();
  int n4 = (n1 * n2 / 2) * n3;
  // Allocate the vector
  NumericVector output = no_init(n4);
  NumericVector el(n3);
  // Loop and fill
  int index = 0;
  for (int i = 0; i < n1; i++){
    List l1 = list[i];
    for(int j = i+1; j < n2; j++){
      el = l1[j];
      std::copy(el.begin(), el.end(), output.begin() + index);
      index += n3;
    }
  }
  return output;
}

// [[Rcpp::export]]
NumericMatrix cppuppertosym(NumericMatrix m){
  NumericMatrix m1 = m;
  int n = m.rows();
  for(int i=0; i < (n-1);i++){
    for(int j=(i+1); j < n; j++){
      m1(j,i) = m1(i,j);
    }
  }
  return(m1);
}

// [[Rcpp::export]]
double cppnrm2(NumericVector x){
  const Eigen::Map<Eigen::VectorXd> v(as<Eigen::Map<Eigen::VectorXd> >(x));
  double y = v.norm();
  return(y);
}

// [[Rcpp::export]]
NumericVector cppmatvec(NumericMatrix m, NumericVector v){
  const Eigen::Map<Eigen::MatrixXd> tm(as<Eigen::Map<Eigen::MatrixXd> >(m));
  const Eigen::Map<Eigen::VectorXd> tv(as<Eigen::Map<Eigen::VectorXd> >(v));
  Eigen::VectorXd prod = tm*tv;
  SEXP s = Rcpp::wrap(prod);
  NumericVector w(s);
  return(w);
}

// [[Rcpp::export]]
NumericVector cppsparsematvec(Eigen::SparseMatrix<double>& m, NumericVector v){
  const Eigen::Map<Eigen::VectorXd> tv(as<Eigen::Map<Eigen::VectorXd> >(v));
  Eigen::VectorXd y = m*tv;
  SEXP s = Rcpp::wrap(y);
  NumericVector w(s);
  return(w);
}

// [[Rcpp::export]]
NumericMatrix cppmatinv(Eigen::MatrixXd m){
  Eigen::MatrixXd sol = m.inverse();
  SEXP s = Rcpp::wrap(sol);
  NumericMatrix w(s);
  return(s);
}

// [[Rcpp::export]]
NumericVector cppmatvecold(NumericMatrix m, NumericVector v){
  int nRow = m.rows();
  int nCol = m.cols();
  NumericVector ans(nRow);
  double v_j;
  for(int j = 0; j < nCol; j++){
    v_j = v[j];
    for(int i = 0; i < nRow; i++){
      ans[i] += m(i,j) * v_j;
    }
  }
  return(ans);
}

// [[Rcpp::export]]
NumericVector cppl2prox(NumericVector x, double lam){
  if(lam==0){
    return(x);
  }
  else{
    double nrm = cppnrm2(x);
    double t = 1-(lam/nrm);
    double s = 0;
    if(t > 0){
      s = t;
    }
    if(s==0){
      return(0 * x);
    }
    else{
      return(s * x);
    }
  }
}



// [[Rcpp::export]]
List cppsubju(NumericVector a, NumericVector subj, int n){
  List l(n-1);
  for(int i = 0; i < n-1; i++){
    List li(n);
    LogicalVector mask1 = subj==(i+1);
    NumericVector a1 = a[mask1];
    for(int j = (i+1); j < n; j++){
      LogicalVector mask2 = subj==(j+1);
      NumericVector a2 = a[mask2];
      NumericVector ty = a1 - a2;
      li[j] = ty;
    }
    l[i] = li;
  }
  return(l);
}

// [[Rcpp::export]]
List cppfeatu(NumericVector a, NumericVector feat, int p){
  List l(p-1);
  for(int i = 0; i < p-1; i++){
    List li(p);
    LogicalVector mask1 = feat==(i+1);
    NumericVector a1 = a[mask1];
    for(int j = (i+1); j < p; j++){
      LogicalVector mask2 = feat==(j+1);
      NumericVector a2 = a[mask2];
      NumericVector ty = a1 - a2;
      li[j] = ty;
    }
    l[i] = li;
  }
  return(l);
}



// [[Rcpp::export]]
NumericVector cppaupdate(NumericVector a, List blist, List clist, List deltalist, List etalist, NumericVector subj,NumericVector feat, NumericVector Xty, Eigen::SparseMatrix<double>& XtXrhoDeltai, double rho){
  List blist0 = blist[0];
  NumericVector blist01 = blist0[1];
  int na = a.length();
  int n1 = blist.length();
  int n2 = blist0.length();
  LogicalVector mask1;
  LogicalVector mask2;
  NumericVector d1(na);
  for(int i =0; i < n1; i++){
    List blisti = blist[i];
    List deltalisti = deltalist[i];
    mask1 = subj==(i+1);
    for(int j = (i+1); j < n2; j++){
      NumericVector bj = blisti[j];
      NumericVector dj = deltalisti[j];
      NumericVector d1s(na);
      mask2 = subj==(j+1);
      NumericVector x = bj - dj;
      NumericVector nx = -x;
      d1s[mask1] = x;
      d1s[mask2] = nx;
      d1 = d1 + d1s;
    }
  }
  
  List clist0 = clist[0];
  NumericVector clist01 = clist0[1];
  n1 = clist.length();
  n2 = clist0.length();
  LogicalVector mask3;
  LogicalVector mask4;
  NumericVector d2(na);
  for(int i = 0; i < n1; i++){
    List clisti = clist[i];
    List etalisti = etalist[i];
    mask3 = feat==(i+1);
    for(int j = (i+1); j < n2; j++){
      NumericVector cj = clisti[j];
      NumericVector ej = etalisti[j];
      NumericVector d2s(na);
      mask4 = feat==(j+1);
      NumericVector x = cj - ej;
      NumericVector nx = -x;
      d2s[mask3] = x;
      d2s[mask4] = nx;
      d1 = d1 + d2s;
    }
  }
  NumericVector d = d1 + d2;
  NumericVector z = Xty + rho * d;
  NumericVector out = cppsparsematvec(XtXrhoDeltai, z);
  return(out);
}

// [[Rcpp::export]]
List cppbupdate(NumericVector a, List blist, List deltalist, NumericVector subj, List lamlist){
  List blist0 = blist[0];
  NumericVector blist01 = blist0[1];
  int n1 = blist.length();
  int n2 = blist0.length();
  // int n3 = blist01.length();
  List bnewlist(n1);
  LogicalVector mask1;
  LogicalVector mask2;
  NumericVector a1;
  NumericVector a2;
  NumericVector dj;
  double lamlistij;
  // double ty;
  for(int i = 0; i < n1; i++){
    List deltalisti = deltalist[i];
    List bnewlisti(n2);
    List lamlisti = lamlist[i];
    mask1 = subj==(i+1);
    a1 = a[mask1];
    for(int j = (i+1); j < n2; j++){
      mask2 = subj==(j+1);
      a2 = a[mask2];
      dj = deltalisti[j];
      lamlistij = lamlisti[j];
      NumericVector ty = a1 - a2 + dj;
      if(lamlistij==0) bnewlisti[j] = ty;
      else bnewlisti[j] = cppl2prox(ty, lamlistij);
      NumericVector bnewlistij = cppl2prox(ty, lamlistij);
      bnewlisti[j] = bnewlistij;
    }
    bnewlist[i] = bnewlisti;
  }
  return(bnewlist);
}

// [[Rcpp::export]]
List cppcupdate(NumericVector a, List clist, List etalist, NumericVector feat, List lamlist){
  List clist0 = clist[0];
  NumericVector clist01 = clist0[1];
  int n1 = clist.length();
  int n2 = clist0.length();
  List cnewlist(n1);
  LogicalVector mask1;
  LogicalVector mask2;
  NumericVector a1;
  NumericVector a2;
  NumericVector dj;
  double lamlistij;
  // double ty;
  for(int i = 0; i < n1; i++){
    List etalisti = etalist[i];
    List cnewlisti(n2);
    List lamlisti = lamlist[i];
    mask1 = feat==(i+1);
    a1 = a[mask1];
    for(int j = (i+1); j < n2; j++){
      mask2 = feat==(j+1);
      a2 = a[mask2];
      dj = etalisti[j];
      lamlistij = lamlisti[j];
      NumericVector ty = a1 - a2 + dj;
      if(lamlistij==0) cnewlisti[j] = ty;
      else cnewlisti[j] = cppl2prox(ty, lamlistij);
    }
    cnewlist[i] = cnewlisti;
  }
  return(cnewlist);
}

// [[Rcpp::export]]
List cppdupdate(NumericVector a, List blist, List deltalist,NumericVector subj){
  List deltalist0 = deltalist[0];
  NumericVector deltalist01 = deltalist0[1];
  int n1 = deltalist.length();
  int n2 = deltalist0.length();
  List deltanewlist(n1);
  LogicalVector mask1;
  LogicalVector mask2;
  NumericVector a1;
  NumericVector a2;
  NumericVector dj;
  NumericVector bj;
  for(int i = 0; i < n1; i++){
    List deltalisti = deltalist[i];
    List blisti = blist[i];
    List deltanewlisti(n2);
    mask1 = subj==(i+1);
    a1 = a[mask1];
    for(int j = (i+1); j < n2; j++){
      mask2 = subj==(j+1);
      a2 = a[mask2];
      dj = deltalisti[j];
      bj = blisti[j];
      NumericVector deltanewlistij = dj + a1 - a2 - bj;
      deltanewlisti[j] = deltanewlistij;
    }
    deltanewlist[i] = deltanewlisti;
  }
  return(deltanewlist);
}

// [[Rcpp::export]]
List cppeupdate(NumericVector a, List clist, List etalist,NumericVector feat){
  List etalist0 = etalist[0];
  NumericVector etalist01 = etalist0[1];
  int n1 = etalist.length();
  int n2 = etalist0.length();
  List etanewlist(n1);
  LogicalVector mask1;
  LogicalVector mask2;
  NumericVector a1;
  NumericVector a2;
  NumericVector ej;
  NumericVector cj;
  for(int i = 0; i < n1; i++){
    List etalisti = etalist[i];
    List clisti = clist[i];
    List etanewlisti(n2);
    mask1 = feat==(i+1);
    for(int j = (i+1); j < n2; j++){
      mask2 = feat==(j+1);
      a1 = a[mask1];
      a2 = a[mask2];
      ej = etalisti[j];
      cj = clisti[j];
      NumericVector etanewlistij = ej + a1 - a2 - cj;
      etanewlisti[j] = etanewlistij;
    }
    etanewlist[i] = etanewlisti;
  }
  return(etanewlist);
}


// [[Rcpp::export]]
double cppobj(NumericVector a, NumericVector subj, NumericVector feat, double m1, double m2, double lambda, NumericVector w1, NumericVector w2,Eigen::SparseMatrix<double>& X,NumericVector y){
  int n = y.length();
  int n1 = w1.length();
  NumericVector d1(n1);
  int n2 = w2.length();
  NumericVector d2(n2);
  LogicalVector mask1;
  LogicalVector mask2;
  int k=0;
  for(int i = 0; i < (m1-1); i++){
    mask1 = subj==(i+1);
    for(int j = (i+1); j < m1; j++){
      if(w1[k]!=0){
        mask2 = subj==(j+1);
        NumericVector x = a[mask1] - a[mask2];
        d1[k] = cppnrm2(x);
        k++;        
      }
    }
  }
  k= 0;
  for(int i = 0; i < (m2-1); i++){
    mask1 = feat==(i+1);
    for(int j = (i+1); j < m2; j++){
      if(w2[k]!=0){
        mask2 = feat==(j+1);
        NumericVector x = a[mask1] - a[mask2];
        d2[k] = cppnrm2(x);
        k++;        
      }
    }
  }
  double l1 = 0;
  for(int i = 0; i < n1; i++){
    if(w1[i]!=0) l1 += w1[i] * d1[i];
  }
  double l2 = 0;
  for(int i = 0; i < n2; i++){
    if(w2[i]!=0) l2 += w2[i] * d2[i];
  }
  double pen = lambda * (l1 + l2);
  
  NumericVector yhat = cppsparsematvec(X, a);
  double rs = cppnrm2(y - yhat);
  double mse = (rs*rs) / (2*n);
  return(mse + pen);
}


// [[Rcpp::export]]
List cppadmm(NumericVector a, List blist, List clist, List deltalist, List etalist, NumericVector subj, NumericVector feat,int n, int p, Eigen::SparseMatrix<double>& X,NumericVector y, double lambda, NumericVector wsubj, NumericVector wfeat,  List wlrsubj, List wlrfeat, NumericVector Xty, Eigen::SparseMatrix<double>& XtXrhoDeltai, double rho,  int niter, double tolrel, double tolabs, bool trace, bool loud){
  NumericVector losses(niter), olda, oldu, oldz, u, z;
  //NumericVector u1, u2, z1, z2, r, s;
  //double e1, e2, sn1, sn2;
  double loss = 0;
  double oldloss = 0;
  int it;
  for(it=0; it < niter; it++){
    olda = a;
    oldu = u;
    a = cppaupdate(a,blist,clist,deltalist,etalist,subj,feat,Xty,XtXrhoDeltai,rho);
    blist = cppbupdate(a,blist,deltalist,subj,wlrsubj);
    clist = cppcupdate(a,clist,etalist,feat,wlrfeat);
    deltalist = cppdupdate(a,blist,deltalist,subj);
    etalist = cppeupdate(a,clist,etalist,feat);
    
    // u1 = cppunlist(cppsubju(a, subj, n));
    // u2 = cppunlist(cppfeatu(a, feat, p));
    // u = cppconcat(u1, u2); 
    // z2 = cppconcat(cppunlist(deltalist),cppunlist(etalist)); 
    // r = u - z1;
    // s = rho * (oldu - u);
    // e1 = cppnrm2(u);
    // e2 = cppnrm2(z2/rho);
    // sn1 = sqrt_double(u.length());
    // sn2 = sqrt_double(z2.length());
    
    if(trace){
      oldloss = loss;
      loss =  cppobj(a,subj, feat, n, p, lambda,wsubj,wfeat,X,y); 
      losses[it] = loss;
      if (it > 0) {
        // if ((cppnrm2(r) < sn1 * tolabs + tolrel * e1) & (cppnrm2(s) < sn2 * tolabs + tolrel * e2)){
        //if(fabs(loss-oldloss)/oldloss < tolrel){
        //  if(loud) Rprintf("reached convergence\n");
        //  break;
        //}
        if(fabs(loss-oldloss)/oldloss < tolrel){
          if(loud) Rprintf("reached convergence\n");
          break;
        }
      }
    }
    else{
      oldz = z;
      z = cppconcat(cppunlist(blist),cppunlist(clist)); 
      if (it > 0) {
        if(fabs(cppnrm2(z-oldz))/cppnrm2(oldz) < tolrel){
          if(loud) Rprintf("reached convergence\n");
          break;
        }
      }
      losses[it] = fabs(cppnrm2(z-oldz))/cppnrm2(oldz);
    } 
    if(it == (niter - 1)) {
      if(loud) Rprintf("reached max iterations\n");
    }
  }
  NumericVector outlosses(it);
  for(int i = 0; i < it; i++){
    outlosses[i] = losses[i];
  }
  List l(2);
  l[0] = a;
  l[1] = outlosses;
  return(l);
}

