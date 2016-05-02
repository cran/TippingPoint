#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// Modified version of Welch’s t-tests using Rcpp.

// // [[Rcpp::export]]

NumericVector Welch_cpp(NumericVector meanYt_mis, NumericVector meanYc_mis, NumericVector treat, NumericVector Yobs) {

  int Nt = sum(treat==1);
  int Nc = sum(treat==0);
  int L=  meanYt_mis.size();

  int Ntmis = sum(is_na(Yobs) & (treat==1));
  int Ntobs = Nt - Ntmis;

  NumericVector YobsT=Yobs[!(is_na(Yobs)) & (treat==1)];
  NumericVector YobsC=Yobs[!(is_na(Yobs)) & (treat==0)];

  double meanYobsT = mean(YobsT);
  int Ncmis = sum(is_na(Yobs) & (treat==0));
  int Ncobs = Nc - Ncmis;
  double meanYobsC = mean(YobsC);
  NumericVector p_val(L);
  for (int i = 0; i < L; i++) {
    double d = (meanYobsT*Ntobs + meanYt_mis(i)*Ntmis)/Nt - (meanYobsC*Ncobs + meanYc_mis(i)*Ncmis)/Nc;
    double Vart = (var(YobsT)*(Ntobs-1) + Ntobs*Ntmis/Nt*(pow(meanYobsT - meanYt_mis[i],2.0)))/(Ntobs);
    double Varc = (var(YobsC)*(Ncobs-1) + Ncobs*Ncmis/Nc*(pow(meanYobsC - meanYc_mis[i],2.0)))/(Ncobs);
    double Welch_DF = pow((Vart/Nt + Varc/Nc),2.0)/(pow((Vart/Nt),2.0)/(Ntobs)+pow((Varc/Nc),2.0)/(Ncobs));
    double t_stat = d/sqrt(Vart/Nt + Varc/Nc);
    double pt_stat = ::Rf_pt((-abs(t_stat)), Welch_DF,1,0);
    p_val[i] = 2*pt_stat;
  }
  return(p_val);
}

