using namespace std;

double Kbar(NumericVector x){
  int n = x.size();
  NumericVector y (n);
  for(int i = 0; i < n; i++){
    if(abs(x[i]) <= 2){
      y[i] = (3.0 / 16.0) * x[i] ^ 2 * pow(min(1,x[i]+1),3) - (9.0 / 16.0) *  pow(x[i], 2) * min(1,x[i]+1) -
      (3.0 / 16.0) * pow(x[i],2) * pow(max(-1,x[i]-1),3) + (9.0 / 16.0) * pow(x[i],2) * max(-1,x[i]-1) -
      (9.0 / 32.0) * x[i] * pow(min(1,x[i]+1),4) + (9.0 / 16.0) * x[i] * pow(min(1,x[i]+1), 2) +
      (9.0 / 32.0) * x[i] * pow(max(-1,x[i]-1), 4) - (9.0 / 16.0) * x[i] * pow(max(-1,x[i]-1), 2) +
      (9.0 / 80.0) * pow(min(1,x[i]+1), 5) - (3.0 / 8.0) * pow(min(1,x[i]+1), 3) + (9.0 / 16.0) * min(1,x[i]+1) -
      (9.0 / 80.0) * pow(max(-1,x[i]-1), 5) + (3.0 / 8.0) * pow(max(-1,x[i]-1), 3) - (9.0 / 16.0) * max(-1,x[i]-1) ;
    }
  }
  return y;
}

//m [[Rcpp:export]]
double fNorm(NumericVector x, double h){
 int n = x.size();
 double par = 0.0;
 dobule constant = 1/(pow(n,2)*h);
 for (int i = 0; i < n; i++){
   for (int j = 0; j < n; j++){
     if( abs((x[i]-x[j])/h) <= 2 ){
       par = par + Kbar((x[i]-x[j])/h);
     }
   }
 }
 double sum = par/constant;
 return sum;
}
