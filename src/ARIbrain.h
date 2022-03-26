#ifndef ARIBRAIN_H
#define ARIBRAIN_H

int Find(int, std::vector<int>&);
void Union(int, int, std::vector<int>&, std::vector<int>&, std::vector<int>&);
int getCategory(double, double, double, int);
int findConcentration(Rcpp::NumericVector&, double, int, double, int);
Rcpp::IntegerVector findDiscoveries(Rcpp::IntegerVector&, Rcpp::NumericVector&, double, int, double, int, int);

#endif
