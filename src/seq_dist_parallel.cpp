// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppParallel.h>
using namespace RcppParallel;

IntegerVector seq_to_idx(CharacterVector x, NumericMatrix mtx);

struct PairDist : public Worker {

  const RMatrix<double> mtx;
  const RMatrix<int> seqs_mtx;
  const RVector<double> weight;

  // output matrix to write to
  RMatrix<double> ret;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  PairDist(const NumericMatrix mtx,
           const IntegerMatrix seqs_mtx,
           const NumericVector weight,
           NumericMatrix ret) : mtx(mtx), seqs_mtx(seqs_mtx), weight(weight), ret(ret) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j <= i; j++) {

        RMatrix<int>::Row x = seqs_mtx.row(i);
        RMatrix<int>::Row y = seqs_mtx.row(j);
        double this_score = 0;
        for (int k = 0; k < x.size(); k++) {
          this_score = this_score + mtx(x[k], y[k]) * weight[k];
        }
        ret(i, j) = this_score;
        ret(j, i) = this_score;
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix parallel_pairwise_score(List seqs, NumericMatrix mtx, NumericVector weight) {
  StringVector names = seqs.names();
  int nrow = seqs.size();
  if (nrow == 0) {
    throw std::range_error("No input sequence");
  }
  CharacterVector first_row = seqs[0];
  int ncol = first_row.size();

  // fill in this input mtx
  IntegerMatrix seqs_mtx(nrow, ncol);
  for(int i = 0; i < nrow; i++) {
    IntegerVector this_idx = seq_to_idx(seqs[i], mtx);
    if (this_idx.size() != ncol) {
      throw std::range_error("Sequences differ in length");
    }
    for (int j = 0; j < ncol; j++) {
      seqs_mtx(i, j) = this_idx[j];
    }
  }

  // allocate the matrix we will return
  NumericMatrix ret(nrow, nrow);

  // create the worker
  PairDist pairDist(mtx, seqs_mtx, weight, ret);

  // call it with parallelFor
  parallelFor(0, nrow, pairDist);

  colnames(ret) = names;
  rownames(ret) = names;
  return ret;
}
