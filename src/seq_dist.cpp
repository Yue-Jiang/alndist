#include <Rcpp.h>
#include <map>
using namespace Rcpp;

//' Helper function to convert a character vector to an integer vector of positions in substitution matrix.
//' @param x Sequence as a character vector.
//' @param mtx Substitution matrix as a numeric matrix.
//'
//' @return A integer vector converting input sequence to 0-based index in col/row names of mtx.
//'
//' @export
// [[Rcpp::export]]
IntegerVector seq_to_idx(CharacterVector x, NumericMatrix mtx) {
  List mtx_names = mtx.attr("dimnames");
  CharacterVector rownames = mtx_names[0];
  CharacterVector colnames = mtx_names[1];
  std::map<std::string, int> seq_map;
  int n = rownames.size();
  for (int i = 0; i < n; i++) {
    std::string this_name = as<std::string>(rownames[i]);
    seq_map[this_name] = i;
  }
  int l = x.size();
  IntegerVector ret(l);
  for (int i = 0; i < l; i++) {
    std::string this_letter = as<std::string>(x[i]);
    std::transform(this_letter.begin(), this_letter.end(), this_letter.begin(), toupper);
    ret[i] = seq_map.at(this_letter);
  }
  return ret;
}

//' Calculate scores for two aligned sequences given a substitution matrix.
//' Optionaly weight for each position in the alignment.
//' @param x Sequence 1 as a character vector.
//' @param y Sequence 2 as a character vector. Same length as sequence 1.
//' @param mtx Substitution matrix as a numeric matrix.
//' @param weight A numeric vector specifying weights for each position in sequences 1 and 2.
//' Same length as sequences 1 and 2.
//'
//' @return A single weighted score.
//'
//' @export
// [[Rcpp::export]]
double two_seq_score(CharacterVector x, CharacterVector y,
                     NumericMatrix mtx, NumericVector weight) {
  if (x.size() != y.size()) {
    throw std::range_error("Sequences differ in length");
  }
  if (x.size() != weight.size()) {
    throw std::range_error("Weight not of the same length as sequences");
  }
  IntegerVector xx = seq_to_idx(x, mtx);
  IntegerVector yy = seq_to_idx(y, mtx);
  double ret;
  for (int i = 0; i < x.size(); i++) {
    ret = ret + mtx(xx[i], yy[i]) * weight[i];
  }
  return ret;
}

//' Calculate pairwise scores for aligned sequences given a substitution matrix.
//' Optionaly weight for each position in the alignment.
//' @param seqs A list of character vectors, each represents a sequence in the alignment.
//' @param mtx Substitution matrix as a numeric matrix.
//' @param weight A numeric vector specifying weights for each position in aligned sequences.
//' Same length as sequences 1 and 2.
//'
//' @return A NumericMatrix of pairwise scores.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix pairwise_score(List seqs, NumericMatrix mtx, NumericVector weight) {
  StringVector names = seqs.names();
  int nseqs = seqs.size();
  std::vector< IntegerVector > seqs_idx;
  for(int i = 0; i < nseqs; i++) {
    IntegerVector this_idx = seq_to_idx(seqs[i], mtx);
    seqs_idx.push_back(this_idx);
  }
  NumericMatrix ret(nseqs, nseqs);
  for (int i = 0; i < nseqs; i++) {
    for (int j = i; j < nseqs; j++) {
      IntegerVector x = seqs_idx[i];
      IntegerVector y = seqs_idx[j];
      if (x.size() != y.size()) {
        throw std::range_error("Sequences differ in length");
      }
      if (x.size() != weight.size()) {
        throw std::range_error("Weight not of the same length as sequences");
      }
      double this_score = 0;
      for (int k = 0; k < x.size(); k++) {
        this_score = this_score + mtx(x[k], y[k]) * weight[k];
      }
      ret(i, j) = this_score;
      ret(j, i) = this_score;
    }
  }
  colnames(ret) = names;
  rownames(ret) = names;
  return ret;
}
