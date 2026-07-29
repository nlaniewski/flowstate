// https://github.com/i-cyto/flowCoreUtils (SamGG ; Samuel Granjeaud)

#include <Rcpp.h>
#include <fstream>
#include <stdexcept>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fcs_rcpp_read_data(
    const std::string& file_path,
    long byte_offset,
    long n_row,
    long n_par,
    bool swap
) {
  // validate dimensions before allocating
  if (n_row <= 0) stop("n_row must be positive");
  if (n_par <= 0) stop("n_par must be positive");

  // use long throughout to avoid int overflow on large files
  // (e.g. 5M events x 50 params = 250M values, overflows int32)
  long n_vals = n_row * n_par;

  // open file in binary mode; ifstream closes itself on scope exit (RAII)
  std::ifstream con(file_path, std::ios::binary);
  if (!con.is_open())
    stop("Cannot open file: " + file_path);

  // seek to data segment offset
  con.seekg(byte_offset, std::ios::beg);
  if (con.fail())
    stop("Failed to seek to offset " + std::to_string(byte_offset));

  // read n_vals float32 values into a flat buffer
  std::vector<float> buf(n_vals);
  con.read(reinterpret_cast<char*>(buf.data()), n_vals * sizeof(float));

  // detect partial reads (truncated file or wrong offset)
  if (con.gcount() != static_cast<std::streamsize>(n_vals * sizeof(float)))
    stop("Short read — file may be truncated or byte_offset/n_row/n_par are incorrect. "
           "Expected " + std::to_string(n_vals * sizeof(float)) + " bytes, "
           "got " + std::to_string(con.gcount()));

  // byte-swap if file endian differs from host (big-endian FCS on little-endian host)
  if (swap) {
    for (long k = 0; k < n_vals; k++) {
      char* p = reinterpret_cast<char*>(&buf[k]);
      std::swap(p[0], p[3]);
      std::swap(p[1], p[2]);
    }
  }

  // fill NumericMatrix from row-major float buffer
  // Rcpp matrices are column-major, so we must index (row, col) explicitly —
  // a flat std::copy would transpose the data incorrectly
  NumericMatrix data_mat(n_row, n_par);
  double* out = data_mat.begin();
  const float* in = buf.data();

  for (long row = 0; row < n_row; row++) {
    const float* row_start = in + row * n_par;
    for (long col = 0; col < n_par; col++) {
      // data_mat(row, col) in column-major storage = out[col * n_row + row]
      out[col * n_row + row] = static_cast<double>(row_start[col]);
    }
  }

  return data_mat;
}
