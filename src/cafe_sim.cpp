#include <Rcpp.h>
#include <set>

using namespace Rcpp;
using namespace std;


int wrap_coord(int x, int n) {
  while (x < 0) x = x + n;
  while (x >= n) x = x - n;
  return x;
}

// Convert cell index to row position assuming column-major order
int index_to_row(int index, int n_rows) {
  return index % n_rows;
}

// Convert cell index to column position assuming column-major order of elements
int index_to_col(int index, int n_rows) {
  return index / n_rows;
}

// Convert row,col position to cell index assuming column-major order of elements
int rowcol_to_index(int r, int c, int n_rows) {
  return c * n_rows + r;
}

// Get the neighbour indices for a cell identified by row,col (assuming
// column-major order of elements). If diagonal is false the neighbourhood is
// the four nearest orthogonal cells; if true it is the eight nearest cells.
//
std::vector<int> get_rowcol_neighbour_indices(const int r, const int c,
                                              const int n_rows, const int n_cols,
                                              const bool diagonal = false,
                                              const bool wrap = false) {
  std::vector<int> nbrs;

  int rnbr, cnbr;

  // north
  rnbr = r - 1;
  if (wrap && rnbr < 0) {
    rnbr = wrap_coord(rnbr, n_rows);
  }
  if (rnbr >= 0) nbrs.push_back( rowcol_to_index(rnbr, c, n_rows));

  // north-east
  if (diagonal) {
    rnbr = r - 1;
    cnbr = c + 1;

    if (wrap) {
      if (rnbr < 0) rnbr = wrap_coord(rnbr, n_rows);
      if (cnbr >= n_cols) cnbr = wrap_coord(cnbr, n_cols);
    }

    if (rnbr >= 0 && cnbr < n_cols) nbrs.push_back( rowcol_to_index(rnbr, cnbr, n_rows));
  }

  // east
  cnbr = c + 1;
  if (wrap && cnbr >= n_cols) {
    cnbr = wrap_coord(cnbr, n_cols);
  }

  if (cnbr < n_cols) nbrs.push_back( rowcol_to_index(r, cnbr, n_rows));

  // south-east
  if (diagonal) {
    rnbr = r + 1;
    cnbr = c + 1;

    if (wrap) {
      if (rnbr >= n_rows) rnbr = wrap_coord(rnbr, n_rows);
      if (cnbr >= n_cols) cnbr = wrap_coord(cnbr, n_cols);
    }

    if (rnbr < n_rows && cnbr < n_cols) nbrs.push_back( rowcol_to_index(rnbr, cnbr, n_rows));
  }

  // south
  rnbr = r + 1;
  if (wrap && rnbr >= n_rows) {
    rnbr = wrap_coord(rnbr, n_rows);
  }
  if (rnbr < n_rows) nbrs.push_back( rowcol_to_index(rnbr, c, n_rows));

  // south-west
  if (diagonal) {
    rnbr = r + 1;
    cnbr = c - 1;

    if (wrap) {
      if (rnbr >= n_rows) rnbr = wrap_coord(rnbr, n_rows);
      if (cnbr < 0) cnbr = wrap_coord(cnbr, n_cols);
    }

    if (rnbr < n_rows && cnbr >= 0) nbrs.push_back( rowcol_to_index(rnbr, cnbr, n_rows));
  }

  // west
  cnbr = c - 1;
  if (wrap && cnbr < 0) {
    cnbr = wrap_coord(cnbr, n_cols);
  }
  if (cnbr >= 0) nbrs.push_back( rowcol_to_index(r, cnbr, n_rows) );

  // north-west
  if (diagonal) {
    rnbr = r - 1;
    cnbr = c - 1;

    if (wrap) {
      if (rnbr < 0) rnbr = wrap_coord(rnbr, n_rows);
      if (cnbr < 0) cnbr = wrap_coord(cnbr, n_cols);
    }

    if (rnbr >= 0 && cnbr >= 0) nbrs.push_back( rowcol_to_index(rnbr, cnbr, n_rows));
  }

  return nbrs;
}


// Get the neighbour indices for a cell identified by its index (assuming
// column-major order of elements). If diagonal is false the neighbourhood is
// the four nearest orthogonal cells; if true it is the eight nearest cells.
//
std::vector<int> get_neighbour_indices(const int index,
                                       const int n_rows, const int n_cols,
                                       const bool diagonal = false,
                                       const bool wrap = false) {
  int r = index_to_row(index, n_rows);
  int c = index_to_col(index, n_rows);
  return get_rowcol_neighbour_indices(r, c, n_rows, n_cols, diagonal, wrap);
}


// [[Rcpp::export]]
IntegerVector TESTER(int n) {
  IntegerVector y(n);

  const int n_rows = 500;
  const int n_cols = 500;

  for (int i = 0; i < n; i++) {
    int icell = floor(R::runif(0, n_rows * n_cols));
    y[i] = icell;
  }

  return y;
}


// Simulate a fire (if ignition is possible) and return a list with elements:
//   - IntegerMatrix where non-zero values indicate burnt cells
//   - number of cells burnt
//   - index of first ignition cell
//
//' @export
// [[Rcpp::export]]
List doFire(const IntegerMatrix& tsf,
            const Function fn_prob_tsf,
            const bool diagonal = false,
            const int max_ignition_attempts = 10) {

  // Constants for cell state
  const int UNBURNT = 0;
  const int BURNING = 1;

  // Landscape size
  const int n_rows = tsf.rows();
  const int n_cols = tsf.cols();

  // Matrix to hold cell fire state (implicitly filled with 0 = UNBURNT)
  IntegerMatrix landscape(n_rows, n_cols);
  int n_cells_burnt = 0;

  // Attempt to ignite a randomly selected cell
  int ignition_index = -1;

  for (int n_attempts = 0; n_attempts < max_ignition_attempts; n_attempts++) {
    int icell = floor(R::runif(0, n_rows * n_cols));
    double p = Rcpp::as<double>( fn_prob_tsf(tsf[icell]) );

    if (R::runif(0, 1) < p) {
      ignition_index = icell;
      break;
    }
  }

  if (ignition_index >= 0) {
    landscape[ignition_index] = BURNING;
    n_cells_burnt = 1;

    // Initial spread cells are neighbours of the first burning cell
    std::set<int> spreadcells;
    std::set<int>::iterator it_spread;

    std::vector<int> nbrs = get_neighbour_indices(ignition_index, n_rows, n_cols, diagonal, false);

    for (std::vector<int>::iterator itnbr = nbrs.begin(); itnbr != nbrs.end(); ++itnbr) {
      spreadcells.insert(*itnbr);
    }

    while (!spreadcells.empty()) {
      // Select a spread cell at random and remove it from the spread list.
      // NOTE: This implies each spread cell only gets one shot unless it is added
      // once more to the spread list as the neighbour of a new burning cell.
      //
      it_spread = spreadcells.begin();
      int k = floor(R::runif(0, spreadcells.size()));
      for (; k > 0; k--) it_spread++ ;

      int spread_index = *it_spread;
      spreadcells.erase(*it_spread);

      double prob_burn = Rcpp::as<double>( fn_prob_tsf(tsf[spread_index]) );
      double p = R::runif(0, 1);
      if (p < prob_burn) {
        landscape[spread_index] = BURNING;
        n_cells_burnt++ ;

        std::vector<int> nbrs = get_neighbour_indices(spread_index, n_rows, n_cols, diagonal, false);

        for (std::vector<int>::iterator itnbr = nbrs.begin(); itnbr != nbrs.end(); ++itnbr) {
          if (landscape[*itnbr] == UNBURNT) spreadcells.insert(*itnbr);
        }
      }
    }
  }

  return List::create(Named("landscape") = landscape,
                      Named("ncells") = n_cells_burnt,
                      Named("ignition_index") = ignition_index);
}


//' @export
// [[Rcpp::export]]
IntegerMatrix cafeSim(const IntegerMatrix initial_tsf,
                      const int n_times,
                      const Function fn_prob_tsf,
                      const bool diagonal = false,
                      const int max_ignition_attempts = 10,
                      const bool display_progress = true) {

  IntegerMatrix tsf = clone(initial_tsf);

  double progress_time_step = n_times / 10;
  int last_progress_time = 0;

  for (int itime = 0; itime < n_times; itime++) {
    if (itime - last_progress_time > progress_time_step) {
      last_progress_time = itime;
      printf("=");
    }

    List fire_results = doFire(tsf, fn_prob_tsf, diagonal, max_ignition_attempts);
    IntegerMatrix landscape = fire_results["landscape"];

    // Update TSF for burnt and unburnt cells
    for (int index = 0; index < tsf.size(); index++) {
      if (landscape[index] > 0) tsf[index] = 0;
      else tsf[index] = tsf[index] + 1;
    }
  }

  printf("\n");

  return tsf;
}


