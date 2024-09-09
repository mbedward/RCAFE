#include <Rcpp.h>
#include <set>

using namespace Rcpp;
using namespace std;


// Constants for regime types
const int REGIME_TYPE_WILDFIRE = 1;
const int REGIME_TYPE_PRESCRIBED_FIRE = 2;

// Accessor functions for R to use to set active bindings to the constants
// [[Rcpp::export]]
int get_wildfire_type(){ return REGIME_TYPE_WILDFIRE; }
int get_prescribed_fire_type() { return REGIME_TYPE_PRESCRIBED_FIRE; }


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


List do_wildfire(const IntegerMatrix& tsf, const List& regime) {
  // Constants for cell state
  const int UNBURNT = 0;
  const int BURNING = 1;

  // Landscape size
  const int n_rows = tsf.rows();
  const int n_cols = tsf.cols();

  // Matrix to hold cell fire state (implicitly filled with 0 = UNBURNT)
  IntegerMatrix landscape(n_rows, n_cols);
  int n_cells_burnt = 0;

  // Regime parameters
  Function fn_prob_tsf = regime["fn_prob_tsf"];
  int max_ignition_attempts = regime["max_ignition_attempts"];
  bool diagonal = regime["diagonal"];

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


List do_prescribed_fire(const IntegerMatrix& tsf, const List& regime) {
  // Constants for cell state
  const int UNBURNT = 0;
  const int BURNING = 1;

  // Landscape size
  const int n_rows = tsf.rows();
  const int n_cols = tsf.cols();

  // Matrix to hold cell fire state (implicitly filled with 0 = UNBURNT)
  IntegerMatrix landscape(n_rows, n_cols);
  int n_cells_burnt = 0;

  // Regime parameters
  Function fn_prob_tsf = regime["fn_prob_tsf"];
  double target_prop_landscape = regime["prop_landscape"];
  int target_num_cells = round(target_prop_landscape * n_rows * n_cols);

  // Candidate cell indices (zero-based)
  IntegerVector cell_indices( Rcpp::sample(n_rows * n_cols, target_num_cells, false, R_NilValue, false) );

  for (int i = 0; i < target_num_cells; i++) {
    int cell_index = cell_indices(i);
    double prob_burn = Rcpp::as<double>( fn_prob_tsf( tsf[cell_index] ) );
    double p = R::runif(0, 1);

    if (p < prob_burn) {
      landscape[cell_index] = BURNING;
      n_cells_burnt++ ;
    }
  }

  return List::create(Named("landscape") = landscape,
                      Named("ncells") = n_cells_burnt);
}


// Simulate a fire (if ignition is possible) and return a list with elements:
//   - IntegerMatrix where non-zero values indicate burnt cells
//   - number of cells burnt
//   - index of first ignition cell
//
//' @export
// [[Rcpp::export]]
List doFire(const IntegerMatrix& tsf, const List& regime) {
  // Get the fire regime parameters
  int regime_type = regime["regime_type"];  // 1 is wildfire, 2 is prescribed fire
  if (regime_type < 1 || regime_type > 2) stop("Unknown regime type constant: %d", regime_type);

  String regime_name = regime["name"];

  List fire_res = NULL;

  if (regime_type == REGIME_TYPE_WILDFIRE) {
    fire_res = do_wildfire(tsf, regime);

  } else if (regime_type == REGIME_TYPE_PRESCRIBED_FIRE) {
    fire_res = do_prescribed_fire(tsf, regime);

  } else {
    // shouldn't be here because of earlier check, but just in case...
    stop("Unknown regime type constant %d", regime_type);
  }

  return fire_res;
}
