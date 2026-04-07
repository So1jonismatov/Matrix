#include "s21_matrix.h"

int s21_determinant(matrix_t *A, double *result) {
  if (!A || !result) return S21_RESULT_ERROR;
  if (A->rows != A->columns) return S21_RESULT_ERROR;

  int n = A->rows;
  double det = 0.0;

  if (n == 1) {
    det = A->matrix[0][0];
  } else if (n == 2) {
    det = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else {
    for (int col = 0; col < n; col++) {
      // Build minor for element (0, col)
      matrix_t minor;
      s21_create_matrix(n - 1, n - 1, &minor);

      int mi = 0;
      for (int i = 1; i < n; i++) {  // skip row 0
        int mj = 0;
        for (int j = 0; j < n; j++) {
          if (j == col) continue;  // skip current col
          minor.matrix[mi][mj] = A->matrix[i][j];
          mj++;
        }
        mi++;
      }

      double minor_det = 0.0;
      s21_determinant(&minor, &minor_det);
      double sign = (col % 2 == 0) ? 1.0 : -1.0;
      det += sign * A->matrix[0][col] * minor_det;

      s21_remove_matrix(&minor);
    }
  }

  *result = det;
  return S21_RESULT_OK;
}
