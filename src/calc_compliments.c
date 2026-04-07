#include "s21_matrix.h"
// buni tugirla negr

int calc_minor(matrix_t *A, matrix_t *result) {
  if (!A || !result) return S21_RESULT_ERROR;
  if (A->rows != A->columns || A->rows < 2) return S21_RESULT_ERROR;

  int n = A->rows;
  s21_create_matrix(n, n, result);

  matrix_t temp;
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      s21_create_matrix(n - 1, n - 1, &temp);

      int mi = 0;
      for (int i = 0; i < n; i++) {
        if (i == row) continue;
        int mj = 0;
        for (int j = 0; j < n; j++) {
          if (j == col) continue;
          temp.matrix[mi][mj] = A->matrix[i][j];
          mj++;
        }
        mi++;
      }

      double det_minor;
      s21_determinant(&temp, &det_minor);
      result->matrix[row][col] = det_minor;

      s21_remove_matrix(&temp);  // free the submatrix
    }
  }

  return S21_RESULT_OK;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (!A || !result) return S21_RESULT_ERROR;
  if (A->rows != A->columns || A->rows < 2) return S21_RESULT_ERROR;

  matrix_t minors;
  if (calc_minor(A, &minors) != S21_RESULT_OK) {
    return S21_RESULT_ERROR;
  }

  if (s21_create_matrix(A->rows, A->columns, result) != S21_RESULT_OK) {
    s21_remove_matrix(&minors);
    return S21_RESULT_ERROR;
  }

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
      result->matrix[i][j] = minors.matrix[i][j] * sign;
    }
  }

  s21_remove_matrix(&minors);
  return S21_RESULT_OK;
}
