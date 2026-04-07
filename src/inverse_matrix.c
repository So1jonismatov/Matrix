#include <math.h>
#include <stdlib.h>

#include "s21_matrix.h"

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL) return S21_RESULT_ERROR;

  if (A->rows != A->columns) return S21_RESULT_ERROR;

  double det = 0.0;
  if (s21_determinant(A, &det) != S21_RESULT_OK) return S21_RESULT_ERROR;

  if (fabs(det) < 1e-7) return S21_RESULT_ERROR;

  // 3. 1x1 special case: inverse is just 1/element
  if (A->rows == 1) {
    s21_create_matrix(1, 1, result);
    result->matrix[0][0] = 1.0 / A->matrix[0][0];
    return S21_RESULT_OK;
  }

  // 4. Cofactor matrix
  matrix_t cofactors;
  if (s21_calc_complements(A, &cofactors) != S21_RESULT_OK)
    return S21_RESULT_ERROR;

  // 5. Transpose -> adjugate
  matrix_t adj;
  if (s21_transpose(&cofactors, &adj) != S21_RESULT_OK) {
    s21_remove_matrix(&cofactors);
    return S21_RESULT_ERROR;
  }

  // 6. Multiply by 1/det
  s21_mult_number(&adj, 1.0 / det, result);

  // 7. Cleanup
  s21_remove_matrix(&cofactors);
  s21_remove_matrix(&adj);

  return S21_RESULT_OK;
}
