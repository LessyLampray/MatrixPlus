#ifndef CPP1_S21_MATRIXPLUS_SRC_TEST_H_
#define CPP1_S21_MATRIXPLUS_SRC_TEST_H_

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

#include "s21_matrix_oop.h"

class S21MatrixTest : public testing::Test {
 protected:
  S21Matrix *matrix_1x1;
  S21Matrix *matrix_2x3;
  S21Matrix *matrix_5x5;
  S21Matrix *matrix_12x21;
  S21Matrix *matrix_21x21;

  void SetUp();
  void TearDown();
  void FillMatrixWithRandomDouble(S21Matrix &matrix);
};

#endif  // CPP1_S21_MATRIXPLUS_SRC_TEST_H_