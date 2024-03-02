#ifndef CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_
#define CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

class S21Matrix {
 private:
  // Attributes
  int rows_, cols_;
  double** matrix_;
  //  Privat function
  void CreateMatrix();
  void DeleteMatrix();
  void CopyMatrix(const S21Matrix& other);
  void MatrixMinors(const S21Matrix& other, int rows_i, int cols_j) const;
  double MinorsDeterminants() const;

 public:
  // Constructors & Destructor
  S21Matrix();                            // Default constructor
  S21Matrix(int rows, int cols);          // Constructor with parameters
  S21Matrix(const S21Matrix& other);      // Copy
  S21Matrix(S21Matrix&& other) noexcept;  // Move
  ~S21Matrix();                           // Destructor

  // Get/Set
  int GetRows() const;
  void SetRows(int rows);
  int GetCols() const;
  void SetCols(int cols);

  // Function
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose();
  S21Matrix CalcComplements() const;
  double Determinant() const;
  S21Matrix InverseMatrix();

  // Operators
  bool operator==(const S21Matrix& other);
  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(double x) const;
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix& operator=(S21Matrix&& other);
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(double x);
  double& operator()(int rows, int cols);
  double operator()(int rows, int cols) const;
};

#endif  // CPP1_S21_MATRIXPLUS_SRC_S21_MATRIX_OOP_H_