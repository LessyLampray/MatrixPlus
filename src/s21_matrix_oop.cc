#include "s21_matrix_oop.h"

// Constructors & destructor

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  CreateMatrix();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_), matrix_(nullptr) {
  CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix::~S21Matrix() { DeleteMatrix(); }

// Privat function

void S21Matrix::CreateMatrix() {
  if (rows_ < 0 || cols_ < 0) {
    throw std::invalid_argument("CreateMatrix: incorrect values");
  }
  matrix_ = new double* [rows_] {};
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]{};
  }
}

void S21Matrix::DeleteMatrix() {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
    matrix_ = nullptr;
  }
}

void S21Matrix::CopyMatrix(const S21Matrix& other) {
  if (rows_ != 0 && cols_ != 0) {
    CreateMatrix();
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::MatrixMinors(const S21Matrix& other, int rows_i,
                             int cols_j) const {
  int row_plus = 0;
  for (int i = 0; i < other.rows_; ++i) {
    int cols_plus = 0;
    if (i == rows_i) {
      row_plus = 1;
    }
    for (int j = 0; j < other.cols_; ++j) {
      if (j == cols_j) {
        cols_plus = 1;
      }
      other.matrix_[i][j] = matrix_[i + row_plus][j + cols_plus];
    }
  }
}

double S21Matrix::MinorsDeterminants() const {
  double determ = 0.0;
  S21Matrix tmp_matrix(rows_ - 1, cols_ - 1);
  for (int x = 0; x < cols_; ++x) {
    MatrixMinors(tmp_matrix, 0, x);
    double tmp = tmp_matrix.Determinant();
    determ += pow(-1, x + 2) * matrix_[0][x] * tmp;
  }
  return determ;
}

// Get/Set

int S21Matrix::GetRows() const { return rows_; }

void S21Matrix::SetRows(int new_rows) {
  if (new_rows <= 0) {
    throw std::invalid_argument("SetRows: incorrect values");
  }
  if (new_rows != rows_) {
    S21Matrix copy(*this);
    DeleteMatrix();
    rows_ = new_rows;
    cols_ = copy.cols_;
    CreateMatrix();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (i < copy.rows_) {
          matrix_[i][j] = copy.matrix_[i][j];
        }
      }
    }
  }
}

int S21Matrix::GetCols() const { return cols_; }

void S21Matrix::SetCols(int new_cols) {
  if (new_cols <= 0) {
    throw std::invalid_argument("SetCols: incorrect values");
  }
  if (new_cols != cols_) {
    S21Matrix copy(*this);
    DeleteMatrix();
    rows_ = copy.rows_;
    cols_ = new_cols;
    CreateMatrix();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (j < copy.cols_) {
          matrix_[i][j] = copy.matrix_[i][j];
        }
      }
    }
  }
}

// Function

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("EqMatrix: different size of matrix");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 0.0000001) {
        return false;
      }
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("SumMatrix: different size of matrix");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("SubMatrix: different size of matrix");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument("MulMatrix: incorrect matrix");
  }
  S21Matrix swap(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < other.rows_; k++) {
        swap.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = swap;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix new_matrix(cols_, rows_);
  for (int i = 0; i < cols_; i++) {
    for (int j = 0; j < rows_; j++) {
      new_matrix.matrix_[i][j] = matrix_[j][i];
    }
  }
  return new_matrix;
}

S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::range_error("CalcComplements: Matrix must be square");
  }
  S21Matrix result(rows_, cols_);
  if (rows_ == 1) {
    result.matrix_[0][0] = matrix_[0][0];
  }
  if (rows_ > 1) {
    double determ = 0.0;
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        S21Matrix tmp(rows_ - 1, cols_ - 1);
        MatrixMinors(tmp, i, j);
        determ = tmp.Determinant();
        result.matrix_[i][j] = pow(-1, i + j) * determ;
      }
    }
  }
  return result;
}

double S21Matrix::Determinant() const {
  double result = 0.0;
  if (rows_ != cols_) {
    throw std::range_error("Determinant: Matrix must be square");
  }
  if (rows_ == 1 && cols_ == 1) {
    result = matrix_[0][0];
  } else if (rows_ == 2 && cols_ == 2) {
    result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else if (rows_ > 2 && cols_ > 2) {
    result = MinorsDeterminants();
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  double determ = Determinant();
  if (fabs(determ) < 0.0000001) {
    throw std::invalid_argument("InverseMatrix: Determinant = 0");
  }
  S21Matrix res_matrix(rows_, cols_);
  if (rows_ == 1) {
    res_matrix.matrix_[0][0] = 1 / determ;
  } else {
    S21Matrix alg_dop = CalcComplements().Transpose();
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        res_matrix.matrix_[i][j] = 1 / determ * alg_dop.matrix_[i][j];
      }
    }
  }
  return res_matrix;
}

// Operators

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(double x) const {
  S21Matrix res(*this);
  res.MulNumber(x);
  return res;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    DeleteMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    CopyMatrix(other);
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (this != &other) {
    DeleteMatrix();
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(double x) {
  MulNumber(x);
  return *this;
}

double& S21Matrix::operator()(int rows, int cols) {
  if (rows < 0 || rows >= rows_ || cols < 0 || cols >= cols_) {
    throw std::out_of_range("operator(): Out of range");
  }
  return matrix_[rows][cols];
}

double S21Matrix::operator()(int rows, int cols) const {
  if (rows < 0 || rows >= rows_ || cols < 0 || cols >= cols_) {
    throw std::out_of_range("operator(): Out of range");
  }
  return matrix_[rows][cols];
}