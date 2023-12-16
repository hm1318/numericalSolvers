#pragma once
#ifndef CLASS_MATRIX
#define CLASS_MATRIX

#include "vector.h"
#include <iostream>
#include <iomanip>
#include <cmath>

class Vector; // forward declaration of Vector class for matrix-vector multiplication
class Matrix { 

private:
	double* matrixElement;
	int matrixElementCount{ 0 }; // number of elements
	int rows{ 0 };
	int columns{ 0 };

public:
	// parameterized constructor
	Matrix(int rows, int columns) : matrixElementCount(rows* columns), rows(rows), columns(columns)
	{
		matrixElement = new double[matrixElementCount];
	}

	// destructor
	~Matrix() {
		delete[] matrixElement;
	}

	// user-defined copy constructor
	Matrix(const Matrix& sourceMatrix):rows(sourceMatrix.rows), columns(sourceMatrix.columns)  {
		matrixElementCount = sourceMatrix.matrixElementCount;
		matrixElement = new double[matrixElementCount]; // new mmRatrix pointer on to heap

		for (int i{ 0 }; i < matrixElementCount; i++)
		{
			*(matrixElement + i) = *(sourceMatrix.matrixElement + i);
		}

	}

	// move constructor
	Matrix(Matrix&& other) noexcept
		: matrixElement(other.matrixElement)
		, matrixElementCount(other.matrixElementCount) {
	
		other.matrixElement = nullptr;
	}

	// Positional overload for element retrieval
	double operator() (int, int) const;

	// Operator overloads for matrix-arithmetic
	friend Vector operator* (const Matrix&, const Vector&);  // allows matrix-vector multiplication
	friend Matrix operator* (const Matrix&, const double scalar); // allows element-wise multiplication of a matrix with a real scalar value
	friend Matrix operator+ (const Matrix&, const Matrix&);  // element-wise addition of compatible matrices
	friend Matrix operator- (const Matrix&, const Matrix&);  // element-wise subtraction of compatible matrices
	friend ostream& operator<< (ostream&, const Matrix&); // stream-operator for displaying matrices

	// matrix-specific functions
	void fill(double); // populates the calling matrix with parameter
	void set(const int, const int, const double); // set a position to have a specific value
	void diagonal(double x, int di, bool validity);
	void tridiagonal(double, double, double);
	void Identity(); 
	void zeros(); 
	
	int getSize(int dimension) const; // 0 for x-size, 1 for y-size
	friend Matrix horzCat(const Matrix& left, const Matrix& right); 
	friend Matrix vertCat(const Matrix& top, const Matrix& bottom); 
	
};

void Matrix::set(const int row, const int column, const double value) {
	*(matrixElement + row * getSize(1) + column) = value;
}

double Matrix::operator() (int r, int c) const {

	if (r < rows && c < columns) {
		return *(matrixElement + r * columns + c);
	}
	else {
		cout << "Element doesn't exist." << endl;
		exit(EXIT_FAILURE);
	}
}

int Matrix::getSize(int dimension) const {
	switch (dimension) {
	default:
		cout << "Dimension not/incorrectly specified. Number of rows returned." << endl;
		return rows;
	case 0:
		return rows;
	case 1:
		return columns;
	
	}
}

Matrix operator* (const Matrix& M, const double multiplier) {

	Matrix product(M);

	for (int i{ 0 }; i < M.matrixElementCount; i++) {
		*(product.matrixElement + i) *= multiplier;
	}
	return product;
}

Vector operator* (const Matrix& M, const Vector& V) // matrix of size N x M, vector of size M x 1 
{

	Vector pVec(V.vecElementCount); // product vector with size M x 1
	pVec.fill(0);

	if (M.columns == pVec.getSize()) // validating the size compatibility
	{
		for (int r{ 0 }; r < M.rows; r++) // iterating through row and column
		{
			for (int c{ 0 }; c < M.columns; c++) {
					*(pVec.vecValue + r) += M(r, c) * V[c];	
				
			}
		}
	}
	return pVec;
}

ostream& operator<< (ostream& out, const Matrix& a) { // overload for operator display
	int i{ 0 };

	for (int j{ 0 }; j < a.matrixElementCount; j++) {
		if (i == 0 && j == 0) {
			out << right << setw(8) << *(a.matrixElement) << "  "; // first element
		}
		else if (j % a.columns == 0) {
			++i;
			out << right << setw(8) << endl << *(a.matrixElement + i * a.columns + (j - i * a.columns)) << "  ";
		}
		else {
			out << right << setw(8) << *(a.matrixElement + i * a.columns + (j - i * a.columns)) << "  ";
		}
	}
	out << endl << endl;
	return out;
}

Matrix operator+ (const Matrix& a, const Matrix& b) {

	int r{ a.rows };
	int c{ a.columns };

	Matrix C(r, c);
	if (a.rows == b.rows && a.columns == b.columns) {  //TODO ' exception handling for incompatible matrices'
		for (int i{ 0 }; i < r; i++) {
			for (int j{ 0 }; j < c; j++) {
				*(C.matrixElement + i * c + j) = a(i, j) + b(i, j);  //TODO: how are the elements stored in memory, depending on order is this the most efficient calling routine?
			}
		}
	}
	return C;
}

Matrix operator- (const Matrix& a, const Matrix& b) {
	return a + (b * -1);
}

void Matrix::diagonal(double x, const int di, bool validity) // valid for square matrices (rows = columns)  //TODO: definition of what variables are mathematically
{
	if (validity) {
		if (columns != rows) {
			throw logic_error("Must be a square matrix");
		}

		// when filling non-main diagonals, the effective size reduces hence
		const int restriction{ rows - abs(di) };
		if (di < 0) {

			int i{ abs(di) };
			int j{ 0 };
			while (j < restriction && i < rows)
			{
				*(matrixElement + i * columns + j) = x;
				++i;
				++j;
			}
		}
		if (di > 0) {
			int i{ 0 };
			int j{ abs(di) };
			while (i < restriction && j < columns)
			{
				*(matrixElement + i * columns + j) = x;
				++i;
				++j;
			}
		}
		if (di == 0) {
			int i{ 0 };
			int j{ 0 };
			while (j < columns && i < rows) {
				*(matrixElement + i * columns + j) = x;
				++j;
				++i;
			}
		}
	}
	else {
		zeros();
		try {
			diagonal(x, di, true);
		}
		catch (const logic_error& e) {
			cerr << e.what();
		}
	}
}

void Matrix::Identity()
{
	try {
		diagonal(1.0, 0, false);
	}
	catch (const logic_error& e) {
		cerr << e.what();
	}
}

void Matrix::zeros()
{
	fill(0);
}

void Matrix::fill(double x) {

	for (int i{ 0 }; i < matrixElementCount; i++) {
		*(matrixElement + i) = x;
	}
}

void Matrix::tridiagonal(double a, double b, double c) // equivalent to a banded matrix B[a, b, c]
{
	bool isZero{ false }; // diagonal should only initially accept zero matrices
	zeros();
	isZero = true;
	try {
		// filling the main diagonal
		diagonal(b, 0, isZero);

		// filling the upper diagonal
		diagonal(c, 1, isZero);

		// filling the lower diagonal
		diagonal(a, -1, isZero);
	}
	catch (const logic_error& e) {
		cerr << e.what();
	}
}

Matrix horzCat(const Matrix& left, const Matrix& right) {

	// size compatibility validation
	if (left.rows != right.rows) {
		throw logic_error("Incompatible matrix sizes for concatenation");
	}

	// dimensions of the final matrix F
	int rowF{ left.rows };
	int columnF{ left.columns + right.columns };

	Matrix F(rowF, columnF);

	// recall that the Matrix class is storing values in a 1D array

	// populating the result row by row

	for (int r{ 0 }; r < rowF; r++) {
		for (int c{ 0 }; c < columnF; c++) {
			if (c < left.columns) { // populates the first part of the row with the left matrix elements first
				*(F.matrixElement + r * columnF + c) = left(r, c);
			}
			else {
				/* right matrix does not have as many columns as the concatenated matrix so needs compensation with the current
				* for loop configuration
				*/
				*(F.matrixElement + r * columnF + c) = right(r, (c - left.columns));
			}
		}
	}
	return F;
}

Matrix vertCat(const Matrix& top, const Matrix& bottom) {

	// size compatibility validation
	if (top.columns != bottom.columns) {
		throw logic_error("Incompatible matrix sizes for concatenation");
	}

	// dimensions of final matrix G
	int rowG{ top.rows + bottom.rows };
	int columnG{ top.columns };

	// result of concatenation
	Matrix G(rowG, columnG);

	for (int r{ 0 }; r < rowG; r++) {
		for (int c{ 0 }; c < columnG; c++) {
			if (r < top.rows) { // populating the top part of G first
				*(G.matrixElement + r * columnG + c) = top(r, c);
			}
			else {
				/*the bottom matrix does not have as many columns as the concatenated matrix
				* so we need to compensate the indices
				*/
				*(G.matrixElement + r * columnG + c) = bottom((r - top.rows), c);
			}
		}
	}
	return G;
}



#endif