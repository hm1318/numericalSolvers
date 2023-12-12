#pragma once
#ifndef CLASS_MATRIX
#define CLASS_MATRIX

#include "vector.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
class Vector;
class Matrix { // only dense matrices, functionality for sparse matrices not supported

private:
	int* matrixElement; 
	int matrixElementCount; 
	int dimensionSize{ 0 }; 

public:
	// parameterized constructor
	explicit Matrix(int size) 
	{
		dimensionSize = size; 
		matrixElementCount = pow(size, 2); 
		matrixElement = new double[matrixElementCount]; 
	}

	// user-defined destructor
	~Matrix() {
		delete[] matrixElement;
	}

	// user-defined copy constructor
	Matrix(const Matrix& sourceMatrix) {
		
		dimensionSize = sourceMatrix.dimensionSize; 
		matrixElementCount = sourceMatrix.matrixElementCount;

		matrixElement = new double[matrixElementCount]; // new matrix pointer on to heap

		for (int i{0}; i < matrixElementCount; i++)
		{
			*(matrixElement + i) = *(sourceMatrix.matrixElement + i);
		}
		
	}

	// move constructor
	Matrix(Matrix&& other) noexcept
		: matrixElement(nullptr)
		, matrixElementCount(0)
		, dimensionSize(0)
	{
		matrixElement = other.matrixElement;
		matrixElementCount = other.matrixElementCount;
		dimensionSize = other.dimensionSize;

		other.matrixElement = nullptr;
	}

	// Positional overload for element retrieval
	int operator() (int, int) const; 

	// Operator overloads for matrix-arithmetic
	friend Vector operator* (const Matrix&, const Vector&);  // allows matrix-vector multiplication
	friend void operator* (Matrix&, const int multiplier); // allows matrix-scalar multiplication
	friend Matrix operator+ (const Matrix&, const Matrix&);  // element-wise addition of compatible matrices
	friend Matrix operator- (const Matrix&, const Matrix&);  // element-wise subtraction of compatible matrices
	friend ostream& operator<< (ostream&, const Matrix&); // stream-operator for displaying matrices

	// matrix-specific functions
	void fill(double x); // populates the calling matrix with x
	void horzCat(Matrix, Matrix); // horizontally concatenates two matrices
	void vertCat(Matrix, Matrix); // vertically concatenates two matrices
	void diagonal(Matrix&, double); // populates the leading diagonal of the calling matrix with x

};
int Matrix::operator() (int r, int c) const {
	if (r < dimensionSize && c < dimensionSize) {
		return *(matrixElement + r * dimensionSize + c);
	}
	throw logic_error("Element doesn't exist. Indices must be in correct range.");
}

Vector operator*(const Matrix& M, const Vector& V) // matrix of size N x M, vector of size M x 1 
{
	Vector pVec(V.vecElementCount); // product vector with size M x 1
	pVec.fill(0);

	if (M.dimensionSize == pVec.getSize()) // validating the size compatibility
	{
		for (int r{ 0 }; r < M.dimensionSize; r++) // iterating through row and column
		{
			for (int c{ 0 }; c < M.dimensionSize; c++) 
			{
				*(pVec.vecValue + r) += M(r, c)*V[c];
				
			}
		}
	}
	return pVec;
}

void operator*(Matrix& A, const int multiplier)
{
	*A.matrixElement *= multiplier;

}

ostream& operator<< (ostream& out, const Matrix& a) { // over
	int i{ 0 };

	for (int j{ 0 }; j < a.matrixElementCount; j++) {
		if (i == 0 && j == 0) {
			out << *(a.matrixElement) << "  "; // first element
		}
		else if (j % a.dimensionSize == 0) {
			++i;
			out << endl << *(a.matrixElement + i * a.dimensionSize + (j - i * a.dimensionSize)) << "  ";
		}
		else {
			out << *(a.matrixElement + i * a.dimensionSize + (j - i * a.dimensionSize)) << "  ";
		}
	}
	out << endl;
	return out;
}

Matrix operator+ (const Matrix& a, const Matrix& b) {

	int sq_size{ a.dimensionSize };

	Matrix c(a.dimensionSize);

	for (int i{ 0 }; i < sq_size; i++) {
		for (int j{ 0 }; j < sq_size; j++) {
			*(c.matrixElement + i * sq_size + j) = a(i, j) + b(i, j);
		}
	}
	return c;
}

Matrix operator- (const Matrix& a, const Matrix& b)
{
	int sq_size{ a.dimensionSize };

	Matrix c(a.dimensionSize);

	for (int i{ 0 }; i < sq_size; i++) {
		for (int j{ 0 }; j < sq_size; j++) {
			*(c.matrixElement + i * sq_size + j) = a(i, j) - b(i, j);
		}
	}
	return c;
}

void Matrix::fill(int x) {

	for (int i{ 0 }; i < matrixElementCount; i++) {
		*(matrixElement + i) = x;
	}
}
void Matrix::Identity(Matrix& I)
{
	int j{0};
	for(int i{0}; i < dimensionSize; i++)
	{
		*(matrixElement + i*dimensionSize + j) = 1;
		++j;
	}
		
}




#endif