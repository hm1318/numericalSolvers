#pragma once
#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;
class Matrix;
class Vector {
private:
	int* vecValue; 
	int vecElementCount; 

public:

	// parameterized constructor
	Vector(int size) {
		vecElementCount = size;
		cout << __FUNCSIG__ << endl;
		vecValue = new int[vecElementCount]; // "new" because array would be destroyed at end of constructor
	}

	// destructor
	~Vector()
	{
		delete[] vecValue; // remove from heap to prevent memory leak
	}

	// user defined copy constructor
	Vector(const Vector& vec) : Vector(vec.vecElementCount) {
		vecElementCount = vec.vecElementCount ;
		vecValue = new int[vec.vecElementCount];
	}

	// user defined move constructor
	Vector(Vector&& other)
		:vecValue{ other.vecValue },
		vecElementCount{ other.vecElementCount }
	{
		other.vecValue = nullptr;
	}

	// user defined move assignment
	Vector& operator= (Vector&& other) noexcept
	{
		 if (&other == this)
		 {
			 return *this;
		 }
		delete vecValue;

		vecValue = new int[other.vecElementCount];
		*vecValue = *other.vecValue;

		return *this;
	}

	// get the value at a certain index
	int get(int index) const;
	int operator[] (int index) const { // "const" enforces immutability to the compiler of the calling object
		return get(index);
	}

	void set(int index, int value); // set an index of the
	void fill(int value);
	void debug();
	void rAssign(); // random assignment of values to the vector
	int getSize();

	Vector& operator=(const Vector& a) {
		cout << "Copy assignment called" << endl;
		if (this == &a) // if object is already the same as argument then nothing needs to be done
		{
			return *this;
		}

		int* new_vecValue = new int[a.vecElementCount]; // allocate on to heap
		memcpy(new_vecValue, a.vecValue, vecElementCount); // copy over data
		delete[] vecValue; // deallocate 
		vecValue = new_vecValue; // re-assign pointer

		cout << "Copy assignment finished " << endl;
	}
	friend Vector operator* (const Matrix&, const Vector&);

};


void Vector::debug() {

	// sequentially display the elements of the vector from index 0

	for (int i{ 0 }; i < vecElementCount; i++)
	{
		std::cout << *(vecValue + i) << std::endl;
	}
	
}

int Vector::getSize() {
	return vecElementCount;
}
void Vector::rAssign() {
	srand(time(0));

	for (int j{ 0 }; j < vecElementCount; j++)
	{
		*(vecValue + j) = rand() % 1000;
	}
}
void Vector::fill(int value)
{
	for (int i{0}; i < vecElementCount; i++)
	{
		*(vecValue + i) = value ;
	}
}

int Vector::get(int index) const {
	return *(vecValue + index);
}

void Vector::set(int index, int value) {
	*(vecValue + index) = value;
}

#endif
