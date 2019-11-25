#pragma once
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <math.h>
#include <stdexcept>

class Vector2 {
public:
	float x;
	float y;

	//Constants
	static const Vector2 zero;
	static const Vector2 one;

	//Operators

	float magnitude() const {
		return sqrtf(powf(x, 2) + powf(y, 2));
	}

	Vector2 unit() const {
		Vector2 unit(*this);
		unit /= magnitude();
		return unit;
	}

	float ScalarProduct(Vector2 v) {
		return x * v.x + y + v.y;
	}

	//Rotate PI/2
	Vector2 perpendicularPositive() const {
		return Vector2(-y, x);
	}

	//Rotate -PI/2
	Vector2 perpendicularNegative() const {
		return Vector2(y, -x);
	}

	Vector2& operator+=(const Vector2& vec) {
		x += vec.x;
		y += vec.y;
		return *this;
	}

	Vector2& operator-=(const Vector2& vec) {
		x -= vec.x;
		y -= vec.y;
		return *this;
	}

	Vector2& operator*=(float k) {
		x *= k;
		y *= k;
		return *this;
	}

	Vector2& operator/=(float k) {
		x /= k;
		y /= k;
		return *this;
	}

	Vector2 operator+(const Vector2& vec) {
		return Vector2(x + vec.x, y + vec.y);
	}

	Vector2 operator-(const Vector2& vec) {
		return Vector2(x - vec.x, y - vec.y);
	}

	Vector2 operator*(float k) {
		return Vector2(x*k, y*k);
	}

	Vector2 operator/(float k) {
		return Vector2(x / k, y / k);
	}

	Vector2() : x(0), y(0) {}
	Vector2(float _x, float _y) : x(_x), y(_y) {}
	Vector2(const Vector2& vector) : x(vector.x), y(vector.y) {}
};
const Vector2 Vector2::zero = Vector2(0, 0);
const Vector2 Vector2::one = Vector2(1, 1);

typedef Vector2* PVector2;

template<typename T>
class Matrix {
private:
	T* pMatrix;
	int dimX, dimY;

public:
	//y = row x = column
	int getDimensionX() const { return dimX; }
	int getDimensionY() const { return dimY; }
public:
	//(row, column)
	Matrix(int row, int column) : dimX(column), dimY(row) {
		pMatrix = new T[dimX*dimY];
		const int && szByte = dimX * dimY * sizeof(T);
	}

	Matrix(const Matrix& matrix) : dimX(matrix.dimX), dimY(matrix.dimY) {
		pMatrix = new T[dimX*dimY];

		const int && szByte = dimX * dimY * sizeof(T);
		memcpy(pMatrix, matrix.pMatrix, szByte);
	}

	~Matrix() {
		delete[] pMatrix;
	}

	T& At(int row, int column) const {
		if ((column < 0 && column >= dimX) || (row < 0 && row >= dimY))
			throw std::out_of_range("Matrix indexing");

		return *(pMatrix + (column + dimX * row));
	}

	Matrix operator+(const Matrix& mat) const {
		if (dimX != mat.dimX || dimY != mat.dimY)
			throw std::out_of_range("Matrixes are not having same size");

		Matrix r(dimY, dimX);
		for (int x = 0; x < dimX; ++x)
		{
			for (int y = 0; y < dimY; ++y)
			{
				int && c = (x + dimX * y);
				*(r.pMatrix + c) = *(pMatrix + c) + *(mat.pMatrix + c);
			}
		}
		return r;
	}

	Matrix operator-(const Matrix& mat) const {
		if (dimX != mat.dimX || dimY != mat.dimY)
			throw std::length_error("Matrixes are not having same size");

		Matrix r(dimY, dimX);
		for (int x = 0; x < dimX; ++x)
		{
			for (int y = 0; y < dimY; ++y)
			{
				int && c = (x + dimX * y);
				*(r.pMatrix + c) = *(pMatrix + c) - *(mat.pMatrix + c);
			}
		}
		return r;
	}

	Matrix operator*(const Matrix& mat) const {
		if (dimX != mat.dimY)
			throw std::length_error("Not matching size");

		Matrix r(dimY, mat.dimX);
		for (int row = 0; row < dimY; ++row) {
			for (int column = 0; column < mat.dimX; ++column) {
				T val = 0;
				for (int i = 0; i < dimX; ++i) {
					val += At(row, i)*mat.At(i, column);
				}
				r.At(row, column) = val;
			}
		}

		return r;
	}

	Matrix operator*(const T& k) const {
		Matrix r(dimY, dimX);
		for (int x = 0; x < dimX; ++x)
		{
			for (int y = 0; y < dimY; ++y)
			{
				int && c = (x + dimX * y);
				*(r.pMatrix + c) = *(pMatrix + c) *k;
			}
		}
		return r;
	}


	Matrix& operator=(const Matrix& mat) {
		delete[] pMatrix;

		dimX = mat.dimX;
		dimY = mat.dimY;
		pMatrix = new T[dimX*dimY];

		const int && szByte = dimX * dimY * sizeof(T);
		memcpy(pMatrix, mat.pMatrix, szByte);

		return *this;
	}

	Matrix& operator+=(const Matrix& mat) {
		if (dimX != mat.dimX || dimY != mat.dimY)
			throw std::out_of_range("Matrixes are not having same size");

		for (int x = 0; x < dimX; ++x)
		{
			for (int y = 0; y < dimY; ++y)
			{
				int && c = (x + dimX * y);
				*(pMatrix + c) = *(pMatrix + c) + *(mat.pMatrix + c);
			}
		}
		return *this;
	}

	Matrix& operator-=(const Matrix& mat) {
		if (dimX != mat.dimX || dimY != mat.dimY)
			throw std::length_error("Matrixes are not having same size");

		for (int x = 0; x < dimX; ++x)
		{
			for (int y = 0; y < dimY; ++y)
			{
				int && c = (x + dimX * y);
				*(pMatrix + c) = *(pMatrix + c) - *(mat.pMatrix + c);
			}
		}
		return *this;
	}

	void Set(T* arr)
	{
		const int && szByte = dimX * dimY * sizeof(T);
		memcpy(pMatrix, arr, szByte);
	}

	static Matrix MakeScalarMatrix(int size, const T& value) {
		Matrix mat(size, size);
		for (int i = 0; i < size; ++i) {
			mat.At(i, i) = value;
		}
		return mat;
	}


	static Matrix MakeTransposedMatrix(const Matrix& source, const T& value) {
		Matrix mat(source.dimX, source.dimY);

		for (int x = 0; x < source.dimX; ++x) {
			for (int y = 0; y < source.dimY; ++y) {
				mat.At(x, y) = source.At(y, x);
			}
		}
		return mat;
	}
};
