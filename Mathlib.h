#pragma once
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <math.h>
#include <stdexcept>
#include <vector>
#include <iostream>

#ifdef _MATHLIB_USE_DOUBLE
#define _m_d_float double
#define _m_sqrt sqrt
#define _m_pow pow
#else
#define _m_d_float float
#define _m_sqrt sqrtf
#define _m_pow powf
#endif

#ifdef _MATHLIB_USE_LONGINT
#define _m_d_int long int
#else
#define _m_d_int int
#endif

namespace mathlib
{
	namespace functions //Testing
	{
		std::vector<unsigned long> getFactors(unsigned long n) {
			if (n <= 1) {
				return std::vector <unsigned long>(1);
			}
			std::vector<unsigned long> factors;
			unsigned int cfactor = 2;
			while (true) {
				if (n <= cfactor) {
					factors.push_back(n);
					break;
				}
				if (n % cfactor != 0) {
					cfactor++;
				}
				else {
					factors.push_back(cfactor);
					n /= cfactor;
				}
			}
			return factors;
		}

		template<typename T>
		void minmax(T& big, T& small) {
			if (small > big)
			{
				T tmp = big;
				big = small;
				small = tmp;
			}
		}

		int gcd(int a, int b) {
			if (a == b)
				return a;
			minmax(a, b);

			while (b != 0) {
				int r = a % b;
				a = b;
				b = r;
			}
			return a;
		}

		int lcm(int a, int b) {
			return (a * b) / gcd(a, b);
		}

		template <typename T>
		T median(std::vector<T> v, bool isSorted = false) {
			if (!isSorted) {
				std::sort(v.begin(), v.end());
			}
			if (v.size() % 2 == 0) {
				int i = v.size() / 2;
				return (v[i] + v[i + 1]) / 2;
			}
			else {
				return v[((v.size() - 1) / 2) + 1];
			}
		}

		//reflexive
		unsigned long _rfx_factorial(unsigned long n) {
			if (n <= 1)
				return 1;
			return n * _rfx_factorial(n - 1);
		}

		unsigned long factorial(unsigned long n) {
			if (n < 2)
				return 1;
			unsigned long o = 1;
			for (int i = n; i != 1; i--) {
				o *= i;
			}
			return o;
		}

		std::vector<unsigned long> factorialvector(unsigned long n) {
			if (n < 2)
				return std::vector<unsigned long>(1);
			std::vector<unsigned long> o;
			for (int i = n; i != 0; i--) {
				o.push_back(i);
			}
			return o;
		}

		inline unsigned long sigma(unsigned long k) {
			return (k*(k + 1)) / 2;
		}
	}
	namespace numbers
	{
		class Vector2 {
		public:
			_m_d_float x;
			_m_d_float y;

			//Constants
			static const Vector2 zero;
			static const Vector2 one;

			//Operators

			_m_d_float magnitude() const {
				return _m_sqrt(_m_pow(x, 2) + _m_pow(y, 2));
			}

			Vector2 unit() const {
				Vector2 unit(*this);
				unit /= magnitude();
				return unit;
			}

			_m_d_float ScalarProduct(Vector2 v) {
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

			Vector2& operator*=(_m_d_float k) {
				x *= k;
				y *= k;
				return *this;
			}

			Vector2& operator/=(_m_d_float k) {
				x /= k;
				y /= k;
				return *this;
			}

			const Vector2 operator+(const Vector2& vec) const {
				return Vector2(x + vec.x, y + vec.y);
			}

			const Vector2 operator-(const Vector2& vec) const {
				return Vector2(x - vec.x, y - vec.y);
			}

			const Vector2 operator*(float k) const {
				return Vector2(x*k, y*k);
			}

			const Vector2 operator/(float k) const {
				return Vector2(x / k, y / k);
			}

			Vector2() : x(0), y(0) {}
			Vector2(_m_d_float _x, _m_d_float _y) : x(_x), y(_y) {}
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
			}

			Matrix(const Matrix& matrix) : dimX(matrix.dimX), dimY(matrix.dimY) {
				pMatrix = new T[dimX*dimY];

				const int && szByte = dimX * dimY * sizeof(T);
				memcpy(pMatrix, matrix.pMatrix, szByte);
			}

			~Matrix() {
				delete[] pMatrix;
			}

			inline T& At(int row, int column) const {
				if ((column < 0 && column >= dimX) || (row < 0 && row >= dimY))
					throw std::out_of_range("Matrix indexing");

				return *(pMatrix + (column + dimX * row));
			}

			void Transpose()
			{
				T* tMatrix = new T[dimY*dimX];
				for (int x = 0; x < dimX; ++x) {
					for (int y = 0; y < dimY; ++y) {
						*(tMatrix + (y + dimX * x)) = *(pMatrix + (x + dimX * y));
					}
				}
				delete[] pMatrix;
				pMatrix = tMatrix;
				const int newY = dimX;
				dimX = dimY;
				dimY = newY;
			}

			Matrix Transposed() const {
				Matrix t(*this);
				t.Transpose();
				return t;
			}

			const Matrix operator+(const Matrix& mat) const {
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

			const Matrix operator-(const Matrix& mat) const {
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

			const Matrix operator*(const Matrix& mat) const {
				if (dimX != mat.dimY)
					throw std::length_error("Not matching size");

				Matrix r(dimY, mat.dimX);
				for (int row = 0; row < dimY; ++row) {
					for (int column = 0; column < mat.dimX; ++column) {
						T val = 0;
						for (int i = 0; i < dimX; ++i) {
							//val += At(row, i)*mat.At(i, column);
							val += *(pMatrix + (i + dimX * row)) * *(mat.pMatrix + (column + mat.dimX * i));

						}
						//r.At(row, column) = val;
						*(r.pMatrix + (column + r.dimX * row)) = val;
					}
				}

				return r;
			}

			const Matrix operator*(const T& k) const {
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

			Matrix& operator*=(const T& k) {
				for (int x = 0; x < dimX; ++x)
				{
					for (int y = 0; y < dimY; ++y)
					{
						int && c = (x + dimX * y);
						*(pMatrix + c) = *(pMatrix + c) * k;
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

			template<typename Ty>
			friend std::ostream& operator<<(std::ostream& os, const Matrix<Ty>& mat);
		};

		template<typename Ty>
		std::ostream& operator<<(std::ostream& os, const Matrix<Ty>& mat) {
			for (int i = 0; i < mat.dimY; ++i) {
				os << "\n";
				for (int j = 0; j < mat.dimX; ++j) {
					os << *(mat.pMatrix + (j + mat.dimX * i)) << ' ';
				}
				os << "\n";
			}
			return os;
		}

		//Testing
		class RationalNumber {
		public:
			_m_d_int numerator, denominator;

			_m_d_float estimate() const {
				if (denominator == 0)
					throw std::runtime_error("Mathlib divide by zero");
				return (_m_d_float)numerator / denominator;
			}

			void reduct() {
				if (abs(numerator) == 1)
					return;
				const _m_d_int _gcd = functions::gcd(numerator, denominator);
				if (_gcd == 1)
					return;
				numerator /= _gcd;
				denominator /= _gcd;
			}

			const RationalNumber operator*(const RationalNumber& num) const {
				RationalNumber r(numerator*num.numerator, denominator*num.denominator);
				r.reduct();
				return r;
			}

			const RationalNumber operator/(const RationalNumber& num) const {
				RationalNumber r(numerator*num.denominator, denominator*num.numerator);
				r.reduct();
				return r;
			}

			const RationalNumber operator+(const RationalNumber& num) const {
				if (denominator*num.denominator == 0)
					throw std::runtime_error("Mathlib zero in denominator");

				if (denominator == num.denominator)
					return RationalNumber(numerator + num.numerator, denominator);
				const _m_d_int _lcm = functions::lcm(abs(denominator), abs(num.denominator));
				return RationalNumber(numerator*(_lcm/denominator)+ num.numerator * (_lcm / num.denominator), _lcm);
			}

			const RationalNumber operator-(const RationalNumber& num) const {
				if (denominator*num.denominator == 0)
					throw std::runtime_error("Mathlib zero in denominator");

				if (denominator == num.denominator)
					return RationalNumber(numerator - num.numerator, denominator);
				const _m_d_int _lcm = functions::lcm(abs(denominator), abs(num.denominator));
				return RationalNumber(numerator*(_lcm / denominator) - num.numerator * (_lcm / num.denominator), _lcm);
			}

			RationalNumber& operator=(const RationalNumber& num) {
				numerator = num.numerator;
				denominator = num.denominator;
				return *this;
			}

			RationalNumber& operator*=(const RationalNumber& num) {
				numerator = numerator * num.numerator;
				denominator = denominator * num.denominator;
				return *this;
			}

			RationalNumber& operator/=(const RationalNumber& num) {
				numerator = numerator * num.denominator;
				denominator = denominator * num.numerator;
				return *this;
			}

			RationalNumber& operator+=(const RationalNumber& num) {
				if (denominator*num.denominator == 0)
					throw std::runtime_error("Mathlib zero in denominator");

				if (denominator == num.denominator)
					numerator += num.numerator;
				else {
					const _m_d_int _lcm = functions::lcm(abs(denominator), abs(num.denominator));
					numerator = numerator * (_lcm / denominator) + num.numerator * (_lcm / num.denominator);
					denominator = _lcm;
				}
				return *this;
			}

			RationalNumber& operator-=(const RationalNumber& num) {
				if (denominator*num.denominator == 0)
					throw std::runtime_error("Mathlib zero in denominator");
				
				if (denominator == num.denominator)
					numerator -= num.numerator;
				else {
					const _m_d_int _lcm = functions::lcm(abs(denominator), abs(num.denominator));
					numerator = numerator * (_lcm / denominator) - num.numerator * (_lcm / num.denominator);
					denominator = _lcm;
				}
				return *this;
			}

			RationalNumber& operator+=(const _m_d_int& k) {
				numerator += k * denominator;
				return *this;
			}

			RationalNumber& operator-=(const _m_d_int& k) {
				numerator -= k * denominator;
				return *this;
			}

			RationalNumber& operator*=(const _m_d_int& k) {
				numerator *= k;
				return *this;
			}

			RationalNumber& operator/=(const _m_d_int& k) {
				if (k == 0)
					throw std::runtime_error("Mathlib divide by zero");
				if (k < 0)
				{
					denominator *= -k;
					numerator *= -1;
				}
				else
					denominator *= k;
				return *this;
			}

			const RationalNumber operator+(const _m_d_int& k) {
				RationalNumber n(numerator, denominator);
				n.numerator += k * n.denominator;
				return n;
			}

			const RationalNumber operator-(const _m_d_int& k) {
				RationalNumber n(numerator, denominator);
				n.numerator -= k * n.denominator;
				return n;
			}

			const RationalNumber operator*(const _m_d_int& k) {
				RationalNumber n(numerator, denominator);
				n.numerator *= k;
				return n;
			}

			const RationalNumber operator/(const _m_d_int& k) {
				if (k == 0)
					throw std::runtime_error("Mathlib divide by zero");
				RationalNumber n(numerator, denominator);
				if (k < 0)
				{
					n.denominator *= -k;
					n.numerator *= -1;
				}
				else
					n.denominator *= k;
				return n;
			}

			RationalNumber() {
				numerator = 0;
				denominator = 1;
			}

			RationalNumber(_m_d_int number) {
				numerator = number;
				denominator = 1;
			}

			RationalNumber(_m_d_int nu, _m_d_int de) {
				numerator = nu;
				denominator = de;
			}

			friend std::ostream& operator<<(std::ostream& os, const RationalNumber& rn);
		};

		std::ostream& operator<<(std::ostream& os, const RationalNumber& rn) {
			//os << rn.numerator << '/' << rn.denominator;
			os << rn.estimate();
			return os;
		}

		typedef RationalNumber Fraction;

		class RealNumber {

		};
	}
}