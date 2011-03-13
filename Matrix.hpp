/** -*- C++ -*-
 * @file Matrix.hpp
 * @author Charly LERSTEAU
 * @date 2011-03-13
 * 
 * Copyright (c) 2011 Charly LERSTEAU
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#ifndef GEOM_MATRIX_HPP
#define GEOM_MATRIX_HPP

#include "Vector.hpp"
#include <cmath>

namespace geom
{

/**
 * @brief Generic matrix class template.
 *
 * A class template for M x N matrix.
 */
template <int M, int N, class Real = float>
class Matrix
{
public:
	/**
	 * @brief Null matrix constructor.
	 */
	Matrix()
	{
		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				_m[ i ][ j ] = 0;
		}
	}

	/**
	 * @brief Copy constructor.
	 */
	Matrix( const Matrix & matrix )
	{
		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				_m[ i ][ j ] = matrix._m[ i ][ j ];
		}
	}

	/**
	 * @brief Constructor from a C-array of Real.
	 * @param mat A C-array of Real (M*N elements).
	 */
	Matrix( const Real * mat )
	{
		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				_m[ i ][ j ] = mat[ i * N + j ];
		}
	}

	/**
	 * @brief Const row accessor.
	 * @param i 0 is a0x, 1 is a1x, 2 is a2x, ...
	 * @return The i-th row.
	 */
	inline const Real * operator[]( int i ) const
	{
		return _m[ i ];
	}

	/**
	 * @brief Row accessor.
	 * @param i 0 is x, 1 is y, 2 is z, ...
	 * @return The i-th row.
	 */
	inline Real * operator[]( int i )
	{
		return _m[ i ];
	}

	/**
	 * @brief Const cell accessor.
	 * @param i
	 * @param j
	 * @return aij element.
	 */
	inline const Real& operator()( int i, int j ) const
	{
		return _m[ i ][ j ];
	}

	/**
	 * @brief Cell accessor.
	 * @param i
	 * @param j
	 * @return aij element.
	 */
	inline Real& operator()( int i, int j )
	{
		return _m[ i ][ j ];
	}

	/**
	 * @brief Row accessor.
	 * @param i
	 * @return The i-th row.
	 */
	Vector<N, Real> row( int i ) const
	{
		return Vector<N, Real>( _m[ i ] );
	}

	/**
	 * @brief Column accessor.
	 * @param j
	 * @return The j-th column.
	 */
	Vector<M, Real> column( int j ) const
	{
		Vector<M, Real> vect;
		for ( int i = 0; i < M; i++ )
			vect[ i ] = _m[ i ][ j ];
		return vect;
	}

	/**
	 * @brief Replace a row.
	 */
	void setRow( int i, const geom::Vector<N, Real> & v )
	{
		for ( int j = 0; j < N; j++ )
			_m[ i ][ j ] = v[ j ];
	}

	/**
	 * @brief Replace every row.
	 */
	void setRows( const geom::Vector<N, Real> * v )
	{
		for ( int i = 0; i < M; i++ )
		{
			setRow( i, v[ i ] );
		}
	}

	/**
	 * @brief Replace a column.
	 */
	void setColumn( int j, const geom::Vector<N, Real> & v )
	{
		for ( int i = 0; i < M; i++ )
			_m[ i ][ j ] = v[ i ];
	}

	/**
	 * @brief Replace every column.
	 */
	void setColumns( const geom::Vector<N, Real> * v )
	{
		for ( int j = 0; j < N; j++ )
			setColumn( j, v[ j ] );
	}

	/**
	 * @brief Pointer accessor.
	 */
	inline const Real * ptr() const
	{
		return (const Real *)_m;
	}

	/**
	 * @brief Equal to operator.
	 */
	bool operator==( const Matrix & matrix ) const
	{
		bool res = true;

		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				res = res && _m[ i ][ j ] == matrix[ i ][ j ];
		}
		return res;
	}

	/**
	 * @brief Not equal to operator.
	 */
	bool operator!=( const Matrix & matrix ) const
	{
		bool res = false;

		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				res = res || _m[ i ][ j ] != matrix[ i ][ j ];
		}
		return res;
	}

	/**
	 * @brief Addition operator.
	 */
	const Matrix operator+( const Matrix & matrix ) const
	{
		Matrix res;

		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				res[ i ][ j ] = _m[ i ][ j ] + matrix[ i ][ j ];
		}
		return res;
	}

	/**
	 * @brief Addition assignment operator.
	 */
	Matrix & operator+=( const Matrix & matrix )
	{
		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				_m[ i ][ j ] += matrix[ i ][ j ];
		}
		return *this;
	}

	/**
	 * @brief Subtraction operator.
	 */
	const Matrix operator-( const Matrix & matrix ) const
	{
		Matrix res;

		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				res[ i ][ j ] = _m[ i ][ j ] - matrix[ i ][ j ];
		}
		return res;
	}

	/**
	 * @brief Subtraction assignment operator.
	 */
	Matrix & operator-=( const Matrix & matrix )
	{
		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				_m[ i ][ j ] -= matrix[ i ][ j ];
		}
		return *this;
	}

	/**
	 * @brief Minus operator.
	 */
	const Matrix operator-() const
	{
		Matrix res;

		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				res[ i ][ j ] = -_m[ i ][ j ];
		}
		return res;
	}

	/**
	 * @brief Dot product.
	 */
	const Matrix operator*( const Matrix & matrix ) const
	{
		Matrix res;

		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
			{
				res[ i ][ j ] = 0;
				for ( int k = 0; k < N; k++ )
					res[ i ][ j ] += _m[ i ][ k ] * matrix[ k ][ j ];
			}
		}
		return res;
	}

	/**
	 * @brief Multiplication with a scalar.
	 */
	const Matrix operator*( const Real & r ) const
	{
		Matrix res;

		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				res[ i ][ j ] = _m[ i ][ j ] * r;
		}
		return res;
	}

	/**
	 * @brief Multiplication assignment with a scalar.
	 */
	Matrix & operator*=( const Real & r )
	{
		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				_m[ i ][ j ] *= r;
		}
		return *this;
	}

	/**
	 * @brief Division with a scalar.
	 */
	const Matrix operator/( const Real & r ) const
	{
		Matrix res;

		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				res[ i ][ j ] = _m[ i ][ j ] / r;
		}
		return res;
	}

	/**
	 * @brief Division assignment with a scalar.
	 */
	Matrix & operator/=( const Real & r )
	{
		for ( int i = 0; i < M; i++ )
		{
			for ( int j = 0; j < N; j++ )
				_m[ i ][ j ] /= r;
		}
		return *this;
	}

private:
	Real _m[ M ][ N ]; /**< Data. */
};

} // namespace

// Operators

/*
>	         1
>	         .
>	         .
>	         N
>	
>	1 . . N  1
>	.     .  .
>	.     .  .
>	M . . .  M
*/

template <int M, int N, class Real>
const geom::Vector<M, Real> operator*( const geom::Matrix<M, N, Real>& matrix, const geom::Vector<N, Real>& vect )
{
	geom::Vector<M, Real> res;
	for ( int i = 0; i < M; i++ )
		res[ i ] = matrix.row( i ) * vect;
	return res;
}

/*
>	         1 . . N
>	         .     .
>	         .     .
>	         M . . .
>	
>	1 . . M  1 . . N
*/

template <int M, int N, class Real>
const geom::Vector<N, Real> operator*( const geom::Vector<M, Real>& vect, const geom::Matrix<M, N, Real>& matrix )
{
	geom::Vector<N, Real> res;
	for ( int j = 0; j < N; j++ )
		res[ j ] = vect * matrix.column( j );
	return res;
}

// Definitions

namespace geom
{

typedef Matrix<2, 2, float> Matrix2f;
typedef Matrix<3, 3, float> Matrix3f;
typedef Matrix<4, 4, float> Matrix4f;
typedef Matrix<2, 2, double> Matrix2d;
typedef Matrix<3, 3, double> Matrix3d;
typedef Matrix<4, 4, double> Matrix4d;

} // namespace

// Complements

#include <ostream>

/**
 * @brief Displays matrix.
 */
template<int M, int N, class Real>
std::ostream & operator<<( std::ostream & os, const geom::Matrix<M, N, Real> & m )
{
	//os << "(";
	for ( int i = 0; i < M; i++ )
	{
		os << "(";
		for ( int j = 0; j < N; j++ )
			os << m[ i ][ j ] << ( j < N - 1 ? ", " : "" );
		os << ")" << ( i < M - 1 ? ",\n" : "" );
	}
	//os << ")";
	return os;
}

// OpenSceneGraph

#include <osg/Uniform>

/**
 * @brief Converts a geom::Matrix<2, 2> into a osg::Matrix2.
 */
template<class Real>
const geom::Matrix<2, 2, Real> & operator>>( const geom::Matrix<2, 2, Real> & u, osg::Matrix2 & v )
{
	v.set( (const Real *)u.ptr() );
	return u;
}

/**
 * @brief Converts a geom::Matrix<3, 3> into a osg::Matrix3.
 */
template<int M, int N, class Real>
const geom::Matrix<3, 3, Real> & operator>>( const geom::Matrix<3, 3, Real> & u, osg::Matrix3 & v )
{
	v.set( (const Real *)u.ptr() );
	return u;
}

#endif

