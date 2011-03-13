/** -*- C++ -*-
 * @file Vector.hpp
 * @author Charly LERSTEAU
 * @date 2011-03-12
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

#ifndef GEOM_VECTOR_HPP
#define GEOM_VECTOR_HPP

#include <cmath>

namespace geom
{

/**
 * @brief Generic vector/point class template.
 *
 * A class template for N-dimension vectors.
 */
template<int N, class Real = float>
class Vector
{
public:
	/**
	 * @brief Origin constructor.
	 */
	Vector()
	{
		for ( int i = 0; i < N; i++ )
			_v[ i ] = 0.;
		_w = 1.;
	}

	/**
	 * @brief Copy constructor.
	 */
	Vector( const Vector & vect )
	{
		for ( int i = 0; i < N; i++ )
			_v[ i ] = vect[ i ];
		_w = vect.weight();
	}

	/**
	 * @brief Constructor from a C-array of Real.
	 * @param vect A C-array of Real.
	 * @param w The weight (default : 1).
	 */
	Vector( const Real * vect, const Real & w = 1. )
	{
		for ( int i = 0; i < N; i++ )
			_v[ i ] = vect[ i ];
		_w = w;
	}

	/**
	 * @brief Const coordinate accessor.
	 * @param i 0 is x, 1 is y, 2 is z, ...
	 * @return The i-th coordinate.
	 */
	inline const Real & operator[]( int i ) const
	{
		return _v[ i ];
	}

	/**
	 * @brief Coordinate accessor.
	 * @param i 0 is x, 1 is y, 2 is z, ...
	 * @return The i-th coordinate.
	 */
	inline Real & operator[]( int i )
	{
		return _v[ i ];
	}

	/**
	 * @brief Equal to operator.
	 */
	inline bool operator==( const Vector & vect ) const
	{
		bool res = true;

		for ( int i = 0; i < N; i++ )
			res = res && _v[ i ] == vect[ i ];
		return res;
	}

	/**
	 * @brief Not equal to operator.
	 */
	inline bool operator!=( const Vector & vect ) const
	{
		bool res = false;

		for ( int i = 0; i < N; i++ )
			res = res || _v[ i ] != vect[ i ];
		return res;
	}

	/**
	 * @brief Addition operator.
	 */
	inline const Vector operator+( const Vector & vect ) const
	{
		Vector res;

		for ( int i = 0; i < N; i++ )
			res[ i ] = _v[ i ] + vect[ i ];
		return res;
	}

	/**
	 * @brief Addition assignment operator.
	 */
	inline Vector & operator+=( const Vector & vect )
	{
		for ( int i = 0; i < N; i++ )
			_v[ i ] += vect[ i ];
		return *this;
	}

	/**
	 * @brief Subtraction operator.
	 */
	inline const Vector operator-( const Vector & vect ) const
	{
		Vector res;

		for ( int i = 0; i < N; i++ )
			res[ i ] = _v[ i ] - vect[ i ];
		return res;
	}

	/**
	 * @brief Subtraction assignment operator.
	 */
	inline Vector & operator-=( const Vector & vect )
	{
		Vector res;

		for ( int i = 0; i < N; i++ )
			_v[ i ] -= vect[ i ];
		return *this;
	}

	/**
	 * @brief Minus operator.
	 */
	inline const Vector operator-() const
	{
		Vector res;

		for ( int i = 0; i < N; i++ )
			res[ i ] = -_v[ i ];
		return res;
	}

	/**
	 * @brief Dot product.
	 */
	inline Real operator*( const Vector & vect ) const
	{
		Real r = 0.;

		for ( int i = 0; i < N; i++ )
			r += _v[ i ] * vect[ i ];
		return r;
	}

	/**
	 * @brief Multiplication with a scalar.
	 */
	inline const Vector operator*( const Real & r ) const
	{
		Vector vect;

		for ( int i = 0; i < N; i++ )
			vect[ i ] = _v[ i ] * r;
		return vect;
	}

	/**
	 * @brief Multiplication assignment with a scalar.
	 */
	inline Vector & operator*=( const Real & r )
	{
		for ( int i = 0; i < N; i++ )
			_v[ i ] *= r;
		return *this;
	}

	/**
	 * @brief Division with a scalar.
	 */
	inline const Vector operator/( const Real & r ) const
	{
		Vector res;

		for ( int i = 0; i < N; i++ )
			res[ i ] = _v[ i ] / r;
		return res;
	}

	/**
	 * @brief Division assignment with a scalar.
	 */
	inline Vector & operator/=( const Real & r )
	{
		for ( int i = 0; i < N; i++ )
			_v[ i ] /= r;
		return *this;
	}

	/**
	 * @brief Returns the const weight.
	 * @return The weight (const).
	 */
	inline const Real & weight() const { return _w; }

	/**
	 * @brief Returns the weight.
	 * @return The weight.
	 */
	inline Real & weight() { return _w; }

	/**
	 * @brief Returns the length.
	 * @return The length.
	 */
	inline Real length() const { return sqrt( (*this) * (*this) ); }

	/**
	 * @brief Makes an unit vector.
	 */
	inline void normalize()
	{
		Real len;
		int i;

		len = length();
		for ( i = 0; i < N; i++ )
			_v[ i ] /= len;
	}

private:
	Real _v[ N ]; /**< Coordinates data. */
	Real _w;      /**< Weight data. */
};

} // namespace

/**
 * @brief Cross product in 3D.
 */
template<class Real>
inline const geom::Vector<3, Real> operator^( const geom::Vector<3, Real> & u, const geom::Vector<3, Real> & v )
{
	geom::Vector<3, Real> w;

	w[ 0 ] = u[ 1 ] * v[ 2 ] - u[ 2 ] * v[ 1 ];
	w[ 1 ] = u[ 2 ] * v[ 0 ] - u[ 0 ] * v[ 2 ];
	w[ 2 ] = u[ 0 ] * v[ 1 ] - u[ 1 ] * v[ 0 ];
	return w;
}

/**
 * @brief Cross product in 7D.
 */
template<class Real>
inline const geom::Vector<7, Real> operator^( const geom::Vector<7, Real> & u, const geom::Vector<7, Real> & v )
{
	geom::Vector<7, Real> w;

	w[ 0 ] = u[ 1 ] * v[ 3 ] - u[ 3 ] * v[ 1 ] + u[ 2 ] * v[ 6 ] - u[ 6 ] * v[ 2 ] + u[ 4 ] * v[ 5 ] - u[ 5 ] * v[ 4 ];
	w[ 1 ] = u[ 2 ] * v[ 4 ] - u[ 4 ] * v[ 2 ] + u[ 3 ] * v[ 0 ] - u[ 0 ] * v[ 3 ] + u[ 5 ] * v[ 6 ] - u[ 6 ] * v[ 5 ];
	w[ 2 ] = u[ 3 ] * v[ 5 ] - u[ 5 ] * v[ 3 ] + u[ 4 ] * v[ 1 ] - u[ 1 ] * v[ 4 ] + u[ 6 ] * v[ 0 ] - u[ 0 ] * v[ 6 ];
	w[ 3 ] = u[ 4 ] * v[ 6 ] - u[ 6 ] * v[ 4 ] + u[ 5 ] * v[ 2 ] - u[ 2 ] * v[ 5 ] + u[ 0 ] * v[ 1 ] - u[ 1 ] * v[ 0 ];
	w[ 4 ] = u[ 5 ] * v[ 0 ] - u[ 0 ] * v[ 5 ] + u[ 6 ] * v[ 3 ] - u[ 3 ] * v[ 6 ] + u[ 1 ] * v[ 2 ] - u[ 2 ] * v[ 1 ];
	w[ 5 ] = u[ 6 ] * v[ 1 ] - u[ 1 ] * v[ 6 ] + u[ 0 ] * v[ 4 ] - u[ 4 ] * v[ 0 ] + u[ 2 ] * v[ 3 ] - u[ 3 ] * v[ 2 ];
	w[ 6 ] = u[ 0 ] * v[ 2 ] - u[ 2 ] * v[ 0 ] + u[ 1 ] * v[ 5 ] - u[ 5 ] * v[ 1 ] + u[ 3 ] * v[ 4 ] - u[ 4 ] * v[ 3 ];
	return w;
}

// Definitions

namespace geom
{

template <class Real = float>
class Vector2 : public Vector<2, Real>
{
public:
	/**
	 * @brief Constructor from a C-array of Real.
	 */
	Vector2( const Real * vect, const Real & w = 1. ) :
		Vector<2, Real>( vect, w )
	{
	}

	/**
	 * @brief Constructor from two coordinates.
	 */
	Vector2( const Real & x, const Real & y, const Real & w = 1. ) :
		Vector<2, Real>()
	{
		(*this)[ 0 ] = x;
		(*this)[ 1 ] = y;
		(*this).weight() = w;
	}
};

template <class Real = float>
class Vector3 : public Vector<3, Real>
{
public:
	/**
	 * @brief Constructor from a C-array of Real.
	 */
	Vector3( const Real * vect, const Real & w = 1. ) :
		Vector<3, Real>( vect, w )
	{
	}

	/**
	 * @brief Constructor from three coordinates.
	 */
	Vector3( const Real & x, const Real & y, const Real & z, const Real & w = 1. ) :
		Vector<3, Real>()
	{
		(*this)[ 0 ] = x;
		(*this)[ 1 ] = y;
		(*this)[ 2 ] = z;
		(*this).weight() = w;
	}
};

template <class Real = float>
class Vector4 : public Vector<4, Real>
{
public:
	/**
	 * @brief Constructor from a C-array of Real.
	 */
	Vector4( const Real * vect, const Real & w = 1. ) :
		Vector<4, Real>( vect, w )
	{
	}

	/**
	 * @brief Constructor from four coordinates.
	 */
	Vector4( const Real & x, const Real & y, const Real & z, const Real & t, const Real & w = 1. ) :
		Vector<4, Real>()
	{
		(*this)[ 0 ] = x;
		(*this)[ 1 ] = y;
		(*this)[ 2 ] = z;
		(*this)[ 3 ] = t;
		(*this).weight() = w;
	}
};

typedef Vector2<float> Vector2f;
typedef Vector2<double> Vector2d;
typedef Vector3<float> Vector3f;
typedef Vector3<double> Vector3d;
typedef Vector4<float> Vector4f;
typedef Vector4<double> Vector4d;

} // namespace

// Complements

#include <ostream>

/**
 * @brief Displays coordinates.
 */
template<int N, class Real>
std::ostream & operator<<( std::ostream & os, const geom::Vector<N, Real> & v )
{
	os << "(";
	for ( int i = 0; i < N; i++ )
		os << v[ i ] << ", ";
	os << v.weight() << ")";
	return os;
}

// OpenSceneGraph

#include <osg/Vec2>
#include <osg/Vec3>

/**
 * @brief Converts a geom::Vector<2> into a osg::Vec2.
 */
template<class Real>
const geom::Vector<2, Real> & operator>>( const geom::Vector<2, Real> & u, osg::Vec2 & v )
{
	v.set( u[ 0 ], u[ 1 ] );
	return u;
}

/**
 * @brief Converts a geom::Vector<3> into a osg::Vec3.
 */
template<class Real>
const geom::Vector<3, Real> & operator>>( const geom::Vector<3, Real> & u, osg::Vec3 & v )
{
	v.set( u[ 0 ], u[ 1 ], u[ 2 ] );
	return u;
}

#endif

