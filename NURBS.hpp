/** -*- C++ -*-
 * @file NURBS.hpp
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

#ifndef CURVE_NURBS_HPP
#define CURVE_NURBS_HPP

#include "Spline.hpp"

namespace curve
{

/**
 * @brief NURBS curve base class.
 *
 * A class template for NURBS curves.
 */
template <int N, class Real = float>
class NURBS : public Spline<N, Real>
{
public:
	typedef Spline<N, Real> Parent;
	typedef geom::Vector<N, Real> Point;

	/**
	 * @brief Computes C(t).
	 * @param t The parameter t.
	 * @return The computed point.
	 */
	virtual inline Point operator()( const Real& t ) const;

	/**
	 * @brief Computes C(k)(t).
	 * @param t The parameter t.
	 * @param k The order k.
	 * @return The computed point.
	 */
	virtual inline Point derivative( const Real& t, int k = 1 ) const;

protected:
	// Please see below.
	class Matrix;

	/**
	 * @see Algorithm A2.1, page 68, The NURBS Book (Springer 1997).
	 */
	int findSpan( int n, int p, Real u, const std::vector<Real> & U ) const;

	/**
	 * @param N Array of size p+1
	 * @see Algorithm A2.2, page 70, The NURBS Book (Springer 1997).
	 */
	void basisFuns( int i, Real u, int p, const std::vector<Real> & U, Real * N_ ) const;

	/**
	 * @see Algorithm A4.1, page 124, The NURBS Book (Springer 1997).
	 */
	void curvePoint( int n, int p, const std::vector<Real> & U, const std::vector<Point> & Pw, Real u, Point & C ) const;

	/**
	 * @param ders Two-dimensional array of size n+1 x p+1
	 * @see Algorithm A2.3, page 72, The NURBS Book (Springer 1997).
	 */
	void dersBasisFuns( int i, Real u, int p, int n, const std::vector<Real> & U, Matrix & ders ) const;

	/**
	 * @param CK Array of size d+1
	 * @see Algorithm A3.2, page 93, The NURBS Book (Springer 1997).
	 */
	void curveDerivs( int n, int p, const std::vector<Real> & U, const std::vector<Point> & P, Real u, int d, Point * CK ) const;

	/**
	 * @brief Keeps u value in the good interval.
	 * @param u The value to correct.
	 */
	void adjustParameter( Real & u ) const;

	/**
	 * @brief A very small matrix class to pass Real[][] in function arguments.
	 */
	class Matrix
	{
	public:
		Matrix( int m, int n ) : rows( m ), cols( n )
		{
			data = new Real[ m * n * sizeof( Real ) ];
		}
		~Matrix()
		{
			delete [] data;
		}
		Real * operator[] ( int i )
		{
			return &data[ rows * i ];
		}
		int rows, cols;
		Real * data;
	};
};

// -----------------------------------------------------------------------------

template <int N, class Real>
typename NURBS<N, Real>::Point NURBS<N, Real>::operator() ( const Real& t ) const
{
	Point C;
	Real u = t;
	int n, p;

	p = this->getDegree();
	n = this->controlPoints().size() - 1;

	adjustParameter( u );
	curvePoint( n, p, Parent::_knotVector, Parent::_controlPoints, u, C );
	return C;
}

template <int N, class Real>
typename NURBS<N, Real>::Point NURBS<N, Real>::derivative( const Real& t, int d ) const
{
	Point CK[ d+1 ];
	Real u = t;
	int n, p;

	p = this->getDegree();
	n = this->controlPoints().size() - 1;

	adjustParameter( u );
	curveDerivs( n, p, Parent::_knotVector, Parent::_controlPoints, u, d, CK );
	return CK[ d ];
}

template <int N, class Real>
int NURBS<N, Real>::findSpan( int n, int p, Real u, const std::vector<Real> & U ) const
{
	int low, high, mid;

	// Special case
	if ( u <= U.at( p ) ) return p;
	if ( u >= U.at( n+1 ) ) return n;

	// Do binary search
	low = p;
	high = n + 1;
	mid = ( low + high ) / 2;
	while ( u < U.at( mid ) || u >= U.at( mid + 1 ) )
	{
		if ( u < U.at( mid ) )
			high = mid;
		else
			low = mid;
		mid = ( low + high ) / 2;
	}
	return mid;
}

template <int N, class Real>
void NURBS<N, Real>::basisFuns( int i, Real u, int p, const std::vector<Real> & U, Real * N_ ) const
{
	Real left[ p+1 ], right[ p+1 ];
	Real saved, temp;
	int j, r;

	N_[ 0 ] = 1.;
	for ( j = 1; j <= p; j++ )
	{
		left [ j ] = u - U.at( i+1-j );
		right[ j ] = U.at( i+j ) - u;
		saved = 0.;
		for ( r = 0; r < j; r++ )
		{
			temp = N_[ r ] / ( right[ r+1 ] + left[ j-r ] );
			N_[ r ] = saved + right[ r+1 ] * temp;
			saved = left[ j-r ] * temp;
		}
		N_[ j ] = saved;
	}
}

template <int N, class Real>
void NURBS<N, Real>::curvePoint( int n, int p, const std::vector<Real> & U, const std::vector<Point> & Pw, Real u, Point & C ) const
{
	Real N_[ p+1 ];
	Point Cw, Pi;
	int span, j;

	span = findSpan( n, p, u, U );
	basisFuns( span, u, p, U, N_ );
	for ( j = 0; j <= p; j++ )
	{
		Pi = Pw[ span-p+j ];
		Cw += Pi * Pi.weight() * N_[ j ];
	}
	// Divide by weight
	C = Cw / Cw.weight();
}

template <int N, class Real>
void NURBS<N, Real>::dersBasisFuns( int i, Real u, int p, int n, const std::vector<Real> & U, Matrix & ders ) const
{
	Real ndu[ p+1 ][ p+1 ], a[ 2 ][ p+1 ];
	Real left[ p+1 ], right[ p+1 ];
	Real saved, temp, d;
	int j, k, r, rk, pk, j1, j2, s1, s2;

	ndu[ 0 ][ 0 ] = 1.;
	for ( j = 1; j <= p; j++ )
	{
		left[ j ] = u - U.at( i+1-j );
		right[ j ] = U.at( i+j ) - u;
		saved = 0.;
		for ( r = 0; r < j; r++ )
		{
			// Lower triangle
			ndu[ j ][ r ] = right[ r+1 ] + left[ j-r ];
			temp = ndu[ r ][ j-1 ] / ndu[ j ][ r ];

			// Upper triangle
			ndu[ r ][ j ] = saved + right[ r+1 ] * temp;
			saved = left[ j-r ] * temp;
		}
		ndu[ j ][ j ] = saved;
	}

	// Load the basis functions
	for ( j = 0; j <= p; j++ )
		ders[ 0 ][ j ] = ndu[ j ][ p ];

	// This section computes the derivatives
	for ( r = 0; r <= p; r++ )
	{
		// Alternate rows in array a
		s1 = 0;
		s2 = 1;
		a[ 0 ][ 0 ] = 1.;
		for ( k = 1; k <= n; k++ )
		{
			d = 0.;
			rk = r - k;
			pk = p - k;
			if ( r >= k )
			{
				a[ s2 ][ 0 ] = a[ s1 ][ 0 ] / ndu[ pk+1 ][ rk ];
				d = a[ s2 ][ 0 ] * ndu[ rk ][ pk ];
			}
			if ( rk >= -1 )
				j1 = 1;
			else
				j1 = -rk;
			if ( r-1 <= pk )
				j2 = k - 1;
			else
				j2 = p - r;
			for ( j = j1; j <= j2; j++ )
			{
				a[ s2 ][ j ] = ( a[ s1 ][ j ] - a[ s1 ][ j-1 ] ) / ndu[ pk+1 ][ rk+j ];
				d += a[ s2 ][ j ] * ndu[ rk+j ][ pk ];
			}
			if ( r <= pk )
			{
				a[ s2 ][ k ] = -a[ s1 ][ k-1 ] / ndu[ pk+1 ][ r ];
				d += a[ s2 ][ k ] * ndu[ r ][ pk ];
			}
			ders[ k ][ r ] = d;
			// Switch rows
			j = s1;
			s1 = s2;
			s2 = j;
		}
	}

	// Multiply through by the correct factors
	r = p;
	for ( k = 1; k <= n; k++ )
	{
		for ( j = 0; j <= p; j++ )
			ders[ k ][ j ] *= r;
		r *= ( p - k );
	}
}

template <int N, class Real>
void NURBS<N, Real>::curveDerivs( int n, int p, const std::vector<Real> & U, const std::vector<Point> & P, Real u, int d, Point * CK ) const
{
	Matrix nders( d+1, p+1 );
	int j, k, du, span;

	du = std::min( d, p );

	span = findSpan( n, p, u, U );
	dersBasisFuns( span, u, p, du, U, nders );

	for ( k = 0; k <= du; k++ )
	{
		for ( j = 0; j <= p; j++ )
			CK[ k ] += P[ span-p+j ] * nders[ k ][ j ];
	}
}

template <int N, class Real>
void NURBS<N, Real>::adjustParameter( Real & u ) const
{
	Real minKnot, maxKnot, total;

	minKnot = this->knotVector().front();
	maxKnot = this->knotVector().back();

	if ( this->isClamped() )
	{
		u = std::max( u, minKnot );
		u = std::min( u, maxKnot );
	}
	else
	{
		total = fabs( maxKnot - minKnot );
		while ( u < minKnot ) u += total;
		while ( u >= maxKnot ) u -= total;
	}
}

} // namespace

#endif

