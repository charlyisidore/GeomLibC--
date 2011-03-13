/** -*- C++ -*-
 * @file Spline.hpp
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

#ifndef CURVE_SPLINE_HPP
#define CURVE_SPLINE_HPP

#include "Parametric.hpp"
#include <vector>

namespace curve
{

/**
 * @brief Spline curve base class.
 *
 * A class template for spline curves.
 */
template <int N, class Real>
class Spline : public Parametric<N, Real>
{
public:
	typedef geom::Vector<N, Real> Point;

	/**
	 * @brief Constructor for uniform splines with no control points.
	 * @param degree Degree of the curve.
	 */
	Spline( int degree = 3 );

	/**
	 * @brief Constructor for uniform splines.
	 * @param points Control points.
	 * @param degree Degree of the curve.
	 */
	Spline( const std::vector<Point> & points, int degree = 3 );

	/**
	 * @brief Constructor for splines with custom knot vector.
	 * @param points Control points.
	 * @param knots Knot vector.
	 * @param degree Degree of the curve.
	 */
	Spline( const std::vector<Point> & points, const std::vector<Real> & knots, int degree = 3 );

	/**
	 * @brief Copy constructor.
	 */
	Spline( const Spline<N, Real> & curve );

	/**
	 * @brief Returns the degree.
	 */
	int getDegree() const { return _degree; }

	/**
	 * @brief Modifies the degree.
	 */
	void setDegree( int n ) { _degree = n; }

	/**
	 * @brief Checks if the knot vector is uniform.
	 */
	bool isUniform() const { return _uniform; }

	/**
	 * @brief Makes the knot vector uniform or not.
	 */
	void setUniform( bool uniform ) { _uniform = uniform; }

	/**
	 * @brief Checks if the curve is clamped.
	 */
	bool isClamped() const { return _clamped; }

	/**
	 * @brief Makes the curve clamped or not.
	 */
	void setClamped( bool clamped ) { _clamped = clamped; }

	/**
	 * @brief Inserts a control point before specified position.
	 */
	void insertControlPoint( typename std::vector<Point>::iterator position, Point point );

	/**
	 * @brief Appends a control point.
	 */
	void pushControlPoint( Point point );

	/**
	 * @brief Removes a control point.
	 */
	void eraseControlPoint( typename std::vector<Point>::iterator position );

	/**
	 * @brief Edits a control point.
	 */
	void replaceControlPoint( typename std::vector<Point>::iterator position, Point point );

	/**
	 * @brief Returns an array of control points.
	 */
	const std::vector<Point>& controlPoints() const { return _controlPoints; }

	/**
	 * @brief Replaces the array of control points.
	 */
	void setControlPoints( const std::vector<Point>& controlPoints ) { _controlPoints = controlPoints; }

	/**
	 * @brief Returns the knot vector.
	 */
	const std::vector<Real>& knotVector() const { return _knotVector; }

	/**
	 * @brief Computes the total arc length.
	 */
	template <class IntegralType>
	Real length( const IntegralType& integral ) const;

	/**
	 * @brief Computes the total arc length with Simpson integration.
	 */
	Real length() const;

protected:
	std::vector<Point> _controlPoints;
	std::vector<Real> _knotVector;
	int _degree;
	bool _uniform;
	bool _clamped;

	/**
	 * @brief Computes a uniform knot vector.
	 */
	void computeUniformKnotVector();
};

// -----------------------------------------------------------------------------

template <int N, class Real>
Spline<N, Real>::Spline( int degree ) :
	_controlPoints(),
	_knotVector   (),
	_degree ( degree ),
	_uniform( true ),
	_clamped( true )
{
	computeUniformKnotVector();
}

template <int N, class Real>
Spline<N, Real>::Spline( const std::vector<Point> & points, int degree ) :
	_controlPoints( points ),
	_knotVector   (),
	_degree ( degree ),
	_uniform( true ),
	_clamped( true )
{
	computeUniformKnotVector();
}

template <int N, class Real>
Spline<N, Real>::Spline( const std::vector<Point> & points, const std::vector<Real> & knots, int degree ) :
	_controlPoints( points ),
	_knotVector   ( knots ),
	_degree ( degree ),
	_uniform( false ),
	_clamped( true )
{
	computeUniformKnotVector();
}

template <int N, class Real>
Spline<N, Real>::Spline( const Spline<N, Real> & curve ) :
	_controlPoints( curve._controlPoints ),
	_knotVector   ( curve._knotVector ),
	_degree ( curve._degree ),
	_uniform( curve._uniform ),
	_clamped( curve._clamped )
{
}

template <int N, class Real>
void Spline<N, Real>::insertControlPoint( typename std::vector<Point>::iterator position, Point point )
{
	_controlPoints.insert( position, point );
	computeUniformKnotVector();
}

template <int N, class Real>
void Spline<N, Real>::pushControlPoint( Point point )
{
	_controlPoints.push_back( point );
	computeUniformKnotVector();
}

template <int N, class Real>
void Spline<N, Real>::eraseControlPoint( typename std::vector<Point>::iterator position )
{
	_controlPoints.erase( position );
	computeUniformKnotVector();
}

template <int N, class Real>
void Spline<N, Real>::replaceControlPoint( typename std::vector<Point>::iterator position, Point point )
{
	_controlPoints.at( position ) = point;
}

template <int N, class Real>
template <class IntegralType>
Real Spline<N, Real>::length( const IntegralType& integral ) const
{
	Real a, b;
	a = _knotVector.front();
	b = _knotVector.back();
	return Parametric<N, Real>::length( a, b, integral );
}

template <int N, class Real>
Real Spline<N, Real>::length() const
{
	return length( integral::Simpson<Real>() );
}

template <int N, class Real>
void Spline<N, Real>::computeUniformKnotVector()
{
	int i, n, numPoints, numKnots;

	// If non-uniform, let user define its knots.
	if ( !_uniform ) return;

	numPoints = _controlPoints.size();
	numKnots  = numPoints + _degree + 1;
	n = numPoints - _degree;

	_knotVector.clear();

	if ( _clamped )
	{
		for ( i = 0; i <= _degree;  i++ ) _knotVector.push_back( 0. );
		for ( /* */; i < numPoints; i++ ) _knotVector.push_back( (Real)( i - _degree ) / (Real)( n ) );
		for ( /* */; i < numKnots;  i++ ) _knotVector.push_back( 1. );
	}
	else
	{
		for ( i = 0; i < numKnots; i++ ) _knotVector.push_back( (Real)( i ) / (Real)( numKnots - 1 ) );
	}
}

} // namespace

#endif

