/** -*- C++ -*-
 * @file Parametric.hpp
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

#ifndef CURVE_PARAMETRIC_HPP
#define CURVE_PARAMETRIC_HPP

#include "Vector.hpp"
#include "Integral.hpp"
#include "Simpson.hpp"
#include <functional>

namespace curve
{

/**
 * @brief Parametric curve base class.
 *
 * A class template for parametric curves.
 */
template <int N, class Real>
class Parametric : public std::unary_function<Real, geom::Vector<N, Real> >
{
public:
	typedef geom::Vector<N, Real> Point;

	/**
	 * @brief Computes C(t).
	 * @param t The parameter t.
	 * @return The computed point.
	 */
	virtual Point operator()( const Real& t ) const = 0;

	/**
	 * @brief Computes C(k)(t).
	 * @param t The parameter t.
	 * @param k The order k.
	 * @return The computed point.
	 */
	virtual Point derivative( const Real& t, int k = 1 ) const = 0;

	/**
	 * @brief Computes the arc length between a and b.
	 * @param a
	 * @param b
	 * @param integral An Integral object.
	 * @return The length between a and b.
	 */
	template <class IntegralType>
	Real length( const Real& a, const Real& b, const IntegralType& integral ) const;

	/**
	 * @brief Computes the arc length between a and b with Simpson integration.
	 * @param a
	 * @param b
	 * @return The length between a and b.
	 */
	Real length( const Real& a, const Real& b ) const;

private:
	/**
	 * @brief Computes the norm of the speed as a functor (std::unary_function).
	 */
	class Speed : public std::unary_function<Real, Real>
	{
	public:
		Speed( const Parametric<N, Real> * p ) : _p( p ) {}
		Real operator()( const Real& t ) const
		{
			Point v;
			v = _p->derivative( t, 1 );
			return v.length();
		}
	private:
		const Parametric<N, Real> * _p;
	};
};

/**
 * @brief Null curve.
 */
template <int N, class Real>
class Null : public Parametric<N, Real>
{
public:
	virtual geom::Vector<N, Real> operator() ( const Real& ) const
	{
		return geom::Vector<N, Real>();
	}

	virtual geom::Vector<N, Real> derivative( const Real&, int ) const
	{
		return geom::Vector<N, Real>();
	}
};

// -----------------------------------------------------------------------------

template <int N, class Real>
template <class IntegralType>
Real Parametric<N, Real>::length( const Real& a, const Real& b, const IntegralType& integral ) const
{
	Speed norm( this );
	return integral( norm, a, b );
}

template <int N, class Real>
Real Parametric<N, Real>::length( const Real& a, const Real& b ) const
{
	return length( a, b, integral::Simpson<Real>() );
}

} // namespace

#endif

