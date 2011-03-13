/** -*- C++ -*-
 * @file Simpson.hpp
 * @author Charly LERSTEAU
 * @date 2011-02-18
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

#ifndef INTEGRAL_SIMPSON_HPP
#define INTEGRAL_SIMPSON_HPP

#include "Integral.hpp"
#include <cmath>

namespace integral
{

/**
 * @brief Adaptive Simpson integral class.
 */
template <class Real = float>
class Simpson : public Integral<Real>
{
public:
	Simpson( Real accuracy = 1e-6, int max = 5 ) :
		Integral<Real>(), _accuracy( accuracy ), _maxRecursionDepth( max ) {}

	template <class FunctionType>
	Real operator()( const FunctionType& f, const Real& a, const Real& b ) const;

	Real getAccuracy() const                 { return _accuracy; }
	void setAccuracy( const Real& accuracy ) { _accuracy = accuracy; }
	int getMaxRecursionDepth() const         { return _maxRecursionDepth; }
	void setMaxRecursionDepth( int max )     { _maxRecursionDepth = max; }

protected:
	Real _accuracy;
	int _maxRecursionDepth;

	template <class FunctionType>
	Real aux( const FunctionType& f, const Real& a, const Real& b, const Real& eps, const Real& S, const Real& fa, const Real& fb, const Real& fc, const int& bottom ) const;
};

// -----------------------------------------------------------------------------

template <class Real>
template <class FunctionType>
Real Simpson<Real>::operator()( const FunctionType& f, const Real& a, const Real& b ) const
{
	Real c, h, fa, fb, fc, S;

	c = ( a + b ) / 2.;
	h = b - a;
	fa = f( a );
	fb = f( b );
	fc = f( c );
	S = ( h / 6 ) * ( fa + 4 * fc + fb );
	return aux( f, a, b, _accuracy, S, fa, fb, fc, _maxRecursionDepth );
}

template <class Real>
template <class FunctionType>
Real Simpson<Real>::aux( const FunctionType& f, const Real& a, const Real& b, const Real& eps, const Real& S, const Real& fa, const Real& fb, const Real& fc, const int& bottom ) const
{
	Real c, d, e, h, fd, fe, Sleft, Sright, S2;

	c = ( a + b ) / 2.;
	h = b - a;
	d = ( a + c ) / 2.;
	e = ( c + b ) / 2.;
	fd = f( d );
	fe = f( e );
	Sleft = ( h / 12. ) * ( fa + 4 * fd + fc );
	Sright = ( h / 12. ) * ( fc + 4 * fe + fb );
	S2 = Sleft + Sright;

	if ( bottom <= 0 || fabs( S2 - S ) <= 15. * eps )
		return S2 + ( S2 - S ) / 15.;

	return aux( f, a, c, eps / 2., Sleft, fa, fc, fd, bottom - 1 ) +
		aux( f, c, b, eps / 2., Sright, fc, fb, fe, bottom - 1 );
}

} // namespace

#endif
