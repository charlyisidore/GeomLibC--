/** -*- C++ -*-
 * @file Frenet.hpp
 * @author Charly LERSTEAU
 * @date 2011-02-20
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

#ifndef FRAME_FRENET_HPP
#define FRAME_FRENET_HPP

#include "Parametric.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"
#include "Frame.hpp"

namespace frame
{

/**
 * @brief Frenet frame generator class.
 *
 * A class to compute the Frenet frame.
 */
template <class Real = float>
class Frenet : public Curve<3, Real>
{
public:
	/**
	 * @brief Empty constructor.
	 */
	Frenet() :
		Curve<3, Real>()
	{
	}

	/**
	 * @brief Constructor from a curve.
	 * @param curve The axial curve.
	 */
	template <class CurveType>
	Frenet( const CurveType & curve ) :
		Curve<3, Real>( curve )
	{
	}

	/**
	 * @brief Computes the frame.
	 * @param t The parameter t along the curve.
	 * @return The computed frame (Matrix 3x3).
	 */
	geom::Matrix<3, 3, Real> operator() ( const Real& t ) const
	{
		// No curve => Null matrix
		if ( this->getCurve() == 0 )
			return geom::Matrix<3, 3, Real>();

		geom::Matrix<3, 3, Real> r;
		geom::Vector<3, Real> v[ 3 ];
		geom::Vector<3, Real> d, a;

		d = this->getCurve()->derivative( t, 1 );
		a = this->getCurve()->derivative( t, 2 );

		v[ 0 ] = d;
		v[ 1 ] = d ^ ( a ^ d );
		v[ 2 ] = d ^ a;

		v[ 0 ].normalize();
		v[ 1 ].normalize();
		v[ 2 ].normalize();

		r.setColumns( v );
		return r;
	}
};

} // namespace

#endif
