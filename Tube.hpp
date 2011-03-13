/** -*- C++ -*-
 * @file Tube.hpp
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

#ifndef SURFACE_TUBE_HPP
#define SURFACE_TUBE_HPP

#include "Parametric.hpp"
#include "Vector.hpp"
#include "Frame.hpp"
#include <functional>
#include <cmath>

namespace surface
{

/**
 * @brief Tube surface class template.
 *
 * A class template for tube.
 */
template <class Real = float>
class Tube : public std::binary_function<Real, Real, geom::Vector<3, Real> >
{
public:
	/**
	 * @brief Constructor from a curve.
	 * @param curve The axial curve.
	 * @param radius Radius of the tube.
	 */
	template <class CurveType, class FrameType>
	Tube( const CurveType& curve, const FrameType&, Real radius = 1. ) :
		_frame ( new FrameType( curve ) ),
		_radius( radius )
	{
	}

	/**
	 * @brief Deletes the allocated frame.
	 */
	~Tube()
	{
		delete _frame;
	}

	inline const curve::Parametric<3, Real> * getCurve() const { return _frame->getCurve(); }
	inline const frame::Curve<3, Real> * getFrame() const { return _frame; }
	inline Real getRadius() const { return _radius; }

	/**
	 * @param radius The new radius.
	 */
	void setRadius( const Real& radius )
	{
		_radius = radius;
	}

	/**
	 * @brief Computes S(t,u).
	 * @param t The parameter t along the curve, t0 <= t <= tn.
	 * @param u The angle u along the ring, 0 <= u <= 2*pi.
	 * @return The computed point.
	 */
	geom::Vector<3, Real> operator() ( const Real& t, const Real& u ) const
	{
		geom::Matrix<3, 3, Real> mTNB;

		mTNB = (*_frame)( t );

		geom::Vector<3, Real> vP( (*_frame->getCurve())( t ) );
		geom::Vector<3, Real> vN( mTNB[ 1 ] );
		geom::Vector<3, Real> vB( mTNB[ 2 ] );

		return vP + vN * _radius * cos( u ) + vB * _radius * sin( u );
	}

protected:
	const frame::Curve<3, Real> * _frame;
	Real _radius;
};

} // namespace

#endif
