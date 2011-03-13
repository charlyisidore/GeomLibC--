/** -*- C++ -*-
 * @file Frame.hpp
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

#ifndef FRAME_HPP
#define FRAME_HPP

#include <functional>
#include "Parametric.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"

namespace frame
{

/**
 * @brief Curve frame generator base class.
 *
 * A class to compute the various frames based on curves.
 */
template <int N, class Real = float>
class Curve : public std::unary_function<Real, geom::Matrix<N, N, Real> >
{
public:
	/**
	 * @brief Empty constructor.
	 */
	Curve() :
		_curve( 0 )
	{
	}

	/**
	 * @brief Constructor from a curve (which is cloned).
	 * @param curve The axial curve.
	 */
	template <class CurveType>
	Curve( const CurveType & curve ) :
		_curve( new CurveType( curve ) )
	{
	}

	/**
	 * @brief Deletes the cloned curve.
	 */
	~Curve()
	{
		if ( _curve ) delete _curve;
	}

	/**
	 * @brief Returns a pointer to the cloned curve.
	 */
	inline const curve::Parametric<N, Real> * getCurve() const { return _curve; }

	/**
	 * @brief Removes the old curve and clones a new curve.
	 */
	template <class CurveType>
	void setCurve( const CurveType & curve )
	{
		if ( _curve ) delete _curve;
		_curve = new CurveType( curve );
	}

	/**
	 * @brief Computes the frame.
	 * @param t The parameter t along the curve.
	 * @return The computed frame (Matrix NxN).
	 */
	virtual geom::Matrix<N, N, Real> operator() ( const Real& t ) const = 0;

private:
	const curve::Parametric<N, Real> * _curve;
};

} // namespace

#endif
