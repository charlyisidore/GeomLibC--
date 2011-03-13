#include "NURBS.hpp"
#include "Tube.hpp"
#include "Frenet.hpp"
#include <iostream>

int main()
{
	curve::NURBS<3> cur;
	cur.pushControlPoint( geom::Vector3f( 0, 0, 0 ) );
	cur.pushControlPoint( geom::Vector3f( 1, 0, 0 ) );
	cur.pushControlPoint( geom::Vector3f( 1, 1, 0 ) );
	cur.pushControlPoint( geom::Vector3f( 0, 1, 0 ) );
	std::cout << cur( 0.5 ) << std::endl;

	surface::Tube<float> surf( cur, frame::Frenet<float>() );
	std::cout << surf( 0.5, 0 ) << std::endl;
	return 0;
}

