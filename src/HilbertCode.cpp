#include <math.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/hilbert_sort.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <boost/iterator/counting_iterator.hpp>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous_d.h>
#include <iostream>
#include <vector>
#include "HilbertCode.h"

#include <stdio.h>
#include <stdlib.h>  //malloc

typedef CGAL::Cartesian_d<double> Kernel;
typedef Kernel::Point_d                                Point_d;

typedef 
  CGAL::Spatial_sort_traits_adapter_d<Kernel,Point_d*>  Search_traits_d;


void Hilbert_Sort_CGAL(double *x, int dx, int N, double *J)
{
	int i,k;
	double *work = new double[dx];

	std::vector<Point_d> points; 
	for(i=0;i<N;i++)
	{
		for(k=0;k<dx;k++)
		{
			work[k]=x[dx*i+k];
		}
		
		points.push_back(Point_d(dx,work+0, work+dx));
	}

	std::vector<std::ptrdiff_t> indices;
  	indices.reserve(points.size());
  
  	std::copy(boost::counting_iterator<std::ptrdiff_t>(0),
            boost::counting_iterator<std::ptrdiff_t>(points.size()),
            std::back_inserter(indices));
  
  	CGAL::hilbert_sort( indices.begin(),indices.end(),Search_traits_d(&(points[0])) );
  
  	for (i=0;i<N;i++)
  	{
		 J[i]=indices[i];
  	}

	delete [] work;
	work=NULL;
}



















