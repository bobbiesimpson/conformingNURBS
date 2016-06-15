#include <memory>
#include <iostream>

#include "NURBSSurface.h"
#include "NURBSCommon.h"
#include "IElem.h"
#include "Point.h"
#include "IElemIntegrate.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkQuad.h"
#include "vtkPointData.h"

using namespace nurbs::nurbshelper;
using namespace nurbs::elem;

namespace nurbs
{
	
	Point3D NURBSSurface::getCoord( const double s, const double t ) const
	{
		Point4D temp;
		std::vector<UIntVec> indices = getIndices(s,t);
		UIntVec span = getSpan(s,t);
		DoubleVecVec basis = BsplineBasisImpl(s,t,span);
		for(std::size_t i = 0; i < getOrder(S) + 1; ++i)
			for(std::size_t j = 0; j < getOrder(T) + 1; ++j)
				temp += mCPts[ indices[T][j] ][ indices[S][i] ]
					* basis[S][i] * basis[T][j];
		return temp.asCartesian();
	}

	double NURBSSurface::jacDet(const double s, const double t) const
	{
		return cross(tangent(s,t,S), tangent(s,t,T)).length();
	}

	Point3D NURBSSurface::normal(const double s, const double t) const
	{
		Point3D n = cross(tangent(s,t,S), tangent(s,t,T));
		return n / n.length();
	}
	
	void NURBSSurface::printData( std::ostream& ost ) const
	{
		ost << "--------- NURBS SURFACE DATA -------------\n\n";
		ost << "order ( s direction ) = " << getOrder( S ) << "\n";
		ost << "order ( t direction ) = " << getOrder( T ) << "\n";	
		ost << "knot vector ( s direction ) = {";
		std::copy( mKnotVec[ S ].begin(), mKnotVec[ S ].end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mKnotVec[ S ].back() << "}\n";
		ost << "knot vector ( t direction ) = {";
		std::copy( mKnotVec[ T ].begin(), mKnotVec[ T ].end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mKnotVec[ T ].back() << "}\n";
		ost << "unique knot vector ( s direction ) = {";
		std::copy( mUniqueKnotVec[ S ].begin(), mUniqueKnotVec[ S ].end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mUniqueKnotVec[ S ].back() << "}\n";
		ost << "unique knot vector ( t direction ) = {";
		std::copy( mUniqueKnotVec[ T ].begin(), mUniqueKnotVec[ T ].end() - 1,
				   std::ostream_iterator< double >( ost, "," ) );
		ost << mUniqueKnotVec[ T ].back() << "}\n";
		ost << "Control points:\n";
		size_t count = 1;
		for( auto row : mCPts )
		{
			for( auto p : row )
				ost << "(" << count++ << ") " << p << "\n";
		}
		ost << "\n--------------------------------------\n";		
	}

	

	void NURBSSurface::outputVTKFile( const std::string& filename, const uint ngridpts ) const
	{
		const uint seg_n = ngridpts - 1;   // number of cells in each parametric direction
		vtkSmartPointer< vtkUnstructuredGrid > grid = vtkUnstructuredGrid::New();
		vtkSmartPointer< vtkPoints > points = vtkPoints::New();

		vtkSmartPointer<vtkDoubleArray> pw_data = vtkDoubleArray::New();
		pw_data->SetNumberOfComponents(2);
		pw_data->SetName("plane wave data");
		pw_data->SetComponentName(0, "real");
		pw_data->SetComponentName(1, "imag");
		
		uint sample_offset = 0;
		for( IElem ielem( *this ); !ielem.isDone(); ++ielem )
		{
			const Element element = ielem.getCurrentEl();
			//std::cout << element << "\n";
			uint count = 0;
			for( ISamplePt isamplept( element, ngridpts ); !isamplept.isDone(); ++isamplept )
			{
				const ParamPt samplept = isamplept.getCurrentPt();
				const Point3D phys_coord = getCoord( samplept.s, samplept.t );
				points->InsertPoint( sample_offset + count, phys_coord.data() );
				const double k = 50.0;
				std::complex<double> pw = std::exp(std::complex<double>(0.0, k * samplept.s));
				pw_data->InsertComponent(sample_offset + count, 0, std::real(pw));
				pw_data->InsertComponent(sample_offset + count, 1, std::imag(pw));
				++count;
				//std::cout << samplept << "\n";
				//std::cout << getCoord( samplept.s, samplept.t ) << "\n";
									 
			}
			for( uint t = 0; t < seg_n; ++t )
			{
				for( uint s = 0; s < seg_n; ++s )
				{
					vtkSmartPointer< vtkCell > cell = vtkQuad::New();
					cell->GetPointIds()->SetId( 0, sample_offset + t * ngridpts + s );
					cell->GetPointIds()->SetId( 1, sample_offset + t * ngridpts + s + 1 );
					cell->GetPointIds()->SetId( 2, sample_offset + ( t + 1 ) * ngridpts + s + 1 );
					cell->GetPointIds()->SetId( 3, sample_offset + ( t + 1 ) * ngridpts + s );
					grid->InsertNextCell( cell->GetCellType(), cell->GetPointIds() );
				}
			}
			sample_offset += ngridpts * ngridpts;
        }
		grid->SetPoints( points );
		grid->GetPointData()->AddArray(pw_data);
		vtkSmartPointer< vtkXMLUnstructuredGridWriter > writer = vtkXMLUnstructuredGridWriter::New();
		writer->SetFileName( filename.c_str() );
		writer->SetInputData( grid );
		if( !writer->Write() )
			error( "Cannot write vtk file" );
	}
	
	Point3D NURBSSurface::tangent(const double s, const double t, const ParamDir dir) const
	{
		const ParamDir dir2 = (dir == S) ? T : S;
		Point4D temp;
		std::vector<UIntVec> indices = getIndices(s,t);
		UIntVec span = getSpan(s,t);
		DoubleVecVec basis = BsplineBasisImpl(s,t,span);
		DoubleVecVec ders = BsplineBasisDersImpl(s,t,span);
		Point3D a, ader;
		double w = 0.0, wder = 0.0;
		for(std::size_t i = 0; i < getOrder(S) + 1; ++i ) {
			for(std::size_t j = 0; j < getOrder(T) + 1; ++j ) {
				const std::size_t der_i = (dir == S) ? i : j;
				const std::size_t nder_i = (dir == S) ? j : i;
				const Point4D& p = mCPts[ indices[T][j] ][ indices[S][i] ];
				a += p.asUnweighted() * basis[S][i] * basis[T][j];
				ader += p.asUnweighted() * ders[dir][der_i] * basis[dir2][nder_i];
				w += p.getWeight() * basis[S][i] * basis[T][j];
				wder += p.getWeight() * ders[dir][der_i] * basis[dir2][nder_i];
			}
		}
		return getNonRationalDeriv( { a, ader }, { w, wder } );
	}

	DoubleVecVec NURBSSurface::BsplineBasisImpl(const double s, const double t,
												const UIntVec& span) const
	{
		DoubleVec s_basis = getBsplineBasis( s, span[S], mKnotVec[ 0 ], getOrder( S ) );
		DoubleVec t_basis = getBsplineBasis( t, span[T], mKnotVec[ 1 ], getOrder( T ) );
		return {s_basis, t_basis};
	}

	DoubleVecVec NURBSSurface::BsplineBasisDersImpl(const double s,
													const double t,
													const UIntVec& span,
													const DerivOrder order) const
	{
		DoubleVecVec s_basis = getBsplineBasisDers(s, span[S], mKnotVec[ 0 ], getOrder(S), order);
		DoubleVecVec t_basis = getBsplineBasisDers(t, span[T], mKnotVec[ 1 ], getOrder(T), order);
		return {s_basis.at(order), t_basis.at(order)};
	}	
	
	UIntVec NURBSSurface::getSpan(const double s, const double t) const
	{
		return {getKnotSpan(s, mKnotVec[0], getOrder(S)),
				getKnotSpan(t, mKnotVec[1], getOrder(T)) };
	}

	std::vector<UIntVec> NURBSSurface::getIndices(const double s, const double t) const
	{
		return {getBasisFnIndices(s, mKnotVec[0], getOrder(S)),
				getBasisFnIndices(t, mKnotVec[1], getOrder(T)) };
	}
	
	std::ostream& operator<<( std::ostream& ost, const NURBSSurface& s )
	{
		s.printData( ost );
		return ost;
	}
	
	
}
