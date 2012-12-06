#include <QtGui/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/base/Lambda2To1.h"
#include "DGtal/kernel/SquaredEuclideanDistance.h"
#include "DGtal/kernel/CanonicSCellEmbedder.h"
#include "DGtal/topology/DistanceVisitor.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "BasicHPolytopeND.h"
#include "NuConvexSet.h"

template <typename TVector>
bool viewHPolytope( Viewer3D & viewer )
{
  typedef TVector Vector;
  typedef BasicHPolytopeND<TVector> HPolytope;
  typedef typename HPolytope::ClosedHalfSpace HalfSpace;
  
  HPolytope P;
  P.add( HalfSpace( Vector( 1, 1, 1 ), 5 ) );
  P.add( HalfSpace( Vector( -1, -1, -1 ), 5 ) );
  P.add( HalfSpace( Vector( -5, -4, 3 ), 12 ) );
  P.add( HalfSpace( Vector( 5, 4, -3 ), 12 ) );
  trace.info() << P << std::endl;
  
  typedef SpaceND<Vector::dimension> Space;
  typedef HyperRectDomain<Space> Domain;
  typedef typename Space::Point Point;
  Domain domain( Point( -6, -6, -6 ), Point( 6, 6, 6 ) );
  for ( typename Domain::ConstIterator it = domain.begin(),
	  itE = domain.end(); it != itE; ++it )
    if ( P( *it ) ) viewer << *it;
  return true;
}

template <typename DigitalSurface, typename VertexEmbedder >
bool viewNuConvexSet( Viewer3D & viewer,
		      const DigitalSurface & digSurf, 
		      const VertexEmbedder & embedder )
{
  typedef DigitalSurface Graph;
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename DigitalSurface::Vertex Vertex;
  typedef typename KSpace::Space Space;
  typedef typename Space::RealPoint RealPoint;
  typedef typename RealPoint::Coordinate Scalar;
  typedef SquaredEuclideanDistance<RealPoint> SqED;
  typedef Lambda2To1<SqED, RealPoint, RealPoint, Scalar> SqEDToPoint;
  typedef Composer<VertexEmbedder, SqEDToPoint, Scalar> VertexFunctor;
  typedef DistanceVisitor< Graph, VertexFunctor > Visitor;

  typedef NuConvexSet< Space, Visitor, VertexEmbedder, DGtal::int64_t > MyNuConvexSet;
  
  Vertex p = *( digSurf.begin() );
  SqED sqed;
  SqEDToPoint distanceToPoint( sqed, embedder( p ) );
  VertexFunctor vfunctor( embedder, distanceToPoint );
  // DistanceVisitor< Graph, VertexFunctor > visitor( g, p, vfunctor );

  MyNuConvexSet nuConvex( Visitor( digSurf, vfunctor, p ), embedder );
  nuConvex.init( 1, 1, 100 );
  nuConvex.compute( -1.0 );
  return true;
}

void usage( int, char** argv )
{
  std::cerr << "Usage: " << argv[ 0 ] << " <fileName.vol> <minT> <maxT>" << std::endl;
  std::cerr << "\t - displays the boundary of the shape stored in vol file <fileName.vol>." << std::endl;
  std::cerr << "\t - voxel v belongs to the shape iff its value I(v) follows minT <= I(v) <= maxT." << std::endl;
}

int main( int argc, char** argv )
{
  if ( argc < 4 )
    {
      usage( argc, argv );
      return 1;
    }

  QApplication application(argc,argv);
  Viewer3D viewer;
  viewer.show();

  std::string inputFilename = argv[ 1 ];
  unsigned int minThreshold = atoi( argv[ 2 ] );
  unsigned int maxThreshold = atoi( argv[ 3 ] );

  //! [volDistanceTraversal-readVol]
  trace.beginBlock( "Reading vol file into an image." );
  using namespace Z3i;
  typedef ImageSelector < Domain, int>::Type Image;
  Image image = VolReader<Image>::importVol(inputFilename);
  DigitalSet set3d (image.domain());
  SetPredicate<DigitalSet> set3dPredicate( set3d );
  SetFromImage<DigitalSet>::append<Image>(set3d, image, 
                                          minThreshold, maxThreshold);
  trace.endBlock();
  //! [volDistanceTraversal-readVol]

  //! [volDistanceTraversal-KSpace]
  trace.beginBlock( "Construct the Khalimsky space from the image domain." );
  KSpace ks;
  bool space_ok = ks.init( image.domain().lowerBound(), 
                           image.domain().upperBound(), true );
  if (!space_ok)
    {
      trace.error() << "Error in the Khamisky space construction."<<std::endl;
      return 2;
    }
  trace.endBlock();
  //! [volDistanceTraversal-KSpace]

  //! [volDistanceTraversal-SurfelAdjacency]
  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
  MySurfelAdjacency surfAdj( true ); // interior in all directions.
  //! [volDistanceTraversal-SurfelAdjacency]

  //! [volDistanceTraversal-SetUpDigitalSurface]
  trace.beginBlock( "Set up digital surface." );
  typedef LightImplicitDigitalSurface<KSpace, SetPredicate<DigitalSet> > 
    MyDigitalSurfaceContainer;
  typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
  SCell bel = Surfaces<KSpace>::findABel( ks, set3dPredicate, 100000 );
  MyDigitalSurfaceContainer* ptrSurfContainer = 
    new MyDigitalSurfaceContainer( ks, set3dPredicate, surfAdj, bel );
  MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
  trace.endBlock();
  //! [volDistanceTraversal-SetUpDigitalSurface]

  typedef CanonicSCellEmbedder<KSpace> VertexEmbedder;
  VertexEmbedder embedder( ks );
  
  viewNuConvexSet( viewer, digSurf, embedder );

  // typedef SpaceND<3>::Vector Vector;
  // bool ok = viewHPolytope<Vector>( viewer );
  // viewer<< Viewer3D::updateDisplay;
  // application.exec();
  return true ? 0 : 1;
}
