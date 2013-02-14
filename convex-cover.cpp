#include <QtGui/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/base/Lambda2To1.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/kernel/CanonicSCellEmbedder.h"
#include "DGtal/graph/BreadthFirstVisitor.h"
#include "DGtal/graph/DistanceVisitor.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "BasicHPolytopeND.h"
#include "NuConvexSet.h"
#include "DigitalSurface2InnerPointFunctor.h"

using namespace DGtal;

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

template <typename DigitalSurface>
struct SurfelAreaEstimator
{
  typedef typename DigitalSurface::Vertex Vertex;
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename KSpace::Space Space;
  typedef typename Space::RealPoint RealPoint;
  typedef typename Space::RealVector RealVector;
  typedef typename RealVector::Component Scalar;
  SurfelAreaEstimator( const DigitalSurface & digSurf )
    : myDigSurf( digSurf )
  {}
  inline
  Scalar operator()( const Vertex & v, const RealVector & n ) const
  {
    const KSpace & ks = myDigSurf.container().space();
    Dimension k = ks.sOrthDir( v );
    bool direct = ks.sDirect( v, k );
    // normal toward outside. 
    return direct ? n[ k ] : -n[ k ];
  }
  
  const DigitalSurface & myDigSurf;
};

template <typename DigitalSurface, typename VertexEmbedder >
bool viewNuConvexSet( Viewer3D & viewer,
		      const DigitalSurface & digSurf, 
		      const VertexEmbedder & embedder,
                      unsigned int nup,
                      unsigned int nuq,
                      typename DigitalSurface::Vertex p )
{
  typedef DigitalSurface Graph;
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename DigitalSurface::Vertex Vertex;
  typedef typename KSpace::Space Space;
  typedef typename Space::RealPoint RealPoint;
  typedef typename RealPoint::Coordinate Scalar;
  typedef ExactPredicateLpSeparableMetric<Space,2> Distance;
  typedef std::binder1st< Distance > DistanceToPoint; 
  typedef Composer<VertexEmbedder, DistanceToPoint, Scalar> VertexFunctor;
  typedef DistanceVisitor< Graph, VertexFunctor > Visitor;

  typedef NuConvexSet< Space, Visitor, VertexEmbedder, DGtal::int64_t > MyNuConvexSet;
  
  Distance distance;
  DistanceToPoint distanceToPoint = std::bind1st( distance, embedder( p ) );
  VertexFunctor vfunctor( embedder, distanceToPoint );
  // DistanceVisitor< Graph, VertexFunctor > visitor( g, p, vfunctor );
  const KSpace & ks = digSurf.container().space();
  MyNuConvexSet nuConvex( Visitor( digSurf, vfunctor, p ), embedder );
  nuConvex.setExtensionMode( false );
  nuConvex.init( nup, nuq, 400 );
  nuConvex.compute( -1.0 );
  viewer << CustomColors3D( Color::White, Color::White ) 
         << ks.unsigns( p );
  unsigned int nb = 0;
  Color c( random() % 256, random() % 256, random() % 256 );
  for ( typename MyNuConvexSet::ConstIterator it = nuConvex.begin(),
          itE = nuConvex.end(); it != itE; ++it, ++nb )
    viewer << CustomColors3D( Color::Black, c )
           << ks.unsigns( *it );
  std::cerr << "nu-convex has " << nb << " surfels." << std::endl;
  unsigned int nb2 = 0;
  // for ( typename MyNuConvexSet::ConstIterator it = nuConvex.myRejectedVertices.begin(),
  //         itE = nuConvex.myRejectedVertices.end(); it != itE; ++it, ++nb2 )
  //   viewer << CustomColors3D( Color::Black, Color::Blue ) 
  //          << ks.unsigns( *it );
  // std::cerr << "nu-convex has " << nb2 << " rejected surfels." << std::endl;
  SurfelAreaEstimator<DigitalSurface> areaEstimator( digSurf );
  MaximalPlaneSummary<Space> mps;
  nuConvex.summarize( mps, areaEstimator );
  std::cerr << mps << std::endl;
  
  return true;
}

void usage( int, char** argv )
{
  std::cerr << "Usage: " << argv[ 0 ] << " <fileName.vol> <minT> <maxT> <p> <q> <n>" << std::endl;
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
  unsigned int p = argc >= 5 ? atoi( argv[ 4 ] ) : 1;
  unsigned int q = argc >= 6 ? atoi( argv[ 5 ] ) : 1;
  unsigned int n = argc >= 7 ? atoi( argv[ 6 ] ) : 0;
  unsigned int step = argc >= 8 ? atoi( argv[ 7 ] ) : 500;
  //! [volDistanceTraversal-readVol]
  trace.beginBlock( "Reading vol file into an image." );
  using namespace Z3i;
  typedef ImageSelector < Domain, int>::Type Image;
  Image image = VolReader<Image>::importVol(inputFilename);
  DigitalSet set3d (image.domain());
  //  SetPredicate<DigitalSet> set3dPredicate( set3d );
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
  typedef LightImplicitDigitalSurface<KSpace, DigitalSet > 
    MyDigitalSurfaceContainer;
  typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
  SCell bel = Surfaces<KSpace>::findABel( ks, set3d, 100000 );
  MyDigitalSurfaceContainer* ptrSurfContainer = 
    new MyDigitalSurfaceContainer( ks, set3d, surfAdj, bel );
  MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
  trace.endBlock();
  //! [volDistanceTraversal-SetUpDigitalSurface]

  typedef DigitalSurface2InnerPointFunctor<MyDigitalSurface> VertexEmbedder;
  VertexEmbedder embedder( digSurf );
  MyDigitalSurface::ConstIterator it = digSurf.begin();
  MyDigitalSurface::ConstIterator itE = digSurf.end();
  for ( unsigned int i = 0; ( it != itE ) && ( i < n ); ++i )
    ++it;
  while ( it != itE )
    {
      viewNuConvexSet( viewer, digSurf, embedder, p, q, *it );
      for ( unsigned int i = 0; ( it != itE ) && ( i < step ); ++i )
        ++it;
    }
  for ( MyDigitalSurface::ConstIterator it = digSurf.begin(),
          itE = digSurf.end(); it != itE; ++it )
    viewer << CustomColors3D( Color::Black, Color::White ) 
           << ks.unsigns( *it );

  // typedef SpaceND<3>::Vector Vector;
  // bool ok = viewHPolytope<Vector>( viewer );
  viewer<< Viewer3D::updateDisplay;
  application.exec();
  return true ? 0 : 1;
}
