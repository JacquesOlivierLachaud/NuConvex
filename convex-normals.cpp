#include <QtGui/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/base/Lambda2To1.h"
#include "DGtal/kernel/SquaredEuclideanDistance.h"
#include "DGtal/kernel/CanonicSCellEmbedder.h"
#include "DGtal/topology/BreadthFirstVisitor.h"
#include "DGtal/topology/DistanceVisitor.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "BasicHPolytopeND.h"
#include "NuConvexSet.h"
#include "DigitalSurface2InnerPointFunctor.h"
#include "TangentialCover.h"


typedef float Vec3f[ 3 ];

template <typename KSpace>
void outputCellInColor( std::ostream & out, 
			const KSpace & ks, 
			typename KSpace::SCell s, 
                        const Color & color )
{
  typename KSpace::Point x = ks.sKCoords( s );
  bool sign = ks.sSign( s );
  out << "Cell " << x[ 0 ] << ' ' << x[ 1 ] << ' ' << x[ 2 ] << ' '
      << sign 
      << ' ' << ((double)color.red())/255.0
      << ' ' << ((double)color.green())/255.0
      << ' ' << ((double)color.blue())/255.0 << std::endl;

}

template <typename KSpace>
void outputCellInColorWithNormal( std::ostream & out, 
				  const KSpace & ks, 
				  typename KSpace::SCell s, 
                                  const Color & color,
				  typename KSpace::Space::RealVector & n )
{
  typename KSpace::Point x = ks.sKCoords( s );
  bool sign = ks.sSign( s );
  out << "CellN " << x[ 0 ] << ' ' << x[ 1 ] << ' ' << x[ 2 ] << ' '
      << sign
      << ' ' << ((double)color.red())/255.0
      << ' ' << ((double)color.green())/255.0
      << ' ' << ((double)color.blue())/255.0
      << ' ' << n[ 0 ] << ' ' << n[ 1 ] << ' ' << n[ 2 ]
      << endl;
}

template <typename KSpace>
void outputLinelsOfSurfelInColor( std::ostream & out, 
				  const KSpace & ks, 
				  typename KSpace::SCell surfel, 
                                  const Color & color )
{
  for ( typename KSpace::DirIterator q = ks.sDirs( surfel );
	q != 0; ++q )
    {
      outputCellInColor( out, ks, ks.sIncident( surfel, *q, true ), color );
      outputCellInColor( out, ks, ks.sIncident( surfel, *q, false ), color );
    }
}

template <typename KSpace, typename Iterator>
void outputSignedSetInColor( std::ostream & out, 
			     const KSpace & ks, 
			     Iterator b, Iterator e, 
                             const Color & color )
{
  for ( ; b != e; ++b )
    outputCellInColor( out, ks, *b, color );
}

void
fillRandomColors( unsigned int nb, Color* colors )
{
  for ( unsigned int i = 0; i < nb; ++i )
    {
      colors[ i ] = Color( random() % 256, random() % 256, random() % 256 );
    }
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
void outputNuConvexSetNormals( ostream & out,
                               const DigitalSurface & digSurf, 
                               const VertexEmbedder & embedder,
                               unsigned int nup,
                               unsigned int nuq )
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
  
  SqED sqed;
  const KSpace & ks = digSurf.container().space();
  SurfelAreaEstimator<DigitalSurface> areaEstimator( digSurf );
  MaximalPlaneSummary<Space> mps;
  Color color( 150, 150, 180 );
  for ( typename DigitalSurface::ConstIterator 
          it = digSurf.begin(), itE = digSurf.end();
        it != itE; ++it )
    {
      Vertex p = *it;
      SqEDToPoint distanceToPoint( sqed, embedder( p ) );
      VertexFunctor vfunctor( embedder, distanceToPoint );
      MyNuConvexSet nuConvex( Visitor( digSurf, vfunctor, p ), embedder );
      nuConvex.setExtensionMode( false );
      nuConvex.init( nup, nuq, 400 );
      nuConvex.compute( -1.0 );
      nuConvex.summarize( mps, areaEstimator );
      std::cerr << mps << std::endl;
      outputCellInColorWithNormal( out, 
                                   ks, p, color, mps.normal );
    }
}

template < typename DigitalSurface, 
           typename VertexEmbedder >
void outputNuConvexSetNormals2( ostream & out,
                                const DigitalSurface & digSurf, 
                                const VertexEmbedder & embedder,
                                unsigned int nup,
                                unsigned int nuq,
                                unsigned int nbMax )
{
  typedef TangentialCover<DigitalSurface,VertexEmbedder, DGtal::int64_t>
    MyTangentialCover;
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename KSpace::Space::RealVector RealVector;
  typedef typename DigitalSurface::Vertex Vertex;
  const KSpace & ks = digSurf.container().space();
  MyTangentialCover tgtCover;
  tgtCover.init( digSurf, embedder, nup, nuq, 400 );
  tgtCover.compute( nbMax );
  Color color( 150, 150, 180 );
  RealVector normal;
  for ( typename DigitalSurface::ConstIterator 
          it = digSurf.begin(), itE = digSurf.end();
        it != itE; ++it )
    {
      Vertex p = *it;
      tgtCover.getEstimatedNormal( normal, p, MyTangentialCover::SimpleAveraging );
      outputCellInColorWithNormal( out, 
                                   ks, p, color, normal );
    }
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

  std::string inputFilename = argv[ 1 ];
  unsigned int minThreshold = atoi( argv[ 2 ] );
  unsigned int maxThreshold = atoi( argv[ 3 ] );
  unsigned int p = argc >= 5 ? atoi( argv[ 4 ] ) : 1;
  unsigned int q = argc >= 6 ? atoi( argv[ 5 ] ) : 1;
  unsigned int nbMax = argc >= 7 ? atoi( argv[ 6 ] ) : 1;
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

  typedef DigitalSurface2InnerPointFunctor<MyDigitalSurface> VertexEmbedder;
  VertexEmbedder embedder( digSurf );
  ofstream outFile( "titi.txt" );
  outputNuConvexSetNormals2( outFile, digSurf, embedder, p, q, nbMax );
  outFile.close();
  return true ? 0 : 1;
}
