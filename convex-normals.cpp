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
#include "DGtal/topology/DigitalSetBoundary.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "BasicHPolytopeND.h"
#include "NuConvexSet.h"
#include "DigitalSurface2InnerPointFunctor.h"
#include "TangentialCover.h"


static const bool ColorAccordingToNbMP = false;

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


// Good version
template < typename DigitalSurface, 
           typename VertexEmbedder >
void 
outputNuConvexSetNormals( ostream & out,
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
  typedef typename MyTangentialCover::Index Index;
  const KSpace & ks = digSurf.container().space();
  MyTangentialCover tgtCover;
  tgtCover.init( digSurf, embedder, nup, nuq, 400 );
  tgtCover.computeOnce( nbMax );
  trace.info() << std::endl;
  HueShadeColorMap<Index,1> hueShade( 0, nbMax );
  Color color( 150, 150, 180 );
  RealVector normal;
  for ( typename DigitalSurface::ConstIterator 
          it = digSurf.begin(), itE = digSurf.end();
        it != itE; ++it )
    {
      typedef typename MyTangentialCover::MaximalPlaneIndicesConstIterator MPIConstIterator;
      Vertex p = *it;
      tgtCover.getEstimatedNormal( normal, p, MyTangentialCover::SimpleAveraging );
      unsigned int nbMP = 0;
      for ( MPIConstIterator itMPI = tgtCover.begin( p ), itMPIEnd = tgtCover.end( p );
            itMPI != itMPIEnd; ++itMPI )
        ++nbMP;
      Color c = ColorAccordingToNbMP ? hueShade( nbMP ) : color;
      outputCellInColorWithNormal( out, 
                                   ks, p, c, normal );
    }
}

void usage( int, char** argv )
{
  std::cerr << "Usage: " << argv[ 0 ] << " <fileName.vol> <minT> <maxT> <p> <q> <n>" << std::endl;
  std::cerr << "Computes maximal planes and displays normals of the shape stored in vol file <fileName.vol>." << std::endl;
  std::cerr << "\t - voxel v belongs to the shape iff its value I(v) follows minT < I(v) <= maxT." << std::endl;
  std::cerr << "\t - the rational number p/q is the width of maximal planes." << std::endl;
  std::cerr << "\t - the integer number n is the maximal number of maximal planes kept for each surfel of the shape." 
            << " Choosing 1 is a segmentation of the shape boundary into planes." << std::endl;
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
  //! [convex-normals-readVol]
  trace.beginBlock( "Reading vol file into an image." );
  using namespace Z3i;
  typedef ImageSelector < Domain, int>::Type Image;
  Image image = VolReader<Image>::importVol(inputFilename);
  DigitalSet set3d (image.domain());
  SetPredicate<DigitalSet> set3dPredicate( set3d );
  SetFromImage<DigitalSet>::append<Image>(set3d, image, 
                                          minThreshold, maxThreshold);
  trace.endBlock();
  //! [convex-normals-readVol]

  //! [convex-normals-KSpace]
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
  //! [convex-normals-KSpace]

  //! [convex-normals-SurfelAdjacency]
  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
  MySurfelAdjacency surfAdj( false ); // interior in all directions.
  //! [convex-normals-SurfelAdjacency]

  //! [convex-normals-SetUpDigitalSurface]
  trace.beginBlock( "Set up digital surface." );
  typedef DigitalSetBoundary<KSpace, DigitalSet > MyDigitalSurfaceContainer;
  typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
  MyDigitalSurfaceContainer* ptrSurfContainer = 
    new MyDigitalSurfaceContainer( ks, set3d, surfAdj );
  // typedef LightImplicitDigitalSurface<KSpace, SetPredicate<DigitalSet> > 
  //   MyDigitalSurfaceContainer;
  // typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
  // SCell bel = Surfaces<KSpace>::findABel( ks, set3dPredicate, 100000 );
  // MyDigitalSurfaceContainer* ptrSurfContainer = 
  //   new MyDigitalSurfaceContainer( ks, set3dPredicate, surfAdj, bel );
  MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
  trace.endBlock();
  //! [convex-normals-SetUpDigitalSurface]

  typedef DigitalSurface2InnerPointFunctor<MyDigitalSurface> VertexEmbedder;
  VertexEmbedder embedder( digSurf );
  ofstream outFile( "titi.txt" );
  outputNuConvexSetNormals( outFile, digSurf, embedder, p, q, nbMax );
  outFile.close();
  return true ? 0 : 1;
}
