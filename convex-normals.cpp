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
#include "MetricCluster.h"
#include "FuzzyPartition.h"


static const bool ColorAccordingToNbMP = false;
static const bool ColorAccordingToNbClusters = false;
static const bool ColorAccordingToRandom = true;
static const bool FuzzyPartitioning = true;
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
                          unsigned int nbMax,
                          double ratio, double simthreshold )
{
  typedef TangentialCover<DigitalSurface,VertexEmbedder, DGtal::int64_t>
    MyTangentialCover;
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename KSpace::Space::RealVector RealVector;
  typedef typename DigitalSurface::Vertex Vertex;
  typedef typename MyTangentialCover::Index Index;

  static const typename MyTangentialCover::AveragingMode averaging = 
    MyTangentialCover::RadiusAndDistanceAveraging; //SimpleAveraging; // DistanceAveraging;
  // DistanceAveraging: seems best
  // RadiusAndDistanceAveraging: seems also best
  // SimpleAveraging: fair
  // InOutAveraging: same as SimpleAveraging
  // MaxProjectedPlane: not good.
  const KSpace & ks = digSurf.container().space();
  MyTangentialCover tgtCover;
  tgtCover.init( digSurf, embedder, nup, nuq, 400 );
  tgtCover.computeOnce( nbMax );
  trace.info() << std::endl;
  trace.beginBlock( "Sorting planes" );
  tgtCover.sortPlanes();
  tgtCover.purgePlanes( ratio );
  trace.endBlock();

  // Random colors
  std::map< Vertex, Color > mapVertexColor;
  for ( typename DigitalSurface::ConstIterator it = digSurf.begin(), 
          itE = digSurf.end(); it != itE; ++it )
    {
      Color x( random() % 256, random() % 256, random() % 256 );
      mapVertexColor[ *it ] = x;
    }

  if ( FuzzyPartitioning )
    {
      // typedef typename MyTangentialCover::NormalSimilarity Similarity;
      typedef typename MyTangentialCover::NormalAreaSimilarity Similarity;
      typedef FuzzyPartition< DigitalSurface, 
                              Similarity > Partition;
      trace.beginBlock( "Partition surface" );
      // Similarity similarity( tgtCover, averaging );
      Similarity similarity( tgtCover, 0.05, averaging );
      Partition partition( digSurf, similarity );
      unsigned int nbC = partition.partition( simthreshold );
      trace.info() << "- Partition has " << nbC << " components." << std::endl;
      trace.endBlock();
      for ( typename DigitalSurface::ConstIterator it = digSurf.begin(), 
              itE = digSurf.end(); it != itE; ++it )
        {
          Vertex rep1 = partition.representative( *it );
          if ( rep1 != *it )
            mapVertexColor[ *it ] = mapVertexColor[ rep1 ];
        }
    }

  HueShadeColorMap<Index,1> hueShade( 0, nbMax );
  HueShadeColorMap<unsigned int,1> hueShade2( 0, 10 );
  Color color( 150, 150, 180 );
  Color c;
  RealVector normal;
  for ( typename DigitalSurface::ConstIterator 
          it = digSurf.begin(), itE = digSurf.end();
        it != itE; ++it )
    {
      typedef typename MyTangentialCover::MaximalPlaneSummaryIndicesConstIterator MPSIConstIterator;
      Vertex p = *it;
      tgtCover.getEstimatedNormal( normal, p, averaging );
      unsigned int nbMP = 0;
      c = color;
      if ( ColorAccordingToNbMP )
        {
          for ( MPSIConstIterator itMPSI = tgtCover.begin( p ), itMPSIEnd = tgtCover.end( p );
                itMPSI != itMPSIEnd; ++itMPSI )
            ++nbMP;
          c = hueShade( nbMP );
        }
      if ( ColorAccordingToNbClusters )
        {
          typedef SquaredEuclideanDistance<RealVector> SqDistance;
          typedef MetricCluster< RealVector, SqDistance > Cluster;
          typedef std::vector< RealVector > Data;
          Data normals;
          SqDistance sqDistance;
          for ( MPSIConstIterator itMPSI = tgtCover.begin( p ), itMPSIEnd = tgtCover.end( p );
                itMPSI != itMPSIEnd; ++itMPSI )
            normals.push_back( tgtCover.maximalPlaneSummary( *itMPSI ).normal );
          // normals.resize( normals.size() / 2 + 1 );
          Cluster cluster;
          cluster.init( sqDistance, normals.begin(), normals.end() );
          double indexI_1 = cluster.indexSigma();
          trace.info() << "#MP=" << normals.size() << std::endl;
          unsigned int best = 1;
          for ( unsigned int k = 2; k < (normals.size()-1) && ( k < 10 ); ++k )
            {
              cluster.randomClusters( k );
              cluster.Lloyd( 100 );
              double indexI_2 = cluster.indexSigma(); // cluster.indexI( 2 );
              if ( indexI_1 > indexI_2 ) break;
              indexI_1 = indexI_2;
              best = k;
            }
          c = (best == 1) ? Color::Red : hueShade2( best );
        }
      if ( ColorAccordingToRandom )
        {
          c = mapVertexColor[ p ];
        }
      outputCellInColorWithNormal( out, 
                                   ks, p, c, normal );
    }
}

void usage( int, char** argv )
{
  std::cerr << "Usage: " << argv[ 0 ] << " <fileName.vol> <minT> <maxT> <p> <q> <n> <ratio> <simthreshold>" << std::endl;
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
  double ratio = argc >= 8 ? atof( argv[ 7 ] ) : 0.5;
  double simthreshold = argc >= 9 ? atof( argv[ 8 ] ) : 0.8;
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
  outputNuConvexSetNormals( outFile, digSurf, embedder, p, q, nbMax, 
                            ratio, simthreshold );
  outFile.close();
  return true ? 0 : 1;
}
