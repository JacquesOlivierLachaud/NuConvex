#include <QtGui/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/base/Lambda2To1.h"
#include "DGtal/kernel/SquaredEuclideanDistance.h"
#include "DGtal/kernel/CanonicCellEmbedder.h"
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
#include "DGtal/io/colormaps/ColorBrightnessColorMap.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "BasicHPolytopeND.h"
#include "NuConvexSet.h"
#include "DigitalSurface2InnerPointFunctor.h"
#include "TangentialCover.h"
#include "MetricCluster.h"
#include "FuzzyPartition.h"


static const bool ColorAccordingToNbMP = false;
static const bool ColorAccordingToNbClusters = false;
static const bool ColorAccordingToRandom = false; // for seeing the partition
static const bool MarkPartitionBoundaries = true; // for seeing the partition boundaries.
static const bool FuzzyPartitioning = true;
static const bool OutputOFFSurface = true;
static const Color DefaultColor = Color( 255, 200, 100 );
static const Color EdgeColor = Color( 155, 0, 0 );
static const Color SuspectColor = Color( 255, 0, 0 );

typedef float Vec3f[ 3 ];

template <typename TKSpace>
struct Map2SCellEmbedderAdapter {
public:
  typedef TKSpace KSpace;
  typedef typename KSpace::SCell SCell;
  typedef typename KSpace::Space Space;
  typedef typename Space::RealPoint RealPoint;
  typedef SCell Argument;
  typedef RealPoint Value;

  typedef std::map< SCell, RealPoint > Map;

  inline Map2SCellEmbedderAdapter( const Map & aMap )
    : myMap( &aMap )
  {}
  inline Map2SCellEmbedderAdapter( const Map2SCellEmbedderAdapter & other )
    : myMap( other.myMap )
  {}
  inline
  Map2SCellEmbedderAdapter&
  operator=( const Map2SCellEmbedderAdapter & other )
  {
    myMap = other.myMap;
    return *this;
  }
  
  inline const Value & operator()( const Argument & arg ) const
  {
    ASSERT( myMap != 0 );
    typename Map::const_iterator it = myMap->find( arg );
    if ( it != myMap->end() )
      return it->second;
    else
      return myDefaultValue;
  }

private:
  const Map* myMap;
  Value myDefaultValue;
};

template <typename TKSpace>
struct Map2CellEmbedderAdapter {
public:
  typedef TKSpace KSpace;
  typedef typename KSpace::Cell Cell;
  typedef typename KSpace::Space Space;
  typedef typename Space::RealPoint RealPoint;
  typedef Cell Argument;
  typedef RealPoint Value;

  typedef std::map< Cell, RealPoint > Map;

  inline Map2CellEmbedderAdapter( const Map & aMap )
    : myMap( &aMap )
  {}
  inline Map2CellEmbedderAdapter( const Map2CellEmbedderAdapter & other )
    : myMap( other.myMap )
  {}
  inline
  Map2CellEmbedderAdapter&
  operator=( const Map2CellEmbedderAdapter & other )
  {
    myMap = other.myMap;
    return *this;
  }
  
  inline const Value & operator()( const Argument & arg ) const
  {
    ASSERT( myMap != 0 );
    typename Map::const_iterator it = myMap->find( arg );
    if ( it != myMap->end() )
      return it->second;
    else
      return myDefaultValue;
  }

private:
  const Map* myMap;
  Value myDefaultValue;
};



//-----------------------------------------------------------------------------  
template <typename TDigitalSurface, typename SCellEmbedder>
void
exportDigitalSurfaceAs3DCOFF
( std::ostream & out,
  const TDigitalSurface & digSurf,
  const SCellEmbedder & cellEmbedder,
  const SCellEmbedder & cellColor )
{
  BOOST_CONCEPT_ASSERT(( CSCellEmbedder< SCellEmbedder > ));

  typedef TDigitalSurface DigitalSurface;
  typedef typename DigitalSurface::Vertex Vertex;
  typedef typename DigitalSurface::FaceSet FaceSet;
  typedef typename DigitalSurface::VertexRange VertexRange;
  typedef typename DigitalSurface::ConstIterator ConstIterator;
  typedef typename SCellEmbedder::SCell MySCell;
  typedef typename SCellEmbedder::RealPoint RealPoint;
  typedef DGtal::uint64_t Number;

  // Numbers all vertices.
  std::map<Vertex, Number> index;
  Number nbv = 0;
  for ( ConstIterator it = digSurf.begin(), it_end = digSurf.end();
        it != it_end; ++it )
    index[ *it ] = nbv++;
  // Get faces
  // std::cerr << "- " << nbv << " vertices." << std::endl;
  FaceSet faces = digSurf.allClosedFaces();
  // Compute the number of edges and faces.
  Number nbe = 0;
  Number nbf = 0;
  for ( typename FaceSet::const_iterator
          itf = faces.begin(), itf_end = faces.end();
        itf != itf_end; ++itf )
    {
      if ( itf->isClosed() ) 
        { nbe += itf->nbVertices; ++nbf; }
      else
        { nbe += itf->nbVertices - 1; }
    }
  // std::cerr << "- " << nbf << " faces." << std::endl;
  // Outputs OFF header. (strangely COFF does not work with geomview.)
  out << "OFF" << std::endl
      << "# Generated by DGtal::DigitalSurface." << std::endl
      << nbv << " " << nbf << " " << ( nbe / 2 ) << std::endl;
  // Outputs vertex coordinates (the 3 first ones).
  RealPoint p, v;
  for ( ConstIterator it = digSurf.begin(), it_end = digSurf.end();
        it != it_end; ++it )
    {
      MySCell c = *it;
      p = cellEmbedder( c );
      v = cellColor( c );
      //cembedder.embedSCell( *it, p, v );
      // double norm = v.norm();
      // if ( norm != 0.0 ) v /= norm;
      out << p[ 0 ] << " " << p[ 1 ] << " " << p[ 2 ] << " "
          << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << std::endl;
      // double areaD = NumberTraits<Coordinate>::castToDouble(area)*2.0; 
    }
  // Outputs closed faces.
  for ( typename FaceSet::const_iterator
          itf = faces.begin(), itf_end = faces.end();
        itf != itf_end; ++itf )
    {
      if ( itf->isClosed() ) 
        {
          out << itf->nbVertices;
          VertexRange vtcs = digSurf.verticesAroundFace( *itf );
          for ( typename VertexRange::const_iterator
                  itv = vtcs.begin(), itv_end = vtcs.end();
                itv != itv_end; ++itv )
            out << " " << index[ *itv ];
          out << std::endl;
        }
    }
}



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
  typedef typename KSpace::Point Point;
  typedef typename KSpace::Space::RealVector RealVector;
  typedef typename DigitalSurface::Vertex Vertex;
  typedef typename DigitalSurface::Arc Arc;
  typedef typename DigitalSurface::ArcRange ArcRange;
  typedef typename DigitalSurface::SCell SCell;
  typedef typename MyTangentialCover::Index Index;

  static const typename MyTangentialCover::AveragingMode averaging = 
    MyTangentialCover::DistanceAveraging;
  // DistanceAveraging: seems best
  // RadiusAndDistanceAveraging: seems also best (only a problem when some radius is 0).
  // SimpleAveraging: fair
  // InOutAveraging: same as SimpleAveraging
  // MaxProjectedPlane: not good.
  const KSpace & ks = digSurf.container().space();
  MyTangentialCover tgtCover;
  MyTangentialCover::current = &tgtCover;
  tgtCover.init( digSurf, embedder, nup, nuq, 400 );
  tgtCover.computeOnceMemoryLess( nbMax, false, 1000000000 );
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
      if ( ColorAccordingToRandom )
        {
          Color x( random() % 256, random() % 256, random() % 256 );
          mapVertexColor[ *it ] = x;
        }
      else
        {
          mapVertexColor[ *it ] = DefaultColor;
        }
    }

  if ( FuzzyPartitioning )
    {
      // Shade colors for dissimilar edges.
      ColorBrightnessColorMap<double> whiteShade( 0.0, simthreshold, Color::White );

      typedef typename MyTangentialCover::NormalSimilarity Similarity;
      //typedef typename MyTangentialCover::NormalAreaSimilarity Similarity;
      typedef FuzzyPartition< DigitalSurface, 
                              Similarity > Partition;
      Similarity similarity( tgtCover, averaging );
      //Similarity similarity( tgtCover, 0.05, averaging );
      
      // Iterating surface partitionning to select maximal planes of a
      // vertex within the ones having center in the same region.
      std::vector< bool > isVtxSuspect;
      unsigned int nbChanged;
      do {
        trace.beginBlock( "Partition surface" );
        Partition partition( digSurf, similarity );
        unsigned int nbC = partition.partition( simthreshold );
        trace.info() << "- Partition has " << nbC << " components." << std::endl;
        trace.endBlock();
        // Recompute normal.
        trace.beginBlock( "Recompute maximal planes." );
        typedef typename MyTangentialCover::MaximalPlaneSummaryIndices MPSIndices;
        nbChanged = 0;
        long nbMPChanged = 0;
        long nbSuspect = 0;
        isVtxSuspect.clear();
        for ( typename DigitalSurface::ConstIterator it = digSurf.begin(), 
                itE = digSurf.end(); it != itE; ++it )
          {
            Vertex vtx = *it;
            Vertex rep = partition.representative( vtx );
            MPSIndices & mpsi = tgtCover.planeIndices( vtx );
            MPSIndices mpsiTmp;
            for ( typename MPSIndices::const_iterator itMPSI = mpsi.begin(),
                    itMPSIEnd = mpsi.end(); itMPSI != itMPSIEnd; ++itMPSI )
              {
                const Vertex & vtxMP = tgtCover.center( *itMPSI );
                if ( partition.representative( vtxMP ) == rep )
                  mpsiTmp.push_back( *itMPSI );
              }
            isVtxSuspect.push_back( mpsiTmp.size() == 0 );
            if ( ( mpsiTmp.size() != mpsi.size() ) && ( mpsiTmp.size() != 0 ) )
              {
                nbChanged += 1;
                nbMPChanged += mpsi.size() - mpsiTmp.size();
                mpsi.swap( mpsiTmp );
              }
            else if ( mpsiTmp.size() == 0 )
              nbSuspect += 1;
          }
        trace.info() << "- nb vertex changed = " << nbChanged 
                     << ", nb plane changed = " << nbMPChanged 
                     << ", nb suspect = " << nbSuspect << std::endl;
        trace.endBlock();

        // Check if last iteration of loop.
        if ( nbChanged == 0 )
          { // Update colors
            unsigned int i = 0;
            for ( typename DigitalSurface::ConstIterator it = digSurf.begin(), 
                    itE = digSurf.end(); it != itE; ++it, ++i )
              {
                // Take care of color.
                Vertex v1 = *it;
                Vertex rep1 = partition.representative( v1 );
                if ( v1 != rep1 )
                  mapVertexColor[ v1 ] = mapVertexColor[ rep1 ];
                if ( isVtxSuspect[ i ] )
                  mapVertexColor[ v1 ] = SuspectColor;
                // Display separating edges.
                ArcRange outArcs = digSurf.outArcs( v1 );
                for ( typename ArcRange::const_iterator itArcs = outArcs.begin(),
                        itArcsEnd = outArcs.end(); itArcs != itArcsEnd; ++itArcs )
                  {
                    Vertex v2 = digSurf.head( *itArcs );
                    Vertex rep2 = partition.representative( v2 );
                    if ( v1 < v2 )
                      {
                        double s = similarity( v1, v2 );
                        if ( s < simthreshold )
                          {
                            outputCellInColor( out, ks, digSurf.separator( *itArcs ),
                                               whiteShade( s ) );
                            if ( MarkPartitionBoundaries 
                                 && ( rep1 != rep2 )
                                 && ! isVtxSuspect[ i ] )
                              {
                                if ( v1 != rep1 ) mapVertexColor[ v1 ] = EdgeColor;
                                if ( v2 != rep2 ) mapVertexColor[ v2 ] = EdgeColor;
                                // mapVertexColor[ v2 ] = Color::Green;
                                //trace.info() << "- " << rep1 << " | " << partition.representative( v2 ) << std::endl;
                              }
                          }
                      }
                  }
              }
            }
      } while ( nbChanged != 0 ); // end of iteration

    }

  // Exporting as OFF
  if ( OutputOFFSurface )
    {
      trace.beginBlock( "Output OFF surface <reconstruction.off>, <partition.off> and <origin.off>" );
      typedef Map2CellEmbedderAdapter< KSpace > CellEmbedder;
      typedef Map2SCellEmbedderAdapter< KSpace > SCellEmbedder;
      typedef typename CellEmbedder::Map CellMap;
      typedef typename SCellEmbedder::Map SCellMap;
      CanonicCellEmbedder<KSpace> canonicCellEmbedder( ks );
      CellMap mapUVtx2RealPoint;
      SCellMap mapVtx2RealPoint;
      //      double maxDist = (double)nup / (double)nuq;
      for ( typename DigitalSurface::ConstIterator it = digSurf.begin(), 
              itE = digSurf.end(); it != itE; ++it )
        {
          Vertex vtx = *it;
          RealVector n1;
          Point p = embedder( vtx );
          RealVector x( (double) p[ 0 ], (double) p[ 1 ], (double) p[ 2 ] );
          // better at 0.5 on cms
          tgtCover.getEstimatedProjection( n1, vtx, x, 0.5, averaging );
          // Limiting distance does not seem an improvment.
          // double d1 = (n1-x).norm();
          // if ( d1 > maxDist ) n1 = x + (n1-x)*maxDist/d1;
          // x = canonicCellEmbedder( ks.unsigns( vtx ) );
          mapUVtx2RealPoint[ ks.unsigns( vtx ) ] = n1;
          mapVtx2RealPoint[ vtx ] = n1;
        }
      CellEmbedder cellEmbedder( mapUVtx2RealPoint );
      SCellEmbedder scellEmbedder( mapVtx2RealPoint );
      ofstream outFile( "reconstruction.off" );
      digSurf.exportEmbeddedSurfaceAs3DOFF( outFile, cellEmbedder );
      outFile.close();
      ofstream outFile2( "origin.off" );
      digSurf.exportEmbeddedSurfaceAs3DOFF( outFile2, canonicCellEmbedder );
      outFile2.close();
      // Exporting as colored OFF.
      SCellMap mapVtx2Color;
      for ( typename std::map< Vertex, Color >::const_iterator 
              it = mapVertexColor.begin(), itE = mapVertexColor.end();
            it != itE; ++it )
        {
          const Color & color = it->second;
          RealVector vcolor( ((double) color.red()) / 255.0,
                            ((double) color.green()) / 255.0,
                            ((double) color.blue()) / 255.0 );
          mapVtx2Color[ it->first ] = vcolor;
        }
      SCellEmbedder scellColor( mapVtx2Color );
      ofstream outFile3( "partition.off" );
      exportDigitalSurfaceAs3DCOFF( outFile3, digSurf, scellEmbedder, scellColor );
      outFile3.close();
      trace.endBlock();
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
      c = DefaultColor;
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

using namespace Z3i;
typedef ImageSelector < Domain, int>::Type Image;
typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
typedef DigitalSetBoundary<KSpace, DigitalSet > MyDigitalSurfaceContainer;
typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
typedef DigitalSurface2InnerPointFunctor<MyDigitalSurface> VertexEmbedder;

namespace DGtal {
  template<>
  TangentialCover< MyDigitalSurface,VertexEmbedder, DGtal::int64_t>*
  TangentialCover< MyDigitalSurface,VertexEmbedder, DGtal::int64_t>::current = 0;
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
  MySurfelAdjacency surfAdj( false ); // interior in all directions.
  //! [convex-normals-SurfelAdjacency]

  //! [convex-normals-SetUpDigitalSurface]
  trace.beginBlock( "Set up digital surface." );
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

  VertexEmbedder embedder( digSurf );
  ofstream outFile( "titi.txt" );
  outputNuConvexSetNormals( outFile, digSurf, embedder, p, q, nbMax, 
                            ratio, simthreshold );
  outFile.close();
  return true ? 0 : 1;
}
