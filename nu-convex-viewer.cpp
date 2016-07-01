//#include <QtGui/qapplication.h>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/topology/CanonicCellEmbedder.h"
#include "DGtal/topology/CanonicSCellEmbedder.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/DigitalSetBoundary.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/math/MPolynomial.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/implicit/ImplicitPolynomial3Shape.h"
#include "DGtal/shapes/implicit/ImplicitFunctionDiff1LinearCellEmbedder.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/Color.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/colormaps/ColorBrightnessColorMap.h"
#include "DGtal/io/readers/MPolynomialReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "BasicHPolytopeND.h"
#include "NuConvexSet.h"
#include "DigitalSurface2InnerPointFunctor.h"
#include "TangentialCover.h"

using namespace DGtal;

template < typename Viewer,
           typename DigitalSurface,
           typename VertexEmbedder >
void
viewNuConvex( Viewer & viewer,
              const DigitalSurface & digSurf,
              const VertexEmbedder & embedder,
              unsigned int idx,
              unsigned int nup,
              unsigned int nuq )
{
  typedef typename DigitalSurface::KSpace KSpace;
  typedef typename KSpace::Space Space;
  typedef typename KSpace::Point Point;
  typedef typename KSpace::Space::RealVector RealVector;
  typedef typename KSpace::Space::RealPoint RealPoint;
  typedef typename DigitalSurface::Vertex Vertex;
  typedef typename DigitalSurface::Arc Arc;
  typedef typename DigitalSurface::ArcRange ArcRange;
  typedef typename DigitalSurface::SCell SCell;
  typedef SquaredEuclideanDistance<RealPoint> Distance;
  typedef typename Distance::Value Scalar;
  typedef std::binder1st< Distance > DistanceToPoint; 
  typedef functors::Composer<VertexEmbedder, DistanceToPoint, Scalar> VertexFunctor;
  typedef DistanceBreadthFirstVisitor< DigitalSurface, VertexFunctor > Visitor;
  typedef typename Visitor::Node MyNode;
  typedef NuConvexSet< Space, Visitor, 
                       VertexEmbedder, DGtal::int64_t > MyNuConvexSet;

  idx = idx % digSurf.size();
  typename DigitalSurface::ConstIterator it = digSurf.begin();
  for ( ; idx != 0; --idx ) ++it;
  Vertex bel = *it;

  Distance distance;
  DistanceToPoint distanceToPoint = std::bind1st( distance, embedder( bel ) );
  VertexFunctor vfunctor( embedder, distanceToPoint );
  Visitor visitor( digSurf, vfunctor, bel );
  const KSpace & ks = digSurf.container().space();
  MyNuConvexSet nuconvex( visitor, embedder );
  nuconvex.init( nup, nuq, ks.size( 0 ) );
  nuconvex.setExtensionMode( false );
  // nuconvex.setExtensionMode( true );
  nuconvex.compute( -1.0 );

  viewer << CustomColors3D( Color::Black, Color::White )
         << ks.unsigns( bel );
  viewer << CustomColors3D( Color( 128, 128, 128 ), Color( 128, 128, 128 ) );
  for ( typename MyNuConvexSet::ConstIterator it = nuconvex.begin(),
          itE = nuconvex.end(); it != itE; ++it )
    {
      viewer << ks.unsigns( *it );
    }

  MyNode node;
  Visitor visitor2( visitor );
  for ( ; ! visitor.finished(); visitor.expand() )
    node = visitor.current();
  Scalar maxDist = node.second;

  HueShadeColorMap<Scalar,1> hueShade( 0, maxDist );
  visitor2.expand();
  while ( ! visitor2.finished() )
    {
      MyNode n = visitor2.current(); 
      Color c = hueShade( n.second );
      viewer << CustomColors3D( Color::Red, c )
             << ks.unsigns( n.first );
      visitor2.expand();
    }
  viewer << Viewer3D<Space,KSpace>::updateDisplay;
}

/** 
 * Missing parameter error message.
 * 
 * @param param 
 */
void missingParam(std::string param)
{
  trace.error() << " Parameter: "<<param<<" is required..";
  trace.info() <<std::endl;
  exit(1);
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  using namespace Z3i;

  typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
  typedef KSpace::SurfelSet SurfelSet;
  MySurfelAdjacency surfAdj( true ); // interior in all directions.
  

  // parse command line ----------------------------------------------
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("polynomial,P", po::value<std::string>(), "The implicit 3d polynomial (like x^2+2*y^2-z^3)")
    ("lower_x,x",  po::value<double>()->default_value(-1.0), "x-coordinate of lower bound" )
    ("lower_y,y",  po::value<double>()->default_value(-1.0), "y-coordinate of lower bound" )
    ("lower_z,z",  po::value<double>()->default_value(-1.0), "z-coordinate of lower bound" )
    ("upper_x,X",  po::value<double>()->default_value(1.0), "x-coordinate of upper bound" )
    ("upper_y,Y",  po::value<double>()->default_value(1.0), "y-coordinate of upper bound" )
    ("upper_z,Z",  po::value<double>()->default_value(1.0), "z-coordinate of upper bound" )
    ("grid,g",  po::value<double>()->default_value(0.01), "grid step" )
    ("vol,v", po::value<std::string>(), "The .vol file defining the shape implicity as min < value <= max")
    ("min,m", po::value<int>(), "The minimum threshold defining the shape in the .vol file")
    ("max,M", po::value<int>(), "The maximum threshold defining the shape in the .vol file")
    ("idx,i", po::value<int>()->default_value(0), "The index of the starting surfel")
    ("nu_p,p", po::value<int>()->default_value(1), "The integral numerator of the width nu")
    ("nu_q,q", po::value<int>()->default_value(1), "The integral denominator of the width nu")
    ("output,o",   po::value<std::string>()->default_value("surface.xyz"), "Base name of the file containing the surface sampling (format per line is x y z nx ny nz )" )
    ("surface,s",   po::value<std::string>()->default_value("surface.off"), "Base name of the file containing the surface triangulation (off format)" );

  bool parseOK = true;
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);  
  } catch(const std::exception& ex){
    parseOK = false;
    trace.error() << "Error checking program options: "<< ex.what()<< std::endl;
  }
  po::notify(vm);    
  if ( ! parseOK || vm.count("help") || ( argc <= 1 ) )
    {
      trace.info() << "View nu-convex sets either along the digitization of a polynomial surface or along the boundary of a digital shape contained in a .vol file." <<std::endl 
                   << "Basic usage: "<<std::endl
                   << "\t nu-convex-viewer -p <polynomial> [otherOptions]"<<std::endl
                   << "\t nu-convex-viewer -v <filename.vol> [otherOptions]"<<std::endl
                   << general_opt << "\n";
      return 0;
    }

  //Parse options
  if ( ( ! vm.count( "polynomial" ) )
       && ( ! vm.count( "vol" ) ) ) 
    missingParam("--polynomial or --vol");

  unsigned int idx = (unsigned int) vm["idx"].as<int>();
  unsigned int nup = (unsigned int) vm["nu_p"].as<int>();
  unsigned int nuq = (unsigned int) vm["nu_q"].as<int>();

  QApplication application(argc,argv);
  Viewer3D<Space,KSpace> viewer;

  if ( vm.count( "polynomial" ) )
    { // The user has specified a polynomial.
      std::string polynomialName = vm["polynomial"].as<std::string>();
      trace.beginBlock( "Making polynomial surface." );
      typedef Space::RealPoint RealPoint;
      typedef Space::RealVector RealVector;
      typedef RealPoint::Coordinate Scalar;
      typedef MPolynomial<3, Scalar> Polynomial3;
      typedef MPolynomialReader<3, Scalar> Polynomial3Reader;
      typedef ImplicitPolynomial3Shape<Space> ImplicitShape;
      typedef GaussDigitizer<Space,ImplicitShape> DigitalShape; 
      typedef DigitalShape::PointEmbedder DigitalEmbedder;
      typedef SetOfSurfels< KSpace, SurfelSet > MySetOfSurfels;
      typedef DigitalSurface< MySetOfSurfels > MyDigitalSurface;
      typedef DigitalSurface2InnerPointFunctor<MyDigitalSurface> VertexEmbedder;
      Polynomial3 P;
      Polynomial3Reader reader;
      std::string::const_iterator iter 
        = reader.read( P, polynomialName.begin(), polynomialName.end() );
      if ( iter != polynomialName.end() )
        {
          trace.error() << "ERROR when reading polynomial: read only <" 
                        << polynomialName.substr( 0, iter - polynomialName.begin() )
                        << ">, created polynomial is P=" << P << std::endl;
          return 2;
        }
      trace.info() << "P( X_0, X_1, X_2 ) = " << P << std::endl;
      ImplicitShape ishape( P );
      DigitalShape dshape;
      dshape.attach( ishape );
      RealPoint p1( vm["lower_x"].as<double>(),
                    vm["lower_y"].as<double>(),
                    vm["lower_z"].as<double>() );
      RealPoint p2( vm["upper_x"].as<double>(),
                    vm["upper_y"].as<double>(),
                    vm["upper_z"].as<double>() );
      Scalar step = vm["grid"].as<double>();
      dshape.init( RealPoint( p1 ), RealPoint( p2 ), step );
      Domain domain = dshape.getDomain();
      trace.endBlock();
      // Construct the Khalimsky space from the image domain
      KSpace K;
      // NB: it is \b necessary to work with a \b closed cellular space
      // since umbrellas use separators and pivots, which must exist for
      // arbitrary surfels.
      bool space_ok = K.init( domain.lowerBound(), 
                              domain.upperBound(), true // necessary
                              );
      if (!space_ok)
        {
          trace.error() << "ERROR in the Khamisky space construction." << std::endl;
          return 3;
        }
      viewer.setKSpace( K );
      viewer.show(); 

      trace.beginBlock( "Extracting boundary by scanning the space. " );
      MySetOfSurfels theSetOfSurfels( K, surfAdj );
      Surfaces<KSpace>::sMakeBoundary( theSetOfSurfels.surfelSet(),
                                       K, dshape,
                                       domain.lowerBound(),
                                       domain.upperBound() );
      MyDigitalSurface digSurf( theSetOfSurfels );
      trace.info() << "Digital surface has " << digSurf.size() << " surfels."
                   << std::endl;
      trace.endBlock();
      VertexEmbedder embedder( digSurf );

      viewNuConvex( viewer, digSurf, embedder, idx, nup, nuq );
    }
  else if ( vm.count( "vol" ) )
    { // The user has specified a .vol file.
      typedef ImageSelector < Domain, int>::Type Image;
      //typedef DigitalSetBoundary<KSpace, DigitalSet > MyDigitalSurfaceContainer;
      typedef SetOfSurfels< KSpace, SurfelSet > MyDigitalSurfaceContainer;
      typedef DigitalSurface< MyDigitalSurfaceContainer > MyDigitalSurface;
      // typedef ImplicitDigitalSurface<KSpace, DigitalSet > MyDigitalSurfaceContainer;
      // typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
      typedef DigitalSurface2InnerPointFunctor<MyDigitalSurface> VertexEmbedder;
      std::string inputFilename = vm["vol"].as<std::string>();
      unsigned int minThreshold = (unsigned int) vm["min"].as<int>();
      unsigned int maxThreshold = (unsigned int) vm["max"].as<int>();
      trace.beginBlock( "Reading vol file into an image." );
      Image image = VolReader<Image>::importVol(inputFilename);
      DigitalSet set3d (image.domain());
      SetFromImage<DigitalSet>::append<Image>(set3d, image, 
                                              minThreshold, maxThreshold);
      trace.endBlock();
      trace.beginBlock( "Construct the Khalimsky space from the image domain." );
      KSpace K;
      bool space_ok = K.init( image.domain().lowerBound(), 
                               image.domain().upperBound(), true );
      if (!space_ok)
        {
          trace.error() << "Error in the Khamisky space construction."<<std::endl;
          return 2;
        }
      trace.endBlock();
      viewer.setKSpace( K );
      viewer.show(); 
      trace.beginBlock( "Set up digital surface." );
      MyDigitalSurfaceContainer surfContainer( K, surfAdj );
      MyDigitalSurface::Vertex bel = Surfaces<KSpace>::findABel( K, set3d, 100000 );
      Surfaces<KSpace>::trackBoundary( surfContainer.surfelSet(),
                                       K, surfAdj, 
                                       set3d, bel );
      MyDigitalSurface digSurf( surfContainer );
      // MyDigitalSurfaceContainer* ptrSurfContainer = 
      //   new MyDigitalSurfaceContainer( K, set3d, surfAdj, bel );
      // MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
      trace.info() << "Digital surface has " << digSurf.size() << " surfels."
                   << std::endl;
      trace.endBlock();
      VertexEmbedder embedder( digSurf );

      viewNuConvex( viewer, digSurf, embedder, idx, nup, nuq );
    }

  return application.exec();
}
