#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/SquaredEuclideanDistance.h"
#include "MetricCluster.h"

using namespace DGtal;

template <typename Cluster>
void
clustering( Cluster & cluster, const std::string & name )
{
  trace.beginBlock( name );
  for ( unsigned int k = 1; k < 4; ++k )
    {
      cluster.randomClusters( k );
      cluster.Lloyd( 100 );
      trace.info() << "I=" << cluster.indexI( 2 ) 
                   << " Is=" << cluster.indexSigma() 
                   << " " << cluster << std::endl;
      if ( k == 2 )
        {
          bool one = cluster.onlyOne( 3.0 );
          trace.info() << "Only One = " << ( one ? "true" : "false" ) << std::endl;
        }
      trace.info() << " normal_similarity=" << cluster.normalSimilarity()
                   << std::endl;
    }
  trace.endBlock();
}

int main( int argc, char** argv )
{
  typedef Z3i::Space Space;
  typedef Space::RealPoint RealPoint;
  typedef RealPoint RealVector;
  typedef SquaredEuclideanDistance<RealPoint> SqDistance;
  typedef MetricCluster< RealVector, SqDistance > Cluster;
  typedef std::vector< RealVector > Data;

  SqDistance distance;
  {
    Data data;
    data.push_back( RealVector( 1, 0, 0 ) );
    data.push_back( RealVector( 0, 0.2, 0.9 ) );
    data.push_back( RealVector( 0.9, 0.1, 0.1 ) );
    data.push_back( RealVector( 0.85, -0.2, 0.05 ) );
    data.push_back( RealVector( 0.1, -0.2, 0.95 ) );
    data.push_back( RealVector( 0.2, 0.05, 0.91 ) );
    Cluster cluster;
    cluster.init( distance, data.begin(), data.end() );
    clustering( cluster, "Clustering (expect 2)" );
  }
  {
    Data data;
    data.push_back( RealVector( 1, 0, 0 ) );
    data.push_back( RealVector( 1, 0.01, 0.005 ) );
    data.push_back( RealVector( 0.995, -0.01, -0.01 ) );
    data.push_back( RealVector( 0.98, 0.01, 0.05 ) );
    data.push_back( RealVector( 0.99, -0.02, 0.03 ) );
    data.push_back( RealVector( 0.97, -0.01, -0.02 ) );
    data.push_back( RealVector( 0.98, 0.02, 0.01 ) );
    data.push_back( RealVector( 0.99, 0.01, -0.03 ) );
    data.push_back( RealVector( 0.97, 0.01, 0.01 ) );
    Cluster cluster;
    cluster.init( distance, data.begin(), data.end() );
    clustering( cluster, "Clustering (expect 1)" );
  }
  return 0;
}
