/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file MetricCluster.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2013/01/08
 *
 * Header file for module MetricCluster.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(MetricCluster_RECURSES)
#error Recursive header files inclusion detected in MetricCluster.h
#else // defined(MetricCluster_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MetricCluster_RECURSES

#if !defined MetricCluster_h
/** Prevents repeated inclusion of headers. */
#define MetricCluster_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <set>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class MetricCluster
  /**
   * Description of template class 'MetricCluster' <p> \brief Aim:
   * Represents a partition of a list of data vectors into clusters,
   * such that clusters have a centroid and a distance is defined
   * between data vectors.
   */
  template <typename TVector, typename TSquaredDistance>
  class MetricCluster
  {
    // ----------------------- Standard services ------------------------------
  public:

    typedef TVector Vector;
    typedef TSquaredDistance SquaredDistance;
    typedef typename Vector::Component Scalar;

    /**
     * Destructor.
     */
    ~MetricCluster();

    /**
       Default constructor. The object is empty.
    */
    MetricCluster();

    /**
       Resets the object.
    */
    void clear();

    /**
       The number of data.
    */
    unsigned int N() const;

    /**
       The number of clusters.
    */
    unsigned int K() const;

    /**
       Initializes the distance functor and the list of data vectors.
       The initial clustering is only one class (labelled 0).
    */
    template <typename Iterator>
    void init( const SquaredDistance & sqDistance, 
               Iterator itb, Iterator ite );

    /**
       Given valid myClustering and myClusterSizes, computes myClusterCentroids.
    */
    void computeCentroids();

    /**
       Given valid myClusterCentroids, computes myClustering and myClusterSizes.
    */
    void computeClusters();

    /**
       Pick randomly \a k points as centroids of \a k clusters, then
       updates myClustering and myClusterSizes accordingly.
    */ 
    void randomClusters( unsigned int k );

    /**
       @return the energy (sum of squared distances between data and their cluster centroid).
    */
    Scalar energy() const;

    /**
       @return the L1-energy (sum of distances between data and their cluster centroid), E_K in some papers.
    */
    Scalar energyK() const;

    /**
       @return the maximum distance between centroids, D_K in some papers.
    */
    Scalar distanceK() const;

    /**
       @param v any vector
       @param k any valid cluster ( k < K()).
       @return the squared distance between \a v and the centroid of cluster \a k.
    */
    Scalar squaredDistance( const Vector & v, unsigned int k ) const;

    /**
       @param v any vector
       @param k any valid cluster ( k < K()).
       @return the distance between \a v and the centroid of cluster \a k.
    */
    Scalar distance( const Vector & v, unsigned int k ) const;

    /**
       Performs at most \a max iterations of Lloyd algorithm.
    */
    unsigned int Lloyd( unsigned int max );

    double indexI( double p = 2.0 ) const;
    double indexI1( double p = 2.0 ) const;
    double indexSigma() const;

    /**
       @param k any valid cluster ( k < K()).
       @return the maximal distance between a point in cluster \a k and the centroid of cluster \a k.
    */
    Scalar radius( unsigned int k ) const;

    /**
       If indexI() is maximized for two, you may check this to see if there is only one cluster.
       When called, the MetricCluster should have two clusters.
     */
    bool onlyOne( double upperBoundRatio = 3.0 ) const;

    /**
       @return (biased) variance estimation of cluster \a k.
    */
    Scalar variance( unsigned int k ) const;

    /**
       @return the similarity to the sum of normal laws. If similar,
       should be sqrt(2)/4 in 1D, 1/2 in 2D, sqrt(2)/2 in 3D.
    */
    Scalar normalSimilarity() const;

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
  private:
    // ------------------------- Private Datas --------------------------------
  private:

    SquaredDistance mySquaredDistance;
    std::vector< Vector > myData;
    std::vector< unsigned int > myClustering;
    std::vector< Vector > myClusterCentroids;
    std::vector< unsigned int > myClusterSizes;
    Scalar myE1;
    Scalar mySqE1;

    // ------------------------- Hidden services ------------------------------
  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    MetricCluster ( const MetricCluster & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    MetricCluster & operator= ( const MetricCluster & other );

    // ------------------------- Internals ------------------------------------
  private:

    void computeE1();

  }; // end of class MetricCluster


  /**
   * Overloads 'operator<<' for displaying objects of class 'MetricCluster'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'MetricCluster' to write.
   * @return the output stream after the writing.
   */
  template <typename TVector, typename TSquaredDistance>
  std::ostream&
  operator<< ( std::ostream & out, const MetricCluster<TVector,TSquaredDistance> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "MetricCluster.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MetricCluster_h

#undef MetricCluster_RECURSES
#endif // else defined(MetricCluster_RECURSES)
