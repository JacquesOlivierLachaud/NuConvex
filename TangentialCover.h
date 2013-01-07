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
 * @file TangentialCover.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/12/09
 *
 * Header file for module TangentialCover.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(TangentialCover_RECURSES)
#error Recursive header files inclusion detected in TangentialCover.h
#else // defined(TangentialCover_RECURSES)
/** Prevents recursive inclusion of headers. */
#define TangentialCover_RECURSES

#if !defined TangentialCover_h
/** Prevents repeated inclusion of headers. */
#define TangentialCover_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/base/Lambda2To1.h"
#include "DGtal/kernel/SquaredEuclideanDistance.h"
#include "DGtal/topology/DistanceVisitor.h"
#include "MaximalPlaneSummary.h"
#include "NuConvexSet.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class TangentialCover
  /**
     Description of template class 'TangentialCover' <p> \brief Aim:
     Represents the set of maximal planes at a given scale.
   
     @tparam TDigitalSurface the type of digital surface.
     @tparam TVertex2PointFunctor a functor DigitalSurface::Vertex -> DigitalSurface::Point
     @tparam TInternalInteger the type of integer used by the plane recognition algorithm.

     @todo Does not require a digital surface, a graph is
     sufficient. However, a visitor and an embedder should then be
     provided. This version is for test purposes.
   */
  template < typename TDigitalSurface,
             typename TVertex2PointFunctor,
             typename TInternalInteger = DGtal::int64_t >
  class TangentialCover
  {
    // ----------------------- public types ------------------------------
  public:
    typedef TDigitalSurface DigitalSurface;
    typedef TVertex2PointFunctor Vertex2PointFunctor;
    typedef TInternalInteger InternalInteger;
    typedef DigitalSurface Graph;
    typedef typename DigitalSurface::KSpace KSpace;
    typedef typename DigitalSurface::Vertex Vertex;
    typedef typename KSpace::Space Space;
    typedef MaximalPlaneSummary<Space> MPS;
    typedef DGtal::uint32_t Index;
    typedef typename Space::Point Point;
    typedef typename Space::RealPoint RealPoint;
    typedef typename Space::RealVector RealVector;
    typedef typename RealPoint::Coordinate Scalar;
    typedef SquaredEuclideanDistance<RealPoint> SqED;
    typedef Lambda2To1<SqED, RealPoint, RealPoint, Scalar> SqEDToPoint;
    typedef Composer<Vertex2PointFunctor, SqEDToPoint, Scalar> VertexFunctor;
    typedef DistanceVisitor< Graph, VertexFunctor > Visitor;
    typedef NuConvexSet< Space, Visitor, 
                         Vertex2PointFunctor, DGtal::int64_t > MyNuConvexSet;

    enum AveragingMode { SimpleAveraging, 
			 DistanceAveraging, 
			 RadiusAndDistanceAveraging,
			 InOutAveraging,
			 MaxProjectedPlane };
    /**
       Used for computing maximal plane summary.
    */
    struct SurfelAreaEstimator
    {
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

    /** 
        First maximal plane is the one with greatest area.
    */
    struct MPSAreaComparator {
      inline
      bool operator()( const MPS & mps1, const MPS & mps2 ) const
      {
        return mps1.projectedArea < mps2.projectedArea;
      }
    };

    
    /** 
        The indices of maximal planes are stored in a vector. The
        invariant is that, when non empty, the first index points
        always to the index of the maximal plane which is smallest
        according to the comparator.
    */
    typedef std::vector<Index> MaximalPlaneSummaryIndices;

    typedef std::map<Vertex, MaximalPlaneSummaryIndices> MapVertex2MPSI;
    typedef std::vector<MPS> MaximalPlaneSummaryTable;

    typedef std::set<Index> MaximalPlaneIndices;
    typedef typename MaximalPlaneIndices::const_iterator MaximalPlaneIndicesConstIterator;
    typedef std::map<Vertex, MaximalPlaneIndices> MapVertex2MPI;


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~TangentialCover();

    /**
      Constructor. The object is not valid.
     */
    TangentialCover();

    /**
       @return a const-reference to the associated digital surface.
    */
    const DigitalSurface & surface() const;

    /**
       @param v any vertex of the digital surface.

       @return the digital point associated with \a v.
       @see myVtx2PtFct
    */
    Point coordinates( const Vertex & v ) const;

    /**
       Initializes essential parameters for the tangential cover. Call clear(). 

       @param digSurf the digital surface where the tangential cover
       is computed (aliased).

       @param v2pFct the functor that associates a digital coordinate
       to any vertex (surfel) of the digital surface \a digSurf
       (aliased).

       @param p the numerator of the width given as the irreducible fraction p/q
       @param q the denominator of the width given as the irreducible fraction p/q

       @param diameter the diameter of the digital surface (necessary
       for COBA digital plane algorithm).
    */

    void init( const DigitalSurface & digSurf,
               const Vertex2PointFunctor & v2pFct,
               InternalInteger p, InternalInteger q,
               InternalInteger diameter );

    /**
       Computes all maximal planes on the surface. For each vertex,
       memorizes at most \a nbMaxPerVertex planes (the ones with
       maximal projected area). Calls clear().

       @param nbMaxPerVertex the maximal number of maximal planes
       stored per vertex.

       @param extensionMode extMode when 'true', the nu-convex set
       considers first vertices that are in the current nu-width
       plane, and only after it tries to find a nu-width plane that
       contains also the other potential vertices ; when 'false', the
       nu-convex set considers all potential vertices at the same time
       to find a nu-width plane.
    */
    void compute( unsigned int nbMaxPerVertex,
                  bool extensionMode = false );

    /**
       Computes all maximal planes on the surface. For each vertex,
       memorizes at most \a nbMaxPerVertex planes (the ones with
       maximal projected area). This method takes care that no two
       planes are identical.  Calls clear().

       @param nbMaxPerVertex the maximal number of maximal planes
       stored per vertex.

       @param extensionMode extMode when 'true', the nu-convex set
       considers first vertices that are in the current nu-width
       plane, and only after it tries to find a nu-width plane that
       contains also the other potential vertices ; when 'false', the
       nu-convex set considers all potential vertices at the same time
       to find a nu-width plane.
    */
    void computeOnce( unsigned int nbMaxPerVertex, 
                      bool extensionMode = false );


    /**
       Resets all computations.
    */
    void clear();

    /**
       Combines the normals of all maximal planes containing the vertex. 
       @param n (modified) the normal estimated at [p].
       @param p any vertex.
       @param nd the mode chosen for computing normals.
    */
    void getEstimatedNormal( RealVector & n, Vertex p,
			     AveragingMode nd = SimpleAveraging ) const;

    /**
       Combines the normals of all maximal planes containing the vertex. 
       @param n (modified) the normal estimated at [p].
       @param p any vertex.
       @param coefs the averaging coefficients.
       @see getAveragingCoefficients
    */
    void getEstimatedNormal( RealVector & n, Vertex p,
			     const std::vector<Scalar> & coefs ) const;

    /**
       Computes the averaging coefficients for the given vertex [p]
       according to the choosen averaging mode.

       @param coefs (modified) contains the averaging coefficients (sum is 1.0)
       @param p any vertex.
       @param nd the mode chosen for averaging.
    */
    void getAveragingCoefficients( std::vector<Scalar> & coefs, 
                                   Vertex p, 
                                   AveragingMode nd ) const;

    void displayPlanes( std::ostream & out ) const;

    /**
       @param vtx any vertex of the surface.

       @return a const iterator pointing on the first element of the
       list of maximal planes of vertex \a vtx.

       NB: non-const method because of std::map::operator[].
    */
    MaximalPlaneIndicesConstIterator begin( const Vertex & vtx );
    /**
       @param vtx any vertex of the surface.

       @return a const iterator pointing after the last element of the
       list of maximal planes of vertex \a vtx.

       NB: non-const method because of std::map::operator[].
    */
    MaximalPlaneIndicesConstIterator end( const Vertex & vtx );

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

    /// The associated digital surface onto which the tangential cover
    /// is computed.
    const DigitalSurface* myDigSurf;
    const Vertex2PointFunctor* myVtx2PtFct;
    InternalInteger myP;
    InternalInteger myQ;
    InternalInteger myDiameter;

    MPSAreaComparator myAreaComparator;
    MapVertex2MPSI myMapVtx2MPSI;
    MaximalPlaneSummaryTable myMPSTable;
    unsigned int myNbMaxPerVertex;
    
    /// Stores the map MP index -> covering MP index.
    std::vector< Index > myPlane2CoveringPlane;
    /// Stores for each MP its most centered vertices.
    std::vector< std::list< Vertex > > myMostCenteredVertices;

    /// Stores completely the mapping vtx -> set of its MP.
    MapVertex2MPI myMapVtx2MPI;
    
    
    // ------------------------- Hidden services ------------------------------
  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    TangentialCover ( const TangentialCover & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    TangentialCover & operator= ( const TangentialCover & other );

    // ------------------------- Internals ------------------------------------
  private:

    /**
       Add the maximal plane index \a mp to the list of maximal planes
       covering vertex \a vtx.

       @param vtx a valid vertex.
       @param mp a valid maximal plane index.

       @return 'true' iff the maximal plane was effectively added to
       the list of maximal planes of vertex \a vtx. It is 'false' when
       the maximal plane is too small wrt to the maximal planes
       covering \a vtx.
    */
    bool addMPSToVertex( Vertex vtx, Index mp );

    /**
       @param mp a valid maximal plane index.
       @return the index of the biggest maximal plane covering \a mp.
    */
    Index father( Index mp );

    void setFather( Index mp, Index f );

    void cleanMaximalPlaneIndices( MaximalPlaneIndices & mpi );

  }; // end of class TangentialCover


  /**
   * Overloads 'operator<<' for displaying objects of class 'TangentialCover'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'TangentialCover' to write.
   * @return the output stream after the writing.
   */
  template <typename TDigitalSurface, typename TVertex2PointFunctor,
            typename TInternalInteger >
  std::ostream&
  operator<< ( std::ostream & out, 
               const TangentialCover<TDigitalSurface, TVertex2PointFunctor, TInternalInteger> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "TangentialCover.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined TangentialCover_h

#undef TangentialCover_RECURSES
#endif // else defined(TangentialCover_RECURSES)
