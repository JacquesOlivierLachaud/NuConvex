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
 * @file NuConvexSet.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/12/06
 *
 * Header file for module NuConvexSet.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(NuConvexSet_RECURSES)
#error Recursive header files inclusion detected in NuConvexSet.h
#else // defined(NuConvexSet_RECURSES)
/** Prevents recursive inclusion of headers. */
#define NuConvexSet_RECURSES

#if !defined NuConvexSet_h
/** Prevents repeated inclusion of headers. */
#define NuConvexSet_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
//#define CHORD
#include <iostream>
#include "DGtal/base/Common.h"
#ifdef CHORD
#include "DGtal/geometry/surfaces/ChordGenericNaivePlane.h"
#else
#include "DGtal/geometry/surfaces/COBAGenericNaivePlane.h"
#endif
#include "BasicHPolytopeND.h"
#include "MaximalPlaneSummary.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class NuConvexSet
  /**
     Description of template class 'NuConvexSet' <p> \brief Aim:
     Nu-convex sets aim at mimicking maximal digital straight segments
     on surfaces (in fact, any graph) in Z^3. A nu-convex set has a
     starting vertex, called its center. Then it grows isotropically
     from the center and determines at each step if its width is no
     greater than parameter nu. If at some point, added vertices
     induce a too long width, then problematic vertices are added as
     orthogonal constraints. The nu-convex set continues to grow but
     ignores vertices that are outside the current set of constraints
     (which is by construction a convex set). When no other vertices
     can be added, the nu-convex set is formed.
   */
  template < typename TSpace, typename TVisitor,
	     typename TVertex2PointFunctor, 
	     typename TInternalInteger = DGtal::int64_t >
  class NuConvexSet
  {
    // ----------------------- public types ------------------------------
  public:
    typedef TSpace Space;
    typedef TVisitor Visitor;
    typedef TInternalInteger InternalInteger;
    typedef TVertex2PointFunctor Vertex2PointFunctor;
    typedef typename Visitor::Vertex Vertex;

    BOOST_STATIC_ASSERT(( ConceptUtils::SameType
			   < Vertex, typename Vertex2PointFunctor::Argument >::value ));
#ifdef CHORD
    typedef ChordGenericNaivePlane< typename Space::Point, InternalInteger > GenericNaivePlane;
#else
    typedef COBAGenericNaivePlane< Space, InternalInteger > GenericNaivePlane;
#endif
    typedef typename Vertex2PointFunctor::Value Point;
    typedef typename Point::Coordinate Scalar;
    typedef Point Vector;
    typedef BasicHPolytopeND<Vector> HPolytope;
    typedef typename HPolytope::ClosedHalfSpace ClosedHalfSpace;
    typedef typename Visitor::Node Node;     // vertex
    typedef typename Visitor::Scalar Distance; // distance

    typedef std::vector<Vertex> Container;
    typedef typename Container::const_iterator ConstIterator;
    typedef typename Container::size_type Size;

    class VertexPolytopePredicateAdapter
    {
    public:
      typedef HPolytope PointPredicate;
      typedef typename Vertex2PointFunctor::Argument Vertex;
      typedef typename Vertex2PointFunctor::Value Point;

      inline
      VertexPolytopePredicateAdapter( const PointPredicate & pPred,
                                      const Vertex2PointFunctor & vtx2pt )
        : myPPred( &pPred ), myVtx2PointFct( &vtx2pt )
      {}

      inline
      VertexPolytopePredicateAdapter( const VertexPolytopePredicateAdapter & other )
        : myPPred( other.myPPred ), myVtx2PointFct( other.myVtx2PointFct )
      {}

      inline
      VertexPolytopePredicateAdapter &
      operator=( const VertexPolytopePredicateAdapter & other )
      {
        if ( this != &other )
          {
            myPPred = other.myPPred;
            myVtx2PointFct = other.myVtx2PointFct;
          }
        return *this;
      }
      
      inline
      bool operator()( const Vertex & v ) const
      {
        return (*myPPred)( (*myVtx2PointFct)( v ) ); 
      }

      inline
      void swap( VertexPolytopePredicateAdapter & other )
      {
        std::swap( myPPred, other.myPPred );
        std::swap( myVtx2PointFct, other.myVtx2PointFct );
      }

    private:
      const PointPredicate* myPPred;
      const Vertex2PointFunctor* myVtx2PointFct;
      
    };

  private:
    enum VertexState { Rejected, ExtendedAsIs, Extendable }; 
    typedef std::pair<Node,VertexState> NodeAndState;

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~NuConvexSet();

    /**
       Constructor.
       @param visitor the visitor is duplicated.
       @param vtx2pt the vtx2pt is aliased.
     */
    NuConvexSet( const Visitor & visitor, const Vertex2PointFunctor & vtx2pt );

    /**
       Default extension mode is true.
       
       @param extMode when 'true', the nu-convex set considers first
       vertices that are in the current nu-width plane, and only after
       it tries to find a nu-width plane that contains also the other
       potential vertices ; when 'false', the nu-convex set considers
       all potential vertices at the same time to find a nu-width
       plane.*/ 
    void setExtensionMode( bool extMode );

    /**
       Uses the given Vertex2PointFunctor to return the vertex coordinates.
       @param v any vertex.
    */
    Point coordinates( const Vertex & v ) const;

    /**
       The start vertex is given by the visitor.
       nu = p/q
       @param diameter the diameter is necessary for COBA digital plane algorithm.
    */
    void init( InternalInteger p, InternalInteger q,
	       InternalInteger diameter );

    bool compute( Scalar distanceUpperBound = -1 );

    Size size() const;
    bool empty() const;
    ConstIterator begin() const;
    ConstIterator end() const;

    /// Sorts nu-convex set vertices in the container.
    void sort();

    /// Clears most of "internal" data to occupy as less memory as
    /// possible. Keeps vertices and geometric information. However,
    /// you may not call compute() again to go on with further computation.
    void compact();

    /**
       @param mps the object that stores the geometric summary of the nu-convex set.

       @param vArea an instance of VertexAreaEstimator that is a
       functor (Vertex,RealVector)->Scalar giving the estimated projected area of
       the given vertex under the given angle.

       @tparam VertexAreaEstimator A model of a functor (Vertex,RealVector)->Scalar.
    */
    template <typename VertexAreaEstimator>
    void summarize( MaximalPlaneSummary<Space> & mps,
                    const VertexAreaEstimator & vArea ) const;

    /**
       Exchange the data of 'this' with the data of 'other'.
       @param other (modified) the other instance which is exchanged with 'this'.
    */
    void swap( NuConvexSet & other );


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

    Visitor* myVisitor;
    const Vertex2PointFunctor* myVtx2PointFct;
    GenericNaivePlane myPlane;
    HPolytope myPolytope;
    VertexPolytopePredicateAdapter myVertexPolytopePredicate;
    Container myVertices;
    bool myExtMode;
  public:
    // Container myRejectedVertices;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    NuConvexSet();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    NuConvexSet ( const NuConvexSet & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    NuConvexSet & operator= ( const NuConvexSet & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class NuConvexSet


  /**
   * Overloads 'operator<<' for displaying objects of class 'NuConvexSet'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'NuConvexSet' to write.
   * @return the output stream after the writing.
   */
  template < typename TSpace, typename TVisitor,
	     typename TVertex2PointFunctor, 
	     typename TInternalInteger >
  std::ostream&
  operator<< ( std::ostream & out, 
	       const NuConvexSet< TSpace, TVisitor, TVertex2PointFunctor, TInternalInteger > & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "NuConvexSet.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined NuConvexSet_h

#undef NuConvexSet_RECURSES
#endif // else defined(NuConvexSet_RECURSES)
