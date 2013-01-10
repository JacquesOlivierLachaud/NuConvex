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
 * @file FuzzyPartition.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2013/01/10
 *
 * Header file for module FuzzyPartition.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(FuzzyPartition_RECURSES)
#error Recursive header files inclusion detected in FuzzyPartition.h
#else // defined(FuzzyPartition_RECURSES)
/** Prevents recursive inclusion of headers. */
#define FuzzyPartition_RECURSES

#if !defined FuzzyPartition_h
/** Prevents repeated inclusion of headers. */
#define FuzzyPartition_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <boost/pending/disjoint_sets.hpp>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class FuzzyPartition
  /**
     Description of template class 'FuzzyPartition' <p>
     \brief Aim:

     @tparam the type of graph (model of CUndirectedSimpleGraph).
     @tparam a functor Vertex x Vertex -> [0:1] giving the similarity
     between the two vertices.
   */
  template <typename TGraph, typename TSimilarity>
  class FuzzyPartition
  {
    // ----------------------- Standard services ------------------------------
  public:
    typedef TGraph Graph;
    typedef TSimilarity Similarity;
    typedef typename Graph::Vertex Vertex;
    typedef typename Graph::ConstIterator ConstIterator;
    typedef std::map< Vertex, unsigned int> RankMap;
    typedef std::map< Vertex, Vertex> ParentMap;

    // Adapter to map: Vertex -> unsigned int, used by disjoint_sets.
    // Adapter to map: Vertex -> Vertex, used by disjoint_sets.
    typedef boost::associative_property_map<RankMap> RankAssociativePropertyMap;
    typedef boost::associative_property_map<ParentMap> ParentAssociativePropertyMap;
    typedef boost::disjoint_sets< RankAssociativePropertyMap, ParentAssociativePropertyMap > DisjointSets;

    /**
     * Destructor.
     */
    ~FuzzyPartition();

    /**
       Constructor. Objects are referenced.

       @param graph the graph to partition.

       @param similarity the object defining the similarity between
       two vertices, i.e. a functor Vertex x Vertex -> [0,1] where 1
       is maximum similarity.
    */
    FuzzyPartition( const Graph & graph, const Similarity & similarity );

    /**
       Partition the graph. All vertices such that there is a path of
       similarity greater than threshold between them are in the same
       component.

       @param threshold the minimum similarity (0 means only one component per graph component).

       @return the number of connected components.
    */
    unsigned int partition( double threshold );

    /**
       All vertices in the same class have the same representative
       vertex.

       @param v any valid vertex.
       @return the representative vertex of \a v.
    */
    Vertex representative( const Vertex & v );

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

    const Graph & myGraph;
    const Similarity & mySimilarity;
    RankMap myRankMap;
    RankAssociativePropertyMap myRankAssociativePropertyMap;
    ParentMap myParentMap;
    ParentAssociativePropertyMap myParentAssociativePropertyMap;
    DisjointSets myDisjointSets;

    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    FuzzyPartition();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    FuzzyPartition ( const FuzzyPartition & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    FuzzyPartition & operator= ( const FuzzyPartition & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class FuzzyPartition


  /**
   * Overloads 'operator<<' for displaying objects of class 'FuzzyPartition'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'FuzzyPartition' to write.
   * @return the output stream after the writing.
   */
  template <typename TGraph, typename TSimilarity>
  std::ostream&
  operator<< ( std::ostream & out, const FuzzyPartition<TGraph,TSimilarity> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "FuzzyPartition.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FuzzyPartition_h

#undef FuzzyPartition_RECURSES
#endif // else defined(FuzzyPartition_RECURSES)
