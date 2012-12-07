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
 * @file InteriorSpelDigitalSurfaceEmbedder.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/02/28
 *
 * Header file for module InteriorSpelDigitalSurfaceEmbedder.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(InteriorSpelDigitalSurfaceEmbedder_RECURSES)
#error Recursive header files inclusion detected in InteriorSpelDigitalSurfaceEmbedder.h
#else // defined(InteriorSpelDigitalSurfaceEmbedder_RECURSES)
/** Prevents recursive inclusion of headers. */
#define InteriorSpelDigitalSurfaceEmbedder_RECURSES

#if !defined InteriorSpelDigitalSurfaceEmbedder_h
/** Prevents repeated inclusion of headers. */
#define InteriorSpelDigitalSurfaceEmbedder_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/topology/CCellularGridSpaceND.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class InteriorSpelDigitalSurfaceEmbedder
/**
   Description of class 'InteriorSpelDigitalSurfaceEmbedder' <p>

   \brief Aim: A trivial embedder for digital surfaces, which
   corresponds to the canonic injection of cell centroids into Rn.

   Model of CInteriorSpelDigitalSurfaceEmbedder (and thus of CSCellEmbedder).

   @tparam TDigitalSurface the type of digital surface where the embedder works.
 */
  template <typename TDigitalSurface>
  struct InteriorSpelDigitalSurfaceEmbedder
  {
  public:
    typedef InteriorSpelDigitalSurfaceEmbedder<TDigitalSurface> Self;

    typedef TDigitalSurface Surface;
    typedef typename Surface::KSpace KSpace;
    BOOST_CONCEPT_ASSERT(( CCellularGridSpaceND<KSpace> ));
    typedef typename KSpace::SCell SCell;
    typedef typename KSpace::Space Space;
    typedef typename Space::RealPoint RealPoint;
    typedef SCell Argument;
    typedef RealPoint Value;

    typedef typename Space::Integer Integer;
    typedef typename Space::Point Point;

    // ----------------------- Standard services ------------------------------
  public:
    /**
       Destructor. Nothing special.
    */
    ~InteriorSpelDigitalSurfaceEmbedder();

    /**
       Default constructor. The object is not valid.
    */
    InteriorSpelDigitalSurfaceEmbedder();

    /**
       Constructor from surface. 
    */
    InteriorSpelDigitalSurfaceEmbedder( const Surface & aSurface );

    /**
       Copy constructor.
       @param other the object to clone.
    */
    InteriorSpelDigitalSurfaceEmbedder( const Self & other );

    /**
       Assignment.
       @param other the object to clone.
       @return a reference to 'this'.
    */
    Self & operator=( const Self & other );

    /**
       @return the digital surface.
    */
    const Surface & surface() const;

    /**
       Map a signed cell to its corresponding point in the Euclidean
       space.
       
       @param cell any signed cell in the digital space.
       @return its canconical embedding in the Euclidean space.
    */
    RealPoint embed( const SCell & cell ) const;

    /**
       Map a signed cell to its corresponding point in the Euclidean
       space.
       
       @param cell any signed cell in the digital space.
       @return its canconical embedding in the Euclidean space.
    */
    RealPoint operator()( const SCell & cell ) const;

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
  protected:
    const Surface* mySurface;

    // ------------------------- Private Datas --------------------------------
private:

    // ------------------------- Hidden services ------------------------------
protected:


    // ------------------------- Internals ------------------------------------
private:

}; // end of class InteriorSpelDigitalSurfaceEmbedder


/**
 * Overloads 'operator<<' for displaying objects of class 'InteriorSpelDigitalSurfaceEmbedder'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'InteriorSpelDigitalSurfaceEmbedder' to write.
 * @return the output stream after the writing.
 */
  template <typename TDigitalSurface>
  std::ostream&
  operator<< ( std::ostream & out, const InteriorSpelDigitalSurfaceEmbedder<TDigitalSurface> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/kernel/InteriorSpelDigitalSurfaceEmbedder.ih"


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined InteriorSpelDigitalSurfaceEmbedder_h

#undef InteriorSpelDigitalSurfaceEmbedder_RECURSES
#endif // else defined(InteriorSpelDigitalSurfaceEmbedder_RECURSES)
