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
 * @file DigitalSurface2InnerPointFunctor.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/02/28
 *
 * Header file for module DigitalSurface2InnerPointFunctor.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DigitalSurface2InnerPointFunctor_RECURSES)
#error Recursive header files inclusion detected in DigitalSurface2InnerPointFunctor.h
#else // defined(DigitalSurface2InnerPointFunctor_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DigitalSurface2InnerPointFunctor_RECURSES

#if !defined DigitalSurface2InnerPointFunctor_h
/** Prevents repeated inclusion of headers. */
#define DigitalSurface2InnerPointFunctor_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/topology/CCellularGridSpaceND.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class DigitalSurface2InnerPointFunctor
/**
   Description of class 'DigitalSurface2InnerPointFunctor' <p>

   \brief Aim: A trivial point embedder for digital surfaces, which
   corresponds to the canonic injection of each surfel to its the
   inner point.

   @tparam TDigitalSurface the type of digital surface where the point
   embedder works.
 */
  template <typename TDigitalSurface>
  struct DigitalSurface2InnerPointFunctor
  {
  public:
    typedef DigitalSurface2InnerPointFunctor<TDigitalSurface> Self;

    typedef TDigitalSurface Surface;
    typedef typename Surface::KSpace KSpace;
    BOOST_CONCEPT_ASSERT(( concepts::CCellularGridSpaceND<KSpace> ));
    typedef typename KSpace::SCell SCell;
    typedef typename KSpace::Space Space;
    typedef typename Space::Point Point;
    typedef SCell Argument;
    typedef Point Value;
    typedef typename Space::Integer Integer;

    // ----------------------- Standard services ------------------------------
  public:
    /**
       Destructor. Nothing special.
    */
    ~DigitalSurface2InnerPointFunctor();

    /**
       Default constructor. The object is not valid.
    */
    DigitalSurface2InnerPointFunctor();

    /**
       Constructor from surface. 
    */
    DigitalSurface2InnerPointFunctor( const Surface & aSurface );

    /**
       Copy constructor.
       @param other the object to clone.
    */
    DigitalSurface2InnerPointFunctor( const Self & other );

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
       Map a signed cell to its corresponding inner point in the digital
       space.
       
       @param cell any signed cell in the digital space.
       @return its inner point embedding in the digital space.
    */
    Point embed( const SCell & cell ) const;

    /**
       Map a signed cell to its corresponding inner point in the digital
       space.
       
       @param cell any signed cell in the digital space.
       @return its inner point embedding in the digital space.
    */
    Point operator()( const SCell & cell ) const;

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

}; // end of class DigitalSurface2InnerPointFunctor


/**
 * Overloads 'operator<<' for displaying objects of class 'DigitalSurface2InnerPointFunctor'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'DigitalSurface2InnerPointFunctor' to write.
 * @return the output stream after the writing.
 */
  template <typename TDigitalSurface>
  std::ostream&
  operator<< ( std::ostream & out, const DigitalSurface2InnerPointFunctor<TDigitalSurface> & object );


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DigitalSurface2InnerPointFunctor.ih"


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DigitalSurface2InnerPointFunctor_h

#undef DigitalSurface2InnerPointFunctor_RECURSES
#endif // else defined(DigitalSurface2InnerPointFunctor_RECURSES)
