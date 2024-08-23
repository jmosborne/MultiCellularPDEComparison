/*

Copyright (c) 2005-2024, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef RADIALGROWTHFORCE_HPP_
#define RADIALGROWTHFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"

/**
 * TODO
 *
 * This class works with all off-lattice cell populations.
 */
template<unsigned DIM>
class RadialGrowthForce : public AbstractForce<DIM>
{
private :

    /**
     * Absolute temperature (in Kelvin).
     */
    double mRadialVelocity;

    /**
     * Viscosity of media. We assume that this is measured in units of
     * kg (10 microns)^(-1) h^(-1), and that cell diameters are scaled with
     * a characteristic length of 10 microns.
     */
    double mViscosity;

    /**
     * The Boltzmann constant, which takes the value 4.97033568e-7
     * (10 microns)^2 kg h^(-2) K^(-1).
     */
    static const double msBoltzmannConstant;

    /**
     * Archiving.
     */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mRadialVelocity;
        archive & mViscosity;
    }

public :

    /**
     * Constructor.
     */
    RadialGrowthForce(double radialVelocity = 0.5);

    /**
     * Destructor.
     */
    ~RadialGrowthForce();

    /**
     * Set the ...
     * 
     * @param radialVelocity the velocityof the 
     */
    void SetRadialVelocity(double radialVelocity);

    /**
     * Get the ...
     *
     * @return mRadialVelocity
     */
    double GetRadialVelocity();

    /**
     * Overridden AddForceContribution() method.
     * Note that this method requires cell/node radii to be set.
     *
     * @param rCellPopulation reference to the tissue
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RadialGrowthForce)

#endif /*RADIALGROWTHFORCE_HPP_*/
