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

#include "GrowingDiscGeometryBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
GrowingDiscGeometryBoundaryCondition<DIM>::GrowingDiscGeometryBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                                      c_vector<double, DIM> centre,
                                                                      double radius,
                                                                      double growthRate)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mCentreOfSphere(centre),
      mRadiusOfSphere(radius),
      mGrowthRate(growthRate)
{
    assert(mRadiusOfSphere > 0.0);

    if (DIM == 1)
    {
        EXCEPTION("This boundary condition is not implemented in 1D.");
    }
}

template<unsigned DIM>
const c_vector<double, DIM>& GrowingDiscGeometryBoundaryCondition<DIM>::rGetCentreOfSphere() const
{
    return mCentreOfSphere;
}

template<unsigned DIM>
double GrowingDiscGeometryBoundaryCondition<DIM>::GetRadiusOfSphere() const
{
    return mRadiusOfSphere;
}

template<unsigned DIM>
void GrowingDiscGeometryBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{

    double time = SimulationTime::Instance()->GetTime();

    double current_radius = mRadiusOfSphere + mGrowthRate * time;

    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Find the radial distance between this cell and the surface of the sphere
        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        double radius = norm_2(cell_location - mCentreOfSphere);
        
        // If the cell is ouside the surface of the sphere...
        if (radius > current_radius)
        {
            double max_jiggle = 0.001;

            // ...move the cell back onto the surface of the sphere
            c_vector<double, DIM> location_within_sphere = mCentreOfSphere + (1-max_jiggle*RandomNumberGenerator::Instance()->ranf())*current_radius*(cell_location - mCentreOfSphere)/radius;

            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            p_node->rGetModifiableLocation() = location_within_sphere;
        }
    }
}

template<unsigned DIM>
bool GrowingDiscGeometryBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    // // Iterate over the cell population
    // for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
    //      cell_iter != this->mpCellPopulation->End();
    //      ++cell_iter)
    // {
    //     // Find the radial distance between this cell and the surface of the sphere
    //     c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
    //     double radius = norm_2(cell_location - mCentreOfSphere);

    //     // If the cell is outside the surface of the sphere...
    //     if (radius > mRadiusOfSphere)
    //     {
    //         // ...then the boundary condition is not satisfied
    //         condition_satisfied = false;
    //         break;
    //     }
    // }
    return condition_satisfied;
}

template<unsigned DIM>
void GrowingDiscGeometryBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CentreOfSphere>";
    for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mCentreOfSphere[index] << ",";
    }
    *rParamsFile << mCentreOfSphere[DIM-1] << "</CentreOfSphere>\n";

    *rParamsFile << "\t\t\t<RadiusOfSphere>" << mRadiusOfSphere << "</RadiusOfSphere>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class GrowingDiscGeometryBoundaryCondition<1>;
template class GrowingDiscGeometryBoundaryCondition<2>;
template class GrowingDiscGeometryBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(GrowingDiscGeometryBoundaryCondition)
