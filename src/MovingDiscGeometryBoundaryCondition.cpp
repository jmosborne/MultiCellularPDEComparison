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

#include "MovingDiscGeometryBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

template<unsigned DIM>
MovingDiscGeometryBoundaryCondition<DIM>::MovingDiscGeometryBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                                      c_vector<double, DIM> centre,
                                                                      double radius)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mCentreOfSphere(centre),
      mRadiusOfSphere(radius)
{
    assert(mRadiusOfSphere > 0.0);
    
    if (DIM == 1)
    {
        EXCEPTION("This boundary condition is not implemented in 1D.");
    }
}

template<unsigned DIM>
const c_vector<double, DIM>& MovingDiscGeometryBoundaryCondition<DIM>::rGetCentreOfSphere() const
{
    return mCentreOfSphere;
}

template<unsigned DIM>
double MovingDiscGeometryBoundaryCondition<DIM>::GetRadiusOfSphere() const
{
    return mRadiusOfSphere;
}

template<unsigned DIM>
void MovingDiscGeometryBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    // How far into disc to define cells as on the boundary.
    double edge_distance = 0.1;


    double total_displacement = 0.0;
    unsigned num_boundary_cells = 0;

    // Iterate over the cell population to find what force applied to plane is 
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Find the radial distance between this cell and the surface of the sphere
        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        double radius = norm_2(cell_location - mCentreOfSphere);
        
        // If the cell is ouside the surface of the sphere...
        if (radius > mRadiusOfSphere - edge_distance)
        {
            // Get index of node associated with cell and get pointer to this node
            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            // calculate radial force on boundary
            c_vector<double, DIM> old_node_location;
            old_node_location = rOldLocations.find(p_node)->second;
            double radial_discplacement = radius - norm_2(old_node_location-mCentreOfSphere);

            // Add to total force (to calculate average)
            total_displacement = total_displacement + radial_discplacement;
            num_boundary_cells++;
        }
    }
    
    if (num_boundary_cells > 0)
    {
        mRadiusOfSphere = mRadiusOfSphere + total_displacement/(double)num_boundary_cells;
        //PRINT_VARIABLE(mRadiusOfSphere);
    }

    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Find the radial distance between this cell and the surface of the sphere
        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        double radius = norm_2(cell_location - mCentreOfSphere);
        
        // If the cell is ouside the surface of the sphere...
        if (radius > mRadiusOfSphere)
        {
            double max_jiggle = 0.001;

            // ...move the cell back onto the surface of the sphere
            c_vector<double, DIM> location_within_sphere = mCentreOfSphere + (1-max_jiggle*RandomNumberGenerator::Instance()->ranf())*mRadiusOfSphere*(cell_location - mCentreOfSphere)/radius;

            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            p_node->rGetModifiableLocation() = location_within_sphere;
        }
    }
}

template<unsigned DIM>
bool MovingDiscGeometryBoundaryCondition<DIM>::VerifyBoundaryCondition()
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
void MovingDiscGeometryBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
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
template class MovingDiscGeometryBoundaryCondition<1>;
template class MovingDiscGeometryBoundaryCondition<2>;
template class MovingDiscGeometryBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MovingDiscGeometryBoundaryCondition)
