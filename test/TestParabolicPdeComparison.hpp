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

#ifndef TESTPARABOLICPDECOMPARISION
#define TESTPARABOLICPDECOMPARISION

#include <cxxtest/TestSuite.h>

// This macro prevents errors with GCC 4.8 of form "unable to find numeric literal operator 'operator"" Q'"
// when compiling with -std=gnu++11 (see #2929). \todo: remove when GCC 4.8 is no longer supported.
#define BOOST_MATH_DISABLE_FLOAT128
#include <boost/math/special_functions/bessel.hpp>

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellIdWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "ParabolicGrowingDomainPdeModifier.hpp"

#include "UniformSourceParabolicPde.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "VolumeTrackingModifier.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"
#include "CellLabel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "RadialGrowthForce.hpp"
#include "GrowingDiscGeometryBoundaryCondition.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "Debug.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

static const double M_DT = 0.01;
static const double M_SAMPLING_MULTIPLE = 100;
static const double M_TIME_FOR_SIMULATION = 10;

static const double M_TISSUE_RADIUS = 5; 
static const double M_APOPTOTIC_RADIUS = 2;
static const double M_BOX_HALF_WIDTH = 6;
static const double M_GROWING_BOX_HALF_WIDTH = 11; 
static const double M_BOX_H = 0.25;
static const double M_CONSTANT_UPTAKE = 0.0; 
static const double M_LINEAR_UPTAKE = -0.1; 
static const double M_DIFFUSION_COEFICIENT = 0.5; 
static const double M_DUDT_COEFICIENT = 1.0;

static const double M_BOUNDARY_CONDITION = 1;


class TestParabolicPdeComparison : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(MutableMesh<2,2>* pMesh, std::vector<CellPtr>& rCells, double EquilibriumVolume)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        for (typename TetrahedralMesh<2,2>::NodeIterator node_iter = pMesh->GetNodeIteratorBegin();
                node_iter != pMesh->GetNodeIteratorEnd();
                ++node_iter)
        {
            double initial_condition_at_node = 1.0;
        
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-3.5);
            p_cycle_model->SetQuiescentVolumeFraction(EquilibriumVolume);
            p_cycle_model->SetEquilibriumVolume(0.5*sqrt(3.0));

            p_cycle_model->SetSDuration(0.0);
            p_cycle_model->SetTransitCellG1Duration(3.0);
            p_cycle_model->SetG2Duration(0.0);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("morphogen",initial_condition_at_node); //Specifying initial condition   

            rCells.push_back(p_cell);
        }

    }

    void MakeApoptoicRegion(AbstractCellPopulation<2>& rCellPopulation)
    {
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);
        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            c_vector<double,2>  cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double radius = norm_2(cell_location);
            
            if (radius<M_APOPTOTIC_RADIUS)
            {
                cell_iter->AddCellProperty(p_apoptotic_state);
            }
        }
    }

    boost::shared_ptr<AbstractLinearPde<2,2> > MakePde(AbstractCellPopulation<2>& rCellPopulation, std::string pde_type)
    {
        boost::shared_ptr<AbstractLinearPde<2,2> > p_pde = boost::shared_ptr<AbstractLinearPde<2,2> >();

        if (pde_type.compare("UniformPde")==0)
        {
            ASSIGN_PTR(p_pde, UniformSourceParabolicPde<2>, (M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT));
        }
        else if (pde_type.compare("CellwisePde")==0)
        {
            bool is_volume_scaled = false;
            ASSIGN_PTR(p_pde, CellwiseSourceParabolicPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT, is_volume_scaled));
        }
        else if (pde_type.compare("VolumeScaledCellwisePde")==0)
        {
            bool is_volume_scaled = true;
            ASSIGN_PTR(p_pde, CellwiseSourceParabolicPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT, is_volume_scaled));
        }
        else if (pde_type.compare("AveragedPde")==0)
        {
            bool is_volume_scaled = false;
            ASSIGN_PTR(p_pde, AveragedSourceParabolicPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT, is_volume_scaled));       
        }
        else if (pde_type.compare("VolumeScaledAveragedPde")==0)
        {
            bool is_volume_scaled = true;
            ASSIGN_PTR(p_pde, AveragedSourceParabolicPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT, is_volume_scaled));
        }
        else
        {
            NEVER_REACHED;
        }

        return p_pde;
    }


// Idea: make up a disc mesh using a cell based siulation 

public:
    /* 
     * First test the solutions are the same (as expected) for Uniform uptake (i.e only constant and linear term)
     * on a static tissue
     *  
     * With an apoptotic region 
     * 
     * du/dt = D Grad.(Grad u) + au + b; 
     * 
     */
    void TestParabolicPdes()
    {
        std::string base_type = "Parabolic/";
        
        std::string tissue_types[3] = {"StaticDisc", "GrowingDisc", "ProliferatingDisc"};

        std::string domain_types[3] = {"GrowingDomain", "BoxDomain", "BoxDomainAdvection"};

        std::string pde_types[3][3] = {{"UniformPde", "CellwisePde", "VolumeScaledCellwisePde"},
                                       {"UniformPde", "AveragedPde", "VolumeScaledAveragedPde"},
                                       {"UniformPde", "AveragedPde", "VolumeScaledAveragedPde"}};
        
        for (unsigned tissue_type_index = 0; tissue_type_index != 3; tissue_type_index++)
        {
            std::string tissue_type = tissue_types[tissue_type_index];

            for (unsigned domain_type_index = 0; domain_type_index != 3; domain_type_index++)
            {
                std::string domain_type = domain_types[domain_type_index];

                for (unsigned pde_type_index = 0; pde_type_index != 3; pde_type_index++)
                {
                    std::string pde_type = pde_types[domain_type_index][pde_type_index];

                    std::string output_dir = base_type +tissue_type + "/" + domain_type + "/" + pde_type;
        
                    PRINT_VARIABLE(output_dir);
                    
                    //TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
                    TrianglesMeshReader<2,2> mesh_reader("projects/MultiCellularPDEComparison/test/data/CircularMeshR5");
                    MutableMesh<2,2> mesh;
                    mesh.ConstructFromMeshReader(mesh_reader);
                    //mesh.Scale(M_TISSUE_RADIUS,M_TISSUE_RADIUS);
                    
                    double equilibrium_volume = 100.0;
                    if (tissue_type.compare("ProliferatingDisc")==0)
                    {
                        TRACE("ProliferatingDisc1");
                        equilibrium_volume = 0.8;
                    }
                    std::vector<CellPtr> cells;
                    GenerateCells(&mesh,cells,equilibrium_volume);
                    MeshBasedCellPopulation<2> cell_population(mesh, cells);

                    cell_population.AddCellWriter<CellIdWriter>();
                    cell_population.AddCellWriter<CellMutationStatesWriter>();
                    cell_population.SetWriteVtkAsPoints(true);
                    cell_population.AddPopulationWriter<VoronoiDataWriter>();

                    // bound poulation so finite areas for cell scaling
                    cell_population.SetBoundVoronoiTessellation(true);

                    if (tissue_type.compare("StaticDisc")==0)
                    {
                        TRACE("StaticDisc");
                       // MakeApoptoicRegion(cell_population);
                    }

                    OffLatticeSimulation<2> simulator(cell_population);
                    simulator.SetOutputDirectory(output_dir);
                    simulator.SetDt(M_DT);
                    simulator.SetSamplingTimestepMultiple(M_SAMPLING_MULTIPLE);
                    simulator.SetEndTime(M_TIME_FOR_SIMULATION);
                    
// if (tissue_type.compare("StaticDisc")!=0)
// {
//     TRACE("Not StaticDisc");
//     if (domain_type.compare("GrowingDomain")==0)
//     {
//         TRACE("GrowingDomain");
//         if (pde_type.compare("UniformPde")==0)
//         {
//             TRACE("UniformPde");
            simulator.SetOutputCellVelocities(true);
//         }
//     }
// }

                    MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
                    simulator.AddSimulationModifier(p_modifier);

                    if (tissue_type.compare("GrowingDisc")==0)
                    {
                        TRACE("GrowingDisc");
                        // Use a force to make the domain grow
                        double growth_rate = 0.5;
                        MAKE_PTR_ARGS(RadialGrowthForce<2>, p_growth_force,(growth_rate));
                        simulator.AddForce(p_growth_force);
                    }

                    if (tissue_type.compare("ProliferatingDisc")==0)
                    {
                        TRACE("ProliferatingDisc2");
                        // Create a force law and pass it to the simulation
                        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
                        p_force->SetCutOffLength(1.5);
                        simulator.AddForce(p_force);
                        
                        // Use a force to move the cells with the boundary
                        double growth_rate = 0.5;
                        MAKE_PTR_ARGS(RadialGrowthForce<2>, p_growth_force,(growth_rate));
                        simulator.AddForce(p_growth_force);

                        //Add disc boundary
                        double radius = 5;
                        c_vector<double,3> centre = zero_vector<double>(2);
                        MAKE_PTR_ARGS(GrowingDiscGeometryBoundaryCondition<2>, p_boundary_condition, (&cell_population, centre, radius, growth_rate));
                        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);
                    }

                    // Create PDE 
                    boost::shared_ptr<AbstractLinearPde<2,2> > p_pde = MakePde(cell_population, pde_type);

                    // create boundary condition
                    MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));
    
                    // Create a PDE modifier and set the name of the dependent variable in the PDE
                    if (domain_type.compare("GrowingDomain")==0) //Growing Domain
                    {
                        TRACE("GrowingDomain");
                        bool is_neuman_bcs = false;
                        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs));
                        p_pde_modifier->SetDependentVariableName("morphogen");
                        simulator.AddSimulationModifier(p_pde_modifier);
                    }
                    else // Box domain
                    {
                        TRACE("BoxDomain");
                        assert(domain_type.compare("BoxDomain")==0 || domain_type.compare("BoxDomainAdvection")==0);
                        
                        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
                        double box_half_width = M_BOX_HALF_WIDTH;
                        if (tissue_type.compare("StaticDisc")!=0)
                        {
                            TRACE("Not StaticDisc");
                            box_half_width = M_GROWING_BOX_HALF_WIDTH;
                        }
                        ChastePoint<2> lower(-box_half_width, -box_half_width);
                        ChastePoint<2> upper(box_half_width, box_half_width);
                        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

                        // Create a PDE modifier and set the name of the dependent variable in the PDE
                        bool is_neuman_bcs = false;
                        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs, p_cuboid, M_BOX_H));
                        p_pde_modifier->SetDependentVariableName("morphogen");

                        if (domain_type.compare("BoxDomainAdvection")==0)
                        {
                            TRACE("BoxDomainAdvection - SetSolutionMovingWithCells")
                            // Uses interpolation between cells 
                            p_pde_modifier->SetSolutionMovingWithCells(true);
                        }   

                        // Set the BSC on the elements that don't contain Cells.
                        p_pde_modifier->SetBcsOnBoxBoundary(false);
                        p_pde_modifier->SetBcsOnBoundingSphere(true);

                        simulator.AddSimulationModifier(p_pde_modifier);
                    }

                    // Add data writer to output oxygen to a file for simple comparison
                    boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
                    cell_population.AddCellWriter(p_cell_data_item_writer);

                    simulator.Solve();

                    // Test some simulation statistics
                    // TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 312u); // No birth yet

                    // Test against analytic solution .... From Bessels functions  
                    // for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
                    //         cell_iter != simulator.rGetCellPopulation().End();
                    //         ++cell_iter)
                    // {
                    //     c_vector<double,2>  cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);
                    //     double radius = norm_2(cell_location);
                    
                    //     if (radius<0.001)
                    //     {
                    //         double morphogen = cell_iter->GetCellData()->GetItem("morphogen");

                    //         if (pde_type.compare("VolumeScaledCellwisePde") || pde_type.compare("AveragedPde"))
                    //         {
                    //             // these are the averaged ones with volume not scaled (and cellwise scaled to match)
                    //             TS_ASSERT_DELTA(morphogen, 0.253, 1e-2);
                    //         }
                    //         else
                    //         {
                    //             TS_ASSERT_DELTA(morphogen, 0.5781, 1e-2);
                    //         }
                    //     }
                    // }

                    // Reset for next pde
                    SimulationTime::Instance()->Destroy();
                    SimulationTime::Instance()->SetStartTime(0.0);
                }
            }
        }
    }
    
    // /* 
    //  * Now test the solutions are the same (as expected) for Uniform uptake (i.e only constant and linear term)
    //  * on a growing tissue
    //  *  
    //  * Maybe With an apoptotic region? 
    //  * 
    //  * du/dt = D Grad.(Grad u) + au + b; 
    //  * 
    //  * Tissue growing at drdt = alpha r (i.e growing radially)
    //  * 
    //  */
    // void noTestParabolicPdesOnGrowingDomain()
    // {
    //     std::string base_type = "Parabolic/GrowingDiscTest/";

    //     std::string domain_types[3] = {"GrowingDomain", "BoxDomain"};

    //     std::string pde_types[2][3] = {{"UniformPde", "CellwisePde", "VolumeScaledCellwisePde"},
    //                                    {"UniformPde", "AveragedPde", "VolumeScaledAveragedPde"}};
        
    //     for (unsigned domain_type_index = 0; domain_type_index != 2; domain_type_index++)
    //     {
    //         std::string domain_type = domain_types[domain_type_index];

    //         for (unsigned pde_type_index = 0; pde_type_index != 1; pde_type_index++)
    //         {
    //             std::string pde_type = pde_types[domain_type_index][pde_type_index];

    //             std::string output_dir = base_type + domain_type + "/" + pde_type;
       
    //             PRINT_VARIABLE(output_dir);

    //             //TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
    //             TrianglesMeshReader<2,2> mesh_reader("projects/MultiCellularPDEComparison/test/data/CircularMeshR5");
    //             MutableMesh<2,2> mesh;
    //             mesh.ConstructFromMeshReader(mesh_reader);
    //             //mesh.Scale(M_TISSUE_RADIUS,M_TISSUE_RADIUS);
    //             std::vector<CellPtr> cells;
    //             GenerateCells(mesh,cells,1.0);
    //             MeshBasedCellPopulation<2> cell_population(mesh, cells);

    //             cell_population.AddCellWriter<CellIdWriter>();
    //             cell_population.AddCellWriter<CellMutationStatesWriter>();
    //             cell_population.SetWriteVtkAsPoints(true);
    //             cell_population.AddPopulationWriter<VoronoiDataWriter>();

    //             // bound poulation so finite areas for cell scaling
    //             cell_population.SetBoundVoronoiTessellation(true);

    //             OffLatticeSimulation<2> simulator(cell_population);
    //             simulator.SetOutputDirectory(output_dir);
    //             simulator.SetDt(M_DT);
    //             simulator.SetSamplingTimestepMultiple(M_SAMPLING_MULTIPLE);
    //             simulator.SetEndTime(M_TIME_FOR_SIMULATION);

                    
    //             MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
    //             simulator.AddSimulationModifier(p_modifier);

    //             // Use a force to make the domain grow
    //             MAKE_PTR(RadialGrowthForce<2>, p_growth_force);
    //             simulator.AddForce(p_growth_force);

    //             // Create PDE 
    //             boost::shared_ptr<AbstractLinearPde<2,2> > p_pde = MakePde(cell_population, pde_type);

    //             // create boundary condition
    //             MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));
    
    //             // Create a PDE modifier and set the name of the dependent variable in the PDE
    //             if (domain_type.compare("GrowingDomain")==0) //Growing Domain
    //             {
    //                 bool is_neuman_bcs = false;
    //                 MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs));
    //                 p_pde_modifier->SetDependentVariableName("morphogen");
    //                 simulator.AddSimulationModifier(p_pde_modifier);
    //             }
    //             else // Box domain
    //             {
    //                 assert(domain_type.compare("BoxDomain")==0);
                
    //                 // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
    //                 ChastePoint<2> lower(-10, -10);
    //                 ChastePoint<2> upper(10, 10);
    //                 MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

    //                 // Create a PDE modifier and set the name of the dependent variable in the PDE
    //                 bool is_neuman_bcs = false;
    //                 MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs, p_cuboid, M_BOX_H));
                    
    //                 //p_pde_modifier->SetSolutionMovingWithCells(true);
    //                 p_pde_modifier->SetDependentVariableName("morphogen");

    //                 // Set the BSC on the elements that don't contain Cells.
    //                 p_pde_modifier->SetBcsOnBoxBoundary(false);
    //                 p_pde_modifier->SetBcsOnBoundingSphere(true);

    //                 simulator.AddSimulationModifier(p_pde_modifier);
    //             }

    //             // Add data writer to output oxygen to a file for simple comparison
    //             boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
    //             cell_population.AddCellWriter(p_cell_data_item_writer);

    //             simulator.Solve();

    //             // Reset for next pde
    //             SimulationTime::Instance()->Destroy();
    //             SimulationTime::Instance()->SetStartTime(0.0);
    //         }
    //     }
    // }

};

#endif /*TESTPARABOLICPDECOMPARISION*/
