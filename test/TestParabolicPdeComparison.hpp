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
#include "Debug.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

static const double M_DT = 0.01;
static const double M_TIME_FOR_SIMULATION = 1.0;

static const double M_TISSUE_RADIUS = 5; 
static const double M_APOPTOTIC_RADIUS = 2;
static const double M_BOX_HALF_WIDTH = 6; 
static const double M_BOX_H = 0.1;
static const double M_CONSTANT_UPTAKE = 0.1; 
static const double M_LINEAR_UPTAKE = 0.1; 
static const double M_DIFFUSION_COEFICIENT = 0.5; 
static const double M_DUDT_COEFICIENT = 0.5;

static const double M_BOUNDARY_CONDITION = 0;


class TestEllipticPdeComparison : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double EquilibriumVolume)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        for (unsigned i=0; i<num_cells; i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);
            p_cell->GetCellData()->SetItem("morphogen",0.0); //Specifying initial condition     
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

    boost::shared_ptr<AbstractLinearPde<2,2> > MakePde(AbstractCellPopulation<2>& rCellPopulation, unsigned pde_type)
    {
        boost::shared_ptr<AbstractLinearPde<2,2> > p_pde = boost::shared_ptr<AbstractLinearPde<2,2> >();

        if (pde_type == 0 || pde_type == 3)
        {
            ASSIGN_PTR(p_pde, UniformSourceParabolicPde<2>, (M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT));
        }
        else if (pde_type == 1)
        {
            bool is_volume_scaled = false;
            ASSIGN_PTR(p_pde, CellwiseSourceParabolicPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT, is_volume_scaled));
        }
        else if (pde_type == 2)
        {
            bool is_volume_scaled = true;
            ASSIGN_PTR(p_pde, CellwiseSourceParabolicPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT, is_volume_scaled));
        }
        else if (pde_type == 4)
        {
            bool is_volume_scaled = false;
            ASSIGN_PTR(p_pde, AveragedSourceParabolicPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT, is_volume_scaled));       
        }
        else
        {
            assert(pde_type == 5);
            bool is_volume_scaled = true;
            ASSIGN_PTR(p_pde, AveragedSourceParabolicPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, M_DUDT_COEFICIENT, is_volume_scaled));
        }

        return p_pde;
    }

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
    void noTestParabolicPdes()
    {
        std::string output_dir[6] = {"Parabolic/StaticDisc/GrowingDomain/UniformPde",
                                     "Parabolic/StaticDisc/GrowingDomain/CellwisePde",
                                     "Parabolic/StaticDisc/GrowingDomain/VolumeScaledCellwisePde",
                                     "Parabolic/StaticDisc/BoxDomain/UniformPde",
                                     "Parabolic/StaticDisc/BoxDomain/AveragedPde",
                                     "Parabolic/StaticDisc/BoxDomain/VolumeScaledAveragedPde"};

        for (unsigned pde_type = 0; pde_type != 6; pde_type++)
        {
            PRINT_VARIABLE(output_dir[pde_type]);

            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            mesh.Scale(M_TISSUE_RADIUS,M_TISSUE_RADIUS);
            std::vector<CellPtr> cells;
            GenerateCells(mesh.GetNumNodes(),cells,1.0);
            MeshBasedCellPopulation<2> cell_population(mesh, cells);

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.AddPopulationWriter<VoronoiDataWriter>();

            // bound poulation so finite areas for cell scaling
            cell_population.SetBoundVoronoiTessellation(true);

MakeApoptoicRegion(cell_population);

            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory(output_dir[pde_type]);
            simulator.SetDt(M_DT);
            simulator.SetSamplingTimestepMultiple(10);
            simulator.SetEndTime(M_TIME_FOR_SIMULATION);

            // Create PDE 
            boost::shared_ptr<AbstractLinearPde<2,2> > p_pde = MakePde(cell_population, pde_type);

            // create boundary condition
            MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));
  
            // Create a PDE modifier and set the name of the dependent variable in the PDE
            if (pde_type <3) //Growing Domain
            {
                bool is_neuman_bcs = false;
                MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs));
                p_pde_modifier->SetDependentVariableName("morphogen");
                simulator.AddSimulationModifier(p_pde_modifier);
            }
            else // Box domain
            {
                // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
                ChastePoint<2> lower(-M_BOX_HALF_WIDTH, -M_BOX_HALF_WIDTH);
                ChastePoint<2> upper(M_BOX_HALF_WIDTH, M_BOX_HALF_WIDTH);
                MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

                // Create a PDE modifier and set the name of the dependent variable in the PDE
                bool is_neuman_bcs = false;
                MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs, p_cuboid, M_BOX_H));
                p_pde_modifier->SetDependentVariableName("morphogen");

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
            TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 312u); // No birth yet

            // Test against analytic solution .... From Bessels functions  
            for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
                    cell_iter != simulator.rGetCellPopulation().End();
                    ++cell_iter)
            {
                c_vector<double,2>  cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);
                double radius = norm_2(cell_location);
            
                if (radius<0.001)
                {
                    double morphogen = cell_iter->GetCellData()->GetItem("morphogen");

                    if (pde_type == 2 || pde_type == 4)
                    {
                        // these are the averaged ones with volume not scaled (and cellwise scaled to match)
                        TS_ASSERT_DELTA(morphogen, 0.253, 1e-2);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(morphogen, 0.5781, 1e-2);
                    }
                }
            }

            // Reset for next pde
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }
    
    /* 
     * Now test the solutions are the same (as expected) for Uniform uptake (i.e only constant and linear term)
     * on a growing tissue
     *  
     * Maybe With an apoptotic region? 
     * 
     * du/dt = D Grad.(Grad u) + au + b; 
     * 
     * Tissue growing at drdt = alpha r (i.e growing radially)
     * 
     */
    void TestParabolicPdesOnGrowingDomain()
    {
        std::string output_dir[6] = {"Parabolic/GrowingDisc/GrowingDomain/UniformPde",
                                     "Parabolic/GrowingDisc/GrowingDomain/CellwisePde",
                                     "Parabolic/GrowingDisc/GrowingDomain/VolumeScaledCellwisePde",
                                     "Parabolic/GrowingDisc/BoxDomain/UniformPde",
                                     "Parabolic/GrowingDisc/BoxDomain/AveragedPde",
                                     "Parabolic/GrowingDisc/BoxDomain/VolumeScaledAveragedPde"};

        for (unsigned pde_type = 0; pde_type != 6; pde_type++)
        {
            PRINT_VARIABLE(output_dir[pde_type]);

            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            mesh.Scale(M_TISSUE_RADIUS,M_TISSUE_RADIUS);
            std::vector<CellPtr> cells;
            GenerateCells(mesh.GetNumNodes(),cells,1.0);
            MeshBasedCellPopulation<2> cell_population(mesh, cells);

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.AddPopulationWriter<VoronoiDataWriter>();

            // bound poulation so finite areas for cell scaling
            cell_population.SetBoundVoronoiTessellation(true);

            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory(output_dir[pde_type]);
            simulator.SetDt(M_DT);
            simulator.SetSamplingTimestepMultiple(10);
            simulator.SetEndTime(M_TIME_FOR_SIMULATION);

            // Use a force to make the domain grow
            MAKE_PTR(RadialGrowthForce<2>, p_growth_force);
            simulator.AddForce(p_growth_force);


            // Create PDE 
            boost::shared_ptr<AbstractLinearPde<2,2> > p_pde = MakePde(cell_population, pde_type);

            // create boundary condition
            MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));
  
            // Create a PDE modifier and set the name of the dependent variable in the PDE
            if (pde_type <3) //Growing Domain
            {
                bool is_neuman_bcs = false;
                MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs));
                p_pde_modifier->SetDependentVariableName("morphogen");
                simulator.AddSimulationModifier(p_pde_modifier);
            }
            else // Box domain
            {
                // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
                ChastePoint<2> lower(-15, -15);
                ChastePoint<2> upper(15, 15);
                MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

                // Create a PDE modifier and set the name of the dependent variable in the PDE
                bool is_neuman_bcs = false;
                MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neuman_bcs, p_cuboid, 5*M_BOX_H));
                p_pde_modifier->SetDependentVariableName("morphogen");

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
            TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 312u); // No birth yet

            // Test against analytic solution .... From Bessels functions  
            for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
                    cell_iter != simulator.rGetCellPopulation().End();
                    ++cell_iter)
            {
                c_vector<double,2>  cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);
                double radius = norm_2(cell_location);
            
                if (radius<0.001)
                {
                    double morphogen = cell_iter->GetCellData()->GetItem("morphogen");

                    if (pde_type == 2 || pde_type == 4)
                    {
                        // these are the averaged ones with volume not scaled (and cellwise scaled to match)
                        TS_ASSERT_DELTA(morphogen, 0.253, 1e-2);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(morphogen, 0.5781, 1e-2);
                    }
                }
            }

            // Reset for next pde
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }

};

#endif /*TESTPARABOLICPDECOMPARISION*/
