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

#ifndef TESTELLIPTICPDECOMPARISION
#define TESTELLIPTICPDECOMPARISION

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellIdWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"

#include "UniformSourceEllipticPde.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "AveragedSourceEllipticPde.hpp"
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

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

static const double M_TIME_FOR_SIMULATION = 1.0;


static const double M_TISSUE_RADIUS = 5; 
static const double M_APOPTOTIC_RADIUS = 2;
static const double M_BOX_HALF_WIDTH = 6; 
static const double M_BOX_H = 0.1;
static const double M_CONSTANT_UPTAKE = -0.1; 
static const double M_LINEAR_UPTAKE = -0.1; 
static const double M_DIFFUSION_COEFICIENT = 0.5; 


static const double M_BOUNDARY_CONDITION = 1;

 


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

            rCells.push_back(p_cell);
        }

    }

    void MakeApoptoicRegion(AbstractCellPopulation<2>& rCellPopulation)
    {
        // MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);
        // for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
        //         cell_iter != rCellPopulation.End();
        //         ++cell_iter)
        // {
        //     c_vector<double,2>  cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        //     double radius = norm_2(cell_location);
            
        //     if (radius<M_APOPTOTIC_RADIUS)
        //     {
        //         cell_iter->AddCellProperty(p_apoptotic_state);
        //     }
        // }
    }

public:
    /* 
     * First test the solutions are the same (and match the analytic solution) for Uniform uptake (i.e only constant linear term)
     *
     * Grad.(Grad u) + k = 0;
     * 
     */
    void TestEllipticGrowingDomainUniformPde()
    {
        // // Circular Mesh
        // TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        // TetrahedralMesh<2,2> temp_mesh;
        // temp_mesh.ConstructFromMeshReader(mesh_reader);
        // temp_mesh.Scale(M_TISSUE_RADIUS,M_TISSUE_RADIUS);

        // NodesOnlyMesh<2> mesh;
        // mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // std::vector<CellPtr> cells;
        // GenerateCells(mesh.GetNumNodes(),cells,1.0);

        // NodeBasedCellPopulation<2> cell_population(mesh, cells);
        // cell_population.AddCellWriter<CellIdWriter>();
        // cell_population.AddCellWriter<CellMutationStatesWriter>();

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

        // Make apoptotic region
        MakeApoptoicRegion(cell_population);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Elliptic/CircleConstantUptake/GrowingDomain/UniformPde");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1.0);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add data writer to output oxygen to a file for simple comparison
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("oxygen"));
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
            //double exact_solution = M_BOUNDARY_CONDITION + M_CELL_UPTAKE * 0.25 * (M_TISSUE_RADIUS * M_TISSUE_RADIUS - radius * radius);

            if (radius<0.001)
            {
                double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
                TS_ASSERT_DELTA(oxygen, 0.036710892271287, 1e-2);
            }
        }
    }

    void TestEllipticBoxDomainUniformPde()
    {
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

        // Make apoptotic region
        MakeApoptoicRegion(cell_population);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Elliptic/CircleConstantUptake/BoxDomain/UniformPde");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-M_BOX_HALF_WIDTH, -M_BOX_HALF_WIDTH);
        ChastePoint<2> upper(M_BOX_HALF_WIDTH, M_BOX_HALF_WIDTH);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, M_BOX_H));
        p_pde_modifier->SetDependentVariableName("oxygen");
        
        // Set the BSC on the elements that don't contain Cells.
        p_pde_modifier->SetBcsOnBoxBoundary(false);
        p_pde_modifier->SetBcsOnBoundingSphere(true);

        simulator.AddSimulationModifier(p_pde_modifier);

        // Add data writer to output oxygen to a file for simple comparison
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("oxygen"));
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
            //double exact_solution = M_BOUNDARY_CONDITION + M_CELL_UPTAKE * 0.25 * (M_TISSUE_RADIUS * M_TISSUE_RADIUS - radius * radius);

            if (radius<0.001)
            {
                double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
                TS_ASSERT_DELTA(oxygen, 0.036710892271287, 1e-2);
            }
        }
    }

    /* 
     * First test the solutions are the same (and match the analytic solution) for Cellwise uptake (i.e density dependent linear term)
     */

    void TestEllipticGrowingDomainCellwisePde()
    {
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

        // Make apoptotic region
        MakeApoptoicRegion(cell_population);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Elliptic/CircleConstantUptake/GrowingDomain/CellwisePde");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1.0);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Create PDE and boundary condition objects
        bool is_volume_scaled = false;
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, is_volume_scaled));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add data writer to output oxygen to a file for simple comparison
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("oxygen"));
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
            //double exact_solution = M_BOUNDARY_CONDITION + M_CELL_UPTAKE * 0.25 * (M_TISSUE_RADIUS * M_TISSUE_RADIUS - radius * radius);

            if (radius<0.001)
            {
                double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
                TS_ASSERT_DELTA(oxygen, 0.036710892271287, 1e-2);
            }
        }
    }

    void TestEllipticGrowingDomainVolumeDependentCellwisePde()
    {
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

        // Make apoptotic region
        MakeApoptoicRegion(cell_population);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Elliptic/CircleConstantUptake/GrowingDomain/VolumeScaledCellwisePde");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1.0);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Create PDE and boundary condition objects
        bool is_volume_scaled = true;
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, is_volume_scaled));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add data writer to output oxygen to a file for simple comparison
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("oxygen"));
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
            //double exact_solution = M_BOUNDARY_CONDITION + M_CELL_UPTAKE * 0.25 * (M_TISSUE_RADIUS * M_TISSUE_RADIUS - radius * radius);

            if (radius<0.001)
            {
                double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
                TS_ASSERT_DELTA(oxygen, 0.036710892271287, 1e-2);
            }
        }
    }


    void TestEllipticBoxDomainAveragedPde()
    {
        // Circular Mesh
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

        // Make apoptotic region
        MakeApoptoicRegion(cell_population);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Elliptic/CircleConstantUptake/BoxDomain/AveragedPde");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Create PDE and boundary condition objects
        bool is_volume_scaled = false;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, is_volume_scaled));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-M_BOX_HALF_WIDTH, -M_BOX_HALF_WIDTH);
        ChastePoint<2> upper(M_BOX_HALF_WIDTH, M_BOX_HALF_WIDTH);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, M_BOX_H));
        p_pde_modifier->SetDependentVariableName("oxygen");
        
        // Set the BSC on the elements that don't contain Cells.
        p_pde_modifier->SetBcsOnBoxBoundary(false);
        p_pde_modifier->SetBcsOnBoundingSphere(true);

        simulator.AddSimulationModifier(p_pde_modifier);

        // Add data writer to output oxygen to a file for simple comparison
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("oxygen"));
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
            //double exact_solution = M_BOUNDARY_CONDITION + M_CELL_UPTAKE * 0.25 * (M_TISSUE_RADIUS * M_TISSUE_RADIUS - radius * radius);

            if (radius<0.001)
            {
                double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
                TS_ASSERT_DELTA(oxygen, 0.036710892271287, 1e-2);
            }
        }
    }

void TestEllipticBoxDomainVolumeDependentAveragedPde()
    {
        // Circular Mesh
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

        // Make apoptotic region
        MakeApoptoicRegion(cell_population);

        // bound poulation so finite areas for cell scaling
        cell_population.SetBoundVoronoiTessellation(true);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Elliptic/CircleConstantUptake/BoxDomain/VolumeScaledAveragedPde");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Create PDE and boundary condition objects
        bool is_volume_scaled = false;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, is_volume_scaled));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-M_BOX_HALF_WIDTH, -M_BOX_HALF_WIDTH);
        ChastePoint<2> upper(M_BOX_HALF_WIDTH, M_BOX_HALF_WIDTH);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, M_BOX_H));
        p_pde_modifier->SetDependentVariableName("oxygen");
        
        // Set the BSC on the elements that don't contain Cells.
        p_pde_modifier->SetBcsOnBoxBoundary(false);
        p_pde_modifier->SetBcsOnBoundingSphere(true);

        simulator.AddSimulationModifier(p_pde_modifier);

        // Add data writer to output oxygen to a file for simple comparison
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("oxygen"));
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
            //double exact_solution = M_BOUNDARY_CONDITION + M_CELL_UPTAKE * 0.25 * (M_TISSUE_RADIUS * M_TISSUE_RADIUS - radius * radius);

            if (radius<0.001)
            {
                double oxygen = cell_iter->GetCellData()->GetItem("oxygen");
                TS_ASSERT_DELTA(oxygen, 0.036710892271287, 1e-2);
            }
        }
    }



};

#endif /*TESTELLIPTICPDECOMPARISION*/
