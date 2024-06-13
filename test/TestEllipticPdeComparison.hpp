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
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"

#include "UniformSourceEllipticPde.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "VolumeTrackingModifier.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
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

static const double M_NUM_CELLS_ACROSS = 5; // this is 5^2 initial cells

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

public:

    // Tissue based
    void TestEllipticGrowingDomainPdeModifierWithNodeBasedMonolayer()
    {
        // Circular Mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructFromMeshReader(mesh_reader);
        temp_mesh.Scale(5.0,5.0);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        std::vector<CellPtr> cells;
        GenerateCells(mesh.GetNumNodes(),cells,1.0);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticGrowingMonolayers/Node");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1.0);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        //MAKE_PTR(RepulsionForce<2>, p_force);
        //simulator.AddForce(p_force);

        // Create PDE and boundary condition objects
        //MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -1));
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 912u); // No birth yet

        // Test nothing's changed
        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 1.5, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], sqrt(3.0)/2.0, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(4))->GetCellData()->GetItem("oxygen"), 0.9717, 1e-4);
    }

    void noTestEllipticGrowingDomainPdeModifierWithMeshBasedMonolayer()
    {
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS,M_NUM_CELLS_ACROSS,0);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells,1.0);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticGrowingMonolayers/MeshPoint");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 912u); // No birth yet

        // Test PDE solution
        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 1.5, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], sqrt(3.0)/2.0, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(4))->GetCellData()->GetItem("oxygen"), 0.9717, 1e-4);
    }

    // Box based

    void TestEllipticBoxDomainPdeModifierWithNodeBasedMonolayer()
    {
        // Circular Mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructFromMeshReader(mesh_reader);
        temp_mesh.Scale(5,5);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        std::vector<CellPtr> cells;
        GenerateCells(mesh.GetNumNodes(),cells,1.0);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticBoxMonolayers/Node");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        //MAKE_PTR(RepulsionForce<2>, p_force);
        //simulator.AddForce(p_force);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-1.0));
        //MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -1.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-10.0, -10.0);
        ChastePoint<2> upper(10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 2.0));
        p_pde_modifier->SetDependentVariableName("oxygen");
        
        // Set the BSC on the elements that don't contain Cells.
        p_pde_modifier->SetBcsOnBoxBoundary(false);

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 25u); // No birth yet

        // Test new locations and concentrations (i.e that the results haven't changed since last revision.)
        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 4.0, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.0, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(4))->GetCellData()->GetItem("oxygen"), 0.6786, 1e-4);
    }

    void noTestEllipticBoxDomainPdeModifierWithMeshBasedMonolayer()
    {
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS,M_NUM_CELLS_ACROSS,0);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells,1.0);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EllipticBoxMonolayers/MeshPoint");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test some simulation statistics
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumAllCells(), 25u); // No birth yet

        // Test new locations and concentrations (i.e that the results haven't changed since last revision.)
        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], 4.0, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], 0.0, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(4))->GetCellData()->GetItem("oxygen"), 0.6786, 1e-4);
    }

};

#endif /*TESTELLIPTICPDECOMPARISION*/
