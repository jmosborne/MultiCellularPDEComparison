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

#ifndef TESTGENERATEDISCTISSUE_HPP_
#define TESTGENERATEDISCTISSUE_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <cmath>

#include "CheckpointArchiveTypes.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "GrowingDiscGeometryBoundaryCondition.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "CellIdWriter.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "CellMutationStatesWriter.hpp"
#include "RadialGrowthForce.hpp"
#include "MovingDiscGeometryBoundaryCondition.hpp"

#include "ContactInhibitionCellCycleModel.hpp"
#include "BernoulliTrialCellCycleModel.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "TrianglesMeshWriter.hpp"
#include "Debug.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestGenerateDiscTissue : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void noTestGenerateDisc()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        int num_cells_depth = 12;
        int num_cells_width = 12;
        int num_ghost_width = 0;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, num_ghost_width);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        p_mesh->Translate(-0.5*num_cells_width,-0.45*num_cells_depth*0.5*sqrt(3.0));
        p_mesh->Scale(0.8,0.8);

        /* We now create a vector of cell pointers. */
        std::vector<CellPtr> cells;

        /* We then define the mutation state of the cells we are working with. We will just consider
         * wild type mutations here. */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);

        /* We now create a cell-cycle (only contact inhibited) model for these cells and loop over the
         * nodes of the mesh to create as many elements in the vector of cell pointers as there are
         * in the initial mesh. */
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-2.0*(double)i);
            p_cycle_model->SetQuiescentVolumeFraction(2.0);
            p_cycle_model->SetEquilibriumVolume(0.5*sqrt(3.0));

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();

            cells.push_back(p_cell);
        }

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // bound poulation so finite areas for cell scaling
        //cell_population.SetBoundVoronoiTessellation(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("DiscTisue");
        simulator.SetDt(0.001);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(20.0);

        /* Then, we define the modifier class, which automatically updates the volumes of the cells in `CellData` and passes it to the simulation.*/
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        //Add disc boundary
        double radius = 5;
        c_vector<double,3> centre = zero_vector<double>(2);
        MAKE_PTR_ARGS(GrowingDiscGeometryBoundaryCondition<2>, p_boundary_condition, (&cell_population, centre, radius)); // Circle radius M_RADIUS centre (0,0,0)
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        simulator.Solve();


        //Save mesh for PDEs
        TrianglesMeshWriter<2,2> mesh_writer("TestGenerateDisc", "CircularMeshR5", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void noTestConstantVelocityTest()
    {
        unsigned num_seeds = 2;

        std::string base_output_dir = "Parabolic/VelocityTest/";

        for (unsigned sim_index = 0; sim_index < num_seeds; sim_index++)
        {
            std::cout << " Run number " << sim_index << "... \n" << std::flush;

            // Reseed the random number generator
            RandomNumberGenerator::Instance()->Reseed(100*sim_index);
          
            //Create output directory
            std::stringstream out;
            out << "Run_" << sim_index;
            std::string output_dir = base_output_dir +  out.str();
            PRINT_VARIABLE(output_dir);
            
            //TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
            TrianglesMeshReader<2,2> mesh_reader("projects/MultiCellularPDEComparison/test/data/CircularMeshR5");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            //mesh.Scale(M_TISSUE_RADIUS,M_TISSUE_RADIUS);
            
            double equilibrium_volume = 0.8;
            std::vector<CellPtr> cells;

            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            for (typename TetrahedralMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                    node_iter != mesh.GetNodeIteratorEnd();
                    ++node_iter)
            {
                
                ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
                p_cycle_model->SetDimension(2);
                p_cycle_model->SetBirthTime(-3.5);
                p_cycle_model->SetQuiescentVolumeFraction(equilibrium_volume);
                p_cycle_model->SetEquilibriumVolume(0.5*sqrt(3.0));

                p_cycle_model->SetSDuration(0.0);
                p_cycle_model->SetTransitCellG1Duration(3.0);
                p_cycle_model->SetG2Duration(0.0);

                CellPtr p_cell(new Cell(p_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->InitialiseCellCycleModel();
                
                cells.push_back(p_cell);
            }

            MeshBasedCellPopulation<2> cell_population(mesh, cells);

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.AddPopulationWriter<VoronoiDataWriter>();

            // bound poulation so finite areas for cell scaling
            cell_population.SetBoundVoronoiTessellation(true);

            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory(output_dir);
            simulator.SetDt(0.01);
            simulator.SetSamplingTimestepMultiple(100);
            simulator.SetEndTime(40);
            
            simulator.SetOutputCellVelocities(true);

            MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

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

            simulator.Solve();

            // Reset for next pde
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }

    void TestExponentialVelocityTest()
    {
        unsigned num_seeds = 10;

        std::string base_output_dir = "Parabolic/ExponentialVelocityTest/";

        for (unsigned sim_index = 0; sim_index < num_seeds; sim_index++)
        {
            std::cout << " Run number " << sim_index << "... \n" << std::flush;

            // Reseed the random number generator
            RandomNumberGenerator::Instance()->Reseed(100*sim_index);
          
            //Create output directory
            std::stringstream out;
            out << "Run_" << sim_index;
            std::string output_dir = base_output_dir +  out.str();
            PRINT_VARIABLE(output_dir);
            
            //TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
            TrianglesMeshReader<2,2> mesh_reader("projects/MultiCellularPDEComparison/test/data/CircularMeshR5");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            //mesh.Scale(M_TISSUE_RADIUS,M_TISSUE_RADIUS);
            
      
            MAKE_PTR(WildTypeCellMutationState, p_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);
            
            std::vector<CellPtr> cells;
            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                BernoulliTrialCellCycleModel* p_cycle_model = new BernoulliTrialCellCycleModel();
                p_cycle_model->SetDimension(2);
                p_cycle_model->SetBirthTime(-2.0);
                p_cycle_model->SetDivisionProbability(0.1);
                p_cycle_model->SetMinimumDivisionAge(2.0);

                CellPtr p_cell(new Cell(p_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_transit_type);
                cells.push_back(p_cell);
            }

            // double equilibrium_volume = 0.01;
            // std::vector<CellPtr> cells;

            // MAKE_PTR(WildTypeCellMutationState, p_state);
            // MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            // for (typename TetrahedralMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
            //         node_iter != mesh.GetNodeIteratorEnd();
            //         ++node_iter)
            // {
                
            //     ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            //     p_cycle_model->SetDimension(2);
            //     p_cycle_model->SetBirthTime(-3.5);
            //     p_cycle_model->SetQuiescentVolumeFraction(equilibrium_volume);
            //     p_cycle_model->SetEquilibriumVolume(0.5*sqrt(3.0));

            //     p_cycle_model->SetSDuration(0.0);
            //     p_cycle_model->SetTransitCellG1Duration(10.0);
            //     p_cycle_model->SetG2Duration(0.0);

            //     CellPtr p_cell(new Cell(p_state, p_cycle_model));
            //     p_cell->SetCellProliferativeType(p_transit_type);
            //     p_cell->InitialiseCellCycleModel();
                
            //     cells.push_back(p_cell);
            // }

            MeshBasedCellPopulation<2> cell_population(mesh, cells);

            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellMutationStatesWriter>();
            cell_population.SetWriteVtkAsPoints(true);
            cell_population.AddPopulationWriter<VoronoiDataWriter>();

            // bound poulation so finite areas for cell scaling
            cell_population.SetBoundVoronoiTessellation(true);

            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory(output_dir);
            simulator.SetDt(0.01);
            simulator.SetSamplingTimestepMultiple(100);
            simulator.SetEndTime(25);
            
            simulator.SetOutputCellVelocities(true);

            MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
            simulator.AddSimulationModifier(p_modifier);

            // Create a force law and pass it to the simulation
            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
            p_force->SetCutOffLength(1.5);
            simulator.AddForce(p_force);

            //Add disc boundary
            double radius = 5;
            c_vector<double,3> centre = zero_vector<double>(2);
            MAKE_PTR_ARGS(MovingDiscGeometryBoundaryCondition<2>, p_boundary_condition, (&cell_population, centre, radius));
            simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

            simulator.Solve();

            // Reset for next pde
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }

    }
    

};

#endif /*TESTGENERATEDISCTISSUE_HPP_*/
