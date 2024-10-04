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

#ifndef TESTELLIPTICPDEMESHCONVERGENCE
#define TESTELLIPTICPDEMESHCONVERGENCE

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
#include "Debug.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

static const double M_TIME_FOR_SIMULATION = 1.0;


static const double M_TISSUE_RADIUS = 5; 
static const double M_APOPTOTIC_RADIUS = 2;
static const double M_BOX_HALF_WIDTH = 6; 
static const double M_BOX_H = 0.1;
static const double M_CONSTANT_UPTAKE = 0.0; 
static const double M_LINEAR_UPTAKE = -0.1; 
static const double M_DIFFUSION_COEFICIENT = 1.0; 


static const double M_BOUNDARY_CONDITION = 1;

 


class TestEllipticPdeMeshConvergence : public AbstractCellBasedWithTimingsTestSuite
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
            ASSIGN_PTR(p_pde, UniformSourceEllipticPde<2>, (M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT));
        }
        else if (pde_type.compare("CellwisePde")==0)
        {
            bool is_volume_scaled = false;
            ASSIGN_PTR(p_pde, CellwiseSourceEllipticPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, is_volume_scaled));
        }
        else if (pde_type.compare("VolumeScaledCellwisePde")==0)
        {
            bool is_volume_scaled = true;
            ASSIGN_PTR(p_pde, CellwiseSourceEllipticPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, is_volume_scaled));
        }
        else if (pde_type.compare("AveragedPde")==0)
        {
            bool is_volume_scaled = false;
            ASSIGN_PTR(p_pde, AveragedSourceEllipticPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, is_volume_scaled));       
        }
        else if (pde_type.compare("VolumeScaledAveragedPde")==0)
        {
            bool is_volume_scaled = true;
            ASSIGN_PTR(p_pde, AveragedSourceEllipticPde<2>, (rCellPopulation, M_CONSTANT_UPTAKE, M_LINEAR_UPTAKE, M_DIFFUSION_COEFICIENT, is_volume_scaled));
        }
        else 
        {
            NEVER_REACHED;
        }

        return p_pde;
    }

public:
    /* 
     * First test the solutions are the same (and match the analytic solution) for Uniform uptake (i.e only constant linear term)
     *
     * D Grad.(Grad u) + au + b = 0; 
     * D=1.0, a=-0.1, b=0.0
     * 
     */
    void TestEllipticPdes()
    {
        std::string base_type = "Elliptic";
        
        std::string tissue_types[2] = {"StaticDisc","StaticDiscApoptotic"};

        std::string domain_types[2] = {"GrowingDomain", "BoxDomain"};

        std::string pde_types[2][3] = {{"UniformPde", "CellwisePde", "VolumeScaledCellwisePde"},
                                       {"UniformPde", "AveragedPde", "VolumeScaledAveragedPde"}};
        
        for (unsigned tissue_type_index = 0; tissue_type_index != 2; tissue_type_index++)
        {
            std::string tissue_type = tissue_types[tissue_type_index];
        
            for (unsigned domain_type_index = 0; domain_type_index != 2; domain_type_index++)
            {
                std::string domain_type = domain_types[domain_type_index];

                for (unsigned pde_type_index = 0; pde_type_index != 3; pde_type_index++)
                {
                    std::string pde_type = pde_types[domain_type_index][pde_type_index];

                    std::string output_dir = base_type + "/" + tissue_type + "/" + domain_type + "/" + pde_type;
            
                    PRINT_VARIABLE(output_dir);
            
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


                    if (tissue_type.compare("StaticDiscApoptotic")==0)
                    {
                        TRACE("StaticDiscApoptotic");
                        // Add Apoptotic region 
                        MakeApoptoicRegion(cell_population);
                    }

                    OffLatticeSimulation<2> simulator(cell_population);
                    simulator.SetOutputDirectory(output_dir);
                    simulator.SetDt(1.0);
                    simulator.SetSamplingTimestepMultiple(1.0);
                    simulator.SetEndTime(M_TIME_FOR_SIMULATION);

                    // Create PDE 
                    boost::shared_ptr<AbstractLinearPde<2,2> > p_pde = MakePde(cell_population, pde_type);

                    // create boundary condition
                    MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (M_BOUNDARY_CONDITION));


                    // Create a PDE modifier and set the name of the dependent variable in the PDE
                    if (domain_type.compare("GrowingDomain")==0) //Growing Domain
                    {
                        TRACE("GrowingDomain");
                        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
                        p_pde_modifier->SetDependentVariableName("oxygen");
                        simulator.AddSimulationModifier(p_pde_modifier);
                    }
                    else // Box domain
                    {
                        TRACE("BoxDomain");
                        assert(domain_type.compare("BoxDomain")==0);
                            
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
                    }

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
                    
                        if (radius<0.001)
                        {
                            double oxygen = cell_iter->GetCellData()->GetItem("oxygen");

                            if (tissue_type.compare("StaticDisc")==0)
                            {
                                if ((pde_type.compare("AveragedPde")==0) || (pde_type.compare("VolumeScaledCellwisePde")==0))
                                {
                                    // these are the averaged ones with volume not scaled (and cellwise scaled to match)
                                    TS_ASSERT_DELTA(oxygen, 0.253, 1e-2);
                                }
                                else
                                {
                                    TS_ASSERT_DELTA(oxygen, 0.5781, 1e-2);
                                }
                            }
                            else if (tissue_type.compare("StaticDiscApoptotic")==0)
                            { 
                                if ((pde_type.compare("AveragedPde")==0) || (pde_type.compare("VolumeScaledCellwisePde")==0))
                                {
                                    // these are the averaged ones with volume not scaled (and cellwise scaled to match)
                                    TS_ASSERT_DELTA(oxygen, 0.43, 1e-2);
                                }
                                else if ((pde_type.compare("VolumeScaledAveragedPde")==0) || (pde_type.compare("CellwisePde")==0))
                                
                                {
                                    TS_ASSERT_DELTA(oxygen, 0.734, 1e-2);
                                }
                                else
                                {              
                                    //These are the uniform ones so no difference from non apoptotic case
                                    TS_ASSERT_DELTA(oxygen, 0.5781, 1e-2);
                                }
                            }
                            else
                            {
                                NEVER_REACHED;
                            }
                        }
                    }

                    // Reset for next pde
                    SimulationTime::Instance()->Destroy();
                    SimulationTime::Instance()->SetStartTime(0.0);
                }
            }
        }
    }
};

#endif /*TESTELLIPTICPDEMESHCONVERGENCE*/
