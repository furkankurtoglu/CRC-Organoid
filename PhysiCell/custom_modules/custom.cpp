/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here 

Cell_Definition organoid_cell;
Cell_Definition fibro_cell;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// ---- START -- Default Cell Definitions -- START ---- //
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "default cell"; 
	
	// set default cell cycle model 

	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = NULL; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int glucose_substrate_index = microenvironment.find_density_index( "glucose" );
	int lactate_substrate_index = microenvironment.find_density_index( "lactate" );
	int atp_substrate_index = microenvironment.find_density_index( "atp" );
	int glutamine_substrate_index = microenvironment.find_density_index( "glutamine" );
	

	int start_index = live.find_phase_index( PhysiCell_constants::live );
	int end_index = live.find_phase_index( PhysiCell_constants::live );

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 0.0; 
	
	cell_defaults.phenotype.secretion.uptake_rates[glucose_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[glucose_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[glucose_substrate_index] = 0.0; 
	
	cell_defaults.phenotype.secretion.uptake_rates[lactate_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[lactate_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[lactate_substrate_index] = 0.0; 
	
	cell_defaults.phenotype.secretion.uptake_rates[atp_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[atp_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[atp_substrate_index] = 0.0; 	
	
	cell_defaults.phenotype.secretion.uptake_rates[glutamine_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[glutamine_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[glutamine_substrate_index] = 0.0;
	
	cell_defaults.custom_data.add_variable( "energy", "dimensionless" , parameters.doubles("cell_default_inital_energy") ); 
	cell_defaults.custom_data.add_variable( "energy_creation_rate", "1/min" , parameters.doubles("cell_default_energy_creation_rate") ); 
	cell_defaults.custom_data.add_variable( "energy_use_rate", "1/min" , parameters.doubles("cell_default_energy_use_rate") ); 
	cell_defaults.custom_data.add_variable( "cycle_energy_threshold", "dimensionless" , parameters.doubles("cell_default_cycle_energy_threshold") ); 
	cell_defaults.custom_data.add_variable( "death_energy_threshold", "dimensionless" , parameters.doubles("cell_default_death_energy_threshold") );

	// ---- END -- Default Cell Definitions -- END ---- //
	
	
	
	
	
	// ---- START -- Organoid Cell Definitions -- START ---- //
	
	organoid_cell = cell_defaults; 
	organoid_cell.type = 1; 
	organoid_cell.name = "organoid cell"; 
	
	// make sure the new cell type has its own reference phenotyhpe
	organoid_cell.parameters.pReference_live_phenotype = &( organoid_cell.phenotype ); 
	// Set cell-cell adhesion to 5% of other cells 
	organoid_cell.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "organoid_cell_relative_adhesion" ); // 0.05; 
	
	// Set apoptosis to zero 
	organoid_cell.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "organoid_cell_apoptosis_rate" ); // 0.0; 
	organoid_cell.phenotype.cycle.data.transition_rate(start_index,end_index) *= parameters.doubles( "organoid_cell_relative_cycle_entry_rate" ); // 0.0;
	
	// Set Energy Function
	organoid_cell.functions.update_phenotype = tumor_energy_update_function;
	// ---- END -- Organoid Cell Definitions -- END ---- //
	
	
	
	// ---- START -- Fibroblast Definitions -- START ---- //
	fibro_cell = cell_defaults; 
	fibro_cell.type = 2; 
	fibro_cell.name = "fibroblast"; 
	
	// make sure the new cell type has its own reference phenotyhpe
	
	fibro_cell.parameters.pReference_live_phenotype = &( fibro_cell.phenotype ); 
	// Set cell-cell adhesion to 5% of other cells 
	fibro_cell.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "fibroblast_relative_adhesion" ); // 0.05; 
	
	// Set apoptosis to zero 
	fibro_cell.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "fibroblast_apoptosis_rate" ); // 0.0; 
	fibro_cell.phenotype.cycle.data.transition_rate(start_index,end_index) *= parameters.doubles( "fibroblast_relative_cycle_entry_rate" ); // 0.0;
	// ---- END -- Fibroblast Cell Definitions -- END ---- //	
	
	return; 
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == true )
	{
		std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl; 
		default_microenvironment_options.simulate_2D = false; 
	}	
	Microenvironment* pME = get_default_microenvironment();
	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int glucose_substrate_index = microenvironment.find_density_index( "glucose" );
	int lactate_substrate_index = microenvironment.find_density_index( "lactate" );
	int atp_substrate_index = microenvironment.find_density_index( "atp" );
	int glutamine_substrate_index = microenvironment.find_density_index( "glutamine" );


    if( glucose_substrate_index < 0 )
    {
        pME->add_density( "glucose", "dimensionless" , 0.0 , 0.0 );
        glucose_substrate_index = pME->find_density_index( "glucose" );
        default_microenvironment_options.Dirichlet_condition_vector[glucose_substrate_index] = 0.0;
        default_microenvironment_options.Dirichlet_activation_vector[glucose_substrate_index] = false;
    }
	
    if( lactate_substrate_index < 0 )
    {
        pME->add_density( "lactate", "dimensionless" , 0.0 , 0.0 );
        lactate_substrate_index = pME->find_density_index( "lactate" );
        default_microenvironment_options.Dirichlet_condition_vector[lactate_substrate_index] = 0.0;
        default_microenvironment_options.Dirichlet_activation_vector[lactate_substrate_index] = false;
    }
	
    if( atp_substrate_index < 0 )
    {
        pME->add_density( "atp", "dimensionless" , 0.0 , 0.0 );
        atp_substrate_index = pME->find_density_index( "atp" );
        default_microenvironment_options.Dirichlet_condition_vector[atp_substrate_index] = 0.0;
        default_microenvironment_options.Dirichlet_activation_vector[atp_substrate_index] = false;
    }	
	
    if( glutamine_substrate_index < 0 )
    {
        pME->add_density( "glutamine", "dimensionless" , 0.0 , 0.0 );
        glutamine_substrate_index = pME->find_density_index( "glutamine" );
        default_microenvironment_options.Dirichlet_condition_vector[glutamine_substrate_index] = 0.0;
        default_microenvironment_options.Dirichlet_activation_vector[glutamine_substrate_index] = false;
    }	
	
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 5 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0,0,0,0,0 }; 
	
	
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pCell;

	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	double organoid_distance = parameters.doubles("organoid_distance");
	double initial_tumor_radius =  parameters.doubles("initial_tumor_radius");

	std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,initial_tumor_radius); 
	std::cout << "creating " << positions.size() << " closely-packed tumor cells ... " << std::endl; 
	

	// create organoid
	for( int i=0; i < positions.size(); i++ )
	{
		positions[i][1] += organoid_distance;
		pCell = create_cell(organoid_cell);
		pCell->assign_position( positions[i] );
	}

	return; 
}


std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = false_cell_coloring_cytometry(pCell); 
	
	// if the cell is motile and not dead, paint it black 
	
	if( pCell->phenotype.death.dead == false && 
		pCell->type == 1 )
	{
		 output[0] = "black"; 
		 output[2] = "black"; 	
	}
	
	return output; 
}


void tumor_energy_update_function( Cell* pCell, Phenotype& phenotype , double dt )
{
/* 	static int nE = pCell->custom_data.find_variable_index( "energy" ); 
	static int nA = pCell->custom_data.find_variable_index( "energy_creation_rate" ); 
	static int nB = pCell->custom_data.find_variable_index( "energy_use_rate" ); 
	
	static int nBirth = pCell->custom_data.find_variable_index( "cycle_energy_threshold" );  
	static int nDeath = pCell->custom_data.find_variable_index( "death_energy_threshold" ); 

	static int nAlpha = pCell->custom_data.find_variable_index( "alpha" ); 
	static int nBeta = pCell->custom_data.find_variable_index( "beta" ); 
	static int nGamma = pCell->custom_data.find_variable_index( "gamma" ); 
	static int nRho = pCell->custom_data.find_variable_index( "rho" );
	static int nPhi = pCell->custom_data.find_variable_index( "phi" );
	static int nChi = pCell->custom_data.find_variable_index( "chi" );
	
	
	double O2 = pCell->nearest_density_vector()[nO2]; 	
	double Glucose = pCell->nearest_density_vector()[nGlucose];  
	
	*/
	//pCell->custom_data[nE] += dt*( pCell->custom_data[nA] * pCell->custom_data[nPhi] * (O2) * Glucose * pCell->custom_data[nAlpha] + pCell->custom_data[nChi] * pCell->custom_data[nBeta] * pCell->custom_data[nA] * (Glucose*2*Glucose) - pCell->custom_data[nGamma] * pCell->custom_data[nB] - pCell->custom_data[nRho] * (O2)); 
	
	//int a=0;
	//std::cout<<pCell->custom_data[nE]<<std::endl;
	//std::cin>>a;
	
/* 	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	static int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	static int ncycle_End_i = live.find_phase_index( PhysiCell_constants::live );
	static int ncycle_Start_i = live.find_phase_index( PhysiCell_constants::live ); */
	

	// No Proliferation, Apoptosis, and Necrosis
/* 	phenotype.death.rates[necrosis_model_index] = 0.0; 
	phenotype.cycle.data.transition_rate(ncycle_Start_i,ncycle_End_i) = 0.0;
	phenotype.death.rates[apoptosis_model_index] = 0.0; */

	// die if energy is low 
/* 	if( pCell->custom_data[nE] < pCell->custom_data[nDeath] )
	{
		phenotype.death.rates[apoptosis_model_index] = parameters.doubles("apoptosis_rate"); 
	}

	if( pCell->custom_data[nE] > pCell->custom_data[nBirth])
	{
		phenotype.cycle.data.transition_rate( ncycle_Start_i,ncycle_End_i ) = parameters.doubles("proliferation_rate"); 
		
	}
		 */
		
		
	return;
}

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;
	
}