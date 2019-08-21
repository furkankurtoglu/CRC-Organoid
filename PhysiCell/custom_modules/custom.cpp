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

	cell_defaults.functions.cycle_model = Ki67_basic; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = NULL; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int i_Ki67_negative = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	int i_Ki67_positive = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive );
	
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int glucose_substrate_index = microenvironment.find_density_index( "glucose" );
	int lactate_substrate_index = microenvironment.find_density_index( "lactate" );
	int glutamine_substrate_index = microenvironment.find_density_index( "glutamine" );


	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0.001; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 0.0; 
	
	cell_defaults.phenotype.secretion.uptake_rates[glucose_substrate_index] = 0.001; 
	cell_defaults.phenotype.secretion.secretion_rates[glucose_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[glucose_substrate_index] = 0.0; 
	
	cell_defaults.phenotype.secretion.uptake_rates[lactate_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.secretion_rates[lactate_substrate_index] = 0.01; 
	cell_defaults.phenotype.secretion.saturation_densities[lactate_substrate_index] = 1.0; 
	
	cell_defaults.phenotype.secretion.uptake_rates[glutamine_substrate_index] = 0.001;
	cell_defaults.phenotype.secretion.secretion_rates[glutamine_substrate_index] = 0.0; 
	cell_defaults.phenotype.secretion.saturation_densities[glutamine_substrate_index] = 0.0;
	
	cell_defaults.custom_data.add_variable( "energy", "dimensionless" , parameters.doubles("cell_default_inital_energy") ); 
	cell_defaults.custom_data.add_variable( "cycle_energy_threshold", "dimensionless" , parameters.doubles("cell_default_cycle_energy_threshold") ); 
	cell_defaults.custom_data.add_variable( "death_energy_threshold", "dimensionless" , parameters.doubles("cell_default_death_energy_threshold") );
	cell_defaults.custom_data.add_variable( "G01S_thr", "dimensionless" , parameters.doubles("organoid_cell_G01S_thr") );
	
	
	cell_defaults.custom_data.add_variable( "k_aerobe", "dimensionless" , parameters.doubles("k_aerobe") );
	cell_defaults.custom_data.add_variable( "k_anaerobe", "dimensionless" , parameters.doubles("k_anaerobe") );
	cell_defaults.custom_data.add_variable( "k_glut", "dimensionless" , parameters.doubles("k_glut") );
	cell_defaults.custom_data.add_variable( "k_usage", "dimensionless" , parameters.doubles("k_usage") );
	cell_defaults.custom_data.add_variable( "alpha", "dimensionless" , parameters.doubles("alpha") );
	cell_defaults.custom_data.add_variable( "beta", "dimensionless" , parameters.doubles("beta") );
	cell_defaults.custom_data.add_variable( "gamma", "dimensionless" , parameters.doubles("gamma") );
	// ---- END -- Default Cell Definitions -- END ---- //
	
	
	
	
	
	// ---- START -- Organoid Cell Definitions -- START ---- //
	
	organoid_cell = cell_defaults; 
	organoid_cell.type = 1;
	organoid_cell.name = "organoid cell"; 
	
	// make sure the new cell type has its own reference phenotype
	organoid_cell.parameters.pReference_live_phenotype = &( organoid_cell.phenotype ); 
	// Set cell-cell adhesion to 5% of other cells 
	// organoid_cell.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "organoid_cell_relative_adhesion" ); // 0.05; 
	
	// Set Energy Function
	organoid_cell.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	organoid_cell.functions.update_phenotype = tumor_energy_update_function;
	
	// Set apoptosis to zero 
	organoid_cell.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "organoid_cell_apoptosis_rate" ); // 0.0; 
	organoid_cell.phenotype.cycle.data.transition_rate(i_Ki67_negative,i_Ki67_positive) = 0.0;
	organoid_cell.phenotype.cycle.data.transition_rate(i_Ki67_positive,i_Ki67_negative) = 0.0;
	
	// ---- END -- Organoid Cell Definitions -- END ---- //
	
	
	
	// ---- START -- Fibroblast Definitions -- START ---- //
	fibro_cell = cell_defaults; 
	fibro_cell.type = 2; 
	fibro_cell.name = "fibroblast"; 
	
	// make sure the new cell type has its own reference phenotype
	
	fibro_cell.parameters.pReference_live_phenotype = &( fibro_cell.phenotype ); 
	// Set cell-cell adhesion to 5% of other cells 
	//fibro_cell.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "fibroblast_relative_adhesion" );
	
	// Set apoptosis to zero 
	fibro_cell.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "fibroblast_apoptosis_rate" ); // 0.0; 
	fibro_cell.phenotype.cycle.data.transition_rate(i_Ki67_negative,i_Ki67_positive) = 0.0;
	fibro_cell.phenotype.cycle.data.transition_rate(i_Ki67_positive,i_Ki67_negative) = 0.0;
	
	fibro_cell.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 0.0;
	fibro_cell.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0.0; 
	fibro_cell.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 0.0; 
	
	fibro_cell.phenotype.secretion.uptake_rates[glucose_substrate_index] = 0.0; 
	fibro_cell.phenotype.secretion.secretion_rates[glucose_substrate_index] = 0.0; 
	fibro_cell.phenotype.secretion.saturation_densities[glucose_substrate_index] = 0.0; 
	
	fibro_cell.phenotype.secretion.uptake_rates[lactate_substrate_index] = 0.0; 
	fibro_cell.phenotype.secretion.secretion_rates[lactate_substrate_index] = 0.0; 
	fibro_cell.phenotype.secretion.saturation_densities[lactate_substrate_index] = 0.0; 
	
	fibro_cell.phenotype.secretion.uptake_rates[glutamine_substrate_index] = 0.0;
	fibro_cell.phenotype.secretion.secretion_rates[glutamine_substrate_index] = 0.0; 
	fibro_cell.phenotype.secretion.saturation_densities[glutamine_substrate_index] = 0.0;
	
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
	
	// set tracking internal substrates on
	default_microenvironment_options.track_internalized_substrates_in_each_agent = true;
	

	// ---- START -- Adding Substrates to Microenvironment -- START ---- //	
	
	pME->add_density( "glucose", "concentration" , 1.6e3 , 0.00 );
	int glucose_substrate_index = pME->find_density_index( "glucose" );
	default_microenvironment_options.Dirichlet_condition_vector[glucose_substrate_index] = 10.0;
	default_microenvironment_options.Dirichlet_activation_vector[glucose_substrate_index] = false;
	
	
	pME->add_density( "lactate", "concentration" , 1.6e3 , 0.00 );
	int lactate_substrate_index = pME->find_density_index( "lactate" );
	default_microenvironment_options.Dirichlet_condition_vector[lactate_substrate_index] = 0.0;
	default_microenvironment_options.Dirichlet_activation_vector[lactate_substrate_index] = false;
	
	
	pME->add_density( "glutamine", "concentration" , 1.6e3 , 0.00 );
	int glutamine_substrate_index = pME->find_density_index( "glutamine" );
	default_microenvironment_options.Dirichlet_condition_vector[glutamine_substrate_index] = 10.0;
	default_microenvironment_options.Dirichlet_activation_vector[glutamine_substrate_index] = false;
	
	// ---- END -- Adding Substrates to Microenvironment -- END ---- //
	
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0,10.0,0.0,10.0 }; 
	
	// initialize BioFVM 
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	// ---- START -- Setting Up Fibroblast Surface -- START ---- //
	
	Cell* pCell;

	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	double organoid_distance = parameters.doubles("organoid_distance");
	double initial_tumor_radius =  parameters.doubles("initial_tumor_radius");
	
	std::cout << "creating fibroblasts" << std::endl;
	
	for (int i= -500; i<500; i+=10)
	{
		for (int j= -500; j<500; j+=10)
		{			
			pCell = create_cell(fibro_cell);
			pCell->assign_position(i,-800,j);	
		}
	}
	
	// ---- END -- Setting Up Fibroblast Surface -- END ---- //



	// ---- START -- Setting Up Organoid -- START ---- //
	std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius,initial_tumor_radius); 
	std::cout << "creating " << positions.size() << " closely-packed organoid cells ... " << std::endl; 
	// create organoid
	for( int i=0; i < positions.size(); i++ )
	{
		positions[i][1] += organoid_distance-730;
		pCell = create_cell(organoid_cell);
		pCell->assign_position( positions[i] );
	}
	// ---- END -- Setting Up Organoid -- END ---- //

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
		 output[0] = "red"; 
		 output[2] = "darkred"; 	
	}
	
	if( pCell->phenotype.death.dead == false && 
		pCell->type == 2 )
	{
		 output[0] = "blue"; 
		 output[2] = "darkblue"; 	
	}
	
	return output; 
}


void tumor_energy_update_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	// ---- START -- Tumor Function -- START ---- //
	Microenvironment* pME = get_default_microenvironment();

	static int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	static int glucose_substrate_index = microenvironment.find_density_index( "glucose" );
	static int lactate_substrate_index = microenvironment.find_density_index( "lactate" );
	static int glutamine_substrate_index = microenvironment.find_density_index( "glutamine" );
	
	static int nk_aerobe = pCell->custom_data.find_variable_index( "k_aerobe" );
	static int nk_anaerobe = pCell->custom_data.find_variable_index( "k_anaerobe" );
	static int nk_glut = pCell->custom_data.find_variable_index( "k_glut" );
	static int nk_usage = pCell->custom_data.find_variable_index( "k_usage" );
	static int nalpha = pCell->custom_data.find_variable_index( "alpha" );
	static int nbeta = pCell->custom_data.find_variable_index( "beta" );
	static int ngamma = pCell->custom_data.find_variable_index( "gamma" );
	static int nE = pCell->custom_data.find_variable_index("energy");
	static int nDeath = pCell->custom_data.find_variable_index("death_energy_threshold");
	static int nBirth = pCell->custom_data.find_variable_index("cycle_energy_threshold");
	static int nG01S_thr = pCell->custom_data.find_variable_index("G01S_thr");
	static int i_Ki67_negative = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_negative );
	static int i_Ki67_positive = Ki67_basic.find_phase_index( PhysiCell_constants::Ki67_positive );
	static int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	
	double internal_oxygen = phenotype.molecular.internalized_total_substrates[oxygen_substrate_index];
	double internal_glucose = phenotype.molecular.internalized_total_substrates[glucose_substrate_index];
	double internal_lactate = phenotype.molecular.internalized_total_substrates[lactate_substrate_index];
	double internal_glutamine = phenotype.molecular.internalized_total_substrates[glutamine_substrate_index];

	// Energy function
	pCell->custom_data[nE] += dt*(internal_glucose * internal_oxygen * pCell->custom_data[nk_aerobe] + internal_glucose * pCell->custom_data[nk_anaerobe] + internal_glutamine * pCell->custom_data[nk_glut] - pCell->custom_data[nE] * pCell->custom_data[nk_usage]);
	
	// Metabolite reactions (Consuming & Producing Internalized Metabolites)
	phenotype.molecular.internalized_total_substrates[glucose_substrate_index] += dt*( -1.0 * internal_glucose * internal_oxygen * pCell->custom_data[nk_aerobe] * pCell->custom_data[nalpha] - internal_glucose * pCell->custom_data[nk_anaerobe] * pCell->custom_data[nbeta]);
	phenotype.molecular.internalized_total_substrates[oxygen_substrate_index] += dt*( -1.0 * internal_glucose * internal_oxygen * pCell->custom_data[nk_aerobe] * pCell->custom_data[nalpha]);
	phenotype.molecular.internalized_total_substrates[lactate_substrate_index] += dt*(internal_glucose * pCell->custom_data[nk_anaerobe] * pCell->custom_data[nbeta]);
	phenotype.molecular.internalized_total_substrates[glutamine_substrate_index] += dt*( -1.0 * internal_glutamine * pCell->custom_data[nk_glut] * pCell->custom_data[ngamma]);
	
	//std::cout<< pCell->custom_data[nE] << std::endl;
	//std::cout<< pCell->ID << std::endl;
	
	//std::cout<< phenotype.molecular.internalized_total_substrates[lactate_substrate_index] << std::endl;
	
	if (pCell->ID == 10002)
	{
	//std::cout << "Intracellular Oxygen:"	<< phenotype.molecular.internalized_total_substrates[oxygen_substrate_index] << std::endl;
	//std::cout<< "Intracellular Glucose:"	<< phenotype.molecular.internalized_total_substrates[glucose_substrate_index] << std::endl;
	//std::cout<< "Intracellular Lactate:"	<<phenotype.molecular.internalized_total_substrates[lactate_substrate_index] << std::endl;
	//std::cout<< "Intracellular Glutamine:"	<<phenotype.molecular.internalized_total_substrates[glutamine_substrate_index] << std::endl;
	//std::cout<< "Energy:"	<< pCell->custom_data[nE] << std::endl;
	}
	

	// die if energy is low 
	if( pCell->custom_data[nE] < pCell->custom_data[nDeath] )
	{
	//std::cout<< pCell->custom_data[nE] << std::endl;
	phenotype.death.rates[apoptosis_model_index] = parameters.doubles("apoptosis_rate");
	}
		
	if( pCell->custom_data[nE] > pCell->custom_data[nG01S_thr])
	{
	phenotype.cycle.data.transition_rate( i_Ki67_negative,i_Ki67_positive ) = parameters.doubles("organoid_cell_r01");
    //std::cout<< phenotype.cycle.data.transition_rate( i_Ki67_negative,i_Ki67_positive ) << std::endl;
	}
	
	if( pCell->custom_data[nE] > pCell->custom_data[nBirth])
	{
	phenotype.cycle.data.transition_rate( i_Ki67_positive,i_Ki67_negative ) = parameters.doubles("organoid_cell_r10"); 
	//std::cout<< phenotype.cycle.data.transition_rate( i_Ki67_positive,i_Ki67_negative ) << std::endl;
	}

	//std::cout<< phenotype.molecular.internalized_total_substrates[phenotype.cycle.data.transition_rate( i_Ki67_negative,i_Ki67_positive )] << std::endl;
	//std::cout<< phenotype.molecular.internalized_total_substrates[phenotype.cycle.data.transition_rate( i_Ki67_positive,i_Ki67_negative )] << std::endl;
	
	// ---- END -- Tumor Function -- END ---- //
	return;
}

void metabolic_equation_solver( Cell* pCell,Phenotype& phenotype, double dt )
{
	
	
	
	
	
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