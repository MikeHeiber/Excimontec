// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef OSC_SIM_H
#define OSC_SIM_H

#include "Simulation.h"
#include "Site.h"
#include "Exciton.h"
#include "Parameters.h"
#include "Polaron.h"
#include "Version.h"
#include <algorithm>
#include <numeric>

namespace Excimontec {

	class Site_OSC : public KMC_Lattice::Site {
	public:
		float getEnergy() const { return energy; }
		short getType() const { return type; }
		void setEnergy(const float energy_input) { energy = energy_input; }
		void setType(const short site_type) { type = site_type; }
	private:
		float energy;
		short type = 0; //  type 1 represent donor, type 2 represents acceptor
	};

	class OSC_Sim : public KMC_Lattice::Simulation {
	public:
		// Functions
		OSC_Sim();
		virtual ~OSC_Sim();
		bool init(const Parameters& params, const int id);
		void calculateAllEvents();
		void calculateDOSCorrelation();
		void calculateDOSCorrelation(const double cutoff_radius);
		std::vector<std::pair<double, double>> calculateTransitTimeHist(const std::vector<double>& data, const int counts) const;
		std::vector<double> calculateMobilityData(const std::vector<double>& transit_times) const;
		bool checkFinished() const;
		void createExciton(const bool spin);
		void createExciton(const KMC_Lattice::Coords& coords, const bool spin);
		void createElectron(const KMC_Lattice::Coords& coords);
		void createHole(const KMC_Lattice::Coords& coords);
		bool executeNextEvent();
		void exportEnergies(std::string filename);
		std::vector<std::string> getChargeExtractionMap(const bool charge) const;
		std::vector<std::pair<double, double>> getDOSCorrelationData() const;
		std::vector<double> getDynamicsExcitonEnergies() const;
		std::vector<double> getDynamicsElectronEnergies() const;
		std::vector<double> getDynamicsHoleEnergies() const;
		std::vector<int> getDynamicsTransientSinglets() const;
		std::vector<int> getDynamicsTransientTriplets() const;
		std::vector<int> getDynamicsTransientElectrons() const;
		std::vector<int> getDynamicsTransientHoles() const;
		std::vector<double> getDynamicsTransientTimes() const;
		std::vector<double> getDynamicsExcitonMSDV() const;
		std::vector<double> getDynamicsElectronMSDV() const;
		std::vector<double> getDynamicsHoleMSDV() const;
		std::vector<double> getExcitonDiffusionData() const;
		std::vector<int> getExcitonHopLengthData() const;
		std::vector<double> getExcitonLifetimeData() const;
		double getInternalField() const;
		int getN_excitons_created() const;
		int getN_excitons_created(const short site_type) const;
		int getN_singlet_excitons_dissociated() const;
		int getN_triplet_excitons_dissociated() const;
		int getN_singlet_excitons_recombined() const;
		int getN_triplet_excitons_recombined() const;
		int getN_singlet_singlet_annihilations() const;
		int getN_singlet_triplet_annihilations() const;
		int getN_triplet_triplet_annihilations() const;
		int getN_singlet_polaron_annihilations() const;
		int getN_triplet_polaron_annihilations() const;
		int getN_electrons_created() const;
		int getN_electrons_collected() const;
		int getN_electrons_recombined() const;
		long int getN_events_executed() const;
		int getN_holes_created() const;
		int getN_holes_collected() const;
		int getN_holes_recombined() const;
		int getN_geminate_recombinations() const;
		int getN_bimolecular_recombinations() const;
		int getN_transient_cycles() const;
		std::string getPreviousEventType() const;
		std::vector<float> getSiteEnergies(const short site_type) const;
		float getSiteEnergy(const KMC_Lattice::Coords& coords) const;
		short getSiteType(const KMC_Lattice::Coords& coords) const;
		double getSteadyEquilibrationEnergy() const;
		double getSteadyFermiEnergy() const;
		double getSteadyMobility() const;
		double getSteadyTransportEnergy() const;
		std::vector<int> getToFTransientCounts() const;
		std::vector<double> getToFTransientEnergies() const;
		std::vector<double> getToFTransientTimes() const;
		std::vector<double> getToFTransientVelocities() const;
		std::vector<double> getTransitTimeData() const;
		void outputStatus();
		void reassignSiteEnergies();
	protected:

	private:

		struct ExcitonEventCalcVars {
			int range;
			int dim;
			Exciton::Hop hop_event;
			std::vector<Exciton::Hop> hops_temp;
			Exciton::Dissociation diss_event;
			std::vector<Exciton::Dissociation> dissociations_temp;
			Exciton::Exciton_Annihilation ee_annihilation_event;
			std::vector<Exciton::Exciton_Annihilation> ee_annihilations_temp;
			Exciton::Polaron_Annihilation ep_annihilation_event;
			std::vector<Exciton::Polaron_Annihilation> ep_annihilations_temp;
			std::vector<bool> hops_valid;
			std::vector<bool> dissociations_valid;
			std::vector<bool> ee_annihilations_valid;
			std::vector<bool> ep_annihilations_valid;
			// pre-calculated distances vector that contains the distances to nearby sites used for event execution time calculations
			std::vector<double> distances;
			// pre-calculated isInDissRange and isInFRETRange vectors that contains booleans to indicate whether the nearby sites are within range for the different exciton events to be possible.
			std::vector<bool> isInDissRange;
			std::vector<bool> isInFRETRange;

			ExcitonEventCalcVars() {}

			ExcitonEventCalcVars(OSC_Sim* sim_ptr) {
				range = (int)ceil(((sim_ptr->params.FRET_cutoff > sim_ptr->params.Exciton_dissociation_cutoff) ? (sim_ptr->params.FRET_cutoff) : (sim_ptr->params.Exciton_dissociation_cutoff)) / sim_ptr->lattice.getUnitSize());
				dim = (2 * range + 1);
				hop_event = Exciton::Hop(sim_ptr);
				hops_temp.assign(dim*dim*dim, hop_event);
				diss_event = Exciton::Dissociation(sim_ptr);
				dissociations_temp.assign(dim*dim*dim, diss_event);
				ee_annihilation_event = Exciton::Exciton_Annihilation(sim_ptr);
				ee_annihilations_temp.assign(dim*dim*dim, ee_annihilation_event);
				ep_annihilation_event = Exciton::Polaron_Annihilation(sim_ptr);
				ep_annihilations_temp.assign(dim*dim*dim, ep_annihilation_event);
				hops_valid.assign(dim*dim*dim, false);
				dissociations_valid.assign(dim*dim*dim, false);
				ee_annihilations_valid.assign(dim*dim*dim, false);
				ep_annihilations_valid.assign(dim*dim*dim, false);
				// pre-calculated distances vector that contains the distances to nearby sites used for event execution time calculations
				distances.assign(dim*dim*dim, 0.0);
				// pre-calculated isInDissRange and isInFRETRange vectors that contains booleans to indicate whether the nearby sites are within range for the different exciton events to be possible.
				isInDissRange.assign(dim*dim*dim, false);
				isInFRETRange.assign(dim*dim*dim, false);
				// Initiaize distances, isInDissRange, and isInFRETRange vectors
				for (int i = -range; i <= range; i++) {
					for (int j = -range; j <= range; j++) {
						for (int k = -range; k <= range; k++) {
							int index = (i + range)*dim*dim + (j + range)*dim + (k + range);
							distances[index] = sim_ptr->lattice.getUnitSize()*sqrt((double)(i*i + j * j + k * k));
							if (!((distances[index] - 0.0001) > sim_ptr->params.Exciton_dissociation_cutoff)) {
								isInDissRange[index] = true;
							}
							if (!((distances[index] - 0.0001) > sim_ptr->params.FRET_cutoff)) {
								isInFRETRange[index] = true;
							}
						}
					}
				}
			}
		};
		ExcitonEventCalcVars exciton_event_calc_vars;

		struct PolaronEventCalcVars {
			int range;
			int dim;
			Polaron::Hop hop_event;
			std::vector<Polaron::Hop> hops_temp;
			Polaron::Recombination rec_event;
			std::vector<Polaron::Recombination> recombinations_temp;
			std::vector<bool> hops_valid;
			std::vector<bool> recombinations_valid;
			// pre-calculated distances vector that contains the distances to nearby sites used for event execution time calculations
			std::vector<double> distances;
			std::vector<double> E_deltas;
			// pre-calculated isInRange vector that contains booleans to indicate if the nearby sites are within range for polaron events to be possible.
			std::vector<bool> isInRange;

			PolaronEventCalcVars() {}

			PolaronEventCalcVars(OSC_Sim* sim_ptr) {
				range = (int)ceil(sim_ptr->params.Polaron_hopping_cutoff / sim_ptr->lattice.getUnitSize());
				dim = (2 * range + 1);
				hop_event = Polaron::Hop(sim_ptr);
				hops_temp.assign(dim*dim*dim, hop_event);
				rec_event = Polaron::Recombination(sim_ptr);
				recombinations_temp.assign(dim*dim*dim, rec_event);
				hops_valid.assign(dim*dim*dim, false);
				recombinations_valid.assign(dim*dim*dim, false);
				// pre-calculated distances vector that contains the distances to nearby sites used for event execution time calculations
				distances.assign(dim*dim*dim, 0.0);
				E_deltas.assign(dim*dim*dim, 0.0);
				// pre-calculated isInRange vector that contains booleans to indicate if the nearby sites are within range for polaron events to be possible.
				isInRange.assign(dim*dim*dim, false);
				// Intialize distances and isInRange vectors
				for (int i = -range; i <= range; i++) {
					for (int j = -range; j <= range; j++) {
						for (int k = -range; k <= range; k++) {
							int index = (i + range)*dim*dim + (j + range)*dim + (k + range);
							distances[index] = sim_ptr->lattice.getUnitSize()*sqrt((double)(i*i + j * j + k * k));
							if (!((distances[index] - 0.0001) > sim_ptr->params.Polaron_hopping_cutoff)) {
								isInRange[index] = true;
							}
						}
					}
				}
			}
		};
		PolaronEventCalcVars polaron_event_calc_vars;
		// Input Parameters
		Parameters params;
		// Additional Derived Parameters
		double Transient_start;
		double Transient_end;
		int Transient_pnts_per_decade;
		bool isLightOn;
		double R_exciton_generation_donor;
		double R_exciton_generation_acceptor;
		double Transient_step_size;
		double Transient_creation_time;
		int Transient_index_prev;
		int Transient_singlet_counts_prev;
		int Transient_triplet_counts_prev;
		int Transient_electron_counts_prev;
		int Transient_hole_counts_prev;
		int Coulomb_range;
		double AvgDielectric;
		double Image_interaction_prefactor;
		int N_initial_excitons;
		// Site Data Structure
		std::vector<Site_OSC> sites;
		// Object Data Structures
		std::list<Exciton> excitons;
		std::list<Polaron> electrons;
		std::list<Polaron> holes;
		// Event Data Structures
		std::string previous_event_type = "";
		double previous_event_time = 0;
		std::list<Exciton::Creation> exciton_creation_events;
		std::list<KMC_Lattice::Event*>::const_iterator exciton_creation_it;
		std::list<Exciton::Hop> exciton_hop_events;
		std::list<Exciton::Recombination> exciton_recombination_events;
		std::list<Exciton::Dissociation> exciton_dissociation_events;
		std::list<Exciton::Exciton_Annihilation> exciton_exciton_annihilation_events;
		std::list<Exciton::Polaron_Annihilation> exciton_polaron_annihilation_events;
		std::list<Exciton::Intersystem_Crossing> exciton_intersystem_crossing_events;
		std::list<Polaron::Hop> electron_hop_events;
		std::list<Polaron::Hop> hole_hop_events;
		std::list<Polaron::Recombination> polaron_recombination_events;
		std::list<Polaron::Extraction> electron_extraction_events;
		std::list<Polaron::Extraction> hole_extraction_events;
		// Additional Data Structures
		std::vector<double> Coulomb_table;
		std::vector<double> E_potential;
		std::vector<std::pair<double, double>> DOS_correlation_data;
		std::vector<double> exciton_lifetimes;
		std::vector<double> exciton_diffusion_distances;
		std::vector<int> exciton_hop_distances; // saved in lattice units squared
		std::vector<int> transient_exciton_tags;
		std::vector<int> transient_electron_tags;
		std::vector<int> transient_hole_tags;
		std::vector<int> ToF_positions_prev;
		std::vector<double> transient_exciton_energies_prev;
		std::vector<double> transient_electron_energies_prev;
		std::vector<double> transient_hole_energies_prev;
		std::vector<double> transient_exciton_msdv;
		std::vector<double> transient_electron_msdv;
		std::vector<double> transient_hole_msdv;
		std::vector<int> electron_extraction_data;
		std::vector<int> hole_extraction_data;
		std::vector<double> transient_times;
		std::vector<double> transient_velocities;
		std::vector<double> transient_exciton_energies;
		std::vector<double> transient_electron_energies;
		std::vector<double> transient_hole_energies;
		std::vector<double> transit_times;
		std::vector<int> transient_singlet_counts;
		std::vector<int> transient_triplet_counts;
		std::vector<int> transient_electron_counts;
		std::vector<int> transient_hole_counts;
		double Steady_equilibration_energy_sum = 0.0;
		double Steady_equilibration_time = 0.0;
		double Transport_energy_weighted_sum = 0.0;
		double Transport_energy_sum_of_weights = 0.0;
		// Additional Counters
		int N_donor_sites = 0;
		int N_acceptor_sites = 0;
		int N_excitons_created = 0;
		int N_excitons_created_donor = 0;
		int N_excitons_created_acceptor = 0;
		int N_singlet_excitons_recombined = 0;
		int N_triplet_excitons_recombined = 0;
		int N_singlet_excitons_dissociated = 0;
		int N_triplet_excitons_dissociated = 0;
		int N_singlet_singlet_annihilations = 0;
		int N_singlet_triplet_annihilations = 0;
		int N_triplet_triplet_annihilations = 0;
		int N_singlet_polaron_annihilations = 0;
		int N_triplet_polaron_annihilations = 0;
		int N_exciton_intersystem_crossings = 0;
		int N_exciton_reverse_intersystem_crossings = 0;
		int N_excitons_quenched = 0;
		int N_excitons = 0;
		int N_singlets = 0;
		int N_triplets = 0;
		int N_electrons_created = 0;
		int N_electrons_recombined = 0;
		int N_electrons_collected = 0;
		int N_electrons = 0;
		long int N_events_executed = 0;
		int N_holes_created = 0;
		int N_holes_recombined = 0;
		int N_holes_collected = 0;
		int N_holes = 0;
		int N_geminate_recombinations = 0;
		int N_bimolecular_recombinations = 0;
		int N_electron_surface_recombinations = 0;
		int N_hole_surface_recombinations = 0;
		int N_transient_cycles = 0;
		// Additional Functions
		double calculateCoulomb(const std::list<Polaron>::const_iterator polaron_it, const KMC_Lattice::Coords& coords) const;
		double calculateCoulomb(const bool charge, const KMC_Lattice::Coords& coords) const;
		KMC_Lattice::Coords calculateRandomExcitonCreationCoords();
		void calculateExcitonEvents(Exciton* exciton_ptr);
		void calculateObjectListEvents(const std::vector<KMC_Lattice::Object*>& object_ptr_vec);
		void calculatePolaronEvents(Polaron* polaron_ptr);
		void createCorrelatedDOS(const double correlation_length);
		bool createImportedMorphology();
		void deleteObject(KMC_Lattice::Object* object_ptr);
		// Exciton Event Execution Functions
		bool executeExcitonCreation();
		bool executeExcitonHop(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executeExcitonRecombination(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executeExcitonDissociation(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executeExcitonIntersystemCrossing(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executeExcitonExcitonAnnihilation(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executeExcitonPolaronAnnihilation(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		// General Event Functions
		bool executeObjectHop(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		// Polaron Event Execution Functions
		bool executePolaronHop(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executePolaronRecombination(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executePolaronExtraction(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		KMC_Lattice::Coords generateExciton();
		void generateExciton(const KMC_Lattice::Coords& coords, const bool spin, int tag);
		void generateElectron(const KMC_Lattice::Coords& coords, int tag);
		void generateHole(const KMC_Lattice::Coords& coords, int tag);
		void generateDynamicsExcitons();
		void generateSteadyPolarons();
		void generateToFPolarons();
		std::list<Exciton>::iterator getExcitonIt(const KMC_Lattice::Object* object_ptr);
		std::list<Polaron>::iterator getPolaronIt(const KMC_Lattice::Object* object_ptr);
		bool initializeArchitecture();
		void removeExciton(std::list<Exciton>::iterator exciton_it);
		bool siteContainsHole(const KMC_Lattice::Coords& coords);
		void updateSteadyData();
		void updateTransientData();
	};

}

#endif //OSC_SIM_H
