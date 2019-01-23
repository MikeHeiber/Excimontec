// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef EXCIMONTEC_OSC_SIM_H
#define EXCIMONTEC_OSC_SIM_H

#include "Simulation.h"
#include "Site.h"
#include "Exciton.h"
#include "Parameters.h"
#include "Polaron.h"
#include "Version.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace Excimontec {

	//! \brief This class extends the KMC_Lattice::Simulation class to create a functioning KMC simulation object for organic semiconductor devices.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2017-2019
	class OSC_Sim : public KMC_Lattice::Simulation {
	public:

		//! \brief Constructs an empty simulation object that is uninitialized.
		OSC_Sim();

		//! Default virtual destructor.
		virtual ~OSC_Sim();

		//! \brief Initializes the simulation object so that it is ready to execute a simulation test.
		//! \param params specifies all of the input parameters needed to run the simulation.
		//! \param id defines the desired ID number of the simulation object.
		//! \return true if the initialization is successful.
		//! \return false if there are any errors during initialization.
		bool init(const Parameters& params, const int id);

		//! \brief Calculates the events for all objects in the simulation.
		void calculateAllEvents();

		//! \brief Calculates the transit time probability histogram using the data generated by the time-of-flight charge transport test.
		//! \param data is a vector of transit times for all extracted polarons.
		//! \param counts is the total number of polarons tested.
		//! \returns A pair vector where the first value is the bin position and the second value is the probability.
		std::vector<std::pair<double, double>> calculateTransitTimeHist(const std::vector<double>& data, const int counts) const;

		//! \brief Calculate the mobility data using the transit time data generated by the time-of-flight charge transport test.
		//! \param transit_times is a vector of transit times for all extracted polarons.
		//! \return A vector of calculated mobility values for each extracted polaron.
		std::vector<double> calculateMobilityData(const std::vector<double>& transit_times) const;

		//! \brief Checks whether the simulation test specified by the input parameters is finished.
		//! \return true if the simulation test is finished.
		//! \return false if the simulation test is not yet finished.
		bool checkFinished() const;

		//! \brief Creates an Exciton on a randomly selected unoccupied site.
		//! \param spin specifies the spin state of the created exciton. (true for singlet and false for triplet)
		void createExciton(const bool spin);

		//! \brief Attempts to create an Exciton on the specified lattice site.
		//! An error is generated if the specified lattice site is already occupied or if the coordinates are not in the lattice.
		//! \param coords is the coordinates where the exciton will be created.
		//! \param spin specifies the spin state of the created exciton. (true for singlet and false for triplet)
		void createExciton(const KMC_Lattice::Coords& coords, const bool spin);

		//! \brief Attempts to create an electron (negatively charged Polaron) on the specified lattice site.
		//! An error is generated if the specified lattice site is already occupied or if the coordinates are not in the lattice.
		//! \param coords is the coordinates where the electron will be created.
		void createElectron(const KMC_Lattice::Coords& coords);

		//! \brief Attempts to create a hole (positively charged Polaron) on the specified lattice site.
		//! An error is generated if the specified lattice site is already occupied or if the coordinates are not in the lattice.
		//! \param coords is the coordinates where the hole will be created.
		void createHole(const KMC_Lattice::Coords& coords);

		//! \brief Attempts to execute the next event in the KMC simulation.
		//! An error is generated if an event cannot be executed.
		//! \return true if the next event is executed successfully.
		//! \return false if the next event cannot be executed.
		bool executeNextEvent();

		//! \brief Exports the lattice site energies to a text file.
		//! \param filename is the name of the file that will be created in the working directory.
		void exportEnergies(std::string filename);

		//! \brief Gets the charge extraction map data generated by the time-of-flight charge transport or internal quantum efficiency tests. 
		//! \param charge specifies whether to get electron or hole polaron extraction data.
		//! \return A string vector that can be separately output to a file.
		std::vector<std::string> getChargeExtractionMap(const bool charge) const;

		//! \brief Gets the radial autocorrelation data for the lattice site energies that is generated when using the correlated Gaussian DOS model.
		//! \return A pair vector where the first value is the radial distance and the second is the autocorrelation probability value.
		std::vector<std::pair<double, double>> getDOSCorrelationData() const;

		//! \brief Gets the transient average Exciton energy data generated by the dynamics test.
		//! \returns A vector of data representing how the average exciton energy in units of eV changes over time.
		std::vector<double> getDynamicsExcitonEnergies() const;

		//! \brief Gets the transient average electron (negatively charged Polaron) energy data generated during the dynamics test.
		//! \returns A vector of data representing how the average electron energy in units of eV changes over time.
		std::vector<double> getDynamicsElectronEnergies() const;

		//! \brief Gets the transient average electron (positively charged Polaron) energy data generated during the dynamics test.
		//! \returns A vector of data representing how the average hole energy in units of eV changes over time.
		std::vector<double> getDynamicsHoleEnergies() const;

		//! \brief Gets the transient singlet exciton counts data generated by the dynamics test.
		//! \returns A vector of data representing how the number of singlet excitons in the simulation changes over time.
		std::vector<int> getDynamicsTransientSinglets() const;

		//! \brief Gets the transient triplet exciton counts data generated by the dynamics test.
		//! \returns A vector of data representing how the number of triplet excitons in the simulation changes over time.
		std::vector<int> getDynamicsTransientTriplets() const;

		//! \brief Gets the transient electron (negatively charged Polaron) counts data generated during the dynamics test.
		//! \returns A vector of data representing how the number of electrons in the simulation changes over time.
		std::vector<int> getDynamicsTransientElectrons() const;

		//! \brief Gets the transient hole (positively charged Polaron) counts data generated during the dynamics test.
		//! \returns A vector of data representing how the number of holes in the simulation changes over time.
		std::vector<int> getDynamicsTransientHoles() const;

		//! \brief Gets the transient time data generated during the dynamics test.
		//! \returns A vector of data representing the time points for the transient data.
		std::vector<double> getDynamicsTransientTimes() const;

		//! \brief Gets the transient exciton mean squared displacement velocity data generated during the dynamics test.
		//! The mean squared displacement velocity is equivalent to the diffusion coefficient for excitons. 
		//! \returns A vector of data representing how the exciton mean squared displacement velocity changes over time.
		std::vector<double> getDynamicsExcitonMSDV() const;

		//! \brief Gets the transient electron (negatively charged Polaron) mean squared displacement velocity data generated by the dynamics test.
		//! The mean squared displacement velocity equals the diffusion coefficient at zero-bias conditions
		//! or can be used to calculate the mobility using the applied internal electrical potential.
		//! \returns A vector of data representing how the electron mean squared displacement velocity changes over time.
		std::vector<double> getDynamicsElectronMSDV() const;

		//! \brief Gets the transient hole (positively charged Polaron) mean squared displacement velocity data generated by the dynamics test.
		//! The mean squared displacement velocity equals the diffusion coefficient at zero-bias conditions 
		//! or can be used to calculate the mobility using the applied internal electrical potential.
		//! \returns A vector of data representing how the hole mean squared displacement velocity changes over time.
		std::vector<double> getDynamicsHoleMSDV() const;

		//! \brief Gets the exciton diffusion distance data generated during the exciton diffusion test.
		//! \return A vector of data representing the displacement distance of each exciton tested.
		std::vector<double> getExcitonDiffusionData() const;

		//! \brief Gets the Exciton hop distance data generated during the exciton diffusion test.
		//! /return A vector of data representing the hop distance of all exciton hop events in lattice units.
		std::vector<int> getExcitonHopLengthData() const;

		//! \brief Gets the exciton lifetime data generated during the exciton diffusion test.
		//! \return A vector of data representing the lifetime of each exciton tested.
		std::vector<double> getExcitonLifetimeData() const;

		//! \brief Gets the internal electric field due to applied internal potential.
		//! \return the internal electric field in units of V/cm.
		double getInternalField() const;

		//! \brief Gets the number of excitons that have been created in the simulation.
		//! \return the number of excitons that have been created since the simulation object was initialized.
		int getN_excitons_created() const;

		//! \brief Gets the number of excitons that have been created on the specified site type in the simulation.
		//! \param site_type specifies the site type, where 1 is for donor sites and 2 is for acceptor sites.
		//! \return The number of excitons that have been created on the specified site type since the simulation object was initialized.
		int getN_excitons_created(const short site_type) const;

		//! \brief Gets the number of singlet excitons that have dissociated.
		//! \return The number of singlet excitons that have dissociated since the simulation object was initialized.
		int getN_singlet_excitons_dissociated() const;

		//! \brief Gets the number of triplet excitons that have dissociated.
		//! \return The number of triplet excitons that have dissociated since the simulation object was initialized.
		int getN_triplet_excitons_dissociated() const;

		//! \brief Gets the number of singlet excitons that have recombined.
		//! \return The number of singlet excitons that have recombined since the simulation object was initialized.
		int getN_singlet_excitons_recombined() const;

		//! \brief Gets the number of triplet excitons that have recombined.
		//! \return The number of triplet excitons that have recombined since the simulation object was initialized.
		int getN_triplet_excitons_recombined() const;

		//! \brief Gets the number of singlet-singlet annihilation events that have occurred.
		//! \return The number of singlet-singlet annihilation events that have occurred since the simulation object was initialized.
		int getN_singlet_singlet_annihilations() const;

		//! \brief Gets the number of singlet-triplet annihilation events that have occurred.
		//! \return The number of singlet-triplet annihilation events that have occurred since the simulation object was initialized.
		int getN_singlet_triplet_annihilations() const;

		//! \brief Gets the number of triplet-triplet annihilation events that have occurred.
		//! \return The number of triplet-triplet annihilation events that have occurred since the simulation object was initialized.
		int getN_triplet_triplet_annihilations() const;

		//! \brief Gets the number of singlet-polaron annihilation events that have occurred.
		//! \return The number of singlet-polaron annihilation events that have occurred since the simulation object was initialized.
		int getN_singlet_polaron_annihilations() const;

		//! \brief Gets the number of triplet-polaron annihilation events that have occurred.
		//! \return The number of triplet-polaron annihilation events that have occurred since the simulation object was initialized.
		int getN_triplet_polaron_annihilations() const;

		//! \brief Gets the number of electrons (negatively charged polarons) that have been created in the simulation.
		//! \return The number of electrons that have been created since the simulation object was initialized.
		int getN_electrons_created() const;

		//! \brief Gets the number of electrons (negatively charged polarons) that have been collected (extracted) at electrodes in the simulation.
		//! \return The number of electrons that have been collected (extracted) since the simulation object was initialized.
		int getN_electrons_collected() const;

		//! \brief Gets the number of electrons (negatively charged polarons) that have recombined with holes (positively charged polarons) in the simulation.
		//! \return The number of electron-hole pairs that have recombined since the simulation object was initialized.
		int getN_electrons_recombined() const;

		//! \brief Gets the number of events that have been executed.
		//! \return The number of events that have been executed since the simulation object was initialized.
		long int getN_events_executed() const;

		//! \brief Gets the number of holes (positively charged polarons) that have been created in the simulation.
		//! \return The number of holes that have been created since the simulation object was initialized.
		int getN_holes_created() const;

		//! \brief Gets the number of holes (positively charged polarons) that have been collected (extracted) at electrodes in the simulation.
		//! \return The number of holes that have been collected (extracted) since the simulation object was initialized.
		int getN_holes_collected() const;

		//! \brief Gets the number of holes (positively charged polarons) that have recombined with electrons (negatively charged polarons) in the simulation.
		//! \return The number of electron-hole pairs that have recombined since the simulation object was initialized.
		int getN_holes_recombined() const;

		//! \brief Gets the number of geminate electron-hole recombination events that have occurred in the simulation.
		//! Geminate recombination occurs between an electron and hole originating from the same exciton.
		//! \return The number of geminate electron-hole recombination events that have occurred since the simulation object was initialized.
		int getN_geminate_recombinations() const;

		//! \brief Gets the number of bimolecular electron-hole recombination events that have occurred in the simulation.
		//! Bimolecular recombination occurs between an electron and hole originating from different excitons.
		//! \return The number of bimolecular electron-hole recombination events that have occurred since the simulation object was initialized.
		int getN_bimolecular_recombinations() const;

		//! \brief Gets the number of transient test cycles that have been performed.
		//! During the time-of-flight charge transport test and the dynamics test, 
		//! multiple transient cycles may be used to reach the target number of tested objects.
		//! \return The number of transient test cycles that have been performed since the simulation object was initialized.
		int getN_transient_cycles() const;

		//! \brief Gets the name of event type that was last performed.
		//! \return The event_type name of the previously executed event.
		std::string getPreviousEventType() const;

		//! \brief Gets the energies of the lattice sites with the specified type.
		//! \param site_type specifies the site type, where 1 is for donor sites and 2 is for acceptor sites.
		//! \return A vector of data representing the lattice site energies of sites with the specified type in units of eV.
		std::vector<float> getSiteEnergies(const short site_type) const;

		//! \brief Gets the site energy of the lattice site at the specified coordinates.
		//! An error is generated if the specified lattice site is not in the lattice.
		//! \param coords is the coordinates of the specified lattice site.
		//! \return The energy of the specified site in units of eV.
		//! \return NAN if the coordinates correspond to a site that is not in the lattice.
		float getSiteEnergy(const KMC_Lattice::Coords& coords);

		//! \brief Gets the site type of the lattice site at the specified coordinates.
		//! An error is generated if the specified lattice site is not in the lattice.
		//! \param coords is the coordinates of the specified lattice site.
		//! \return The type of the specified site, where 1 is for a donor site and 2 is for an acceptor site.
		//! \return -1 if the coordinates correspond to a site that is not in the lattice.
		short getSiteType(const KMC_Lattice::Coords& coords);

		//! \brief Gets the average steady state current density calculated during the steady state charge transport test.
		//! \return The average current density in units of mA cm^-2.
		double getSteadyCurrentDensity() const;

		//! \brief Gets the average equilibration energy of the polarons calculated during the steady state charge transport test.
		//! The average equilibration energy corresponds to the average of the density of occupied states at steady state, quasi-equilibrium conditions.
		//! \return The calculated average equilibration energy of the polarons in units of eV.
		double getSteadyEquilibrationEnergy() const;

		//! \brief Gets the average equilibration energy of the polarons calculated during the steady state charge transport test.
		//! The average equilibration energy corresponds to the average of the density of occupied states at steady state, quasi-equilibrium conditions.
		//! This function returns the equilibration energy calculated including the Coulomb potential due to interactions between the polarons.
		//! \return The calculated average equilibration energy of the polarons in units of eV.
		double getSteadyEquilibrationEnergy_Coulomb() const;

		//! \brief Gets the Fermi energy of the system calculated during the steady state charge transport test.
		//! The Fermi energy corresponds to the highest occupied energy state under equilibrium conditions at 0 K.
		//! \return The Fermi energy of the system in units of eV.
		double getSteadyFermiEnergy() const;

		//! \brief Gets the Fermi energy of the system calculated during the steady state charge transport test.
		//! The Fermi energy corresponds to the highest occupied energy state under equilibrium conditions at 0 K.
		//! This Fermi energy calculation includes the Coulomb potential energy due to interactions between the polarons.
		//! \return The Fermi energy of the system in units of eV.
		double getSteadyFermiEnergy_Coulomb() const;

		//! \brief Gets the average steady state charge carrier mobility calculated during the steady state charge transport test.
		//! \return The average charge carrier mobility in units of cm^2 V^-1 s^-1.
		double getSteadyMobility() const;

		//! \brief Gets the average transport energy of the polarons calculated during the steady state charge transport test.
		//! \return The average transport energy of the polarons in units of eV.
		double getSteadyTransportEnergy() const;

		//! \brief Gets the average transport energy of the polarons calculated during the steady state charge transport test.
		//! This returns the transport energy calculated including the Coulomb potential energy due to interactions between the polarons.
		//! \return The average transport energy of the polarons in units of eV.
		double getSteadyTransportEnergy_Coulomb() const;

		//! \brief Gets the transient polaron counts data generated during the time-of-flight charge transport test.
		//! \returns A vector of data representing how the number of polarons in the simulation changes over time.
		std::vector<int> getToFTransientCounts() const;

		//! \brief Gets the transient average polaron energy data generated during the time-of-flight charge transport test.
		//! \returns A vector of data representing how the average polaron energy in units of eV changes over time.
		std::vector<double> getToFTransientEnergies() const;

		//! \brief Gets the transient time data generated during the time-of-flight charge transport test.
		//! \returns A vector of data representing the time points for the transient data.
		std::vector<double> getToFTransientTimes() const;

		//! \brief Gets the transient average polaron velocity data generated during the time-of-flight charge transport test.
		//! The average polaron velocity can be used to calculate the current density and the charge carrier mobility.
		//! \return A vector of data representing how the average polaron velocity in units of cm/s changes over time.
		std::vector<double> getToFTransientVelocities() const;

		//! \brief Gets the transit time data generate during the time-of-flight charge transport test.
		//! \return A vector of data representing the transit time of all extracted polarons.
		std::vector<double> getTransitTimeData() const;

		//! \brief Prints a message to the command line about the current status of the simulation test.
		void outputStatus();

		//! \brief Regenerates the site energies for all sites in the lattice.
		void reassignSiteEnergies();

	protected:

	private:

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
			// precalculated distances vector that contains the distances to nearby sites used for event execution time calculations
			std::vector<double> distances;
			// precalculated isInDissRange and isInFRETRange vectors that contains booleans to indicate whether the nearby sites are within range for the different exciton events to be possible.
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
				// precalculated distances vector that contains the distances to nearby sites used for event execution time calculations
				distances.assign(dim*dim*dim, 0.0);
				// precalculated isInDissRange and isInFRETRange vectors that contains booleans to indicate whether the nearby sites are within range for the different exciton events to be possible.
				isInDissRange.assign(dim*dim*dim, false);
				isInFRETRange.assign(dim*dim*dim, false);
				// Initialize distances, isInDissRange, and isInFRETRange vectors
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
			// precalculated distances vector that contains the distances to nearby sites used for event execution time calculations
			std::vector<double> distances;
			std::vector<double> E_deltas;
			// precalculated isInRange vector that contains booleans to indicate if the nearby sites are within range for polaron events to be possible.
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
				// precalculated distances vector that contains the distances to nearby sites used for event execution time calculations
				distances.assign(dim*dim*dim, 0.0);
				E_deltas.assign(dim*dim*dim, 0.0);
				// precalculated isInRange vector that contains booleans to indicate if the nearby sites are within range for polaron events to be possible.
				isInRange.assign(dim*dim*dim, false);
				// Initialize distances and isInRange vectors
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
		double Steady_equilibration_energy_sum_Coulomb = 0.0;
		double Steady_equilibration_time = 0.0;
		double Steady_Fermi_energy = 0.0;
		double Transport_energy_weighted_sum = 0.0;
		double Transport_energy_weighted_sum_Coulomb = 0.0;
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
		void calculateDOSCorrelation();
		void calculateDOSCorrelation(const double cutoff_radius);
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
		void generateExciton(const KMC_Lattice::Coords& coords, const bool spin, int tag = 0);
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

#endif // EXCIMONTEC_OSC_SIM_H
