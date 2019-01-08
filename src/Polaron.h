// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef POLARON_H
#define POLARON_H

#include "Utils.h"
#include "Object.h"
#include "Event.h"
#include <string>

namespace Excimontec {

	class Polaron : public KMC_Lattice::Object {
	public:

		class Hop : public KMC_Lattice::Event {
		public:

			static const std::string event_type;
			
			Hop() : KMC_Lattice::Event() {}
			
			Hop(KMC_Lattice::Simulation* simulation_ptr) : KMC_Lattice::Event(simulation_ptr) {}
			
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
				// Calculates hopping using the Miller-Abrahams model
				rate_constant = prefactor * exp(-2.0*localization*distance);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta, const double reorganization) {
				// Calculates hopping using the Marcus model
				rate_constant = (prefactor / sqrt(4.0*KMC_Lattice::Pi*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-KMC_Lattice::intpow(reorganization + E_delta, 2) / (4.0*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()));
			}
			std::string getEventType() const { return event_type; }
		private:

		};

		class Recombination : public KMC_Lattice::Event {
		public:

			static const std::string event_type;
			
			Recombination() : KMC_Lattice::Event() {}
			
			Recombination(KMC_Lattice::Simulation* simulation_ptr) : KMC_Lattice::Event(simulation_ptr) {}
			
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
				// Calculates recombination using the Miller-Abrahams model
				rate_constant = prefactor * exp(-2.0*localization*distance);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta, const double reorganization) {
				// Calculates recombination using the Marcus model
				rate_constant = (prefactor / sqrt(4.0*KMC_Lattice::Pi*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-KMC_Lattice::intpow(reorganization + E_delta, 2) / (4.0*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()));
			}
			std::string getEventType() const { return event_type; }

		private:

		};

		class Extraction : public KMC_Lattice::Event {
		public:
			
			static const std::string event_type;
			
			Extraction() : KMC_Lattice::Event() {}
			
			Extraction(KMC_Lattice::Simulation* simulation_ptr) : KMC_Lattice::Event(simulation_ptr) {}
			
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
				// Calculates extraction using the Miller-Abrahams model
				rate_constant = prefactor * exp(-2.0*localization*distance);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}
			
			std::string getEventType() const { return event_type; }

		private:
		};

		static const std::string object_type;

		Polaron(const double time, const int tag_num, const KMC_Lattice::Coords& start_coords, const bool polaron_charge) : KMC_Lattice::Object(time, tag_num, start_coords) { charge = polaron_charge; }
		
		bool getCharge() const { return charge; }
		
		std::string getObjectType() const { return object_type; }

	private:

		bool charge; // false represents negative charge, true represents positive charge

	};

}

#endif // POLARON_H
