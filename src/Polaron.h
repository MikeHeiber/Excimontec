// Copyright (c) 2018 Michael C. Heiber
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

	class Polaron : public Object {
	public:
		static const std::string object_type;
		Polaron(const double time, const int tag_num, const Coords& start_coords, const bool polaron_charge) : Object(time, tag_num, start_coords) { charge = polaron_charge; }
		bool getCharge() const { return charge; }
		std::string getObjectType() const { return object_type; }
	private:
		bool charge; // false represents negative charge, true represents positive charge
	};

	class Polaron_Hop : public Event {
	public:
		static const std::string event_type;
		Polaron_Hop() : Event() {}
		Polaron_Hop(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
			// Calculates hopping using the Miller-Abrahams model
			rate_constant = prefactor * exp(-2.0*localization*distance);
			if (E_delta > 0) {
				rate_constant *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
		}
		void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta, const double reorganization) {
			// Calculates hopping using the Marcus model
			rate_constant = (prefactor / sqrt(4.0*Utils::Pi*reorganization*Utils::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-Utils::intpow(reorganization + E_delta, 2) / (4.0*reorganization*Utils::K_b*sim_ptr->getTemp()));
		}
		std::string getEventType() const { return event_type; }
	private:

	};

	class Polaron_Recombination : public Event {
	public:
		static const std::string event_type;
		Polaron_Recombination() : Event() {}
		Polaron_Recombination(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
			// Calculates recombination using the Miller-Abrahams model
			rate_constant = prefactor * exp(-2.0*localization*distance);
			if (E_delta > 0) {
				rate_constant *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
		}
		void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta, const double reorganization) {
			// Calculates recombination using the Marcus model
			rate_constant = (prefactor / sqrt(4.0*Utils::Pi*reorganization*Utils::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-Utils::intpow(reorganization + E_delta, 2) / (4.0*reorganization*Utils::K_b*sim_ptr->getTemp()));
		}
		std::string getEventType() const { return event_type; }
	private:

	};

	class Polaron_Extraction : public Event {
	public:
		static const std::string event_type;
		Polaron_Extraction() : Event() {}
		Polaron_Extraction(Simulation* simulation_ptr) : Event(simulation_ptr) {}
		void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
			// Calculates extraction using the Miller-Abrahams model
			rate_constant = prefactor * exp(-2.0*localization*distance);
			if (E_delta > 0) {
				rate_constant *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
		}
		std::string getEventType() const { return event_type; }
	private:
	};

}

#endif // POLARON_H
