// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef EXCIMONTEC_POLARON_H
#define EXCIMONTEC_POLARON_H

#include "Utils.h"
#include "Object.h"
#include "Event.h"
#include <string>

namespace Excimontec {

	//! \brief This class extends the KMC_Lattice::Object class to create a polaron object that represents an electron or hole in an organic semiconductor.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2017-2019
	class Polaron : public KMC_Lattice::Object {
	public:

		//! \brief This class extends the KMC_Lattice::Event class to create a polaron hop event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Hop : public KMC_Lattice::Event {
		public:

			//! This static member variable holds the name of the event type, which is "Polaron_Hop".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Hop() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Hop(KMC_Lattice::Simulation* simulation_ptr) : KMC_Lattice::Event(simulation_ptr) {}

			//! \brief Calculates the rate constant for the polaron hop event using the Miller-Abraham hopping model.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the polaron localization property of the starting site.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
				rate_constant = prefactor * exp(-2.0*localization*distance);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}

			//! \brief Calculates the rate constant for the polaron hop event using the Marcus hopping model.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the polaron localization property of the starting site.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			//! \param reorganization is the reorganization energy of the transition.
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta, const double reorganization) {
				rate_constant = (prefactor / sqrt(4.0*KMC_Lattice::Pi*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-KMC_Lattice::intpow(reorganization + E_delta, 2) / (4.0*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()));
			}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Polaron_Hop".
			std::string getEventType() const { return event_type; }

		private:

		};

		//! \brief This class extends the KMC_Lattice::Event class to create a polaron recombination event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Recombination : public KMC_Lattice::Event {
		public:

			//! This static member variable holds the name of the event type, which is "Polaron_Recombination".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Recombination() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Recombination(KMC_Lattice::Simulation* simulation_ptr) : KMC_Lattice::Event(simulation_ptr) {}

			//! \brief Calculates the rate constant for the polaron recombination event using the Miller-Abraham hopping model.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the polaron localization property of the starting site.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
				rate_constant = prefactor * exp(-2.0*localization*distance);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}

			//! \brief Calculates the rate constant for the polaron recombination event using the Marcus hopping model.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the polaron localization property of the starting site.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			//! \param reorganization is the reorganization energy of the transition.
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta, const double reorganization) {
				rate_constant = (prefactor / sqrt(4.0*KMC_Lattice::Pi*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-KMC_Lattice::intpow(reorganization + E_delta, 2) / (4.0*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()));
			}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Polaron_Recombination".
			std::string getEventType() const { return event_type; }

		private:

		};

		//! \brief This class extends the KMC_Lattice::Event class to create a polaron extraction event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Extraction : public KMC_Lattice::Event {
		public:

			//! This static member variable holds the name of the event type, which is "Polaron_Extraction".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Extraction() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Extraction(KMC_Lattice::Simulation* simulation_ptr) : KMC_Lattice::Event(simulation_ptr) {}

			//! \brief Calculates the rate constant for the polaron extraction event using the Miller-Abraham hopping model.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the polaron localization property of the starting site.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
				rate_constant = prefactor * exp(-2.0*localization*distance);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Polaron_Extraction".
			std::string getEventType() const { return event_type; }

		private:
		};

		//! This static member variable holds the name of the object type, which is "Polaron".
		static const std::string object_type;

		//! \brief Constructor that creates and initializes a polaron.
		//! \param time is the simulation time denoting when the polaron was created.
		//! \param tag_num is a unique id number used to distinguish the polaron from other polarons.
		//! \param coords_start is the starting coordinates of the polaron.
		//! \param polaron_charge is the charge state of the polaron. (true for positive and false for negative) 
		Polaron(const double time, const int tag_num, const KMC_Lattice::Coords& coords_start, const bool polaron_charge) : KMC_Lattice::Object(time, tag_num, coords_start) {
			charge = polaron_charge;
		}

		//! \brief Gets the charge state of the polaron.
		//! \returns true if the polaron is a positively charged hole.
		//! \returns false if the polaron is a negatively charged electron.
		bool getCharge() const { return charge; }

		//! \brief Gets the object type string that denotes what type of derived object class this is.
		//! \returns The string "Polaron".
		std::string getObjectType() const { return object_type; }

	private:

		bool charge; // false represents negative charge, true represents positive charge

	};

}

#endif // EXCIMONTEC_POLARON_H
