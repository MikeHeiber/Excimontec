// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef EXCIMONTEC_EXCITON_H
#define EXCIMONTEC_EXCITON_H

#include "Utils.h"
#include "Object.h"
#include "Event.h"
#include "Simulation.h"
#include <string>

namespace Excimontec {

	//! \brief This class extends the KMC_Lattice::Object class to create an exciton object that represents a singlet or triplet exciton in an organic semiconductor.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2017-2019
	class Exciton : public KMC_Lattice::Object {
	public:

		//! \brief This class extends the KMC_Lattice::Event class to create an exciton creation event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Creation : public KMC_Lattice::Event {
		public:
			//! This static member variable holds the name of the event type, which is "Exciton_Creation".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Creation() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Creation(KMC_Lattice::Simulation* simulation_ptr) : KMC_Lattice::Event(simulation_ptr) {}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Exciton_Creation".
			std::string getEventType() const { return event_type; }

		private:
		};

		//! \brief This class extends the KMC_Lattice::Event class to create an exciton hop event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Hop : public KMC_Lattice::Event {
		public:
			//! This static member variable holds the name of the event type, which is "Exciton_Hop".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Hop() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Hop(KMC_Lattice::Simulation* simulation_ptr) : Event(simulation_ptr) {}

			//! \brief Calculates the rate constant for the exciton hop event using the FRET hopping mechanism.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			void calculateRateConstant(const double prefactor, const double distance, const double E_delta) {
				rate_constant = prefactor * KMC_Lattice::intpow(1.0 / distance, 6);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}

			//! \brief Calculates the rate constant for the exciton hop event using the Dexter electron exchange hopping mechanism.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the inverse localization parameter that describes how localized the exciton is.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
				rate_constant = prefactor * exp(-2.0*localization*distance);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Exciton_Hop".
			std::string getEventType() const { return event_type; }

		private:
		};

		//! \brief This class extends the KMC_Lattice::Event class to create an exciton recombination event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Recombination : public KMC_Lattice::Event {
		public:
			//! This static member variable holds the name of the event type, which is "Exciton_Recombination".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Recombination() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Recombination(KMC_Lattice::Simulation* simulation_ptr) : Event(simulation_ptr) {}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Exciton_Recombination".
			std::string getEventType() const { return event_type; }

		private:
		};

		//! \brief This class extends the KMC_Lattice::Event class to create an exciton dissociation event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Dissociation : public KMC_Lattice::Event {
		public:
			//! This static member variable holds the name of the event type, which is "Exciton_Dissociation".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Dissociation() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Dissociation(KMC_Lattice::Simulation* simulation_ptr) : Event(simulation_ptr) {}

			//! \brief Calculates the rate constant for the exciton dissociation event using the Miller-Abrahams polaron hopping mechanism.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the inverse localization parameter that describes how localized the exciton is.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
				rate_constant = prefactor * exp(-2.0*localization*distance);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}

			//! \brief Calculates the rate constant for the exciton dissociation event using the Marcus polaron hopping mechanism.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the inverse localization parameter that describes how localized the exciton is.
			//! \param distance is the distance between the starting site and destination site.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			//! \param reorganization is the reorganization energy for the Marcus electron transfer mechanism.
			void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta, const double reorganization) {
				rate_constant = (prefactor / sqrt(4.0*KMC_Lattice::Pi*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-KMC_Lattice::intpow(reorganization + E_delta, 2) / (4.0*reorganization*KMC_Lattice::K_b*sim_ptr->getTemp()));
			}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Exciton_Dissociation".
			std::string getEventType() const { return event_type; }

		private:
		};

		//! \brief This class extends the KMC_Lattice::Event class to create an exciton intersystem crossing event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Intersystem_Crossing : public KMC_Lattice::Event {
		public:
			//! This static member variable holds the name of the event type, which is "Exciton_Intersystem_Crossing".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Intersystem_Crossing() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Intersystem_Crossing(KMC_Lattice::Simulation* simulation_ptr) : Event(simulation_ptr) {}

			//! \brief Calculates the rate constant for the exciton intersystem crossing event.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param E_delta is the potential energy change that would occur if the event is executed.
			void calculateRateConstant(const double prefactor, const double E_delta) {
				rate_constant = prefactor;
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Exciton_Intersystem_Crossing".
			std::string getEventType() const { return event_type; }

		private:
		};

		//! \brief This class extends the KMC_Lattice::Event class to create an exciton-exciton annihilation event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Exciton_Annihilation : public KMC_Lattice::Event {
		public:
			//! This static member variable holds the name of the event type, which is "Exciton_Exciton_Annihilation".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Exciton_Annihilation() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Exciton_Annihilation(KMC_Lattice::Simulation* simulation_ptr) : Event(simulation_ptr) {}

			//! \brief Calculates the rate constant for the exciton-exciton annihilation event using a FRET hopping mechanism.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param distance is the distance between the starting site and destination site.
			void calculateRateConstant(const double prefactor, const double distance) {
				rate_constant = prefactor * KMC_Lattice::intpow(1.0 / distance, 6);
			}

			//! \brief Calculates the rate constant for the exciton-exciton annihilation event using the Dexter electron exchange hopping mechanism.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the inverse localization parameter that describes how localized the exciton is.
			//! \param distance is the distance between the starting site and destination site.
			void calculateRateConstant(const double prefactor, const double localization, const double distance) {
				rate_constant = prefactor * exp(-2.0*localization*distance);
			}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Exciton_Exciton_Annihilation".
			std::string getEventType() const { return event_type; }

		private:
		};

		//! \brief This class extends the KMC_Lattice::Event class to create an exciton-polaron annihilation event.
		//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
		//! \author Michael C. Heiber
		//! \date 2017-2019
		class Polaron_Annihilation : public KMC_Lattice::Event {
		public:
			//! This static member variable holds the name of the event type, which is "Exciton_Polaron_Annihilation".
			static const std::string event_type;

			//! \brief Constructs an empty event that is uninitialized.
			Polaron_Annihilation() : KMC_Lattice::Event() {}

			//! \brief Constructs and initializes an event.
			//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
			Polaron_Annihilation(KMC_Lattice::Simulation* simulation_ptr) : KMC_Lattice::Event(simulation_ptr) {}

			//! \brief Calculates the rate constant for the exciton-polaron annihilation event using a FRET hopping mechanism.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param distance is the distance between the starting site and destination site.
			void calculateRateConstant(const double prefactor, const double distance) {
				rate_constant = prefactor * KMC_Lattice::intpow(1.0 / distance, 6);
			}

			//! \brief Calculates the rate constant for the exciton-polaron annihilation event using the Dexter electron exchange hopping mechanism.
			//! \param prefactor is the rate constant prefactor for the transition.
			//! \param localization is the inverse localization parameter that describes how localized the exciton is.
			//! \param distance is the distance between the starting site and destination site.
			void calculateRateConstant(const double prefactor, const double localization, const double distance) {
				rate_constant = prefactor * exp(-2.0*localization*distance);
			}

			//! \brief Gets the event type string that denotes what type of derived event class this is.
			//! \returns The string "Exciton_Polaron_Annihilation".
			std::string getEventType() const { return event_type; }

		private:
		};

		//! This static member variable holds the name of the object type, which is "Exciton".
		static const std::string object_type;

		//! \brief Constructor that creates and initializes an exciton.
		//! \param time is the simulation time denoting when the exciton was created.
		//! \param tag_num is a unique id number used to distinguish the exciton from other excitons.
		//! \param coords_start is the starting coordinates of the exciton.
		//! \param exciton_spin is the spin state of the exciton. (true for singlet and false for triplet)
		Exciton(const double time, const int tag_num, const KMC_Lattice::Coords& coords_start, const bool exciton_spin) : KMC_Lattice::Object(time, tag_num, coords_start) {
			spin_state = exciton_spin;
		}

		//! \brief Flips the spin state of the exciton from singlet to triplet or from triplet to singlet.
		void flipSpin() { spin_state = !spin_state; }

		//! \brief Gets the object type string that denotes what type of derived object class this is.
		//! \returns The string "Exciton".
		std::string getObjectType() const { return object_type; }

		//! \brief Gets the current spin state of the exciton.
		//! \returns true if the exciton is in a singlet state.
		//! \returns false if the exciton is in a triplet state.
		bool getSpin() const { return spin_state; }

		//! \brief Sets the spin state of the exciton.
		//! \param spin_state_new indicates what the spin state will be set to.  True for singlet and false for triplet.
		void setSpin(bool spin_state_new) { spin_state = spin_state_new; }

	private:
		bool spin_state; // false represents triplet state, true represents singlet state

	};

}

#endif // EXCIMONTEC_EXCITON_H
