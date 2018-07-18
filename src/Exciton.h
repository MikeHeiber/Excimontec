// Copyright (c) 2018 Michael C. Heiber
// This source file is part of the Excimontec project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The Excimontec project can be found on Github at https://github.com/MikeHeiber/Excimontec

#ifndef EXCITON_H
#define EXCITON_H

#include "Utils.h"
#include "Object.h"
#include "Event.h"
#include "Simulation.h"
#include <string>

namespace Excimontec {
	//! \brief This class extends the Object class to create an exciton object to represent a singlet or triplet exciton in an organic semiconductor.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2018
	class Exciton : public Object {
	public:
		//! This static member variable holds the name of the object, which is "Exciton".
		static const std::string object_type;

		//! \brief Constructor that creates and initializes an exciton.
		//! \param time is the simulation time denoting when the exciton was created.
		//! \param tag_num is a unique id number used to distinguish the exciton from other excitons.
		//! \param coords_start is the Coords struct that represents the starting coordinates of the exciton.
		Exciton(const double time, const int tag_num, const Coords& coords_start) : Object(time, tag_num, coords_start) {}

		//! \brief Flips the spin state of the exciton from singlet to triplet or from triplet to singlet.
		void flipSpin() { spin_state = !spin_state; }

		//! \brief Gets the object type string that denotes what type of Object class this is.
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

	//! \brief This class extends the Event class to create an specific type of exciton event.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2018
	class Exciton_Creation : public Event {
	public:
		//! This static member variable holds the name of the event, which is "Exciton_Creation".
		static const std::string event_type;

		//! \brief Constructs an empty event that is uninitialized.
		Exciton_Creation() : Event() {}

		//! \brief Constructs and initializes an event.
		//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
		Exciton_Creation(Simulation* simulation_ptr) : Event(simulation_ptr) {}

		//! \brief Gets the event type string that denotes what type of Event class this is.
		//! \returns The string "Exciton_Creation".
		std::string getEventType() const { return event_type; }

	private:
	};

	//! \brief This class extends the Event class to create an specific type of exciton event.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2018
	class Exciton_Hop : public Event {
	public:
		//! This static member variable holds the name of the event, which is "Exciton_Hop".
		static const std::string event_type;

		//! \brief Constructs an empty event that is uninitialized.
		Exciton_Hop() : Event() {}

		//! \brief Constructs and initializes an event.
		//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
		Exciton_Hop(Simulation* simulation_ptr) : Event(simulation_ptr) {}

		//! \brief Calculates the rate constant for the exciton hop event using the FRET hopping mechanism.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param distance is the distance between the starting site and destination site.
		//! \param E_delta is the potential energy change that would occur if the event is executed.
		void calculateRateConstant(const double prefactor, const double distance, const double E_delta) {
			rate_constant = prefactor * Utils::intpow(1.0 / distance, 6);
			if (E_delta > 0) {
				rate_constant *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
		}

		//! \brief Calculates the rate constant for the exciton hop event using the Dexter hopping mechanism.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param localization is the inverse localization parameter that describes how localized the exciton is.
		//! \param distance is the distance between the starting site and destination site.
		//! \param E_delta is the potential energy change that would occur if the event is executed.
		void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
			rate_constant = prefactor * exp(-2.0*localization*distance);
			if (E_delta > 0) {
				rate_constant *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
		}

		//! \brief Gets the event type string that denotes what type of Event class this is.
		//! \returns The string "Exciton_Hop".
		std::string getEventType() const { return event_type; }

	private:
	};

	//! \brief This class extends the Event class to create an specific type of exciton event.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2018
	class Exciton_Recombination : public Event {
	public:
		//! This static member variable holds the name of the event, which is "Exciton_Recombination".
		static const std::string event_type;

		//! \brief Constructs an empty event that is uninitialized.
		Exciton_Recombination() : Event() {}

		//! \brief Constructs and initializes an event.
		//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
		Exciton_Recombination(Simulation* simulation_ptr) : Event(simulation_ptr) {}

		//! \brief Gets the event type string that denotes what type of Event class this is.
		//! \returns The string "Exciton_Recombination".
		std::string getEventType() const { return event_type; }

	private:
	};

	//! \brief This class extends the Event class to create an specific type of exciton event.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2018
	class Exciton_Dissociation : public Event {
	public:
		//! This static member variable holds the name of the event, which is "Exciton_Dissociation".
		static const std::string event_type;

		//! \brief Constructs an empty event that is uninitialized.
		Exciton_Dissociation() : Event() {}

		//! \brief Constructs and initializes an event.
		//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
		Exciton_Dissociation(Simulation* simulation_ptr) : Event(simulation_ptr) {}

		//! \brief Calculates the rate constant for the exciton dissociation event using the Miller-Abrahams polaron hopping mechanism.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param localization is the inverse localization parameter that describes how localized the exciton is.
		//! \param distance is the distance between the starting site and destination site.
		//! \param E_delta is the potential energy change that would occur if the event is executed.
		void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta) {
			rate_constant = prefactor * exp(-2.0*localization*distance);
			if (E_delta > 0) {
				rate_constant *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
		}

		//! \brief Calculates the rate constant for the exciton dissociation event using the Marcus polaron hopping mechanism.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param localization is the inverse localization parameter that describes how localized the exciton is.
		//! \param distance is the distance between the starting site and destination site.
		//! \param E_delta is the potential energy change that would occur if the event is executed.
		//! \param reorganization is the reorganization energy for the Marcus electron transfer mechanism.
		void calculateRateConstant(const double prefactor, const double localization, const double distance, const double E_delta, const double reorganization) {
			rate_constant = (prefactor / sqrt(4.0*Utils::Pi*reorganization*Utils::K_b*sim_ptr->getTemp()))*exp(-2.0*localization*distance)*exp(-Utils::intpow(reorganization + E_delta, 2) / (4.0*reorganization*Utils::K_b*sim_ptr->getTemp()));
		}

		//! \brief Gets the event type string that denotes what type of Event class this is.
		//! \returns The string "Exciton_Dissociation".
		std::string getEventType() const { return event_type; }

	private:
	};

	//! \brief This class extends the Event class to create an specific type of exciton event.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2018
	class Exciton_Intersystem_Crossing : public Event {
	public:
		//! This static member variable holds the name of the event, which is "Exciton_Intersystem_Crossing".
		static const std::string event_type;

		//! \brief Constructs an empty event that is uninitialized.
		Exciton_Intersystem_Crossing() : Event() {}

		//! \brief Constructs and initializes an event.
		//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
		Exciton_Intersystem_Crossing(Simulation* simulation_ptr) : Event(simulation_ptr) {}

		//! \brief Calculates the rate constant for the exciton intersystem crossing event.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param E_delta is the potential energy change that would occur if the event is executed.
		void calculateRateConstant(const double prefactor, const double E_delta) {
			rate_constant = prefactor;
			if (E_delta > 0) {
				rate_constant *= exp(-E_delta / (Utils::K_b*sim_ptr->getTemp()));
			}
		}

		//! \brief Gets the event type string that denotes what type of Event class this is.
		//! \returns The string "Exciton_Intersystem_Crossing".
		std::string getEventType() const { return event_type; }

	private:
	};

	//! \brief This class extends the Event class to create an specific type of exciton event.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2018
	class Exciton_Exciton_Annihilation : public Event {
	public:
		//! This static member variable holds the name of the event, which is "Exciton_Exciton_Annihilation".
		static const std::string event_type;

		//! \brief Constructs an empty event that is uninitialized.
		Exciton_Exciton_Annihilation() : Event() {}

		//! \brief Constructs and initializes an event.
		//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
		Exciton_Exciton_Annihilation(Simulation* simulation_ptr) : Event(simulation_ptr) {}

		//! \brief Calculates the rate constant for the exciton-exciton annihilation event using a FRET hopping mechanism.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param distance is the distance between the starting site and destination site.
		void calculateRateConstant(const double prefactor, const double distance) {
			rate_constant = prefactor * Utils::intpow(1.0 / distance, 6);
		}

		//! \brief Calculates the rate constant for the exciton-exciton annihilation event using the Dexter hopping mechanism.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param localization is the inverse localization parameter that describes how localized the exciton is.
		//! \param distance is the distance between the starting site and destination site.
		void calculateRateConstant(const double prefactor, const double localization, const double distance) {
			rate_constant = prefactor * exp(-2.0*localization*distance);
		}

		//! \brief Gets the event type string that denotes what type of Event class this is.
		//! \returns The string "Exciton_Exciton_Annihilation".
		std::string getEventType() const { return event_type; }

	private:
	};

	//! \brief This class extends the Event class to create an specific type of exciton event.
	//! \copyright MIT License.  For more information, see the LICENSE file that accompanies this software package.
	//! \author Michael C. Heiber
	//! \date 2018
	class Exciton_Polaron_Annihilation : public Event {
	public:
		//! This static member variable holds the name of the event, which is "Exciton_Polaron_Annihilation".
		static const std::string event_type;

		//! \brief Constructs an empty event that is uninitialized.
		Exciton_Polaron_Annihilation() : Event() {}

		//! \brief Constructs and initializes an event.
		//! \param simulation_ptr is a pointer to the Simulation object that is associated with the event.
		Exciton_Polaron_Annihilation(Simulation* simulation_ptr) : Event(simulation_ptr) {}

		//! \brief Calculates the rate constant for the exciton-polaron annihilation event using a FRET hopping mechanism.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param distance is the distance between the starting site and destination site.
		void calculateRateConstant(const double prefactor, const double distance) {
			rate_constant = prefactor * Utils::intpow(1.0 / distance, 6);
		}

		//! \brief Calculates the rate constant for the exciton-polaron annihilation event using the Dexter hopping mechanism.
		//! \param prefactor is the rate constant prefactor for the transition.
		//! \param localization is the inverse localization parameter that describes how localized the exciton is.
		//! \param distance is the distance between the starting site and destination site.
		void calculateRateConstant(const double prefactor, const double localization, const double distance) {
			rate_constant = prefactor * exp(-2.0*localization*distance);
		}

		//! \brief Gets the event type string that denotes what type of Event class this is.
		//! \returns The string "Exciton_Polaron_Annihilation".
		std::string getEventType() const { return event_type; }

	private:
	};

}

#endif // EXCITON_H
