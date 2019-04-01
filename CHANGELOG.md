# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0.0-rc.3]- 2019-04-01 - Density of States Integration Bugfix

### Added
- OSC_Sim - Steady_DOS_sampling_count and Steady_DOOS_sampling_counter to keep track of how many times the DOS and DOOS are sampled during the simulation
- test.cpp (SteadyTransportTests) - tests to check the integral of the DOS and DOOS with and without including Coulomb interactions

### Changed
- OSC_Sim (updateSteadyData) - To increment the new DOS and DOOS sampling counters and removed calculation of the DOS because the DOS is does not change during the simulation when not including Coulomb interactions
- OSC_Sim (getSteadyDOS) - Function is not longer const

### Removed

### Fixed
- main.cpp - Output of density of states and density of occupied states data to have the correct units
- OSC_Sim (getSteadyDOOS) - To use the new sampling counter to avoid errors averaging over multiple samplings
- OSC_Sim (getSteadyDOOS_Coulomb) - To use the new sampling counter to avoid errors averaging over multiple samplings
- OSC_Sim (getSteadyDOS) - To calculate the DOS using only one sample at call time because the DOS (without Coulomb) does not change during the simulation
- OSC_Sim (getSteadyDOS_Coulomb) - To use the new sampling counter to avoid errors averaging over multiple samplings


## [v1.0.0-rc.2]- 2019-02-09 - Steady State Charge Transport Test Update

### Added
- .gitignore - Numerous ignore statements to ignore files generated during build and test operations and from Microsoft Visual Studio
- CHANGELOG.md - Notes about all changes in this release
- docs - KMC_Lattice documentation
- docs - Documentation of new functions in the OSC_Sim class
- Exciton - Project name (Excimontec) to header guards
- main.cpp - Output of equilibration energy and transport energy both with and without including the Coulomb potential when running a steady transport test
- main.cpp - Output the current density when running a steady transport test
- main.cpp - Output the density of states and density of occupied states data calculated with and without Coulomb potential adjustments
- OSC_Sim - Project name (Excimontec) to header guards
- OSC_Sim - Public function and member variable documentation
- OSC_Sim (getSiteEnergy) - Checking of the input coordinates validity and generating error if invalid
- OSC_Sim (getSiteType) - Checking of the input coordinates validity and generating error if invalid
- OSC_Sim - include statements for all used components to better document class dependencies
- OSC_Sim (getSteadyCurrentDensity) - New function that returns the average current density from the steady state transport test
- OSC_Sim (getSteadyEquilibrationEnergy_Coulomb) - New function to return the equilibration energy calculated including the Coulomb potential
- OSC_Sim (getSteadyTransportEnergy_Coulomb) - New function to return the transport energy calculated including the Coulomb potential
- OSC_Sim (Steady_equilibration_energy_sum_Coulomb) - New private member variable for calculating the equilibration energy including the Coulomb potential
- OSC_Sim (Transport_energy_weighted_sum_Coulomb) - New private member variable for storing the data used to calculate the transport energy including the Coulomb potential
- OSC_Sim (Steady_hops_per_DOS_sample) - New private member variable for calculating the DOS during the steady transport test
- OSC_Sim (Steady_hops_per_DOOS_sample) - New private member variable for calculating the DOOS during the steady transport test
- OSC_Sim (DOS_bin_size) - New private member variable for calculating the DOOS and DOS during the steady transport test
- OSC_Sim (executePolaronHop) - Calculation of the transport energy including the Coulomb potential
- OSC_Sim (exportEnergies) - New overloaded function allowing the user to output the absolute site energies for electrons or holes
- OSC_Sim (generateSteadyPolarons) - Simpler creation of polarons on random sites when energetic disorder is disabled
- OSC_Sim (generateSteadyPolarons) - Status output about how many polarons are created in the lattice for the steady transport test
- OSC_Sim (generateSteadyPolarons) - Allow creation of polarons on acceptor sites
- OSC_Sim (getSteadyDOOS) - New function for getting density of occupied states data
- OSC_Sim (getSteadyDOOS_Coulomb) - New function for getting density of occupied states data
- OSC_Sim (getSteadyDOS) - New function for getting density of states data
- OSC_Sim (getSteadyDOS_Coulomb) - New function for getting density of states data
- OSC_Sim (getSteadyTransportEnergy) - Function now allows calculation when the sum of weights is negative
- OSC_Sim (getSteadyTransportEnergy_Coulomb) - Function now allows calculation when the sum of weights is negative
- OSC_Sim (outputStatus) - Status output for the steady state charge transport test
- OSC_Sim (updateSteadyData) - Calculation of the equilibration energy including the Coulomb potential
- OSC_Sim (updateSteadyData) - Periodic sampling of the DOOS and DOS
- OSC_Sim (updateSteadyDOS) - New function that uses the input energy of a given site to update the input density of states data
- Parameters - Project name (Excimontec) to header guards
- Parameters - Public function and member variable documentation
- Parameters (Enable_steady_data_output) - New member variable that is set to true by default but could be used in the future to disable DOS data output
- Polaron - Project name (Excimontec) to header guards
- Polaron - Public function and member variable documentation
- README.md - Build instructions link for Windows users
- README.md - Information about new DOS and DOOS data files generated during the steady transport test
- README.md - Examples about what the custom site energies import feature can be used for.
- test.cpp - Added a simple command line status message at the beginning of all test cases
- test.cpp (SteadyTransportTests) - Tests checking the output of the transport and equilibration energies when the steady transport test has not been run
- test.cpp (SteadyTransportTests) - Tests comparing the energies calculated with and without Coulomb interactions and relative positions of the equilibration and transport energies
- test.cpp - (SteadyTransportTests) - Test to check that phase restriction disabling increases the number of available sites for creating the initial polarons in donor-acceptor blends
- test.cpp - (SteadyTransportTests) - Tests to check the peak position of the DOS and DOOS data from the very low field test
- test.cpp (SteadyTransportTests) - Tests to check the transport energy calculated during a medium electric field test condition
- test.cpp (SteadyTransportTests) - Test to check the relative position of the transport energy and the donor HOMO during the medium electric field test
- test.cpp (SteadyTransportTests) - Test to check the absolute position of the transport energy
- test.cpp (SteadyTransportTests) - Test for transport in a random donor-acceptor blend and check for the relative change in transport energy position relative to the neat donor architecture
- test.cpp (SteadyTransportTests) - Test checking that the simulation works when phase restriction in enabled with a random blend
- test.cpp (SteadyTransportTests) - Test checking the magnitude of the current density
- test.cpp (ToFTests) - Test checking the hole extraction map output

### Changed
- KMC_Lattice - KMC_Lattice submodule to v2.1.0-beta.1
- Many files - Copyright statement years to 2017-2019
- docs - Updated docs using Doxygen v1.8.15
- Doxyfile - Project version number to v1.0.0
- Doxyfile - Settings so that KMC_Lattice is included and markdown files are no longer be included in the documentation
- Exciton - Nested derived Exciton event classes into the Exciton class
- Exciton (constructor) - Must now specify the spin state upon construction
- Exciton - Revised documentation for public functions and member variables
- main.cpp - Version string to v1.0.0-rc.2 in preparation for next release
- OSC_Sim - Functions to use new object event class nesting format
- OSC_Sim (Site_OSC) - Store the site type internally as a char instead of short to save memory
- OSC_Sim (generateExciton) - Moved definition of default tag value to the header file function declaration statement
- OSC_Sim (generateExciton) - Implemented the new Exciton constructor where one must specify the spin state
- OSC_Sim - Nested derived Site class (Site_OSC) into the OSC_Sim class as a private class
- OSC_Sim (calculateDOSCorrelation) - Scope of the two functions from public to private
- OSC_Sim (executePolaronHop) - Refactored function and using temporary local variables to make code more readable
- OSC_Sim (executePolaronHop) - Calculation of the transport energy to absolute value that includes the HOMO energy and accounts for donor and acceptor site occupation
- OSC_Sim (executePolaronHop) - Calculation of the transport energy is only performed after the equilibration phase is complete
- OSC_Sim (executePolaronHop) - Calculation of the transport energy is is done using the displacement as the weights instead of the velocity
- OSC_Sim (getChargeExtractionMap) - Refactored code replacing usage of stringstream with addition of substrings
- OSC_Sim (updateSteadyData) - Calculation of equilibration energy to absolute value that includes the HOMO energy and accounts for donor and acceptor site occupation
- Parameters (checkParameters) - Lower the limit for the smallest internal potential that one can during a steady transport simulation to allow very low field simulations
- parameters_default.txt - Default morphology file format to not include the compression specifier suffix 
- Polaron - Nested derived Polaron event classes into the Polaron class
- README.md - Replaced version badges with text links
- README.md - Updated description of steady transport test feature
- README.md - Updated current release status info for KMC_Lattice to v2.1
- test.cpp - All tests to route command line output to a test_log.txt file instead of cluttering the command line making it easier to see the test results
- test.cpp (EnergiesImportTests) - Tests of the new exportEnergies function checking the absolute value of the exported electron and hole site energies
- test.cpp (ExcitonDynamicsTests) - Increased the range of the transient to get a more accurate assessment of the equilibrium energy position
- test.cpp (LoggingTests) - Adjusted parameters to promote additional mechanisms to occur including RISC and exciton recombination
- test.cpp (SteadyTransportTests) - Test of the energy values to compare to absolute energy values including the HOMO energy
- test.cpp (SteadyTransportTests) - Very low field test to make it more accurate and a little bit faster by decreasing the internal potential, lattice size, and Coulomb cutoff radius
- test.cpp (SteadyTransportTests) - Medium field test by reducing the number of iterations to make the test faster
- test.cpp (SteadyTransportTests) - Adjusted parameters of the no disorder mobility test to decrease the test time

### Removed
- docs - Markdown files from the generated documentation
- main.cpp - Output of the Fermi energy during the steady state charge transport test
- OSC_Sim (Steady_Fermi_energy) - private member variable that is no longer used
- OSC_Sim (getSteadyFermiEnergy) - Fermi energy is no longer calculated and DOOS and DOS data is output instead
- test.cpp (SteadyTransportTests) - Tests of the Fermi energy

### Fixed
- CHANGELOG.md - Several spelling mistakes in previous release sections
- main.cpp - Spelling mistake in the results output of steady transport test
- main.cpp - Label error in the results output of the time-of-flight charge transport test, current should have been current density
- main.cpp - Bug bug where status output and logfile reset was only being performed when checking the error status
- OSC_Sim (calculateCoulomb) - Specified input parameter namespaces to avoid Doxygen confusion between the header declaration and source file definition
- OSC_Sim (createExciton) - Specified input parameter namespaces to avoid Doxygen confusion between the header declaration and source file definition
- OSC_Sim (calculatePolaronEvents) - Spelling mistake in error message
- OSC_Sim (reassignSiteEnergies) - Spelling mistake in error message
- OSC_Sim (executePolaronHop) - Transport energy calculation to correctly account for hops across the periodic boundary when calculating the displacement
- OSC_Sim (getSteadyTransportEnergy) - Spelling mistake in the function documentation
- OSC_Sim (getPolaronIt) - Bug that could occur when comparing a list iterator to an iterator from a different list
- OSC_Sim (getSteadyEquilibrationEnergy) - Bug that could cause rounding error due to integer division
- Parameters (checkParameters) - Spelling mistakes in error messages
- Parameters (importParameters) - Spelling mistakes in error messages
- parameters_default.txt - Incorrect units for the exciton hopping and annihilation rate prefactors
- test.cpp - Spelling mistakes in the test comments

## [v1.0.0-rc.1]- 2018-12-11 - Interfacial Energy Shift, Site Energy Import, and Steady State Transport Test Update

### Added
- CHANGELOG.md - New file detailing the changes for each release
- README.md - Link to new Changelog file
- README.md - Copyright statement
- slurm_script.sh - Copyright statement
- Parameters - New class to store all parameters used by the simulation and to contain functions for parsing the parameter file and checking for parameter validity
- Parameters - Default values to parameters used in main
- makefile - Commands to compile and link the Parameters class source files
- makefile - Missing dependency of all Excimontec objects on the KMC_Lattice lib
- test.cpp (ParameterTests) - New test for importing the parameters_default.txt file that checks the importParameters function
- test.cpp (ParameterTests) - New series of tests for how the program handles misspelled boolean values in the parameter file
- parameters_default.txt - Three new parameters (Enable_interfacial_energy_shift, Energy_shift_donor, Energy_shift_acceptor) for the interfacial energy shift model
- Parameters - The three new parameters for the interfacial energy shift feature
- Parameters (importParameters) - Code to read the interfacial energy shift parameters from the parameter file
- Parameters (checkParameters) - Validity checks for the interfacial energy shift parameters
- OSC_Sim (reassignSiteEnergies) - Code to implement the interfacial energy shift model
- test.cpp - The three new parameters for the interfacial energy shift model to the default parameters struct
- test.cpp (InterfacialEnergyShiftTests) - New test function checking the energy shift using a bilayer with and without energetic disorder
- test.cpp (ParameterTests) - Tests to check response of OSC_Sim to initialization with invalid interfacial energy shift parameters
- parameters_default.txt - Two new parameters (Enable_import_energies, Energies_import_filename) for importing site energies from a file
- Parameters - The two new parameters for the site energy import feature
- Parameters (importParameters) - Code to read the site energy import parameters from the parameter file
- Parameters (checkParameters) - Validity checks for the new site energies import parameters
- OSC_Sim (exportEnergies) - New function that creates a site energies file
- OSC_Sim (reassignSiteEnergies) - Code section that imports the energies from the specified file when the site energy import feature is enabled
- New site energies files that have improper format to check for handling of invalid site energies files during testing
- test.cpp - The two new site energies import parameters to the default parameters struct
- test.cpp - (ParameterTests) - Tests to check response of OSC_Sim to initialization with invalid site energies import parameter combinations
- test.cpp (EnergiesImportTests) - New test function checking the export and import of valid energies file and to check how the program handles energies file with improper format or missing data
- Feature allowing users to enable event logging for debugging purposes by adding -enable_logging as a command line argument after the parameter file name
- main.cpp (main) - Code to check command line arguments and enable logging if the -enable_logging argument is passed after the parameter file name
- OSC_Sim - New N_events_executed counter and getN_events_executed function to overwrite the one in the base Simulation class to separately keep track oh how many events have been executed
- OSC_Sim - New private member variable previous_event_type to store the type of event that was just executed and a getPreviousEventType function to retrieve the string for testing purposes
- OSC_Sim (calculateNextEvent) - Code that increments the new N_events_executed and store the previous event type
- test.cpp (LoggingTests) - New test function to check that all events being executed are accurately logged and test the getN_events_executed as well
- test.cpp (ChargeDynamicsTests) - New simple test for an exciton dissociation and charge recombination dynamics simulation
- test.cpp (ExcitonDiffusionTests) - New test to check that increased energetic disorder reduces the diffusion length, which also tests activated exciton hopping functionality
- test.cpp (ExcitonDiffusionTests) - New tests checking the exciton recombination event counters and a test for exciton diffusion in an exponential DOS
- test.cpp (IQETests) - New test with a weakly donating/accepting bilayer to test its impact on the charge separation yield and allow polarons to hop to opposing phase
- test.cpp (ObjectCreationTests) - New test function to test the createExciton, createElectron, and createHole functions of the OSC_Sim class with both valid and invalid input coordinates and site types and to test some of the event counters and the checkFinished function
- test.cpp (ParameterTests) - Several previously untested invalid parameter combinations
- test.cpp (ToFTests) - New hole ToF test using the Marcus hopping model
- test.cpp (ToFTests) - New tests to check that increased energetic disorder reduces the mobility, which also tests activated polaron hopping functionality
- test.cpp (ToFTests) - New test that checks the behavior of the ToF_placement_energy feature
- test.cpp (ToFTests) - New tests that get and analyze the transient energy relaxation data and test various aspects of this behavior
- README - New information about the interfacial energy shift model, the site energy import capability, and the steady transport test
- makefile - New PGI compiler flag enabling output of compiler warnings
- parameters_default.txt - New parameters for the steady transport test: Enable_steady_transport_test, Steady_carrier_density, and N_equilibration_events
- main.cpp (main) - Command line output and results file output for the steady transport test
- OSC_Sim - New public functions for the steady transport test: getSteadyEquilibrationEnergy, getSteadyFermiEnergy, getSteadyMobility, and getSteadyTransportEnergy
- OSC_Sim- New private member variables to store data needed to calculate the final steady transport test results: Steady_equilibration_energy_sum, Steady_equilibration_time, Transport_energy_weighted_sum, and Transport_energy_sum_of_weights
- OSC_Sim - New generateSteadyPolarons function for creating initial polarons a the beginning of the steady transport test
- OSC_Sim - New updateSteadyData function to updating the steady transport test data variables each event iteration
- OSC_Sim (init) - Call to the new generateSteadyPolarons function when the steady transport test is enabled
- OSC_Sim (calculatePolaronEvents) - Code to determine when polarons are trying to cross the z periodic boundary and adjust the potential energy change accordingly
- OSC_Sim (calculatePolaronEvents) - Code to disable calculation of polaron extraction events when the steady transport test is enabled
- OSC_Sim (checkFinished) - Code to check for the test termination conditions for the steady transport test
- OSC_Sim (executeNextEvent) - Call to the new updateSteadyData function before executing each event when the steady transport test is enabled
- OSC_Sim (executePolaronHop) - Code to record the transport energy data when the steady transport test is enabled
- OSC_Sim (updateTransientData) - Code to cast the return of the distance function as an int to prevent compiler warnings
- Parameters - New parameters for the steady transport test: Enable_steady_transport_test, Steady_carrier_density, and N_equilibration_events
- Parameters (checkParameters) - Call to the base class checkParameters function
- Parameters (checkParameters) - New checks for the steady transport test parameter combinations
- Parameters (importParameters) - Code to import new steady transport test parameters from the parameter file
- Parameters (importParameters) - Code to check the status of the input ifstream and throw an exception if there is a problem
- test.cpp - Default parameter values for the new parameters for the steady transport test: Enable_steady_transport_test, Steady_carrier_density, and N_equilibration_events
- test.cpp (ParameterTests) - New tests checking the validity of the parameters when enabling the steady transport test
- test.cpp - New test function SteadyTransportTests that check the output of the steady transport test
- test.cpp (ParameterTests) - New tests for attempting to import parameters using an uninitialized or close ifstream
- OSC_Sim - New overloaded createExciton function that does not specify the creation coordinates and generates them randomly
- OSC_sim (createExciton(Coords)) - Code to check whether the input coordinates are occupied and producing an error if they are
- test.cpp (ChargeDynamicsTests) - New tests to check the transient data
- test.cpp (ExcitonDynamicsTests) - New tests to check exciton relaxation into the tail of a Gaussian DOS, triplet dissociation, and triplet annihilation mechanisms
- test.cpp (IQETests) - New tests for changes in the number of geminate recombination events under different test conditions
- test.cpp (ObjectCreationTests) - New tests checking attempts to create an exciton in a fully occupied lattice
- test.cpp (ParameterTests) - New tests for invalid parameters in the Parameters_Lattice and Parameters_Simulations classes
- test.cpp (SteadyTransportTests) - New tests for attempts to create more initial polarons that there are donor sites
- test.cpp (ToFTests) - New tests for attempts to create more initial polarons that there are donor sites
- test.cpp (ToFTests) - New tests comparing the relaxed mobility and relaxed occupation energy from the transient data against the expected steady state mobility and steady state equilibration energy
- CONTRIBUTING.md - New file with detailed instructions for how new people can contribute to the project
- README.md - Link to the CONTRIBUTING.md file

### Changed
- main.cpp (main) - Refactored code to use the new Parameters class
- OSC_Sim - Refactored all functions to use to use the new Parameters class member
- test.cpp - Refactored all tests to use the new Parameters class
- test.cpp - Renamed ParameterCheckTests to ParameterTests
- OSC_Sim (getSiteEnergy) - Updated function access from private to public so that it can be used during testing
- OSC_Sim (getSiteType) - Updated function access from private to public so that it can be used during testing
- parameters_default.txt - Formatting of interfacial energy shift model parameters
- .travis.yml - Copyright to show correct years, 2017-2018
- .travis.yml - Switched test coverage compiler configuration from GCC 5 to GCC 4.7 to make the test faster
- LICENSE - Copyright to show correct years, 2017-2018
- makefile - Copyright to show correct years, 2017-2018
- OSC_Sim (calculateExcitonCreationCoords) - Refactored code to be more robust and be guaranteed to find an appropriate empty site for creation if one exists
- KMC_Lattice - Submodule to latest version that fixes the chooseNextEvent and removeObject bugs that occurred when an object does not have a valid event
- OSC_Sim (createElectron, createExciton, CreateHole) - Code to catch out_of_range exceptions from the lattice class as the method for checking for invalid input coordinates
- test.cpp (EnergiesImportTests) - Test to use a bilayer so that it tests assignment of site energies to both donor and acceptor type sites
- test.cpp (ExcitonDiffusionTests) - Updated test by increasing N_tests to gather more statistics for checking the numerical accuracy of the lifetime and reducing N_tests when simply checking for relative changes in the diffusion length
- test.cpp (ExcitonDynamicsTests) - Updated test to use reduced Dynamics_transient_end and check how the program handles cutting off the tail end of the transient
- test.cpp (ExcitonDynamicsTests) - Updated test by reducing N_tests to shorten the test time
- Doxyfile - Updated the release number to v1.0.0-rc.1
- KMC_Lattice - Updated to latest version that has additional features needed by the steady transport test
- parameters_default.txt - Updated the version number to v1.0.0-rc.1
- main.cpp (main) - Updated version string to v1.0.0-rc.1
- main.cpp (main) - Updated Parameters object usage to use new format from the updated KMC_Lattice submodule
- Site_OSC - Storage of site energies using pointers to now directly storing a float
- OSC_Sim - Updated all appropriate functions to use float data type for site energies instead of double to save memory
- Parameters - Update base class from the struct used by the old KMC_Lattice submodule version to the new Parameters_Simulation class in the latest KMC_Lattice submodule version
- Parameters (importParameters) - Updated function to use format from new Parameters_Simulation base class in the latest KMC_Lattice version
- test.cpp - Updated functions to use format from new Parameters_Simulation base class in the latest KMC_Lattice version
- test.cpp (InterfacialEnergyShiftTests) - Updated site energies to use float data type instead of double
- test.cpp (ExcitonDynamicsTests) - Increased the lifetime tests tolerance to reduce likelihood of test failure
- main.cpp (main) - Updated results file output to use new getN_singlet_excitons_dissociation and getN_triplet_excitons_dissociated functions instead of the old getN_excitons_dissociated function
- OSC_Sim (calculateExcitonCreationCoords) - Changed function name to calculateRandomExcitonCreationCoords
- OSC_Sim (getN_excitons_dissociated) - Replace the function with getN_singlet_excitons_dissociated and getN_triplet_excitons_dissocated allowing the user to get separate stats
- OSC_Sim (N_excitons_dissociated) - Replace the counter variable with N_singlet_excitons_dissociated and N_triplet_excitons_dissociated to keep separate stats
- OSC_Sim (executeExcitonDissociation) - Changed event counter by counting the singlet and triplet dissociation events separated using the new counter variables
- test.cpp - Reduced default test parameters ToF_pnts_per_decade and Dynamics_pnts_per_decade down to 10 from 20 to have a larger time step to average over
- test.cpp (SteadyTransportTests) - Adjusted N_equilibration_events and N_tests to reduce test time
- README.md - Updated status to note that all features are now implemented and have major testing
- docs - Updated Doxygen documentation

### Removed
- main.cpp (Parameters_main) - the Parameters_main struct because these parameters are now stored on the Parameters class
- main.cpp (importParameters) - the importParameters function because this function is now contained in the Parameters class
- OSC_Sim (Parameters_OPV) - the Parameters_OPV struct because these parameters are now stored in the Parameters class
- OSC_Sim - Many of the parameter variables (replaced by a Parameters class member)
- OSC_Sim (checkParameters) - the checkParameters function because this function is now contained in the Parameters class
- OSC_Sim (createImportedMorphology) - Several unused local variables
- OSC_Sim (getEvents) - The function had no immediate purpose
- Parameters (checkParameters) - Several parameter checks that are now performed by the base Parameters_Simulation class
test.cpp (ParameterTests) - Several tests that are now handled by the KMC_Lattice submodule

### Fixed
- OSC_Sim (calculatePolaronEvents) - Update code so an error is not generated if no valid events are calculated because sometimes polarons can be trapped on sites where there are no valid events
- Parameters (checkParameters) - Bug where the function was not indicating an error when unable to read in the bool value for the Enable_recalc option
- Parameters (checkParameters) - Bug where the check for enabling of both imported energies and interfacial energy shift was implemented wrong
- OSC_Sim (getChargeExtractionMap) - Bug where the function as skipping output of the 0,0 coordinate data

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.5] - 2018-11-01 - Morphology Import and BKL Algorithm Update

### Added
- README.md - More detailed instructions for building and testing Excimontec
- CODEOWNERS.md - New file that designates which users have ownership over which parts of the codebase
- CODE_OF_CONDUCT.md - New file that specifies standards for developer and user interactions
- Excimontec namespace for all source code
- Usage of new KMC_Lattice namespace throughout code (replaces Utils namespace)
- main.cpp (main) - Usage of the new calculateTransitTimeHist function from the OSC_Sim class
- makefile - Commands to the clean target, so that make clean is also run for the KMC_Lattice submodule
- OSC_Sim - morphology files are imported so they are now opened and checked by the OSC_Sim class for easier testing
- Exciton.h - New calculateRateConstant function to the derived Event classes
- OSC_Sim (createImportedMorphology) - Code to allow import of uncompressed morphologies from Ising_OPV v3.2 and v4.0
- OSC_Sim (createImportedMorphology) - Code that parses and checks the morphology file version and handles attempts to import invalid morphology files
- test.cpp (ToFTests) - Electron ToF simulation tests
- test.cpp (ExcitonDynamicsTests) - Triplet exciton lifetime dynamics test
- test.cpp (ExcitonDiffusionTests) - Triplet exciton diffusion test
- test.cpp (IQETests) - Several basic IQE tests for field, temperature, delocalization, and recombination rate dependence of the charge separation yield in a bilayer morphology
- test.cpp (IQETests) - IQE testing for high illumination intensity to check whether higher order exciton-exciton and exciton-polaron annihilations and bimolecular recombination events increase as expected
- test.cpp (ParameterCheckTests) - Many new tests checking for how the OSC_Sim class handles initialization with invalid parameters
- Several new valid and invalid test morphology files from Ising_OPV to be used during morphology file import tests
- test.cpp (MorphologyImportTests) - New test function MorphologyImportTests and added numerous tests of importing valid and invalid sample morphology files to check for correct behavior

### Changed
- docs - Updated Doxygen API documentation
- README.md - Section and text formatting 
- README.md - Name for the transit time histogram file from "ToF_transit_time_dist.txt" to "ToF_transit_time_hist.txt"
- KMC_Lattice - Submodule to latest development branch version that has the BKL algorithm implemented, the Version class, and other bugfixes
- Code formatting for all source files using MS Visual Studio style
- Copyright in files to show correct years, 2017-2018
- main.cpp (main) - Name of output file for transit time data from "ToF_transit_time_dist.txt" to "ToF_transit_time_hist.txt" to specify that it is histogram data
- main.cpp (main) - Usage of calculateTransitTimeDist function to use new calculateTransitTimeHist that was renamed in the OSC_Sim class
- main.cpp (importParameters) - Usage of the new str2bool function from the KMC_Lattice namespace
- OSC_Sim (calculateExcitonEvents, calculatePolaronEvents) - Refactored code to use new BKL algorithm from the KMC_Lattice library
- OSC_Sim - Replaced morphology file ifstream with a morphology filename string
- OSC_Sim (calculateTransitTimeDist) - Function to output vector of pairs to match the format of the rest of the software package
- OSC_Sim (calculateTransitTimeDist) - Renamed function to calculateTransitTimeHist to be more descriptive about the data that it produces
- OSC_Sim (generateDynamicsExcitons) - Reduced frequency of command line status output for dynamics simulations
- OSC_Sim (generateToFPolarons) - Reduced frequency of command line status output for ToF simulations
- OSC_Sim (generateDynamicsExcitons) - Dynamics simulation command line output to print info when starting the first transient cycle
- OSC_Sim (executeExcitonHop, executePolaronHop) - Command line error output to give more details about the event when there is an error
- .travis.yml - Method for running the test executable
- test.cpp (ToFTests) - Refactored test to use updated calculateTransitTimeHist function
- test.cpp (MorphologyImportTests) - Test morphology file names to include their path relative to the root directory so that tests can be run from the root directory
- test.cpp (ExcitonDynamicsTests) - Reduced N_tests to shorten the test calculation time

### Removed
- Usage of Utils namespace throughout code (replaced by KMC_Lattice namespace)
- main.cpp (Parameters_main) - Morphology filename member variable
- main.cpp (main) - Morphology file ifstream
- main.cpp (main) - Opening and checking morphology files for import 
- main.cpp (importParameters) - Usage of the importBooleanParam function that was removed from the KMC_Lattice library
- Exciton.h - removed calculateExecutionTime overloaded definition from all derived Event classes
- Polaron.h - removed calculateExecutionTime overloaded definition from all derived Event classes
- OSC_Sim - Stopped compile warning about undefined pragma lines by commenting the out the pragma lines

### Fixed
- OSC_Sim (checkParameters) - A few parameter check bugs

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.4] - 2018-05-17 - Testing and Continuous Integration Update

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.3] - 2018-02-16

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.2] - 2018-02-08

--------------------------------------------------------------------------------------------------------------------------------

## [v1.0-beta.1] - 2017-10-23
