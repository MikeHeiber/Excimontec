# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

--------------------------------------------------------------------------------------------------------------------------------

## [Unreleased]

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
- parameters_default.txt - Three new parameters (Enable_interfacial_energy_shift, Energy_shift_donor, Energy_shift_acceptor) for the interfacial energy shify model
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
- OSC_Sim - New N_events_executed counter and getN_events_executed function to overwrite the the one in the base Simulation class to separately keep track oh how many events have been executed
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
- KMC_Lattice - Submodule to latest development branch version that fixes the chooseNextEvent and removeObject bugs that occurred when an object does not have a valid event
- OSC_Sim (createElectron, createExciton, CreateHole) - Code to catch out_of_range exceptions from the lattice class as the method for checking for invalid input coordinates
- test.cpp (EnergiesImportTests) - Test to use a bilayer so that it tests assignment of site energies to both donor and acceptor type sites
- test.cpp (ExcitonDiffusionTests) - Updated test by increasing N_tests to gather more statistics for checking the numerical accuracy of the lifetime and reducing N_tests when simply checking for relative changes in the diffusion length
- test.cpp (ExcitonDynamicsTests) - Updated test to use reduced Dynamics_transient_end and check how the program handles cutting off the tail end of the transient
- test.cpp (ExcitonDynamicsTests) - Updated test by reducing N_tests to shorten the test time

### Removed
- main.cpp (Parameters_main) - the Parameters_main struct because these parameters are now stored on the Parameters class
- main.cpp (importParameters) - the importParameters function because this function is now contained in the Parameters class
- OSC_Sim (Parameters_OPV) - the Parameters_OPV struct because these parameters are now stored in the Parameters class
- OSC_Sim - Many of the parameter variables (replaced by a Parameters class member)
- OSC_Sim (checkParameters) - the checkParameters function because this function is now contained in the Parameters class
- OSC_Sim (createImportedMorphology) - Several unused local variables
- OSC_Sim (getEvents) - The function had no immediate purpose

### Fixed
- OSC_Sim (calculatePolaronEvents) - Update code so an error is not generated if no valid events are calculated because sometimes polarons can be trapped on sites where there are no valid events
- Parameters (checkParameters) - Bug where the function was not indicating an error when unable to read in the bool value for the Enable_recalc option
- Parameters (checkParameters) - Bug where the check for enabling of both imported energies and interfacial energy shift was implemented wrong

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
- test.cpp (MorphologyImportTests) - Test morphology file names to include their path relative to the root directory so that tests can be run from the root dir
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