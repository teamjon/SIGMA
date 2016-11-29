-------
C Codes
-------
- grid_basis.c
	-> Used to generate a basis of pulse shapes for a specified r,z,phi range
	-> Currently set to generate a basis for 1/8th/1 segment of the detector
	-> Output binary files stored in /basis_dat_copies/
	-> Output contains xyz, rpz, td, charge pulses
	-> Format defined in ‘pdecomp.h’

- signal_tester_ppc.c
	-> Used to test the underlying code
	-> Enables pulses to be output for a specified interaction pos
	-> If CHATTY enabled, very useful for debugging
	-> Full list of commands available using 'help'

- basis_test.c
	-> Similar to grid_basis.c but outputs ASCII file

- mode2.c / mode2_1evnt.c / mode2_multi.c
	-> Simple code to generate random data in mode 2 format for testing 'read_mode2.c'

- geant_siggen.c
	-> Provides an interface between Geant data and SigGen
	-> Uses 'read_mode2.c' to get info from Geant binary output
	-> Processes events through SigGen
	-> Outputs binary file with Geant info + charge pulses and drift times
	-> Structures for all data type defined in 'geant_siggen.h'
	29/11/16
	-> Adapted to read in SIGMA Geant output using 'read_g4.c'
	-> New data structures for binary files defined in 'geant_siggen.h'

- calc_signal.c
	-> Main script for tracking charge through the detector
	-> Added function 'calc_dt' for calculating drift times
	-> Used charge pulses on 1ns scale in addition to user specified 'start', 'end' pos as a %

- read_mode2.c/read_g4.c
	-> Simple unctions to read a single event from the Geant binary output

- G4_parser.c
	-> Parser to convert the output struct from Geant4 simulations (Marc_G4_Struct) to SigGen input struct (SigGen_4_Struct)

All other codees are SigGen codes, used to calculate the traces

-------
Folders
-------

- Geant_Heather
	-> Contains the binary outputs from Heather's Geant simulations

- Matlab
	-> Contains Matlab scripts used to interpret data and plot results
	-> Also contains all results in structured binary '.mat' format
	-> Structure of each file can be found by looking at the pasers files
-----
Other
-----

- compile.txt
	-> Contains command line instructions for compiling each of the C codes above
	-> Makefile also available

- Marco_fn.h
	-> Header file containing Marco's drift time calculation for better comparison with exp data
	-> Not currently implemented but may be useful to determine dt in same way as exp for test purposes
