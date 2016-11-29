all: gb bt st gs_lbnl gs_sigma g4p

gb: 
	gcc -o grid_basis grid_basis.c point.c cyl_point.c calc_signal.c fields_ppc.c detector_geometry_ppc.c signal_calc_util.c -lm -lreadline

bt: 
	gcc -o basis_test basis_test.c read_basis.c -lm -lreadline

st:
	gcc -o st signal_tester_ppc.c point.c cyl_point.c calc_signal.c fields_ppc.c detector_geometry_ppc.c signal_calc_util.c -lm -lreadline

gs_lbnl:
	gcc -o g_s geant_siggen.c point.c cyl_point.c calc_signal.c fields_ppc.c detector_geometry_ppc.c signal_calc_util.c read_mode2.c -lm -lreadline

gs_sigma:
	gcc -o g_s geant_siggen.c point.c cyl_point.c calc_signal.c fields_ppc.c detector_geometry_ppc.c signal_calc_util.c read_g4.c -lm -lreadline

g4p:
	gcc -o g4p G4_parser.c -lm -lreadline
