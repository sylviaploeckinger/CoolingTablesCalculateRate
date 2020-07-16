# CoolingTablesCalculateRate
Stand-alone C program that reads in the CoolingTables from Ploeckinger & Schaye (in prep.), 
interpolates in 4 dimensions,and returns the total net cooling rate for arbitrary abundance 
ratios, redshifts, densities and temperatures

Individual routines can be directly used and / or adapted to couple to different hydro-dynamic codes.

The cooling tables are read in and initialized in the routine cooling_init_backend. 
The routine cooling_function is a wrapper that loops over internal energies and 
calls the main routine colibre_cooling_rate.

For coupling with a hydro code, the net cooling rate for a resolution element is calcated with:

    cooling_rate = colibre_cooling_rate(log_u_cgs, redshift, n_H_cgs, pow(10., logZZsol), abundance_ratio,
                   n_H_index, d_n_H, met_index, d_met, red_index, d_red, cooling);

where:

log_u_cgs: log10 of the internal energy per mass of the particle in cgs units [erg g-1]
redshift : redshift of the simulation
n_H_cgs  : hydrogen number density in cgs [cm-3]
logZZsol : log10 of the metal mass fraction of the particle over the solar metal mass fraction
abundance_ratio: abundance ratio for the traced elements relative to the abundance ratio assumed in the tables
                 (see routine abundance_ratio_to_solar for details)
indices and deltas for each table dimension (get with get_index_1d) 
cooling  : structure that contains all necessary table values

Run the code:
- requires gsl (https://www.gnu.org/software/gsl/)

1) change the path to the cooling tables in cooling.c
   change gas density and redshift in cooling_function.c (optional)
2) make
3) ./cool
