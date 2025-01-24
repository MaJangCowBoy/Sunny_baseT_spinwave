"""
This code will calculate the goodness of fit.
Each goodness of fit would be feedback into Bayesian optimization algorithm.
* Originally written by: PJ Park
* Modified by: WH Cho, Last modified, 2025-01-20
* compatible with Sunny v0.5.11
"""

# =======================================
# Load Packages
# =======================================
using JSON: JSON
using Dates, JLD2, Statistics, Printf, PythonPlot
using ProgressBars
using Sunny
using  LinearAlgebra, Random, MATLAB, Interpolations
include("model_CoTaS_5var.jl");
include("funcGridData.jl");
include("funcRemNaN.jl");
include("funcBroadening.jl");
include("funcFitSpectrum.jl");
include("funcVisualize.jl");
include("define_variables.jl");

J1  =  1.10970389;
j2  =  0.50299860;  J2  =  j2 * J1;
j3  =  0.05112092;  J3  =  j3 * J1; 
jc1 =  1.27511611;  Jc1 = jc1 * J1;
jc2 = -0.19292577;  Jc2 = jc2 * J1;

vA = jc1 + 1.0 * jc2;  vB = jc1 - 0.5 * jc2;  vC = jc2;
j3_ideal   = 0.5 * (1 - (vA*vB-0.5*vA*vC+2*vB*vC)/sqrt(vA*vA-2*vA*vB+4*vB*vB) );
  # j3 is calculated from jc1 and jc2, which are input parameters.
  # this is for the magnetic ordering vector matching with the experimental data.
  # which is 1/3 a*.

include("Ei_dE.jl");  # Data extracted from Ei = 20 meV / Chopper = 150 Hz / Sample size = 80 mm
Ei = Ei[:];  dE = dE[:];  dE_4SEASONS = linear_interpolation(Ei, dE; extrapolation_bc=Line());


# =======================================
# Main loop for the calculation
# =======================================

seed = Dates.format(now(),"yyyy-mm-dd HH:MM:SS");  seed = parse(Int64,seed[end-1:end]);
rng = MersenneTwister(seed);

Tpara = 8;  kT = Tpara * Sunny.meV_per_K;
formfactors = [FormFactor("Co2"; g_lande = 2)];

dim_small = (3,3,1);  sys_small, crystl = CoTaS_5var(dim_small, J1, j2, j3, jc1, jc2; rng);
randomize_spins!(sys_small);  for _ in 1:10  minimize_energy!(sys_small; maxiters = 3000);  end
print_wrapped_intensities(sys_small)

dim = (30, 30, 20);
sys1, crystl = CoTaS_5var(dim, J1, j2, j3, jc1, jc2; rng);
sys2, crystl = CoTaS_5var(dim, J1, j2, j3, jc1, jc2; rng);

for x in axes(sys1.dipoles,1), y in axes(sys1.dipoles,2), z in axes(sys1.dipoles,3), b in axes(sys1.dipoles,4)
  idx = mod1(x,3);  idy = mod1(y,3);
  sys1.dipoles[x,y,z,b] = +1 * sys_small.dipoles[idx,idy,1,b];
  sys2.dipoles[x,y,z,b] = -1 * sys_small.dipoles[idx,idy,1,b];
end;  minimize_energy!(sys1);  minimize_energy!(sys2);

damping = 0.1;  dt = 0.02;  langevin = Langevin(dt; kT, damping);
ntherm = 4000;  # randomize_spins!(sys);


for i in ProgressBar(1:ntherm)  step!(sys1, langevin);  step!(sys2, langevin);  end

ωmax = 20.1;  Nω = 201;
sc1 = dynamical_correlations(sys1; dt=2*dt, nω=Nω, ωmax=ωmax, apply_g=true);
sc2 = dynamical_correlations(sys2; dt=2*dt, nω=Nω, ωmax=ωmax, apply_g=true);
add_sample!(sc1, sys1);  add_sample!(sc2, sys2);
# sys1 = nothing;  sys2 = nothing;

# =======================================
# Calculate scat. x-sec for each data
# and apply (Q,E) resolution convolution
# =======================================

formula1 = intensity_formula(sc1, :perp; formfactors = formfactors, kT = kT);
formula2 = intensity_formula(sc2, :perp; formfactors = formfactors, kT = kT);
Elist = available_energies(sc1); 

sqw02_1 = intensities_interpolated(sc1, qPathData02, formula1; interpolation = :linear)

sqw02_2 = intensities_interpolated(sc2, qPathData02, formula2; interpolation = :linear)

if res_on == true # unfortunately, QE-broadening is not working well for spaghetti plot now.
  br_sqw02_1 = E_broadening(Elist[4:end], sqw02_1[:,:,:,4:end], dE_4SEASONS);
  br_sqw02_2 = E_broadening(Elist[4:end], sqw02_2[:,:,:,4:end], dE_4SEASONS);
end

sqw03_06_1 = intensities_interpolated(sc1, qpts03, formula1; interpolation = :linear)
sqw07_10_1 = intensities_interpolated(sc1, qpts07, formula1; interpolation = :linear)

sqw03_06_2 = intensities_interpolated(sc2, qpts03, formula2; interpolation = :linear)
sqw07_10_2 = intensities_interpolated(sc2, qpts07, formula2; interpolation = :linear)

if res_on == true
  br_sqw03_06_1 = QE_broaden_and_Sum(sqw03_06_1, xGrd03, yGrd03, Elist, dE_4SEASONS, 0.020, 0.030; Qsc = 1.0);
  br_sqw07_10_1 = QE_broaden_and_Sum(sqw07_10_1, xGrd07, yGrd07, Elist, dE_4SEASONS, 0.015, 0.040; Qsc = 1.0);

  br_sqw03_06_2 = QE_broaden_and_Sum(sqw03_06_2, xGrd03, yGrd03, Elist, dE_4SEASONS, 0.020, 0.030; Qsc = 1.0);
  br_sqw07_10_2 = QE_broaden_and_Sum(sqw07_10_2, xGrd07, yGrd07, Elist, dE_4SEASONS, 0.015, 0.040; Qsc = 1.0);

  # (cf) Version 1 : I used following parameter for broadening.
  # 03-05 : 0.0150, 0.0400  # 06-08 : 0.0200, 0.0300
  # 09    : 0.0150, 0.0750  # 10    : 0.0087, 0.0750
  # (cf) Version 2 : I used following parameter for broadening.
  # 03-05 : 0.0800, 0.0100  # 06-08 : 0.0600, 0.0150
  # 09    : 0.0150, 0.0750  # 10    : 0.0087, 0.0750
end

# =======================================
# Data is average into 1D or 2D data ...
# =======================================

eBot02 = 0.60:0.20:20.00;  eTop02 = 0.80:0.20:20.20;
FitData02_1 = Average_Out_Spectra(br_sqw02_1, Elist[4:end], eBot02[:], eTop02[:]; type = "spaghettiPlot_v3"); 
FitData02_2 = Average_Out_Spectra(br_sqw02_2, Elist[4:end], eBot02[:], eTop02[:]; type = "spaghettiPlot_v3"); 
# ? This is for new version of spaghetti plot
# ! # FitData02 = Average_Out_Spectra(sqw02,     Elist, eBot02[:], eTop02[:]; type = "spaghettiPlot");
# ! # --> this is for old version of spaghetti plot
FitData02 = (FitData02_1 + FitData02_2) / 2;

FitData03_1 = Average_Out_Spectra(br_sqw03_06_1, Elist, eSum03[1], eSum03[2]);
FitData03_2 = Average_Out_Spectra(br_sqw03_06_2, Elist, eSum03[1], eSum03[2]);
FitData04_1 = Average_Out_Spectra(br_sqw03_06_1, Elist, eSum04[1], eSum04[2]);
FitData04_2 = Average_Out_Spectra(br_sqw03_06_2, Elist, eSum04[1], eSum04[2]);
FitData05_1 = Average_Out_Spectra(br_sqw03_06_1, Elist, eSum05[1], eSum05[2]);
FitData05_2 = Average_Out_Spectra(br_sqw03_06_2, Elist, eSum05[1], eSum05[2]);
FitData06_1 = Average_Out_Spectra(br_sqw03_06_1, Elist, eSum06[1], eSum06[2]);
FitData06_2 = Average_Out_Spectra(br_sqw03_06_2, Elist, eSum06[1], eSum06[2]);
FitData07_1 = Average_Out_Spectra(br_sqw07_10_1, Elist, eSum07[1], eSum07[2]);
FitData07_2 = Average_Out_Spectra(br_sqw07_10_2, Elist, eSum07[1], eSum07[2]);
FitData08_1 = Average_Out_Spectra(br_sqw07_10_1, Elist, eSum08[1], eSum08[2]);
FitData08_2 = Average_Out_Spectra(br_sqw07_10_2, Elist, eSum08[1], eSum08[2]);
FitData09_1 = Average_Out_Spectra(br_sqw07_10_1, Elist, eSum09[1], eSum09[2]);
FitData09_2 = Average_Out_Spectra(br_sqw07_10_2, Elist, eSum09[1], eSum09[2]);
FitData10_1 = Average_Out_Spectra(br_sqw07_10_1, Elist, eSum10[1], eSum10[2]);
FitData10_2 = Average_Out_Spectra(br_sqw07_10_2, Elist, eSum10[1], eSum10[2]);

FitData03 = (FitData03_1 + FitData03_2) / 2;
FitData04 = (FitData04_1 + FitData04_2) / 2;
FitData05 = (FitData05_1 + FitData05_2) / 2;
FitData06 = (FitData06_1 + FitData06_2) / 2;
FitData07 = (FitData07_1 + FitData07_2) / 2;
FitData08 = (FitData08_1 + FitData08_2) / 2;
FitData09 = (FitData09_1 + FitData09_2) / 2;
FitData10 = (FitData10_1 + FitData10_2) / 2;

# =======================================
# Save data for the future use, plotting etc,
# =======================================

using HDF5 
fname = "ordered_phase_scattering.h5"
fid = h5open(fname, "w")
write(fid,"J1", J1);  write(fid,"j2", j2);  write(fid,"j3", j3);  write(fid,"jc1", jc1);  write(fid,"jc2", jc2);
write(fid, "FitData02", FitData02);
write(fid, "FitData03", FitData03);  write(fid, "FitData04", FitData04);  write(fid, "FitData05", FitData05);
write(fid, "FitData06", FitData06);  write(fid, "FitData07", FitData07);  write(fid, "FitData08", FitData08);
close(fid)