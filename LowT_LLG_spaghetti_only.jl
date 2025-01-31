using JSON: JSON
using Dates, JLD2, Statistics, Printf, PythonPlot
using ProgressBars, Rotations
using Sunny
using  LinearAlgebra, Random, MATLAB, Interpolations
include("model_CoTaS_5var.jl");
include("funcGridData.jl");
include("funcRemNaN.jl");
include("funcBroadening.jl");
include("funcFitSpectrum.jl");
include("funcVisualize.jl");
include("define_variables.jl");
include("function_bundle.jl");

J1  =  1.10970389;
j2  =  0.50299860;  J2  =  j2 * J1;
j3  =  0.02112092;  J3  =  j3 * J1; 
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
sys_small = system_initialize(sys_small, "3Q", J1);
print_wrapped_intensities(sys_small);
# randomize_spins!(sys_small);  
# for _ in 1:10  minimize_energy!(sys_small; maxiters = 3000);  end

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

# =======================================
# Save data for the future use, plotting etc,
# =======================================

using HDF5 
fname = "spaghetti_plot.h5"
fid = h5open(fname, "w")
write(fid,"J1", J1);  write(fid,"j2", j2);  write(fid,"j3", j3);  write(fid,"jc1", jc1);  write(fid,"jc2", jc2);
write(fid, "FitData02", FitData02);
close(fid)