using JSON: JSON
using Dates, JLD2, Statistics, Printf
using ProgressBars, Rotations
using Sunny
using  LinearAlgebra, Random, MATLAB, Interpolations
include("model_CoTaS_5var.jl");
include("funcGridData.jl");
include("funcRemNaN.jl");
include("funcBroadening.jl");
include("funcFitSpectrum.jl");
include("funcVisualize.jl");
include("function_bundle.jl");

J1  =  1.1955542291322174;
j2  =  0.3002342874679325;  J2  =  j2 * J1;
j3  =  0.0327544014014257;  J3  =  j3 * J1; 
jc1 =  1.3575164668571524;  Jc1 = jc1 * J1;
jc2 = -0.0396314440925254;  Jc2 = jc2 * J1;

vA = jc1 + 1.0 * jc2;  vB = jc1 - 0.5 * jc2;  vC = jc2;
j3_ideal   = 0.5 * (1 - (vA*vB-0.5*vA*vC+2*vB*vC)/sqrt(vA*vA-2*vA*vB+4*vB*vB) );
# j3 = j3_ideal

  # j3 is calculated from jc1 and jc2, which are input parameters.
  # this is for the magnetic ordering vector matching with the experimental data.
  # which is 1/3 a*.

include("Ei_dE.jl");  # Data extracted from Ei = 20 meV / Chopper = 150 Hz / Sample size = 80 mm
Ei = Ei[:];  dE = 0.3 * dE[:];  dE_4SEASONS = linear_interpolation(Ei, dE; extrapolation_bc=Line());

seed = Dates.format(now(),"yyyy-mm-dd HH:MM:SS");  seed = parse(Int64,seed[end-1:end]);
rng = MersenneTwister(seed);

Tpara = 3;  kT = Tpara * Sunny.meV_per_K;

dim_small = (3,3,1);
sys_small, crystl = CoTaS_5var(dim_small, J1, j2, j3, jc1, jc2; rng);
sys_small = system_initialize(sys_small, "3Q", J1);
#@ sys_small, crystl = CoTaS_5var(dim_small, J1, j2, j3, jc1, jc2; rng, b1 = 0.0);
#@ sys_small = system_initialize(sys_small, "1Q_1", J1);
print_wrapped_intensities(sys_small);

dim = (30, 30, 2);  sys, crystl = CoTaS_5var(dim, J1, j2, j3, jc1, jc2; rng);
#@ dim = (30, 30, 2);  sys, crystl = CoTaS_5var(dim, J1, j2, j3, jc1, jc2; rng, b1 = 0.0);

damping = 0.1;  dt = 0.02;  ntherm = 4000;  langevin = Langevin(dt; kT, damping);
  
for x in axes(sys.dipoles,1), y in axes(sys.dipoles,2), 
    z in axes(sys.dipoles,3), b in axes(sys.dipoles,4)
  idx = mod1(x,3);  idy = mod1(y,3);
  sys.dipoles[x,y,z,b] = +1 * sys_small.dipoles[idx,idy,1,b];
end;  minimize_energy!(sys);

for i in ProgressBar(1:100)  step!(sys, langevin);  end;  print_wrapped_intensities(sys);

ωmax = 3.0;  nω = 61;
sc = dynamical_correlations(sys; dt, nω, ωmax);  
add_sample!(sc, sys);

shift  = [ 0.0, 0.0, 1.0];
formfactors = [FormFactor("Co2"; g_lande = 2)];
formula = intensity_formula(sc, :perp; formfactors = formfactors, kT = kT);

N = 240
xGrd = [ x for x in range(-1,1,N)];  yGrd = [ y for y in range(-1,1,N)];

#?
basis1 = [ 1.0, 0.0, 0.0];  basis2 = [-0.5, 1.0, 0.0];
qptsL001 = [ shift + x * basis1 + y * basis2  for x in xGrd, y in yGrd ];

sqwL001A = intensities_interpolated(sc, qptsL001, formula; interpolation = :linear);
Elist = available_energies(sc);
# br_sqwL001A = QE_broaden_and_Sum(sqwL001A, xGrd, yGrd, Elist, dE_4SEASONS, 0.020, 0.030; Qsc = 1.0);
br_sqwL001A = Q_broadening(xGrd, yGrd, sqwL001A, 0.020, 0.026);
br_sqwL001A_2 = E_broadening(Elist, br_sqwL001A, dE_4SEASONS);

#?
basis1 = [ 0.0, 1.0, 0.0];  basis2 = [ 1.0,-0.5, 0.0];  
qptsL001 = [ shift + x * basis1 + y * basis2  for x in xGrd, y in yGrd ];

sqwL001B = intensities_interpolated(sc, qptsL001, formula; interpolation = :linear);
Elist = available_energies(sc);
# br_sqwL001B = QE_broaden_and_Sum(sqwL001B, xGrd, yGrd, Elist, dE_4SEASONS, 0.020, 0.030; Qsc = 1.0);
br_sqwL001B = Q_broadening(xGrd, yGrd, sqwL001B, 0.020, 0.026);
br_sqwL001B_2 = E_broadening(Elist, br_sqwL001B, dE_4SEASONS);

#?
basis1 = [-1.0, 1.0, 0.0];  basis2 = [ 0.5, 0.5, 0.0];  
qptsL001 = [ shift + x * basis1 + y * basis2  for x in xGrd, y in yGrd ];

sqwL001C = intensities_interpolated(sc, qptsL001, formula; interpolation = :linear);
Elist = available_energies(sc);
# br_sqwL001C = QE_broaden_and_Sum(sqwL001C, xGrd, yGrd, Elist, dE_4SEASONS, 0.020, 0.030; Qsc = 1.0);
br_sqwL001C = Q_broadening(xGrd, yGrd, sqwL001C, 0.020, 0.026);
br_sqwL001C_2 = E_broadening(Elist, br_sqwL001C, dE_4SEASONS);

#?
diff, idx = findmin(abs.(Elist .- 1.0));

data_2DcutA = br_sqwL001A_2[:,:,idx];  data_2DcutD = reverse(data_2DcutA, dims = 1);
data_2DcutB = br_sqwL001B_2[:,:,idx];  data_2DcutE = reverse(data_2DcutB, dims = 1);
data_2DcutC = br_sqwL001C_2[:,:,idx];  data_2DcutF = reverse(data_2DcutC, dims = 1);
data2Dcut = (data_2DcutA + data_2DcutB + data_2DcutC + data_2DcutD + data_2DcutE + data_2DcutF) / 6.0;

#?
data2Dcut_flip = reverse(data2Dcut, dims = 2);  data2Dcut = (data2Dcut + data2Dcut_flip) / 2.0;

using CairoMakie
figure = CairoMakie.Figure();  ax = CairoMakie.Axis(figure[1,1]; aspect = 2/√3);
heatmap!(ax, xGrd, yGrd, data2Dcut; colormap = :viridis);  save("2Dcut_prime.png", figure);


#?
using HDF5, Printf
filename = @sprintf("constEcut_2p2_3Q_same_param.h5");
h5open(filename, "w") do file
  write(file, "xGrd", xGrd);  write(file, "yGrd", yGrd);  write(file, "data2Dcut", data2Dcut);
end;