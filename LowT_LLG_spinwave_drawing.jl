"""
This code will calculate the goodness of fit.
Each goodness of fit would be feedback into Bayesian optimization algorithm.
* Originally written by: PJ Park
* Modified by: WH Cho, Last modified, 2025-01-20
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

if isdir("FitCheck") == false  mkdir("FitCheck");  end

envVar = JSON.parsefile("param.json");
res_on = get(envVar,"res_on", false);  
const_on = get(envVar,"const_on", false);
vis_on = get(envVar,"vis_on", false);
DirOfBOdotJSON = get(envVar,"DirOfBOdotJSON", "");

DirOfExpData  = get(envVar,"DirOfExpData", "");   
DirOfGridData = get(envVar,"DirOfGridData", "");  
DirOfResData  = get(envVar,"DirOfResData", "");   
DirOfJsonData = get(envVar,"DirOfJsonData", "");

if DirOfExpData[end]  != "/"  DirOfExpData  = DirOfExpData*"/"   end
if DirOfGridData[end] != "/"  DirOfGridData = DirOfGridData*"/"  end
if DirOfResData[end]  != "/"  DirOfResData  = DirOfResData*"/"   end
if DirOfJsonData[end] != "/"  DirOfJsonData = DirOfJsonData*"/"  end

NameOfExpData  = get(envVar,"NameOfExpData", "");
NameOfGridData = get(envVar,"NameOfGridData", "");
NameOfResData  = get(envVar,"NameOfResData", "");

# params = JSON.parsefile(DirOfBOdotJSON*"BO.json")
# idx = get(params, "ID",    0);  J1  = get(params, "J1",  0.0);  
# j2  = get(params, "J2",  0.0);  jc1 = get(params, "Jc1", 0.0);
# jc2 = get(params, "Jc2", 0.0);  j3  = get(params, "J3",  0.0);

J1  =  1.10970389;
j2  =  0.50299860;  J2  =  j2 * J1;
j3  =  0.05112092;  J3  =  j3 * J1; 
jc1 =  1.27511611;  Jc1 = jc1 * J1;
jc2 = -0.19292577;  Jc2 = jc2 * J1;

vA = jc1 + 1.0 * jc2;  vB = jc1 - 0.5 * jc2;  vC = jc2;
j3   = 0.5 * (1 - (vA*vB-0.5*vA*vC+2*vB*vC)/sqrt(vA*vA-2*vA*vB+4*vB*vB) );
  # j3 is calculated from jc1 and jc2, which are input parameters.
  # this is for the magnetic ordering vector matching with the experimental data.
  # which is 1/3 a*.

idx = 10000;

println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("Execute [", PROGRAM_FILE,"] with param_id = $(idx)")

# =======================================
# Experimental data is loaded here.
# =======================================
AbsExpData  = MatFile(DirOfExpData*NameOfExpData);
AbsGridData = MatFile(DirOfGridData*NameOfGridData);
AbsResData  = MatFile(DirOfResData*NameOfResData);

ExpData02 = get_variable(AbsExpData, "ExpData02");
path02, qGrd02, qItv02, qPtN02, eBot02, eTop02, us02, vs02, ws02, qWid02 = 
  loadGridData(AbsGridData,"02"; type = "spaghettiPlot_v3");

# ! # path02, qGrd02, qItv02, eBot02, eTop02 = loadGridData(AbsGridData,"02"; type = "spaghettiPlot");
# ! # --> this is for old version of spaghetti plot

ExpData03 = get_variable(AbsExpData, "ExpData03");
xBas03, yBas03, zBas03, xGrd03, yGrd03, zSum03, eSum03, HKL03 = loadGridData(AbsGridData,"03");

ExpData04 = get_variable(AbsExpData, "ExpData04");
xBas04, yBas04, zBas04, xGrd04, yGrd04, zSum04, eSum04, HKL04 = loadGridData(AbsGridData,"04");

ExpData05 = get_variable(AbsExpData, "ExpData05");
xBas05, yBas05, zBas05, xGrd05, yGrd05, zSum05, eSum05, HKL05 = loadGridData(AbsGridData,"05");

ExpData06 = get_variable(AbsExpData, "ExpData06");
xBas06, yBas06, zBas06, xGrd06, yGrd06, zSum06, eSum06, HKL06 = loadGridData(AbsGridData,"06");

ExpData07 = get_variable(AbsExpData, "ExpData07");
xBas07, yBas07, zBas07, xGrd07, yGrd07, zSum07, eSum07, HKL07 = loadGridData(AbsGridData,"07");

ExpData08 = get_variable(AbsExpData, "ExpData08");
xBas08, yBas08, zBas08, xGrd08, yGrd08, zSum08, eSum08, HKL08 = loadGridData(AbsGridData,"08");

ExpData09 = get_variable(AbsExpData, "ExpData09");
xBas09, yBas09, zBas09, xGrd09, yGrd09, zSum09, eSum09, HKL09 = loadGridData(AbsGridData,"09");

ExpData10 = get_variable(AbsExpData, "ExpData10");
xBas10, yBas10, zBas10, xGrd10, yGrd10, zSum10, eSum10, HKL10 = loadGridData(AbsGridData,"10");

# EGrd = get_variable(AbsResData, "E");    EGrd = EGrd[:];
# dEGrd = get_variable(AbsResData, "dE");  dEGrd = dEGrd[:];
# dE_4SEASONS = linear_interpolation(EGrd, dEGrd; extrapolation_bc=Line());

include("Ei_dE.jl");  # Data extracted from Ei = 20 meV / Chopper = 150 Hz / Sample size = 80 mm
Ei = Ei[:];  dE = dE[:];  dE_4SEASONS = linear_interpolation(Ei, dE; extrapolation_bc=Line());

println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("Experiment data is loaded with param_id = $(idx)")

# =======================================
# Main loop for the calculation
# =======================================
println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("Main loop is started with param_id = $(idx)")

seed = Dates.format(now(),"yyyy-mm-dd HH:MM:SS");  seed = parse(Int64,seed[end-1:end]);
rng = MersenneTwister(seed);

Tpara = 8;  kT = Tpara * Sunny.meV_per_K;
formfactors = [FormFactor("Co2"; g_lande = 2)];

dim_small = (3,3,1);
sys_small, crystl = CoTaS_5var(dim_small, J1, j2, j3, jc1, jc2; rng);
randomize_spins!(sys_small);
for _ in 1:10  minimize_energy!(sys_small; maxiters = 3000);  end
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
sys1 = nothing;  sys2 = nothing;

# =======================================
# Q-grid is generated here.
# Each data has different Q-grid.
# =======================================
qPathData02 = makeSweepPoints(path02,qPtN02,vs02,ws02,qWid02); 
# ? This is for new version of spaghetti plot
# ! # qPathData02 = makeSweepPath(crystl,path02,qGrd02,qItv02);
# ! # --> this is for old version of spaghetti plot

qpts03 = makeGridPoint(xBas03,yBas03,xGrd03,yGrd03,zBas03,zSum03,HKL03);
qpts04 = makeGridPoint(xBas04,yBas04,xGrd04,yGrd04,zBas04,zSum04,HKL04);
qpts05 = makeGridPoint(xBas05,yBas05,xGrd05,yGrd05,zBas05,zSum05,HKL05);
qpts06 = makeGridPoint(xBas06,yBas06,xGrd06,yGrd06,zBas06,zSum06,HKL06);
qpts07 = makeGridPoint(xBas07,yBas07,xGrd07,yGrd07,zBas07,zSum07,HKL07);
qpts08 = makeGridPoint(xBas08,yBas08,xGrd08,yGrd08,zBas08,zSum08,HKL08);
qpts09 = makeGridPoint(xBas09,yBas09,xGrd09,yGrd09,zBas09,zSum09,HKL09);
qpts10 = makeGridPoint(xBas10,yBas10,xGrd10,yGrd10,zBas10,zSum10,HKL10);

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

println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("Basic SQW calculation is done with param_id = $(idx)")

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

println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("(Q,E) broadening is done with param_id = $(idx)")

# =======================================
# Save data for the future use, plotting etc,
# =======================================

using HDF5 
fname = "ordered_phase_scattering.h5"
fid = h5open(fname, "w")
write(fid, "FitData02", FitData02);
write(fid, "FitData03", FitData03);  write(fid, "FitData04", FitData04);  write(fid, "FitData05", FitData05);
write(fid, "FitData06", FitData06);  write(fid, "FitData07", FitData07);  write(fid, "FitData08", FitData08);
close(fid)

# ====================================================================================== #
# =============== Remaining part would be unnecessary for current use ================== # 
# ====================================================================================== # 

"""
# =======================================
# Apply data coverage
# =======================================
FitData02 = removeNaN(ExpData02, FitData02);
FitData03 = removeNaN(ExpData03, FitData03);
FitData04 = removeNaN(ExpData04, FitData04);
FitData05 = removeNaN(ExpData05, FitData05);
FitData06 = removeNaN(ExpData06, FitData06);
FitData07 = removeNaN(ExpData07, FitData07);
FitData08 = removeNaN(ExpData08, FitData08);
FitData09 = removeNaN(ExpData09, FitData09);
FitData10 = removeNaN(ExpData10, FitData10);

println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("Data coverage is reflected with param_id = $(idx)")

# =======================================
# Calculate chi-squared value
# =======================================

S₁ = zeros(10,1);  S₂ = zeros(10,1);  χ² = zeros(10,1);

S₁[2] , S₂[2] , χ²[2] , ErrMap02 = fitSpectrum(ExpData02,FitData02,const_on);
S₁[3] , S₂[3] , χ²[3] , ErrMap03 = fitSpectrum(ExpData03,FitData03,const_on);
S₁[4] , S₂[4] , χ²[4] , ErrMap04 = fitSpectrum(ExpData04,FitData04,const_on);
S₁[5] , S₂[5] , χ²[5] , ErrMap05 = fitSpectrum(ExpData05,FitData05,const_on);
S₁[6] , S₂[6] , χ²[6] , ErrMap06 = fitSpectrum(ExpData06,FitData06,const_on);
S₁[7] , S₂[7] , χ²[7] , ErrMap07 = fitSpectrum(ExpData07,FitData07,const_on);
S₁[8] , S₂[8] , χ²[8] , ErrMap08 = fitSpectrum(ExpData08,FitData08,const_on);
S₁[9] , S₂[9] , χ²[9] , ErrMap09 = fitSpectrum(ExpData09,FitData09,const_on);
S₁[10], S₂[10], χ²[10], ErrMap10 = fitSpectrum(ExpData10,FitData10,const_on);
χ²sum = sum(χ²[2:end]);

println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("Chi square is obtained with param_id = $(idx)")

if vis_on == true
  visualizeFit(ExpData02, FitData02, ErrMap02,  2, idx; type = "spaghettiPlot");
  visualizeFit(ExpData03, FitData03, ErrMap03,  3, idx);
  visualizeFit(ExpData04, FitData04, ErrMap04,  4, idx);
  visualizeFit(ExpData05, FitData05, ErrMap05,  5, idx);
  visualizeFit(ExpData06, FitData06, ErrMap06,  6, idx);
  visualizeFit(ExpData07, FitData07, ErrMap07,  7, idx);
  visualizeFit(ExpData08, FitData08, ErrMap08,  8, idx);
  visualizeFit(ExpData09, FitData09, ErrMap09,  9, idx);
  visualizeFit(ExpData10, FitData10, ErrMap10, 10, idx);
end

filename = @sprintf("FitData_%d.jld2",idx);
jldsave(filename; FitData02, FitData03, FitData04, FitData05, 
       FitData06, FitData07, FitData08, FitData09, FitData10);

println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("Fitted result visualization is done for id = $(idx)")

# =======================================
# Save the data for Bayesian optimization
# =======================================
println("")
println(Dates.format(now(), "yy-mm-dd HH:MM"))
println("Sunny calculation finished for param_id = $(idx)")
println("Chi-sqaure value for param_id = $(idx) is $(χ²sum)")
println("")

open(DirOfJsonData*"idx_$idx.json", "w") do f
  JSON.print(f, params)
end

open(DirOfJsonData*"loss_$idx.csv", "w") do file
  write(file, string(-χ²sum))
end
"""