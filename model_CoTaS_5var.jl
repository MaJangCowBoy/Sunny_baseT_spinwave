function CoTaS_5var(dim, J1, j2, j3, jc1, jc2; b1 = 0.06, rng = MersenneTwister(1234), Kz = -0.001)

  # a = b = 5.726;  c = 11.878;
  # alp = bet = 90;  gam = 120;
  # lattice_param = lattice_vectors(a,b,c,alp,bet,gam);
  # spgr = 182;
  # latCon = [a, b, c];  latAngle = float([alp, bet, gam]);
  
  # types = ["Co", "Ta", "Ta", "S"];
  # positions = [ [1/3, 2/3, 1/4], [ 0, 0, 0], [ 1/3, 2/3, 0.9991], [ 2/3, 0.0017, 0.132], ];
  crystl = Crystal("CoTaS.cif"; symprec = 1e-3);
  # crystl = Crystal(lattice_param, positions, spgr, setting = "2", types = types);
  magSub = subcrystal(crystl,"Co");
  magS = 1.5;
  randomN = Int64(round(100*rand(rng)));
  sys = System(magSub, dim, [SpinInfo(1, S=magS, g=2)], :dipole; seed=randomN); 
  # option --> :dipole, :heisenberg, :ising

  B1  = b1 * J1;
  J2  = j2 * J1;
  J3  = j3 * J1;
  Jc1 = jc1 * J1;
  Jc2 = jc2 * J1;
  
  if abs(B1) < 10 * eps()  set_exchange!(sys, J1,  Bond(1, 1, [1, 0, 0]));
  else  set_pair_coupling!(sys, (Si,Sj) -> J1*(Si'*Sj) + B1*(Si'*Sj)^2, Bond(1, 1, [1, 0, 0]));
  end
  # set_exchange!(sys, J1,  Bond(1, 1, [1, 0, 0]));
  set_exchange!(sys, Jc1, Bond(1, 2, [0, 0, 0]));
  set_exchange!(sys, Jc2, Bond(1, 2, [1, 1, 0]));
  set_exchange!(sys, J2,  Bond(1, 1, [1, 2, 0]));
  set_exchange!(sys, J3,  Bond(1, 1, [2, 0, 0]));
  set_onsite_coupling!(sys, S -> Kz*S[3]^2, 1);
  
  return sys, crystl;

end