function system_initialize(sys::System, keyword::String, J1::Float64)
        
  if keyword == "1Q_1"
    k = [ 1/3, 0, 0];  ψ = π/6;  ϕ = 120 * π/180;  S0 = [0, 0, 1.5];
    for x in axes(sys.dipoles,1), y in axes(sys.dipoles,2), z in axes(sys.dipoles,3), b in axes(sys.dipoles,4)
      θ = 2π * dot(k,[x,y,z]) + ϕ * Float64(isodd(b));
      sys.dipoles[x,y,z,b] = RotZ(ψ) * RotX(θ) * S0;
    end;
  elseif keyword == "1Q_2"
    k = [ 0,1/3, 0];  ψ = π/2;  ϕ = -120 * π/180;  S0 = [0, 0, 1.5];
    for x in axes(sys.dipoles,1), y in axes(sys.dipoles,2), z in axes(sys.dipoles,3), b in axes(sys.dipoles,4)
      θ = 2π * dot(k,[x,y,z]) + ϕ * Float64(isodd(b));
      sys.dipoles[x,y,z,b] = RotZ(ψ) * RotX(θ) * S0;
    end;
  elseif keyword == "1Q_3"
    k = [-1/3,1/3, 0];  ψ = 5π/6;  ϕ = -120 * π/180;  S0 = [0, 0, 1.5];
    for x in axes(sys.dipoles,1), y in axes(sys.dipoles,2), z in axes(sys.dipoles,3), b in axes(sys.dipoles,4)
      θ = 2π * dot(k,[x,y,z]) + ϕ * Float64(isodd(b));
      sys.dipoles[x,y,z,b] = RotZ(ψ) * RotX(θ) * S0;
    end;
  elseif keyword == "3Q"
    S1 = RotZ(π/3) * RotX(2π/3) * [0, 0, 1.5];
    S2 = RotZ(π/3) * RotX(4π/3) * [0, 0, 1.5];
    S3 = RotZ(π/3) * RotX(6π/3) * [0, 0, 1.5];
    for z in axes(sys.dipoles,3)
      sys.dipoles[1,3,z,1] = S2;  sys.dipoles[2,3,z,1] = S3;  sys.dipoles[3,3,z,1] = S3;
      sys.dipoles[1,2,z,1] = S2;  sys.dipoles[2,2,z,1] = S1;  sys.dipoles[3,2,z,1] = S2;
      sys.dipoles[1,1,z,1] = S1;  sys.dipoles[2,1,z,1] = S1;  sys.dipoles[3,1,z,1] = S3;
    
      sys.dipoles[1,3,z,2] = -S2;  sys.dipoles[2,3,z,2] = -S3;  sys.dipoles[3,3,z,2] = -S2;
      sys.dipoles[1,2,z,2] = -S1;  sys.dipoles[2,2,z,2] = -S1;  sys.dipoles[3,2,z,2] = -S2;
      sys.dipoles[1,1,z,2] = -S1;  sys.dipoles[2,1,z,2] = -S3;  sys.dipoles[3,1,z,2] = -S3;
    end
  else
    error("Invalid keyword.")
  end
      
  dt = 0.1/(J1 * 1.5 * 1.5);  damping = 0.1;  
  langevin = Langevin(dt; damping, kT = 0.0);
  langevin.kT = 0.1 * meV_per_K;  for _ in 1:100000  step!(sys, langevin)  end
  langevin.kT = 0.0;              for _ in 1:100000  step!(sys, langevin)  end
  
  for _ in 1:20  minimize_energy!(sys; maxiters = 3000);  end

  return sys;
end