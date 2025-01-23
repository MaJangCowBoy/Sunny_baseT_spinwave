function gaussian(x, σ)
  return exp(-x^2/(2*σ^2))/(σ*sqrt(2*π))
end

lorentzian_broadening(x, η) = η/(π*(x^2 + η^2))
lorentzian_broadening(η) = x -> lorentzian_broadening(x,η)


function E_broadening(ωvals, is, Intp)

        dims = size(is);  isE = zeros(dims);

        for (ωi, ω) in enumerate(ωvals)
    for (ω₀i, ω₀) in enumerate(ωvals)
      FWHM = Intp(ω₀);  η = FWHM/2.0;
      for qi in CartesianIndices(dims[1:end-1])
        isE[qi,ωi] += is[qi,ω₀i] * lorentzian_broadening(ω-ω₀, η);
      end
    end
  end

  return isE;
end

function Q_broadening(qxvals, qyvals, is, σx, σy; type = "2Dcut")

  dims = size(is);  isQ = zeros(dims);

  if type == "2Dcut"

    # for qi in CartesianIndices(dims[1:end-2])
    #   for qj in CartesianIndices(dims[1:end-2])
    #     for qkω in CartesianIndices(dims[end-1:end])
    #       isQ[qi,qkω] += 
    #         lorentzian_broadening(qxvals[qi[1]]-qxvals[qj[1]], σx) *
    #         lorentzian_broadening(qyvals[qi[2]]-qyvals[qj[2]], σy) *
    #         is[qj,qkω];
    #     end
    #   end
    # end
    for itX = 1:dims[1], itY = 1:dims[2]
      Xstt = max(1,itX-10);  Xend = min(dims[1],itX+10);
      Ystt = max(1,itY-10);  Yend = min(dims[2],itY+10);
      for id1 = Xstt:Xend, id2 = Ystt:Yend
        isQ[id1,id2,:,:] += exp(-(qxvals[itX]-qxvals[id1])^2/(2*σx^2)) *
                            exp(-(qyvals[itY]-qyvals[id2])^2/(2*σy^2)) * 
                            is[itX,itY,:,:];
      end
    end

  elseif type == "1Dcut"
      
    # for qi in 1:dims[1]
    #   for qj in 1:dims[1]
    #     for qqω in CartesianIndices(dims[2:end])
    #       isQ[qi,qqω] += 
    #         lorentzian_broadening(qxvals[qi]-qxvals[qj], σx) *
    #         is[qj,qqω];
    #     end
    #   end
    # end
    for itX = 1:dims[1]
      Xstt = max(1,itX-10);  Xend = min(dims[1],itX+10);
      for id1 = Xstt:Xend
        isQ[id1,:,:,:] += exp(-(qxvals[itX]-qxvals[id1])^2/(2*σx^2)) * is[itX,:,:,:];
      end
    end

  else 
    error("Q_broadening: 2Dcut/1Dcut is only implemented");
  end

  return isQ;
end

"""
QE_broaden_and_Sum, function summary.

Qsc: Q-broadening scaling factor, default = 1.0
     Basically, we need to scale the Q-broadening factor respect to given |Q|, 
     since broadening is get severe when we go to higher |Q|.
     In this case, we just to make sure that the broadening is not too severe.
     Conservatively, we usually set Qsc = 2.0 in our calculation.
     I think in-between 2.0 - 3.0 would be appropriate.
"""

function QE_broaden_and_Sum(is, qxvals, qyvals, ωvals, Intp, σx, σy; Qon = true, Eon = true, Qsc = 1.0)

  if Eon == true && Qon == false
    isE  = E_broadening(ωvals, is, Intp);
    return isE;
  elseif Eon == false && Qon == true
    isQ = Q_broadening(qxvals, qyvals, is, Qsc*σx, Qsc*σy);
    return isQ;
  elseif Eon == true && Qon == true
    isE  = E_broadening(ωvals, is, Intp);
    isQE = Q_broadening(qxvals, qyvals, isE, Qsc*σx, Qsc*σy);
    return isQE;
  else
    error("QE_broaden_and_Sum: this type of Qon/Eon combination is not implemented yet");
    return 0;
  end

end

function Average_Out_Spectra(is, ωvals, ωBot, ωTop; type = "constEcut")

  if type == "spaghettiPlot" # note this option is not checked for errors

    if length(ωBot) != length(ωTop)
      error("Average_Out_Spectra: ωBot and ωTop must have the same length");
    end

    d2d = zeros(size(is,1),length(ωBot));
    
    for (ωi, (ωb, ωt)) in enumerate(zip(ωBot,ωTop))
      _, IdxBot = findmin(x->abs(x-ωb), ωvals);
      _, IdxTop = findmin(x->abs(x-ωt), ωvals);

      d2d[:,ωi] = mean(is[:,IdxBot:IdxTop], dims = 2);
    end
  
  elseif type == "spaghettiPlot_v3"

    if length(ωBot) != length(ωTop)
      error("Average_Out_Spectra: ωBot and ωTop must have the same length");
    end

    d2d_tmp = zeros(size(is,1),length(ωvals));
    d2d     = zeros(size(is,1),length(ωBot));

    for x in 1:size(is,1), ω in 1:length(ωvals)
      d2d_tmp[x,ω] = mean(is[x,:,:,ω]);
    end

    for (ωi, (ωb, ωt)) in enumerate(zip(ωBot,ωTop))
      _, IdxBot = findmin(x->abs(x-ωb), ωvals);
      _, IdxTop = findmin(x->abs(x-ωt), ωvals);

      d2d[:,ωi] = mean(d2d_tmp[:,IdxBot:IdxTop], dims = 2);
    end

  elseif type == "constEcut"

    _, IdxBot = findmin(x->abs(x-ωBot), ωvals);
    _, IdxTop = findmin(x->abs(x-ωTop), ωvals);

          d2d = mean(is[:,:,:,IdxBot:IdxTop], dims = 4);
          d2d = dropdims(d2d, dims = 4);
    d2d = mean(d2d[:,:,:], dims = 3);
          d2d = dropdims(d2d, dims = 3);

  elseif type == "1Dcut"

    _, IdxBot = findmin(x->abs(x-ωBot), ωvals);
    _, IdxTop = findmin(x->abs(x-ωTop), ωvals);

          d2d = mean(is[:,:,:,IdxBot:IdxTop], dims = 4);
          d2d = dropdims(d2d, dims = 4);
    d2d = mean(d2d[:,:,:], dims = 3);
          d2d = dropdims(d2d, dims = 3);
    d2d = mean(d2d[:,:], dims = 2);
    d2d = dropdims(d2d, dims = 2);

  end

  return d2d;
end
      