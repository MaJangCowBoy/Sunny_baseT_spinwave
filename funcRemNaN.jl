function removeNaN(ExpData, FitData)

  if size(FitData) != size(ExpData)

    FitData = permutedims(FitData, (2,1))

    if size(FitData) != size(ExpData) # note that this option is not checked for errors.
      error("The size of the experimental data and the fitted data are not the same.")
    end

  end
  
  tt = findall(x -> isnan(x) == true, ExpData)
  for pos in tt
    FitData[pos] = NaN;
  end  

  return FitData;
end