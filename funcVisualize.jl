function visualizeFit(ExpData, FitData, ErrMap, N, idx; type = "constEcut")

  if type == "spaghettiPlot"
    ExpData = ExpData';  FitData = FitData';  ErrMap = ErrMap';
  end

  for (it,(x,y)) in enumerate(zip(ExpData,FitData))
    if isnan(x) == true && isnan(y) == true
      ErrMap[it] = NaN;
    end
  end

  figure(figsize=(15,3))
  subplot(1,3,1)
  pcolor(ExpData)
  subplot(1,3,2)
  pcolor(FitData)
  subplot(1,3,3)
  pcolor(ErrMap)

  savefig(pwd()*"/FitCheck/"*"fit_$(N)_$(idx).png")

end

function VisualizeOneDimFit(xGrdA,xGrdB,FitData1,FitData2,FitData3,FitData4,FitData5,FitData6,idx)
  figure(figsize=(15,6))
  subplot(2,3,1);  plot(xGrdA,FitData1);
  subplot(2,3,2);  plot(xGrdA,FitData2);
  subplot(2,3,3);  plot(xGrdA,FitData3);
  subplot(2,3,4);  plot(xGrdB,FitData4);
  subplot(2,3,5);  plot(xGrdB,FitData5);
  subplot(2,3,6);  plot(xGrdB,FitData6);
  savefig(pwd()*"/FitCheck/"*"oneDim_$(idx).png");
end