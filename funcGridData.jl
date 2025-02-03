function loadGridData(MatFile,idx; type = "constEcut")
  if type == "spaghettiPlot"

    path = read(MatFile, "paths"*idx);
    qGrd = read(MatFile, "qGrid"*idx);
    qItv = read(MatFile, "qItvl"*idx);
    eBot = read(MatFile, "eBots"*idx);
    eTop = read(MatFile, "eTops"*idx);
    return path, qGrd, qItv, eBot, eTop;
  
  elseif type == "spaghettiPlot_v3"

    path = read(MatFile, "paths"*idx);
    qGrd = read(MatFile, "qGrid"*idx);
    qItv = read(MatFile, "qItvl"*idx);
    eBot = read(MatFile, "eBots"*idx);
    eTop = read(MatFile, "eTops"*idx);
    qPtN = read(MatFile, "qPtsN"*idx);
    us   = read(MatFile, "us"*idx);
    vs   = read(MatFile, "vs"*idx);
    ws   = read(MatFile, "ws"*idx);
    qWid = read(MatFile, "qWid"*idx);
    return path, qGrd, qItv, qPtN, eBot, eTop, us, vs, ws, qWid;

  elseif type == "constEcut"

    xBas = read(MatFile, "xBasis"*idx);
    yBas = read(MatFile, "yBasis"*idx);
    zBas = read(MatFile, "zBasis"*idx);
    xGrd = read(MatFile, "xGrids"*idx);
    yGrd = read(MatFile, "yGrids"*idx);
    zSum = read(MatFile,   "zSum"*idx);
    eSum = read(MatFile,   "eSum"*idx);
    HKL  = read(MatFile,    "HKL"*idx);
    return xBas, yBas, zBas, xGrd, yGrd, zSum, eSum, HKL;

  else

    error("type not defined");
    return 0;
    
  end
end

function makeGridPoint(xBas,yBas,xGrd,yGrd,zBas,zSum,HKL)
  
  A = [xBas; yBas; zBas];  A = Array(transpose(A));
  uGrid = xGrd;  vGrid = yGrd;
  wWidt = range(zSum[1], zSum[2], 9);
  coeffData = [[a,b,w] for a in uGrid, b in vGrid, w in wWidt];
  qGridData = [A*q for q in coeffData];

  qGridData = dropdims(qGridData, dims = 3);
  qGridData = dropdims(qGridData, dims = 1);
  HKL = dropdims(HKL, dims = 1);

  for elem in qGridData
    elem[:] = elem[:] + HKL[:];
  end

  return qGridData;
end

function makeOneDimGrid(xBas,xGrd,yBas,ySum,zBas,zSum,HKL)
  
  A = [xBas; yBas; zBas];  A = Array(transpose(A));
  uGrid = xGrd;  
  vWidt = range(ySum[1], ySum[2], 9);
  wWidt = range(zSum[1], zSum[2], 9);
  coeffData = [[a,b,w] for a in uGrid, b in vWidt, w in wWidt];
  qGridData = [A*q for q in coeffData];

  HKL = dropdims(HKL, dims = 1);

  for elem in qGridData
    elem[:] = elem[:] + HKL[:];
  end

  return qGridData;
end

"""
  this makes spaghetti plot path,
  via using internal command of Su(n)ny program.
  due to incompatibility with Horace, additional factor is needed.
  factor = 0.98 is recommended, but it depends on your specific path.
"""
function makeSweepPath(crystl,path,qGrd,qItv; factor = 0.98)
  
  len = length(qGrd);
  density = 1/qItv;

  points2return = [];
  
  for it = 1:len
    point2start = path[it,:];
    point2end   = path[it+1,:];
    path2sweep  = [point2start, point2end]; 
    points2add, _ = reciprocal_space_path(crystl, path2sweep, density*factor);
    points2return = vcat(points2return,points2add);
  end

  return points2return;
end

"""
This function is for a specific case of spaghetti plot.
--> It is still being developed.
--> It may not work well, but this way is more realistic than the previous one.
--> qWidth is considered and averaging is done.
"""

function makeSweepPoints(path,qPtN,vs,ws,qWid; Nsample = 9)

  len_a⁺ = 1.2618;  len_c⁺ = 0.5291;

  abc = [ +1.0 +1/2 0.0 ;
           0.0 √3/2 0.0 ;
           0.0  0.0 1.0 ];

  qWidY = qWid[1]; qWidZ = qWid[2];
  
  PathNum = sum(qPtN);
  points2return = 
    Vector{Float64}[zeros(3) for _ in (1:PathNum), _ in 1:Nsample, _ in 1:Nsample];

  coeff_vs = [];
  for elem in eachrow(vs)
    elemV = [elem[1],elem[2],elem[3]];
    val2add = [-qWidY/(len_a⁺*norm(abc*elemV)) , qWidY/(len_a⁺*norm(abc*elemV))];
    push!(coeff_vs,val2add);
  end
  coeff_vs = float.(coeff_vs);
  
  coeff_ws = [];
  for elem in eachrow(ws)
    elemV = [elem[1],elem[2],elem[3]];
    val2add = [-qWidZ/(len_c⁺*norm(abc*elemV)) , qWidZ/(len_c⁺*norm(abc*elemV))];
    push!(coeff_ws,val2add);
  end
  coeff_ws = float.(coeff_ws);

  for (it,elem) in enumerate(qPtN)
    stt = sum(qPtN[1:it-1]);
    pathStt = [path[it,1]  , path[it,2]  , path[it,3]  ];  
    pathEnd = [path[it+1,1], path[it+1,2], path[it+1,3]];
    PathCenters = range(pathStt,pathEnd,Int(elem+1));
    for x in 1:Int(elem)
      pts = Int(stt) + x;
      bot_vs = coeff_vs[it][1];  top_vs = coeff_vs[it][2];  
      vSample = range(bot_vs,top_vs,Nsample);
      bot_ws = coeff_ws[it][1];  top_ws = coeff_ws[it][2];
      wSample = range(bot_ws,top_ws,Nsample);
      for (jt,y) in enumerate(vSample), (kt,z) in enumerate(wSample)
        tmp = 0.5 * ( PathCenters[x] + PathCenters[x+1] ) + 
                       [vs[it,1], vs[it,2], vs[it,3]] * y + 
                       [ws[it,1], ws[it,2], ws[it,3]] * z;
        points2return[pts,jt,kt] = tmp;
      end
    end
  end

  return points2return;
end