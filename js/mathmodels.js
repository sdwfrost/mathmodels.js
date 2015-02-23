"use strict";

var mathmodels = (typeof exports === "undefined")?(function mathmodels() {}):(exports);
if(typeof global !== "undefined") { global.mathmodels = mathmodels; }

mathmodels.version = "0.01";

mathmodels.si = function(y0,p,t0,tf,nsteps){
  var f = function(x,y){
    return [-p[0]*y[0]*y[1],p[0]*y[0]*y[1]];
  };
  var sol = numeric.dopri(t0,tf,y0,f,1e-6,2000);
  var ix =  Array.apply(0, Array(nsteps)).map(function(e,i) { return t0+tf*(i/(nsteps-1)); });
  var iy=sol.at(ix);
  var out=iy.map(function(e,i){return {"t":ix[i],"S":e[0],"I":e[1]};});
  return out;
}


mathmodels.levins = function(y0,p,t0,tf,nsteps){
  var f = function(x,y){
    return [p[0]*y[0]*(1-y[0])-p[1]*y[0]];
  };
  var sol = numeric.dopri(t0,tf,y0,f,1e-6,2000);
  var ix =  Array.apply(0, Array(nsteps)).map(function(e,i) { return t0+tf*(i/(nsteps-1)); });
  var iy=sol.at(ix);
  var out=iy.map(function(e,i){return {"t":ix[i],"N":e[0]};});
  return out;
}


mathmodels.predprey = function(y0,p,t0,tf,nsteps){
  var f = function(x,y){
    return [p[0]*y[0]-p[1]*y[0]*y[1],p[2]*y[0]*y[1]-p[3]*y[1]];
  };
  var sol = numeric.dopri(t0,tf,y0,f,1e-6,2000);
  var ix =  Array.apply(0, Array(nsteps)).map(function(e,i) { return t0+tf*(i/(nsteps-1)); });
  var iy=sol.at(ix);
  var out=iy.map(function(e,i){return {"t":ix[i],"X":e[0],"Y":e[1]};});
  return out;
}

mathmodels.lvcompetition = function(y0,p,t0,tf,nsteps){
  var f = function(x,y){
    return [p[0]*y[0]*(1-((y[0]+p[1]*y[1])/p[2])),p[3]*y[1]*(1-((y[1]+p[4]*y[0])/p[5]))];
  };
  var sol = numeric.dopri(t0,tf,y0,f,1e-6,2000);
  var ix =  Array.apply(0, Array(nsteps)).map(function(e,i) { return t0+tf*(i/(nsteps-1)); });
  var iy=sol.at(ix);
  var out=iy.map(function(e,i){return {"t":ix[i],"X1":e[0],"X2":e[1]};});
  return out;
}

mathmodels.logisticmap = function(y0,p,nsteps){
  var y=new Array(nsteps);
  y[0]=y0;
  for(var i=1;i<nsteps;i++){
    y[i]=p[0]*y[i-1]*(1-y[i-1]);
  }
  return y;
}
