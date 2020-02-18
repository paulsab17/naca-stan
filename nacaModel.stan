//Stan model for NACA system


functions{
  real[] NACAModelODE(real t,real[] y,real[] parms,real[] rdata,int[] idata){
    //y is the current state of the system. It is broken up into different vars
    //y[1] = [NACA] N-acetylcysteine amide
    //y[2] = [DTP] 4,4-aldrithiol
    //y[3] = [TP] 4-thiopyridine (-)

    real kChemEff = parms[1];
    real kBleachEff = parms[2];

    real dydt[3];

    //NACA
    dydt[1] = -kChemEff*y[1]*y[2];
    //DTP
    dydt[2] = -kChemEff*y[1]*y[2];
    //TP
    dydt[3] = kChemEff*y[1]*y[2] - kBleachEff*y[3];

    return dydt;
  }

  real[] NACAModel(real t0,real[] t,real[] init,real kChemEff,
      real kBleachEff,real[] rdata,int[] idata){

      real temp[1,3];
      real y[3];
      real p[2];
      
      p[1] = kChemEff;
      p[2] = kBleachEff;

      y = rep_array(0,3);

      if(t0 >= t[1]){
        return init;
      }
      print("Init is: ",init);
      print("parms is: ",p);
      print("times is: ",t);
      print("t0 is: ",t0);
      temp = integrate_ode_bdf(NACAModelODE, init, t0, t, p, rdata, idata);
      print("temp is",temp);

      y = to_array_1d(temp);

      //returns the value of [NACA],[DTP],and [TP]
      print("y is: ",y);
      return y;
  }

  matrix NACAModelVals(real[] time,real logkChemInt,
      real logkBleachInt,real mChem,real mBleach,real delay,
      real[] rdata,int[] idata,real initCys,real initDTP,real gdn){

    real init[3];
    matrix[size(time),3] result;
    int nt;
    real t0;
    real kChemEff;
    real kBleachEff;

    nt = size(time);
    t0 = delay; //start t=0 with some mixing delay

    init = rep_array(0,3);
    init[1] = initCys;  //initial NACA concentration
    init[2] = initDTP;  //initial DTP concentration
    init[3] = 0;        //initial TP concentration
    
    
    kChemEff = pow(10,logkChemInt + mChem * gdn);
    kBleachEff = pow(10,logkBleachInt + mBleach * gdn);
    
    print("logkChemInt: ",logkChemInt);
    print("kChemEff: ",kChemEff);
    
    for(i in 1:nt){      //go through times and calculate for each
      if(time[i]>delay){ //if this time is after delay
        //update initial conditions to what they will be at t=t[i]
        init = NACAModel(t0,time[i:i],init,kChemEff,kBleachEff,rdata,idata);
        //update starting time to t=t[i]
        t0 = time[i];
      }
      //add current state to result matrix
      for(j in 1:3){
        result[i,j] = init[j];
      }
    }
    return result; //return value for each species at each time point
  }

}


data{
  int < lower = 1 > N; // number of data points
  vector<lower=0>[N] time; // Predictor
  vector<lower=0>[N] gdn; // Predictor
  vector<lower=0>[N] absorbance; // Outcome


  int <lower=0> numTrials; //number of individual trials
  int trialStarts[numTrials+1]; //the N each trial starts at, plus the final N+1

  real < lower = 0 > initCys;
  real < lower = 0 > initDTP;

  real prior_logkChemInt;
  real priorSD_logkChemInt;
  real prior_logkBleachInt;
  real priorSD_logkBleachInt;
  real prior_mChem;
  real priorSD_mChem;
  real prior_mBleach;
  real priorSD_mBleach;
  real prior_ext;
  real priorSD_ext;
  real prior_delay;
  real priorSD_delay;
  
  int nDat_gen; //number of simulated data points to generate
}


transformed data {
  real rdata[0];
  int idata[0];
}


parameters{
  real logkChemInt; //base 10 log of kChemInt
  real logkBleachInt; //base 10 log of kBleach
  real mChem;
  real mBleach;
  real < lower = 0 > ext; //extinction coefficient
  real delay;
  real < lower = 0 > sigma; //error SD
  
  //print("Params logkChem:",logkChemInt);
}


transformed parameters {
  vector[N] predictedAbs; //predicted absorbance for each data point
  
  
  for(i in 1:numTrials){
    int start = trialStarts[numTrials];
    int last = trialStarts[numTrials+1]-1;
    int numTimes = last - start + 1;
    vector[numTimes] temp;

    matrix[numTimes,3] concs; //concentrations of each point
    
    real forceLogkChemInt = -3;
    
    print("NACAModelVals times: ",to_array_1d(time[start:last]));
    print("NMV params: kC ",logkChemInt," kB ",logkBleachInt," mC ",mChem," mB ",mBleach);
    concs = NACAModelVals(to_array_1d(time[start:last]),forceLogkChemInt,logkBleachInt,mChem,
      mBleach,delay,rdata,idata,initCys,initDTP,gdn[start]);
      
    temp = concs[,3]; //just extract concentration of TP
    temp = temp * ext; //convert to the absorbances
    for(j in start:last){
      predictedAbs[j] = temp[j-start+1];
    }
  }
}



model {
  print("Model prior: ",prior_logkChemInt);
  print("Model value: ",logkChemInt);
  
  target += normal_lpdf(logkChemInt | prior_logkChemInt,priorSD_logkChemInt);
  target += normal_lpdf(logkBleachInt | prior_logkBleachInt,priorSD_logkBleachInt);
  target += normal_lpdf(mChem | prior_mChem,priorSD_mChem);
  target += normal_lpdf(mBleach | prior_mBleach,priorSD_mBleach);
  target += normal_lpdf(ext | prior_ext,priorSD_ext);
  target += normal_lpdf(delay | prior_delay,priorSD_delay);
  
  target += normal_lpdf(absorbance | predictedAbs,sigma);
  
  /*
  logkChemInt ~ normal(prior_logkChemInt,priorSD_logkChemInt);
  logkBleachInt ~ normal(prior_logkBleachInt,priorSD_logkBleachInt);
  mChem ~ normal(prior_mChem,priorSD_mChem);
  mBleach ~ normal(prior_mBleach,priorSD_mBleach);
  ext ~ normal(prior_ext,priorSD_ext);
  delay ~ normal(prior_delay,priorSD_delay);
  sigma ~ cauchy(0,1);
  
  absorbance ~ normal(predictedAbs,sigma);
  */
}


generated quantities {
  real tMax = 1800;
  real times_gen[nDat_gen];
  matrix[nDat_gen,3] temp;
  real pred_gdn_0[nDat_gen];
  real pred_gdn_3[nDat_gen];
  real pred_gdn_6[nDat_gen];
  
  real forceLogkChemInt = -3;
  real forceLogkBleachInt = -5;
  
  print("GQ logkChemInt: ",logkChemInt);
  
  for(i in 1:nDat_gen){
    times_gen[i] = (exp2(i)/exp2(nDat_gen))*tMax;
  }
  
  temp = NACAModelVals(times_gen,forceLogkChemInt,forceLogkBleachInt,mChem,
    mBleach,delay,rdata,idata,initCys,initDTP,0);
  for(i in 1:nDat_gen){
    pred_gdn_0[i] = temp[i,3];
  }

  
  temp = NACAModelVals(times_gen,forceLogkChemInt,forceLogkBleachInt,mChem,
    mBleach,delay,rdata,idata,initCys,initDTP,3);
  for(i in 1:nDat_gen){
    pred_gdn_3[i] = temp[i,3];
  }

  temp = NACAModelVals(times_gen,forceLogkChemInt,forceLogkBleachInt,mChem,
    mBleach,delay,rdata,idata,initCys,initDTP,6);
  for(i in 1:nDat_gen){
    pred_gdn_6[i] = temp[i,3];
  }
}

