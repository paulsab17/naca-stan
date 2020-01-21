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

    real parms[2];
    parms[1] = kChemEff;
    parms[2] = kBleachEff;

    //create array of state values initialized all to 0
    real x[3];
    x[1] = 0;
    x[2] = 0;
    x[3] = 0;

    if(t0=t[1]){ //t is an array but only t[1] is meaningful, current time
      x = init;
    }else{
      temp = integrate_ode_rk45(NACAModelODE, init, t0,
        t, parms, rdata, idata);
      x = to_array_1d(temp);
    }

    //returns the value of [NACA],[DTP],and [TP]
    return x;
  }

  matrix NACAModelVals(real[] time,real logkChemInt,
      real logkBleachInt,real mChem,real mBleach,real delay,
      real[] rdata,int[] idata,real initCys,real initDTP,real gdn){

    real init[3];
    init[1] = initCys;  //initial NACA concentration
    init[2] = initDTP;  //initial DTP concentration
    init[3] = 0;        //initial TP concentration

    real kChemEff = pow(10,logkChemInt)+mChem*gdn;
    real kBleachEff = pow(10,logkBleachInt)+mBleach*gdn;

    int nt = size(time);
    matrix[nt,3] result;

    int t0 = delay; //start t=0 with some mixing delay

    for(i in 1:nt){ //go through times and calculate for each
      if(time[i]<=t0){ //if this time is before mixDelay, just return init
        for(j in 1:3) result [i,j] = init[j];
      }else{
        //update initial conditions to what they will be at t=t[i]
        init = NACAModel(t0,time[i:i],init,kChemEff,kBleachEff,rdata,idata);
        //update starting time to t=t[i]
        t0 = time[i];
        //add to result matrix
        for(j in 1:3) result [i,j] = init[j];
      }
    }
    return result; //return value for each species at each time point
  }

}

data{
  int < lower = 1 > N; // number of data points
  vector[N] <lower=0> time; // Predictor
  vector[N] <lower=0> gdn; // Predictor
  vector[N] <lower=0> absorbance; // Outcome

  int <lower=0> numTrials; //number of individual trials
  int trialStarts[numTrials]; //the N each trial starts at

  real <lower=0> initCys;
  real <lower=0> initDTP;

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
  real <lower=0> ext; //extinction coefficient
  real delay;
  real <lower=0> sigma; //error SD
}


transformed parameters {
  vector[0] predictedAbs; //predicted absorbance for each data point
  for(i in 1:numTrials){
    int start = trialStarts[numTrials];
    int last;
    if(i<numTrials) last = trialStarts[numTrials+1]-1 else last = N;
    int numTimes = last - start + 1;

    matrix[numTimes,3] concs; //concentrations of each point
    concs = NACAModelVals(time[start:last],logkChemInt,logkBleachInt,mChem,
      mBleach,delay,rdata,idata,initCys,initDTP,gdn[start]);
    vector[numTimes] temp;
    temp = concs[,3]; //just extract concentration of TP
    temp = temp * ext; //convert to the absorbances
    predictedAbs = append_row(predictedAbs,temp);
  }
  int checkSize = size(predictedAbs);
}


model {
  logkChemInt ~ normal(prior_logkChemInt,priorSD_logkChemInt);
  logkBleachInt ~ normal(prior_logkBleachInt,priorSD_logkBleachInt);
  mChem ~ normal(prior_mChem,priorSD_mChem);
  mBleach ~ normal(prior_mBleach,priorSD_mBleach);
  ext ~ normal(prior_ext,priorSD_ext);
  delay ~ normal(prior_delay,priorSD_delay);
  sigma ~ cauchy(0,1);

  absorbance ~ normal(predictedAbs,sigma);
}
generated quantities {
  real tMax = 1800;
  int nDat = 100;
  vector[nDat] times;
  for(i in 1:nDat){
    times[i] = (exp2(i)/exp2(nDat))*tMax;
  }
  matrix[nDat,2] pred_gdn_0;
  matrix[nDat,3] temp;
  temp = NACAModelVals(times,logkChemInt,logkBleachInt,mChem,
    mBleach,delay,rdata,idata,initCys,initDTP,0);
  for(i in 1:nDat){
    pred_gdn_0[i,1] = times[i];
    pred_gdn_0[i,2] = temp[i,3];
  }

  matrix[nDat,2] pred_gdn_3;
  temp = NACAModelVals(times,logkChemInt,logkBleachInt,mChem,
    mBleach,delay,rdata,idata,initCys,initDTP,3);
  for(i in 1:nDat){
    pred_gdn_3[i,1] = times[i];
    pred_gdn_3[i,2] = temp[i,3];
  }

  matrix[nDat,2] pred_gdn_6;
  temp = NACAModelVals(times,logkChemInt,logkBleachInt,mChem,
    mBleach,delay,rdata,idata,initCys,initDTP,6);
  for(i in 1:nDat){
    pred_gdn_6[i,1] = times[i];
    pred_gdn_6[i,2] = temp[i,3];
  }
}
