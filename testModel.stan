//Stan model to test what's going wrong

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

    if(t0 == t[1]){
      return init;
    }

    temp = integrate_ode_rk45(NACAModelODE, init, t0, t, p, rdata, idata);

    y = to_array_1d(temp);

    //returns the value of [NACA],[DTP],and [TP]
    return y;
  }
  
  matrix NACAModelVals(real[] time,real logkChemInt,real logkBleachInt,real mChem,real mBleach,real delay,real[] rdata,int[] idata,real initCys,real initDTP,real gdn){

    real init[3];
    matrix[20,3] result;
    int nt;
    real t0;
    real kChemEff;
    real kBleachEff;

    nt = 20;
    t0 = delay; //start t=0 with some mixing delay

    init = rep_array(0,3);
    init[1] = initCys;  //initial NACA concentration
    init[2] = initDTP;  //initial DTP concentration
    init[3] = 0;        //initial TP concentration

    kChemEff = logkChemInt + mChem + gdn;
    kBleachEff = logkBleachInt + mBleach + gdn;

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
  int trialStarts[numTrials]; //the N each trial starts at


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
}

transformed parameters {
  vector[N] predictedAbs; //predicted absorbance for each data point
  for(i in 1:numTrials-1){
    int start = trialStarts[i];
    int last = trialStarts[i+1];
    int numTimes = last - start + 1;

    matrix[numTimes,3] concs; //concentrations of each point

    for(j in 1:numTimes){
      concs[j,1] = 0;
      concs[j,2] = 0;
      concs[j,3] = 0;
    }
    /*
    concs = NACAModelVals(time[start:last],logkChemInt,logkBleachInt,mChem,
      mBleach,delay,rdata,idata,initCys,initDTP,gdn[start]);
      */

    for(k in 1:numTimes){
      predictedAbs[k+start-1] = concs[k,3];
    }
  }
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
}
