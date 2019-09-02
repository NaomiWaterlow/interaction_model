#include <Rcpp.h>
using namespace Rcpp;

// function to calculate the transmission of RSV
NumericVector transmission_calc_R2 (NumericVector y, NumericMatrix Contact_Structure, int num_grps, float bR)
{ //initialise output vector
  NumericVector transmission_R(num_grps*13);
  //for the age group being infected
  for(int l=0;(num_grps)>l;++l){
    //set up as 0 initially
    transmission_R[l]=0;
    //for each infecting age group
    for (int k=0;(num_grps)>k;++k){
      //wotk out relative infection
      float temp = bR * ( y[1+k*13] * Contact_Structure(l, k) + y[5+k*13] * Contact_Structure(l, k) +
                          y[8+k*13] * Contact_Structure(l, k));
      //add to previous infections of age group
      transmission_R[l]= transmission_R[l]+temp;
    };
  }; return(transmission_R);
}

// function to calculate the transmission of influenza
NumericVector transmission_calc_I2 (NumericVector y, NumericMatrix Contact_Structure, int num_grps, float bI)
{ //initialise output vector
  NumericVector transmission_I(num_grps*13);
  //for the age group being infected
  for(int a=0;(num_grps)>a;++a){
    //set up as 0 initially
    transmission_I[a]=0;
    //for each infecting age group
    for (int b=0;(num_grps)>b;++b){
      //wotk out relative infection
      float temp = bI * ( y[4+b*13] * Contact_Structure(a, b) + y[5+b*13] * Contact_Structure(a, b) +
                          y[6+b*13] * Contact_Structure(a, b));
      //add to previous infections of age group
      transmission_I[a]= transmission_I[a]+temp;
    };
  }; return(transmission_I);
}


//[[Rcpp::export]]
List derivatives (double t, NumericVector y, List parms)
{
  //read in variables
  int num_grps = parms["num_grps"];
  float bR_l = parms["bR"];
  float bR_t = exp(bR_l);
  float bI_l = parms["bI"];
  float bI = exp(bI_l);
  float sigma_both = parms["l_sig"];
  float gammaR_t =  parms["gammaR"];
  float gammaI = parms["gammaI"];
  float rho =  parms["l_rho"];
  float seedI = parms["seedI"];
  float seedR = parms["seedR"];
  float trickleI_t = parms["trickleI"];
  NumericVector RSV_Sus = parms["RSV_Sus"];
  NumericMatrix Contact_Structure = parms["Contact_Structure"];
  NumericVector ydot(num_grps*13);
  
  float trickleI = 0;
  float bR = 0;
  float gammaR = 0;
  
  // start seeding influenza at time 'seedI'
  if (t<seedI){
trickleI = 0;
  } else {trickleI = trickleI_t;}
  
  // start RSV transmission and recovery at time 'seedR'
  if (t<seedR){
    bR = 0;
    gammaR = 0;
  } else {bR = bR_t;
    gammaR = gammaR_t;}
  
// calculate the transmission matrices  
NumericVector transmission_R = transmission_calc_R2 (y, Contact_Structure, num_grps, bR);
NumericVector transmission_I = transmission_calc_I2 (y, Contact_Structure, num_grps, bI);

// per age group
  for(int i = 0; (num_grps) > i; ++i){
    
    //SS
    ydot[0+i*13] = ( - (RSV_Sus[i]*transmission_R[i]*y[0+i*13])
                    - (transmission_I[i]*y[0+i*13])
                    - trickleI );

    //IS
    ydot[1+i*13]= ( (RSV_Sus[i]*transmission_R[i]*y[0+i*13])
                  - (transmission_I[i]*sigma_both*y[1+i*13])
                  - (gammaR*y[1+i*13]) );

    //PS
    ydot[2+i*13] = ( (gammaR*y[1+i*13])
                  - (rho*y[2+i*13])
                  - (sigma_both*transmission_I[i]*y[2+i*13]) );

    //RS
    ydot[3+i*13] = ( (rho*y[2+i*13])
                  - (transmission_I[i]*y[3+i*13]) );

    //SI
    ydot[4+i*13] = ( (transmission_I[i]*y[0+i*13])
                  - (RSV_Sus[i]*transmission_R[i]*sigma_both*y[4+i*13])
                  - (gammaI *y[4+i*13]) 
                  + trickleI );
    
    //II
    ydot[5+i*13] = ( (transmission_I[i]*y[1+i*13]*sigma_both)
                  + (RSV_Sus[i]*transmission_R[i]*sigma_both*y[4+i*13])
                  - (gammaR*y[5+i*13])
                  - (gammaI*y[5+i*13]) );

    //P/RI
    ydot[6+i*13] = ( (transmission_I[i]*(sigma_both*y[2+i*13]+y[3+i*13]) )
                  + (gammaR*y[5+i*13])
                  - (gammaI*y[6+i*13]) );

    //SP
    ydot[7+i*13] = ( (gammaI*y[4+i*13])
                  - (rho*y[7+i*13])
                  - (RSV_Sus[i]*sigma_both*transmission_R[i]*y[7+i*13]) );
                  
    //IP/R
    ydot[8+i*13] = ((RSV_Sus[i]*transmission_R[i]*(sigma_both*y[7+i*13]+y[9+i*13]))
                  + (gammaI*y[5+i*13])
                  - (gammaR*y[8+i*13]) );

    //SR
    ydot[9+i*13] = ( (rho*y[7+i*13])
                   - (RSV_Sus[i]*transmission_R[i]*y[9+i*13]) );
  
         
    //RR
    ydot[10+i*13] = ((gammaR*y[8+i*13])
                   + (gammaI*y[6+i*13]) );

    //R cases
    ydot[11+i*13] = ((RSV_Sus[i]*transmission_R[i]*y[9+i*13])
                   + (RSV_Sus[i]*sigma_both*transmission_R[i]*y[7+i*13])
                   + (RSV_Sus[i]*transmission_R[i]*sigma_both*y[4+i*13])
                   + (RSV_Sus[i]*transmission_R[i]*y[0+i*13]) );

    //I cases
    ydot[12+i*13] = ((transmission_I[i]*y[3+i*13])
                   + (transmission_I[i]*sigma_both*y[2+i*13])
                   + (transmission_I[i]*y[1+i*13]*sigma_both)
                   + (transmission_I[i]*y[0+i*13]) ); 
                  
                  };
  // return the list
  return List::create(_["yout"] = ydot) ;
}


