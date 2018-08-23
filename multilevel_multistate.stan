// Adapted from https://github.com/stan-dev/example-models/blob/master/BPA/Ch.07/cjs_c_c.stan
// Further adapted to make more efficient and easy to read
// This models is derived from section 12.3 of "Stan Modeling Language
// User's Guide and Reference Manual"

/* States:  
* 1 took A 
* 2 took B 
* 3 took C 
* 4 route unknown

* Observed route: 
* 1 A 
* 2 B 
* 3 C 
* 4 not seen
*/



// Every fish starts at 1
// user-defined functions.  Likelihood is the same for every observation.  These functions are for convenience.
functions {
  
  int last_capture(int[] y_i) {
	for (k_rev in 0:(size(y_i) - 1)) {
	  int k;
	  k = size(y_i) - k_rev;
	  if (y_i[k])
		return k;
	}
	return 0;
  }

// function to calculate probability that a fish goes uncaptured at any timestep
  vector prob_uncaptured(int T, vector detect, vector surv) {
	vector[T] chi; // vector equal to the number of timesteps is called chi
	chi[T] = 1.0; // the value of chi at time T is equal to 1 (we definitely wont' observe the fish again)
	for (t in 1:(T - 1)) { // for all time steps except the last one
	  int t_curr; 
	  int t_next;
	  t_curr = T - t; // current timestep is all the timesteps minus the ones that haven't happened
	  t_next = t_curr + 1; // next timestep is current timestep plus 1
	  chi[t_curr] = (1 - surv[t_curr]) // and chi of current timestep equals the probability of death in the current timesetp
		+ surv[t_curr] // plus the probability of surving the current timestep (on log scale)
		* (1 - detect[t_next]) // multiplied by the probability of being missed at the next timestep
		* chi[t_next]; //multiplied by chi of the next timestep.
	}
	return chi;
  }

// function to calculate individual survival probability
  real foo(int[] y_i, vector surv, vector detect, vector chi,
				int first, int last, real route_p) {
	real tar = 0;
	for (t in (first + 1):last) {
	  tar = tar + log(route_p);
	  tar = tar + log(surv[t-1]);
	  if(y_i[t])
		tar = tar + log(detect[t]);
	  else
		tar = tar + log1m(detect[t]);	  	
	}
	tar = tar + log(chi[last]);
	return(tar);
  }

  // need a second one that does not include route and chi for calculation of
  // fish where the route is unknown
  real foo2(int[] y_i, vector surv, vector detect,
				int first, int last) {
	real tar = 0;
	for (t in (first + 1):last) {
	  tar = tar + log(surv[t-1]);
	  if(y_i[t])
		tar = tar + log(detect[t]);
	  else
		tar = tar + log1m(detect[t]);	  	
	}
	return(tar);
  }

}


// data block: have to declara a variable before you can use it.  covariates go in the data block.
data {
  int<lower=1> n_rec; // n unique receivers
  int<lower=1> n_reach; 
  vector[n_reach] reachKm;  // reachKm will have the same length that n_reach elements
  int<lower=1> hyper_group[n_reach];
  // 2012
  int<lower=0> nSac12; //n ind. A
  vector[nSac12] FL_Sac12;
  int<lower=0> nYolo12; //n ind. Y
  vector[nYolo12] FL_Yolo12;
  int<lower=0> nObsSac12; //n ind. A
  int<lower=0> nObsYolo12; //n ind. Y
  int<lower=0,upper=1> Sac12[nSac12, nObsSac12]; //matrix of integer values between 0 and 1 with nrows equal to length of nSac12, and ncols equal to nObsSac12.
  int<lower=0,upper=1> Yolo12[nYolo12, nObsYolo12];
  int<lower=1> routeSac12[nObsSac12 , 3];
  int<lower=1> reachSac12[(nObsSac12-1) , 3];
  int<lower=1> routeYolo12[nObsYolo12];
  int<lower=1> reachYolo12[nObsYolo12-1];
  int<lower=1> group12[nSac12];
  // 2013
  int<lower=0> nSac13; //n ind. A
  vector[nSac13] FL_Sac13;
  int<lower=0> nYolo13; //n ind. Y
  vector[nYolo13] FL_Yolo13;
  int<lower=0> nObsSac13; //n ind. A
  int<lower=0> nObsYolo13; //n ind. Y
  int<lower=0,upper=1> Sac13[nSac13, nObsSac13];
  int<lower=0,upper=1> Yolo13[nYolo13, nObsYolo13];
  int<lower=1> routeSac13[nObsSac13 , 3];
  int<lower=1> reachSac13[(nObsSac13-1) , 3];
  int<lower=1> routeYolo13[nObsYolo13];
  int<lower=1> reachYolo13[nObsYolo13-1];
  int<lower=1> group13[nSac13];
}

// for every encounter history, need to find the last observation
transformed data {
  int<lower=1> lastSac12[nSac12];
  int<lower=1> lastYolo12[nYolo12];
  int<lower=1> lastSac13[nSac13];
  int<lower=1> lastYolo13[nYolo13];
  for (i in 1:nSac12)
	lastSac12[i] = last_capture(Sac12[i]);
  for (i in 1:nSac13)
	lastSac13[i] = last_capture(Sac13[i]);
  for (i in 1:nYolo12)
	lastYolo12[i] = last_capture(Yolo12[i]);
  for (i in 1:nYolo13)
	lastYolo13[i] = last_capture(Yolo13[i]);
}

parameters {
  vector<lower=0,upper=1>[8] base_survival;
  real<lower=0> sigma_base_survival;
  vector<lower=0,upper=1>[n_reach] survival;
  real beta_reachKm;
  real beta_FL;
  vector<lower=0, upper=1>[n_rec] detect;
  simplex[3] route_p12;
  simplex[3] route_p13;
}

transformed parameters {

  // Sacramento
  vector[nObsSac12] detSac12[3];
  vector[nObsSac13] detSac13[3];
  vector[nObsSac12-1] surSac12[3];
  vector[nObsSac13-1] surSac13[3];
  vector[nSac12] FL_offsetSac12 = beta_FL*FL_Sac12;
  vector[nSac13] FL_offsetSac13 = beta_FL*FL_Sac13;
  vector[nObsSac12] chiSac12[3];
  vector[nObsSac13] chiSac13[3];

  // Yolo
  vector[nObsYolo12] detYolo12 = detect[routeYolo12];
  vector[nObsYolo13] detYolo13 = detect[routeYolo13];
  vector[nObsYolo12-1] surYolo12 = survival[reachYolo12];
  vector[nObsYolo13-1] surYolo13 = survival[reachYolo13];
  vector[nYolo12] FL_offsetYolo12 = beta_FL*FL_Yolo12;
  vector[nYolo13] FL_offsetYolo13 = beta_FL*FL_Yolo13;
  vector[nObsYolo12] chiYolo12 = prob_uncaptured(nObsYolo12, detYolo12, surYolo12);
  vector[nObsYolo13] chiYolo13 = prob_uncaptured(nObsYolo13, detYolo13, surYolo13);

  // Loop over the 3 routes
  for(i in 1:3){
	detSac12[i] = detect[routeSac12[,i]];
	detSac13[i] = detect[routeSac13[,i]];
	surSac12[i] = survival[reachSac12[, i]];
	surSac13[i] = survival[reachSac13[, i]];

	// force detection and survival for phantoms
	/* if(i != 1){ */
	/*   surSac13[i][11] = 1.0; */
	/*   detSac13[i][12] = 0.0; */
	/* } */
	/* if(i != 2){ */
	/*   surSac12[i][5] = 1.0; */
	/*   detSac12[i][6] = 0.0; */
	/* } */

	chiSac12[i] = prob_uncaptured(nObsSac12, detSac12[i], surSac12[i]);
	chiSac13[i] = prob_uncaptured(nObsSac13, detSac13[i], surSac13[i]);
  }
  
}

model {
  // Priors
  target += normal_lpdf(sigma_base_survival | 0, 0.1);
  
  // Effect of reach length on survival
  target += normal_lpdf(survival | base_survival[hyper_group]
						+ beta_reachKm * reachKm, sigma_base_survival);

  target += beta_lpdf(base_survival | 2, 1);
  target += normal_lpdf(beta_reachKm | 0, 0.2);
  target += normal_lpdf(beta_FL | 0, 0.02); // based on a centered value of 120mm, and SD of 10mm
  target += beta_lpdf(detect | 0.5, 0.5);
  //target += beta_lpdf(survival | 0.5, 0.5);
  
  // if detected, is just the correct route p()
  // if undetected, needs to be chi[1:3] + log(route_p)
  // Do it really dumb first, then clean up
  for(i in 1:nSac12){
	// If not uncertain about the route, just increment the log_posterior normally
	if(group12[i] != 4)
	  target += foo(Sac12[i], surSac12[group12[i]] + FL_offsetSac12[i], detSac12[group12[i]],
					//chiSac12[group12[i]], precalc by group; old way
					prob_uncaptured(nObsSac12, detSac12[group12[i]], surSac12[group12[i]] + FL_offsetSac12[i]),
					
					1, lastSac12[i], route_p12[group12[i]]); // have to change how Chi is calc here
	if(group12[i] == 4){
	  //Assumes the fish went undetected before split - check this
	  target += foo2(Sac12[i], surSac12[1], detSac12[1],
					1, lastSac12[i]);
	  // Since the fish is not able to take all three routes,
	  // this is a multually exclusive p
	  target += log(
	    prob_uncaptured(nObsSac12, detSac12[1], surSac12[1] + FL_offsetSac12[i])[lastSac12[i]] * route_p12[1] +
					prob_uncaptured(nObsSac12, detSac12[2], surSac12[2] + FL_offsetSac12[i])[lastSac12[i]] * route_p12[2] +
					prob_uncaptured(nObsSac12, detSac12[3], surSac12[3] + FL_offsetSac12[i])[lastSac12[i]] * route_p12[3]); 
	}
  }

  for(i in 1:nSac13){
	// If not uncertain about the route, just increment the log_posterior normally
	if(group13[i] != 4)
	  target += foo(Sac13[i], surSac13[group13[i]] + FL_offsetSac13[i], detSac13[group13[i]],
					// chiSac13[group13[i]], // precal by group; old way
					prob_uncaptured(nObsSac13, detSac13[group13[i]], surSac13[group13[i]] + FL_offsetSac13[i]), // chi
					1, lastSac13[i], route_p13[group13[i]]);
	if(group13[i] == 4){
	  //Assumes the fish went undetected before split - check this
	  target += foo2(Sac13[i], surSac13[1], detSac13[1],
					1, lastSac13[i]);
	  // Since the fish is not able to take all three routes,
	  // this is a multually exclusive p
	  target += log(prob_uncaptured(nObsSac13, detSac13[1], surSac13[1] + FL_offsetSac13[i])[lastSac13[i]] * route_p13[1] +
					prob_uncaptured(nObsSac13, detSac13[2], surSac13[2] + FL_offsetSac13[i])[lastSac13[i]] * route_p13[2] +
					prob_uncaptured(nObsSac13, detSac13[3], surSac13[3] + FL_offsetSac13[i])[lastSac13[i]] * route_p13[3]); 
	}
  }
  
  // Yolo is much easier
  for(i in 1:nYolo12)
	target += foo(Yolo12[i], surYolo12 + FL_offsetYolo12[i], detYolo12,
				 // chiYolo12, 
				  prob_uncaptured(nObsYolo12, detYolo12, surYolo12 + FL_offsetYolo12[i]),
				  1, lastYolo12[i], 1.0);

  for(i in 1:nYolo13)
	target += foo(Yolo13[i], surYolo13 + FL_offsetYolo13[i], detYolo13,
				  // chiYolo13, 
				  prob_uncaptured(nObsYolo13, detYolo13, surYolo13 + FL_offsetYolo13[i]),
				  1, lastYolo13[i], 1.0);
  
}

generated quantities{
  vector[nObsSac12-1] pred_survA12 = exp(cumulative_sum(log(surSac12[1])));
  vector[nObsSac12-1] pred_survB12 = exp(cumulative_sum(log(surSac12[2])));
  vector[nObsSac12-1] pred_survC12 = exp(cumulative_sum(log(surSac12[3])));
  vector[nObsSac13-1] pred_survA13 = exp(cumulative_sum(log(surSac13[1])));
  vector[nObsSac13-1] pred_survB13 = exp(cumulative_sum(log(surSac13[2])));
  vector[nObsSac13-1] pred_survC13 = exp(cumulative_sum(log(surSac13[3])));
  vector[nObsYolo12-1] pred_survYolo12 = exp(cumulative_sum(log(surYolo12)));
  vector[nObsYolo13-1] pred_survYolo13 = exp(cumulative_sum(log(surYolo13)));

  // p(surv) delta = p(a) * p(take a) + ...
  
  vector[nObsSac12-1] pred_survDelta12 = pred_survA12 * route_p12[1] +
	pred_survB12 * route_p12[2] +
	pred_survC12 * route_p12[3];

  vector[nObsSac13-1] pred_survDelta13 = pred_survA13 * route_p13[1] +
	pred_survB13 * route_p13[2] +
	pred_survC13 * route_p13[3];

  vector[nSac12 + nSac13 + nYolo12 + nYolo13] log_lik =
	rep_vector(0.0, nSac12 + nSac13 + nYolo12 + nYolo13);

  for (i in 1:nSac12){
	
	if(group12[i] != 4)
	  log_lik[i] = log_lik[i] +
		foo(Sac12[i], surSac12[group12[i]], detSac12[group12[i]],
			chiSac12[group12[i]], 1, lastSac12[i], route_p12[group12[i]]);
	if(group12[i] == 4){
	  //Assumes the fish went undetected before split - check this
	  log_lik[i] = log_lik[i] +
		foo2(Sac12[i], surSac12[1], detSac12[1],
			 1, lastSac12[i]);
	  // Since the fish is not able to take all three routes,
	  // this is a multually exclusive p
	  log_lik[i] = log_lik[i] +
		log(chiSac12[1][lastSac12[i]] * route_p12[1] +
			chiSac12[2][lastSac12[i]] * route_p12[2] +
			chiSac12[3][lastSac12[i]] * route_p12[3]); 
	}
  }
  for(i in 1:nSac13){
	// If not uncertain about the route, just increment the log_posterior normally
	if(group13[i] != 4)
	  log_lik[i + nSac12] = log_lik[i + nSac12] +
		foo(Sac13[i], surSac13[group13[i]], detSac13[group13[i]],
			chiSac13[group13[i]], 1, lastSac13[i], route_p13[group13[i]]);
	if(group13[i] == 4){
	  //Assumes the fish went undetected before split - check this
	  log_lik[i+nSac12] = log_lik[i+nSac12] +
		foo2(Sac13[i], surSac13[1], detSac13[1],
			 1, lastSac13[i]);
	  // Since the fish is not able to take all three routes,
	  // this is a multually exclusive p
	  log_lik[i+nSac12] = log_lik[i+nSac12] +
		log(chiSac13[1][lastSac13[i]] * route_p13[1] +
			chiSac13[2][lastSac13[i]] * route_p13[2] +
			chiSac13[3][lastSac13[i]] * route_p13[3]); 
	}
  }
  
  // Yolo is much easier
  for(i in 1:nYolo12)
	log_lik[i+nSac12+nSac13] = log_lik[i+nSac12+nSac13] +
	  foo(Yolo12[i], surYolo12, detYolo12,
		  chiYolo12, 1, lastYolo12[i], 1.0);
  
  for(i in 1:nYolo13)
	log_lik[i+nSac12+nSac13+nYolo12] = log_lik[i+nSac12+nSac13+nYolo12] +
	  foo(Yolo13[i], surYolo13, detYolo13,
		  chiYolo13, 1, lastYolo13[i], 1.0);
  
}
