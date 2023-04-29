# Program that calculates the probabilities of rejecting null hypotheses through the three route described in the supplementary materials.
# User can loop over the ten design paraemters find the optimal design.
# tn, tp represent the true response rate for biomarker negatives (A-) and biomarker positives (A+).
# Returns probabilities of rejecting null hypotheses R1, R2 and R3, probability of early termination, PET, and expected sample size EN.

exact_pava = function(k1n,k1p,N1n,N1p,kep,Nep,kn,kp,Nn,Np,tn,tp){
  # PEN is the probability of early stopping of A- only.
  R1 = R2 = R3 = PET = EN = pen = 0
  for (x1n in 0:N1n) {
    y1p.all = k1n*(N1n+N1p)/N1n - x1n
    y1p.p.upper = y1p.all
    y1p.p.lower =  k1p*(N1n+N1p)/N1p - x1n
    PET = PET + dbinom(x1n,N1n,tn) * pbinom(min(k1p-1,y1p.p.lower-1),N1p,tp)
    EN = EN + (N1n+N1p)*dbinom(x1n,N1n,tn) * pbinom(min(k1p-1,y1p.p.lower-1),N1p,tp)
    if(x1n < k1n){ # Enrichment only
      for(x1p in min(k1p,ceiling(y1p.p.lower)):N1p){
        pen = pen + dbinom(x1n,N1n,tn) * dbinom(x1p,N1p,tp) 
        xep = kep - x1p
        R3 = R3 + dbinom(x1n,N1n,tn) * dbinom(x1p,N1p,tp) * (1-pbinom((xep-1),Nep-N1p,tp))
      }
    }
    else{ #Chance of passing both prior to PAVA
      if(ceiling(y1p.p.lower)<= floor(y1p.p.upper)){
        for (x1p in ceiling(y1p.p.lower):floor(y1p.p.upper)) {#PAVA drag down A- to go enrichment Route 3 part 1
          pen = pen + dbinom(x1n, N1n, tn) * dbinom(x1p,N1p,tp) 
          R3 = R3 + dbinom(x1n, N1n, tn) * dbinom(x1p,N1p,tp) * (1-pbinom(kep-x1p-1,Nep-N1p,tp))
        }
      }
      for(xn in min(x1n,k1n):Nn){ #If both goes to second stage
        yp.all = kn*(Nn+Np)/Nn - xn
        yp.p.upper = yp.all
        yp.p.lower =kp*(Nn+Np)/Np - xn
        if(xn >=kn){
          for(x1p in ceiling(y1p.all):N1p){ 
            R1 = R1 + dbinom(x1n,N1n,tn) * dbinom(x1p,N1p,tp) * dbinom(xn-x1n,Nn-N1n,tn) * (1-pbinom(yp.all-x1p,Np-N1p,tp)) #Route 1, done
            if(yp.p.lower < yp.p.upper){ #Chances where 
              R2 = R2 + dbinom(x1n,N1n,tn) * dbinom(x1p,N1p,tp) * dbinom(xn-x1n,Nn-N1n,tn) * (pbinom(yp.p.upper-x1p,Np-N1p,tp)-pbinom(yp.p.lower-x1p,Np-N1p,tp)) #Route 2 part 1
            }
          }
        }else{
          for(x1p in min(k1p,ceiling(y1p.p.lower)):N1p){
            R2 = R2 + dbinom(x1n,N1n,tn) * dbinom(x1p,N1p,tp) * dbinom((xn-x1n),(Nn-N1n),tn)*(1-pbinom((kp-x1p-1),Np-N1p,tp)) 
          }
        }
        
      }
    }
  }
  EN = (N1n+N1p) *PET + (N1n+Nep) * pen + (1-PET - pen) *(Nn+Np)
  return(c(R1,R2,R3,PET,EN))
}

