pava <- function (x, wt=rep(1,length(x)))
{
  n <- length(x)
  if (n <= 1) return (list(estim=x,levelsets = 1))
  if (any(is.na(x)) || any(is.na(wt))) {
    stop ("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)  # Find adjacent violators
    if (!(any(viol))) break
    
    i <- min( (1:(n-1))[viol])        # Pool first pair of violators
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i+1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl]*wt[ilvl]) / sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  list( estim = x, levelsets = lvlsets)
}

# This program implements the isotonic design and the Parashar et al., 2016 design in one function;
# If pava=T, isotonic design is simulated, if pava=F, parashar design is simulated.
parashar = function(k1n,k1p,N1n,N1p,kep,Nep,kn,kp,Nn,Np,tn,tp,pava=T){
  decp = 0
  decn = 0
  pet = 0
  y1p = rbinom(1,N1p,tp) #stage1 B+ responses
  y1n = rbinom(1,N1n,tn) #stage1 B- responses
  if(pava){ #Keep the same rate
    rrp = y1p/N1p
    rrn = y1n/N1n
    resrate = pava(c(rrn,rrp),c(N1n,N1p))$estim
    y1bn = N1n * resrate[1]
    y1bp = N1p * resrate[2]
  }else{ #Keep the same rate
    y1bp=y1p
    y1bn =y1n
  }
  if(y1bn >= kn){ # Route 1 Up
    pet = 1
    decp = 1
    decn = 1
    return(c(pet=pet,decp = decp, decn= decn, N = N1n+N1p))
  }
  else if (y1bn < k1n){ #Route 3
    if(y1bp>=kep){ #Route 3 up
      pet = 1
      decp = 1
      decn = 0
      return(c(pet=pet,decp = decp, decn= decn, N = N1n+N1p))
    }else if(y1bp < k1p){ #Right
      pet = 1
      decp = 0
      decn = 0
      return(c(pet=pet,decp = decp, decn= decn, N = N1n+N1p))
    }else{
      y2bep = rbinom(1,(Nep-N1p),tp) + y1p
      if(y2bep  >= kep){ # Route 3
        decp = 1
        decn = 0
        return(c(pet=pet,decp = decp, decn= decn, N = N1n+Nep))
      }
      else{ #No go 2
        decp = 0
        decn = 0
        return(c(pet=pet,decp = decp, decn= decn, N = N1n+Nep))
      }
    }
  }
  else{
    N2p = Np-N1p
    N2n = Nn - N1n
    y2p = rbinom(1,N2p,tp) + y1p
    y2n = rbinom(1,N2n,tn) + y1n
    if(pava){ #Keep the same rate
      rn = y2n/Nn
      rp = y2p/Np
      resrate1 = pava(c(rn,rp),c(Nn,Np))$estim
      y2bn = Nn * resrate1[1]
      y2bp = Np * resrate1[2]
    }else{ #Keep the same rate
      y2bp=y2p
      y2bn =y2n
    }
    if (y2bn >= kn){
      decp = 1
      decn = 1
      return(c(pet=pet,decp = decp, decn= decn, N = Nn+Np))
    }
    else{
      if (y2bp >= kp){
        decp = 1
        decn = 0
        return(c(pet=pet,decp = decp, decn= decn, N = Nn+Np))
      }else{
        decn = 0
        decp = 0
        return(c(pet=pet,decp = decp, decn= decn, N = Nn+Np))
      }
    }
  }
}




