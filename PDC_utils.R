### from Andreas' GitHub pdc, https://github.com/brandmaier/pdc
codebook <-
    function(x,m=3,t=1, use.fast=T, normalized=T)
    {
        if (!is.numeric(x)) {
            stop("X must be of class numeric!")
        }
        
        if (t*(m-1) >= length(x)) {
            return(NA);
        }
        
        # branch to fast C implementation, if possible
        
        if (use.fast) {
            if (m==2) {
                distribution <- ( .Call("fast_codebook2",x,t) )	 	
            } else if (m==3) {
                distribution <-( .Call("fast_codebook3",x,t) )
            } else if (m==4) {
                distribution <-( .Call("fast_codebook4", x,t))
            } else if (m==5) {
                distribution <-( .Call("fast_codebook5", x,t))
            } else if (m==6) {
                distribution <-( .Call("fast_codebook6", x,t))
            } else if (m==7) {
                distribution <- ( .Call("fast_codebook7", x,t))
            }
        } 
        
        if (!use.fast | m>7) {
            
            
            # generic codeword function
            codeword.func <- codeword;
            
            distribution <- rep.int(0, factorial(m))
            to <- (length(x)-t*(m-1))
            for(i in 1:to)
            {
                
                data <- x[seq.int(i,i+t*(m-1),t)]
                
                # skip data containing NA
                if	(any(is.na(data))) { next; }
                
                # calculate permutation index
                number <- codeword.func(data, m);
                
                distribution[number] <- distribution[number] + 1
            }
            
        }
        
        if (normalized) {	
            return(distribution/sum(distribution))
        } else {
            return(distribution)
        }
    }

entropy <-
    function(dist)
    {
        dist <- dist[dist!=0]
        if (length(dist)==1) return (0);
        return(-sum(dist*log(dist))/log(length(dist)));
        
    }

codebook.entropy <-
    function(data, m, t=1)
    {
        
        ent <- entropy(codebook(data,m, t))
        
        return( ent);
    }
