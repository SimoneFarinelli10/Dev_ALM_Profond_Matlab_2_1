function S = HAC_Estimation(y, maxLag)
 
[T, D] = size(y);
 mu = sum(y)/T;

    function Omega = autocov(hh)
       Omega = zeros(D, D);
       for t = hh+1:T
            Omega = Omega + (y(t,:) - mu)' * (y(t-hh, :) - mu);  
       end    
       Omega = 1/T * Omega;
    end  

    function k = kernel(arg)
      %k = 1 - floor(arg); %Truncated
      k = 1 - arg; %Bartlett
    end    

 S = autocov(0);
 for h = 1:maxLag
    runningOmega = autocov(h); 
    S =  S + kernel(h/(maxLag+1)) * (runningOmega + runningOmega');
 end    


end