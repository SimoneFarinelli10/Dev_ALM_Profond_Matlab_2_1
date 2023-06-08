function out = relW(w)
%relW computes weights wherever possible
    den = sum(w);
    if den == 0.
       out = w;
    else
      if isnan(den)
         out = zeros(size(w)); %to be corrected; introduced to allow for zero non insurance liabilities
      else   
       out = w/den;
      end 
    end
end   