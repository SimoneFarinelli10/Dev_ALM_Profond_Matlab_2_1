function out = cellfunVectorOutput(fun, cellArray)
%overloading of cellfun function for fun with horizontal numeric vector output
%different name to avoid error messages
  N = length(cellArray);
  T = size(fun(cellArray{1}),2);
  out = zeros(N, T);
  for catNum=1:N
      out(catNum,:) = fun(cellArray{catNum});
  end
end

