function y=vec2square(y)
% length(y) = N^2

N=round(sqrt(length(y)));

y=reshape(y,N,N);

end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
