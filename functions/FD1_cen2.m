function yd = FD1_cen2(y,h)
% finite difference of y along columns, with uniform spacing h

yd=(y([2:end,1],:)-y([end,1:end-1],:))/(2*h);

end

% Copyright (C) 2014  Gudmundur Adalsteinsson
% See file LICENCE for licence and warranty details
