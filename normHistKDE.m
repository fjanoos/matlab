% @name: normHistKDE
% @desc : normalized KDE and histogram of samples
% @params : Y - the data (1-d), xi - locations for kde (optional)
%
% @version: 1.0
% @date : 2013.09.19
% @author: fj 
% @notes: none

function [f, n, x] = normHistKDE(Y, xi)


if exist('xi','var')
    [f_,x_,u_] = ksdensity(Y,xi);    
else
    [f_,x_,u_] = ksdensity(Y);    
end
[n_,x__] = hist(Y,x_);

dx = x_(2) - x_(1);
f = f_/(sum(f_)*dx);
n = n_/ (sum(n_)*dx); 
x = x_;

    
