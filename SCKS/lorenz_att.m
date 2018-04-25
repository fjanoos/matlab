function [x1] = lorenz_att(x0,y0,z0,varargin)
%This function will generate the Lorenz process with prescribed time steps
%and also with specified initial conditions. If no parameters are specified
%then beta = 8/3 sigma = 10 ro = 28 h = 0.0005 steps = 10000. Thus the way to use this
%function will be [x y z] = lorenz_att(1,0,0). The variable input arguments
%should be mentioned in the order they have mentioned previously. For
%example if we want to initialize a Lorenz system with initial condition
%beta = 3 sigma = 20 ro = 90 h = 0.001 then the function can be used in the
%following way [x y z] = lorenz_att(1,0,0,3,20,90,0.001)
%The lorenz equation are as follows
% x` = sigma(y - x)
% y` = x(ro - z) - y
% z` = xy - beta*z

tot_arg = size(varargin,2) ;
switch tot_arg
    case 0
        %This is the default case
        beta = 8/3;
        sigma = 10;
        ro = 28;
        h = 0.0005;
        steps = 10000;
    case 1
        beta = cell2mat(varargin(1,1));
        sigma = 10;
        ro = 28;
        h = 0.0005;
        steps = 10000;
    case 2
         beta = cell2mat(varargin(1,1));
         sigma = cell2mat(varargin(1,2));
         ro = 28;
         h = 0.0005;
         steps = 10000;
    case 3
        beta = cell2mat(varargin(1,1));
        sigma = cell2mat(varargin(1,2));
        ro = cell2mat(varargin(1,3));
        h = 0.0005;
        steps = 10000;
    case 4
        beta = cell2mat(varargin(1,1));
        sigma = cell2mat(varargin(1,2));
        ro = cell2mat(varargin(1,3));
        h = cell2mat(varargin(1,4));
        steps = 10000;
    case 5
        beta = cell2mat(varargin(1,1));
        sigma = cell2mat(varargin(1,2));
        ro = cell2mat(varargin(1,3));
        h = cell2mat(varargin(1,4));
        steps = cell2mat(varargin(1,5));
end
%The variable h is the step step that will be used for the numerical
%integration. The numerical integrator is Runge-Kutta 4th order method.
noise_flag = 0;
[x y z] = rk4('lorenz',x0,y0,z0,beta,sigma,ro,h,steps,noise_flag);

x1 = [x(2:end);y(2:end);z(2:end)];
