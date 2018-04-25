function [varargout] = rk4(type,varargin)
%This function is the numerical ode solver using the Runge Kutta 4th order
%method.
name_of_system = type;
switch (name_of_system)
    case 'lorenz'
        %The output variables are initialized for speedy execution 
        x_fin = zeros(1,cell2mat(varargin(1,8)));
        y_fin = zeros(1,cell2mat(varargin(1,8)));
        z_fin = zeros(1,cell2mat(varargin(1,8)));
        
        x_fin(1) = cell2mat(varargin(1,1));
        y_fin(1) = cell2mat(varargin(1,2));
        z_fin(1) = cell2mat(varargin(1,3));
        %Initializing the parameters of the Lorenz process
        beta = cell2mat(varargin(1,4));
        sigma = cell2mat(varargin(1,5));
        ro = cell2mat(varargin(1,6));
        h = cell2mat(varargin(1,7));%This is the time step parameter
        iteration_counter = cell2mat(varargin(1,8));
        if (size(varargin,2) == 8)
            noise_flag = 0;
        elseif(size(varargin,2) == 9)
            noise_flag = cell2mat(varargin(1,9));
        end
        
        %The numerical integration starts from here
        for count = 2:iteration_counter
            k1 = lor_eq(x_fin(count-1),y_fin(count-1),z_fin(count-1),beta,sigma,ro);
            k2 = lor_eq(x_fin(count-1)+0.5*h*k1(1),y_fin(count-1)+0.5*h*k1(2),z_fin(count-1)+0.5*h*k1(3),beta,sigma,ro);
            k3 = lor_eq(x_fin(count-1)+0.5*h*k2(1),y_fin(count-1)+0.5*h*k2(2),z_fin(count-1)+0.5*h*k2(3),beta,sigma,ro);
            k4 = lor_eq(x_fin(count-1)+h*k3(1),y_fin(count-1)+h*k3(2),z_fin(count-1)+h*k3(3),beta,sigma,ro);
            ans_new = [x_fin(count-1);y_fin(count-1);z_fin(count-1)] + h.*(k1 + 2*k2 + 2*k3 + k4)./6;
            if (noise_flag == 1)
                ans_new = ans_new + randn(3,1);
            end
            
            x_fin(count) = ans_new(1);
            y_fin(count) = ans_new(2);
            z_fin(count) = ans_new(3);
        end
        varargout(1) = {x_fin};
        varargout(2) = {y_fin};
        varargout(3) = {z_fin};
%     otherwise
%         varargout(1) = {1};
%         varargout(2) = {2};
%         varargout(3) = {3};
end

        