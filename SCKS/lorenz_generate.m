% functions for generating data from the Lorenz attractor model
% see Havlicek2011 and lorenz_att.m

clear all; close all;
addpath ('C:\Users\fjanoos\LocalData\my-data\research\fMRI\tools\spm8')

T = 1000;

D_x = 3; % state dimension
D_y = 1; % measurement dimension
N_params = 3;

% linear model
% zz = rand(3,3); zz = 0.5*(zz - zz'); 
% gen.params = expm( zz) ; % an orthogonal matrix
% gen.params = gen.params(:);
%f = @(x,u,theta, params)( reshape(params, D_x, D_x)*x );

% the lorenz state 3d function - x is a column vector
% params contains beta, sigma, ro
f = @(x,u,theta, params)( lorenz_att(x(1),x(2),x(3), theta(1), theta(2), theta(3) , 0.005,2) );
%f = @(x,u,theta, params)( lorenz_att(x(1),x(2),x(3), params(1), params(2), params(3) , 0.005,2) );
% lorenz observation model                        
g = @(x,u, theta, params)   ( x(1) + x(2) + x(3) );

% -- generate the data --
gen.params = [8/3;10;28]; %params for generating the lorenz model
gen.Q = 1*eye(D_x); %state noise
gen.R = 0.1*eye(D_y); %observation noise
gen.V = 0*eye(D_x); % innovation variance

gen.x = zeros( D_x, T); % the original states 
gen.y = zeros( D_y, T); % outputs
gen.q = gen.Q*randn( D_x, T);
gen.r = gen.R*randn( D_y,T);
gen.u = gen.V*randn( D_x, T); % innovations to the states

gen.x(:,1) =  [0.0, 0.1, 0.1]'; 
for t = 1:T-1
    gen.x(:,t+1) = f(gen.x(:,t), [], gen.params, gen.params) + gen.q(:,t) ;
    gen.y(:,t) = g(gen.x(:,t), [], gen.params, gen.params) + gen.r(:,t);
end
gen.y(:,T) = g(gen.x(:,T), [], gen.params, gen.params)  + gen.r(:,T);
    
plot3(gen.x(1,:), gen.x(2,:),gen.x(3,:),'-')    
 
%% -- do inversion here --
 
data.y = gen.y;
%data.x_true = X;

% model.n_params = length(gen.params);
% model.params = gen.params;
model.theta = gen.params;
model.x= [10,2,7]';

model.f = f;
model.g = g;

model.R = 0.1;

model.n_x = 3;
model.P_x = eye(3);
model.Q = eye(3)*0.1;

model.n_u = 0;
%model.P_u = [];

model.n_theta = length(gen.params);
model.P_theta = 0.1*eye(model.n_theta );
model.W =  0.1*eye(model.n_theta );

model.method = 2;  
model.smooth_flg = 1; %smooth
model.estimate_noise = 1; %no noise
model.it_max = 1;

out = SCKS( model, data );

close all; 
figure(1); 
for ii = 1:D_x
    subplot(D_x,1,ii); 
    plot( gen.x(ii,:), '--k');  hold all;
    plot( out.x_f(ii,:), '-r' ); 
    plot( out.x_s(ii,:), '-b' ); 
end

figure(2); 
plot( data.y, '--k' ); hold all; 
plot( out.y_f, '-r' ); hold all; 
plot( out.y_s, '-b' ); 

figure(3); 
for ii = 1:N_params
    subplot(N_params,1,ii); 
    plot( repmat(gen.params(ii), T,1), '--k');  hold all;
    plot( out.theta_f(ii,:), '-r' ); 
    plot( out.theta_s(ii,:), '-b' ); 
end

figure(4); clf;
plot3( gen.x(1,:), gen.x(2,:), gen.x(3,:), '--k'); hold all;
plot3( out.x_f(1,:), out.x_f(2,:), out.x_f(3,:), '-r');
plot3( out.x_s(1,:), out.x_s(2,:), out.x_s(3,:), '-b');


% test 