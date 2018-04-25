function [out] = SCKS( model, data )
% Square-root Cubature Smoother solves the augmented system of Havlicek
% Ensure spm toolbox is in path !
%__________________________________________________________________________
% -- State Model --
% x_t = f(x_t, theta_t, u_t, t) + q_t ~ Q
% -- Measurement Model --
% y_t = g(x_t, theta_t, u_t, t) + r_t ~ R
% -- Parameter Model --
% theta_t = theta_{t-1} + w  ~ W
% -- Innovations / Disturbance Model --
% u_t = u_{t-1} + v_t ~ V
%
% x_t : latent states
% u_t : total innovations
% v_t : unobserved part of the innovations (var V)
% e_t : state noise (var Q)
% y_t : measurements
% r_t : measurement noise (var R)
% theta_t : (unknown) parameters of the system
% w_t : noise in the parameters (var W)
% xa = [x, u, theta, q, v, w, r] (dim N) Augmented state
% References:
% [1] Havlicek2011
% [2] Arasaratnam2008
% [3] Arasaratnam2011
%
% NOTE - everything needs to be colum vectors (col-idx is time)
%__________________________________________________________________________
% Input fields
% -- data --
% .y - observations (p x T )
% .inputs = any other parameters that don't need to be estimated. either
%           ninputs x 1 or ninputs x T vector (e.g. observed inputs)
%
%  -- model --
% .f = @f(x,u, theta, inputs, t)  (state, innovation, inputs, t)
% .g = @g(x,u, theta, inputs, t)(state, innovation, inputs, t)  #deprecated#
%
% .parallel_g   = Parallel version of g - takes matrices x, u, theta and column
%               vector for inputs
%               textprogressbar should be insite this function
% .parallel_flag = whether to use parallel g.
% .P_x = prior Cov of states x (his xP)
% .P_u = prior Cov of inputs u (his uP)
% .P_theta = prior Cov of params noise w  (his .wP)
%
% .Q = variance of q (he used W)
% .V = variance of v
% .W = variance of w (he used Q)
% .R = variance of r (he used M.V)
%
% .lb_theta = parameter bounds ([lower, upper] )
% .ub_theta
%
% .n_x = states dimension
% .n_u = innvoations dimension
% .n_y = y-dimension
% .n_theta = theta dimension
% .n_inputs = size of known inputs
%
% Initialization values
% .x - states (initial)
% .u - innovations (initial)
% .theta = initial starting point
%
% Optimization parameters
% .lambda_w - robins munroe factor between 0 & 1  (   line 430  )
% .lambda_q - annealing factor between 0 & 1          ( Havlicek2011)
% .tol = convergence tolerance
% .it_max = max number of iterations
% .smooth_flg = 1  => smooth , else only filter
% .smooth_est = 0  => run the 'g' function again on smoothed estimates.
% .method = 1 => arasaratam's (faster, smaller) default
%           2=> havlicek's method (more accurate)
% .estimate_noise = 0/1 (estimate Q,V,W,R) matrices. off is default
% .callback(t, x, u, theta, y_hat) = a general purpose hook called after each filtering and
% smoothing step.
%__________________________________________________________________________
% Returns
%  model.F    = negative log-likelihood
%__________________________________________________________________________


% run some checks x
%--------------------------------------------------------------------------
M  = model;
f  = fcnchk(M.f,'x','u','theta', 'inputs');
if isfield(M, 'parallel_flag') & M.parallel_flag == 1
    parallel_flag = 1;
    parallel_g = M.parallel_g;
else
    parallel_flag = 0;
    g  = fcnchk(M.g,'x','u','theta', 'inputs');
end
callback = @(t, x,u,theta, inputs, y_hat)(false);  % a null callback
try
callback   = fcnchk(M.callback, 't', 'x','u','theta', 'inputs', 'y_hat');
catch err
    display('no callback defined');
end
%% INITIALISATION:
% =========================================================================
y     = data.y;            % observations
% if size(y,1)>size(y,2)     % GARBAGE LOGIC ...
%     y = y';
% end
try
    ny = M.n_y;
    y = y(1:ny,:);
catch
    ny = size(y,1);
end
T    = size(y,2);            % number of time points
try
    x = M.x(:);
    if isfield( M, 'n_x' )
        nx = M.n_x;
        x = x(1:nx);
    else
        nx = size(x,1);
    end
catch e
    nx = M.n_x;
    x = zeros( nx, 1);
    display('no initial state given')
end
% the process noise annealing multiplier
try lambda_q = M.lambda_q; catch e, lambda_q = 5e-2; end
lambda_q_min = 1e-8;
try
    u = M.u(:);
    if isfield( M, 'n_u' )
        nu = M.n_u;
        u = u(1:nu);
    else
        nu = size(u,1);
    end
catch e
    nu = M.n_u;
    u = zeros( nu, 1);
    display('no initial input given')
end
try
    theta = M.theta(:);
    if isfield( M, 'n_theta' )
        ntheta = M.n_theta;
        theta = theta(1:ntheta);
    else
        ntheta = size(theta,1);
    end
catch e
    ntheta = M.n_theta;
    theta = zeros( ntheta, 1);
    display('no initial parameter given')
end
lb_theta = -Inf(size(theta));
ub_theta = +Inf(size(theta));
try % parameter constrains
    lb_theta  = M.lb_theta;
    ub_theta  = M.ub_theta;
catch e
end;
% the robbins-munroe multiplier
try lambda_w = M.lambda_w; catch e, lambda_w = 5e-2; end
lambda_w_min = 1e-8;
try
    inputs = data.inputs; %inputs that don't need to be estimated
    n_inputs = M.n_inputs;
    if size(inputs,2) ~= 1 && size(inputs,2) ~=T
        error('misspecified inputs');
    end
catch e
    inputs = zeros(n_inputs,1);
end
smooth_flg = 0; %should we do smoothing or just filtering
try smooth_flg = M.smooth_flg ; catch e, end;
smooth_est = 0; %should we compute the smoothed y_s ? can be done in callback !
try smooth_est = M.smooth_est ; catch e, end;
% state, input and parameter covariances (square-roots)
%--------------------------------------------------------------------------
sR      = cell(1,1); % sqrt of r cov (obs noise)
sQ      = cell(1,1); % sqrt of q cov (state noise)
sV      = cell(1,1); %sqrt of v cov (input noise)
sW      = cell(1,1); %sqrt of w cov
estimate_noise = 0;
try
    if M.estimate_noise == 1
        sR      = cell(1,T);  % contrast this with T estimates !
        sQ      = cell(1,T);
        sV      = cell(1,T);
        sW      = cell(1,T);
        display('Estimating noise covariances');
        estimate_noise = 1;
    end
catch e
    display('not estimating noise covariances');
end
% replace with chol
sX      = sparse(real(mysqrt(M.P_x)));   %sqrt of x cov (state)
[sR{:}] = deal(sparse(mysqrt((M.R))));
[sQ{:}] = deal(sparse(mysqrt((M.Q))));
sU  = sparse([]); %sqrt of u cov (inputs)
% only bother if there are inputs
if nu > 0
    sU      = sparse(real(mysqrt(M.P_u))); % prior of the input variance u
    [sV{:}] = deal(sparse(mysqrt((M.V))));
end
sTheta = sparse([]);
% only bother if need to estimate params
if ntheta > 0
    sTheta = sparse(real(mysqrt(M.P_theta)));
    [sW{:}] = deal(sparse(real(mysqrt(M.W))));
end
% augmented state for (x,u,\theta)- mean and cov:
% in havlicek's method the full-aug state also includes (q,v,w)
%--------------------------------------------------------------------------
na = nx+nu+ntheta;
xa_f = zeros(na,T);
xa_s = zeros(na,T);
xa_f(:,1)      = [x; u; theta];  % should be an na x T matrix - forward pass
xa_s(:,1)      = [x; u; theta];  % backward pass - store x_{t+1|t} during forward pass

clear x u theta

sA      = cell(1,T);   %square root of the augmented state cov matrix
[sA{:}] = deal(sparse(na,na));
sA{1}   = blkdiag(sX,sU,sTheta);

% get vector indices for components of full-augmented state vector
xa_mask   = [ones(1,nx),   ones(1,nu)*2,ones(1,ntheta)*3, ...
             ones(1,nx)*4 ,ones(1,nu)*5,ones(1,ntheta)*6, ...
             ones(1,ny)*7];
x_ind    = find(xa_mask==1);
u_ind    = find(xa_mask==2);
theta_ind    = find(xa_mask==3);
q_ind    = find(xa_mask==4);
v_ind    = find(xa_mask==5);
w_ind    = find(xa_mask==6);
r_ind    = find(xa_mask==7);
% Precalculate cubature points:
%--------------------------------------------------------------------------
% the extra factor of 2 is to account for na noise variables
Nf = na;
try
    if M.method == 2 %Havlicek's method
        Nf = na + na +ny;
        display 'using Havlicek method';
    else
        display 'using Arasaratnam method';
    end
catch e
    display 'using Arasaratnam method';
end
nCub       = 2*Nf;                         % number of cubature points
cub_pts = sqrt(Nf)*[eye(Nf) -eye(Nf)];   % cubature points array
%% Main Iteration Scheme
%__________________________________________________________________________
try    max_its = M.it_max;  catch e,    max_its = 50;   end
try    tol = M.tol;         catch e,    tol = 1e-6;     end
Xa = cell(1,T); % the propagated values of augmented state cubature points
y_f = zeros(ny,T); % the predicted function values
resid_f = zeros(ny,T); % the residuals for prediction
y_s = zeros(ny,T); % the function values after smoothing - overfit potential
resid_s = zeros(ny,T); %residuals after smoothing
% loop variables - preallocated to save time
Xa_t = zeros(na,nCub);      %save augmented states xt = f(x_t-1) for each cub pt
Y_t = zeros(ny,nCub);       %predicted obs  yt = g(x_t-1)for each cub pt
Stm1 =[]; St =[];
x_tm1=[];  x_tp=[];  x_t=[];
u_tm1=[]; u_tp=[]; u_t=[];
theta_tm1=[]; theta_tp=[]; theta_t=[];
q_tm1=zeros(nx,1);  q_t=zeros(nx,1); q_tp1=zeros(nx,1);
v_tm1=zeros(nu,1); v_t=zeros(nu,1); v_tp1=zeros(nu,1);
w_tm1=zeros(ntheta,1); w_t=zeros(ntheta,1); w_tp1=zeros(ntheta,1);
r_tm1=zeros(ny,1); r_t=zeros(ny,1); r_tp1=zeros(ny,1);
xF_tm1 =[]; xF_t =[]; xF_tp1 =[];
cpts =[]; XF_t =[];
inputs_t=[];inputs_tm1=[]; inputs_tp1=[];
loglik_max = -Inf; %log likehoods
loglik_fwd = 0; %forward pass log likelihood
loglik_bck = 0; %bakcwad pass log likelihood
if parallel_flag
    x_ttm1 = zeros( nx, nCub); x_t = zeros( nx, nCub);
    u_ttm1 = zeros( nu, nCub); u_t = zeros( nu, nCub);
    theta_ttm1 = zeros( ntheta, nCub); theta_t = zeros( ntheta, nCub);
else
    x_ttm1 = zeros( nx, 1); x_t = zeros( nx, 1);
    u_ttm1 = zeros( nu, 1); u_t = zeros( nu, 1);
    theta_ttm1 = zeros( ntheta, 1); theta_t = zeros( ntheta, 1);
end
t = 1;
try inputs_ = inputs(:, t); catch e,  inputs_ = inputs;      end
callback(t, xa_f(x_ind,t), xa_f(u_ind,t), xa_f(theta_ind,t), inputs_, y_f(:,t) );
for it = 1: max_its %outer loop of forward backward passes
    display(sprintf(' F-B pass %d ... ', it));
    tic;
    loglik_fwd    = 0;
    %% Forward pass (filtering):
    %==================================================================
    textprogressbar(' Filtering ');
    for t = 2:T
        ni_tm1 = min(t-1, length(sQ) ); % t-1 if estimating noise, else 1
        ni_t = min(t, length(sR) ); %noise index for t if estimating noise, else 1
        % S_{t-1|t-1}
        if Nf > na % Havlicek's method - aug with noise also.
            Stm1 = sparse( blkdiag( sA{t-1} , sQ{ni_tm1}, sV{ni_tm1}, sW{ni_tm1}, sR{ni_tm1} ) );
        else % Arasaratam's method
            Stm1 = sA{t-1};
        end
        % cubature points
        cpts = Stm1*cub_pts;
        %inputs at time t-1 and t
        try inputs_tm1 = inputs(:, t-1); catch e,  inputs_tm1 = inputs;      end
        try  inputs_t = inputs(:, t);    catch e,  inputs_t = inputs;    end
        % -- Prediction step --
        for m = 1:nCub
            % the cubature ponts for the original x, u, theta, q, v, w,r
            x_tm1 = xa_f(x_ind,t-1) + cpts(x_ind,m);
            u_tm1 = xa_f(u_ind,t-1) + cpts(u_ind,m);
            theta_tm1 = xa_f(theta_ind,t-1) + cpts(theta_ind,m);
            if ntheta > 0
                theta_tm1 = min( max( theta_tm1, lb_theta), ub_theta) ;% constraints
            end
            if (Nf > na)  % only if we are adding noise to the aug state vector
                q_tm1 = cpts(q_ind,m);
                v_tm1 = cpts(v_ind,m);
                w_m1 = cpts(w_ind,m);
                r_tm1 = cpts(r_ind,m);
            end
            % propagation of cubature points through nonlinear function:
            x_tp = f(x_tm1, u_tm1, theta_tm1, inputs_tm1, t-1) + q_tm1;
            u_tp = u_tm1 + v_tm1;
            theta_tp = theta_tm1 + w_tm1 ; % identity transforms
            % augmented state - only for (x,u,theta)
            Xa_t(:,m) = [x_tp;  u_tp; theta_tp];

            if Nf >  na
                textprogressbar(double( (m+t*nCub)/(T*nCub)*100.0 ));
                % if adding noise to the aug state
                % reuse cubature points - Havlicek's method -
                Y_t(:,m) = g(x_tp, u_tp, theta_tp, inputs_t, t) + r_t;
            end
        end
        %store x_{t|t-1}  for backward pass
        xa_s(:,t-1) = mean(Xa_t, 2);
        %meanshifted propagated cubature points X_{t|t-1}
        % X*_{k+1|k} of Arasaratnam 2011
        Xa{t} = ( Xa_t - repmat( xa_s(:,t-1) ,1, nCub) )/sqrt(nCub);
        % S_{t|t-1} = sqrt of P_{t|t-1}
        if Nf > na      %Havlicek's method
          % sA{t} = Tria([Xa{t}]);  % keep only (x,u,theta) components.
          % not needed !!
        else           %Arasaratam's method
            % S_{x_t+1|t}
            S_tp1t = Tria([Xa{t}, blkdiag(sQ{ni_tm1}, sV{ni_tm1}, sW{ni_tm1}) ] );
            % generate new cubature points for y from S_{t|t-1} - which
            % incorporate state noise.
            cpts = S_tp1t*cub_pts;   %eqn 13 of Arasaratam2011 called X_{k+1|k}
            if parallel_flag   % a parallel call to g
                %display (['t = ', num2str(t)])
                for m = 1:nCub
                    x_ttm1(:,m)     = xa_s(x_ind,t-1) + cpts(x_ind,m);
                    u_ttm1(:,m)     = xa_s(u_ind,t-1) + cpts(u_ind,m);
                    theta_ttm1(:,m) = xa_s(theta_ind,t-1) + cpts(theta_ind,m);
                    if ntheta > 0
                        theta_ttm1(:,m) = min( max( theta_ttm1(:,m), lb_theta), ub_theta) ;  % constraints
                    end
                end
                Y_t = parallel_g(x_ttm1, u_ttm1, theta_ttm1, inputs_t, t);
            else
                for m = 1:nCub
                    textprogressbar(double( (m+t*nCub)/(T*nCub)*100.0 ));
                    x_ttm1     = xa_s(x_ind,t-1) + cpts(x_ind,m);
                    u_ttm1     = xa_s(u_ind,t-1) + cpts(u_ind,m);
                    theta_ttm1 = xa_s(theta_ind,t-1) + cpts(theta_ind,m);
                    if ntheta > 0
                        theta_ttm1 = min( max( theta_ttm1, lb_theta), ub_theta) ;  % constraints
                    end
                    Y_t(:,m) = g(x_ttm1, u_ttm1, theta_ttm1, inputs_t, t);
                end
            end
        end
        %-- update step --
        y_f(:,t) = mean(Y_t,2); %y_{t|t-1}
        Y_t = (Y_t - repmat( y_f(:,t), 1, nCub))/sqrt(nCub); %mean shifted cubature points
        resid_f(:,t) = y(:,t) - y_f(:,t);
        if Nf > na              %Havlicek's method
            Syy = Tria(Y_t);        %S_{y_t|t-1}
            Pxy = Xa{t}*Y_t';       %P_{xy_t|t-1}
            Kt  = (Pxy/Syy')/Syy;   %Kalman gain
            sA{t} = Tria( Xa{t} - Kt*Y_t ); % S_{x_t|t}
        else                   %Arasaratam's method
            Xa_tp1t = cpts/sqrt(nCub);  %these are cubature points with state-noise (eqn 13)
            TT = Tria( [Y_t sR{ni_t} ; ...
                        Xa_tp1t, zeros(na, ny)] );
            Syy = TT(1:ny,1:ny);             % sqrt of P_{yt|t-1}
            Sxy = TT(ny+1:end, 1:ny);        % sqrt of P_{xt,yt|t-1}
            Kt = Sxy/Syy;                   % Kalman gain
            sA{t} = TT(ny+1:end, ny+1:end);  % S_{x_t|t}
        end
        % filtered state x_{t|t}
        xa_f(:,t) = xa_s(:,t-1) + Kt*resid_f(:,t);
        % apply constraints !
        if ntheta > 0
            theta_f = min( max( xa_f(theta_ind,t), lb_theta), ub_theta) ;  % constraints
            if any( xa_f(theta_ind,t) ~=  theta_f )
                warning('bad solution - violates contraints')
            end
            xa_f(theta_ind,t) =  theta_f;
        end
        if estimate_noise == 1
            % update the noise variances
            % Robbins-Munroe instead of a Kalman update of W_t|t
            if ntheta > 0
                cor_fac = Kt(theta_ind, :)*resid_f(:,t);  % kalman only for theta
                new_W = mysqrt( sW{ni_tm1}*sW{ni_tm1}'*(1-lambda_w) ...
                                        + lambda_w*(cor_fac*cor_fac') );
                % diagonalize this - decorrelate parameters !
                sW{ni_t} = diag( diag(new_W) );
                % drop the correction factor - annealing step
                lambda_w = max( 1 /(1/lambda_w +0.001), lambda_w_min ) ;
            end
            % annealing update of Q_t|t
            cor_fac = Kt(x_ind, :)*resid_f(:,t);
            new_Q = mysqrt( sQ{ni_tm1}*sQ{ni_tm1}'*(1-lambda_q) ...
                                            + lambda_q*cor_fac*cor_fac' );

            sQ{ni_t} = diag( diag(new_Q) );  %diagonalize
            lambda_q = max( 1 /(1/lambda_q +0.001), lambda_q_min ) ;
            % what about input noise V_t|t?? - should we anneal ??
% %             cor_fac = Kt(u_ind, :)*resid_f(:,t);
% %             new_V = spm_sqrtm( sV{ni_tm1}*sV{ni_tm1}'*(1-lambda_q) ...
% %                                         + lambda_q*cor_fac*cor_fac' );
% %             sV{ni_t} = diag( diag(new_V) );  %diagonalize

            % what about observation noise R_t|t
        end
        callback(t, xa_f(x_ind,t), xa_f(u_ind,t), xa_f(theta_ind,t), inputs_t, y_f(:,t) );
        % compute the log-likihood of this observation
        r_t = Syy\resid_f(:,t);
        loglik_fwd = loglik_fwd - 2*log(det(Syy)) - r_t'*r_t;
        save SCKS_state.mat
    end % forward pass t
    textprogressbar(' Filtering ... Done');
    if loglik_fwd > loglik_max
        loglik_max = loglik_fwd;
    else
        warning ('your minimization has failed !');
    end
    %% Backward pass (smoothing):
    %==================================================================
    if smooth_flg
        textprogressbar(' Smoothing ');
        xa_s(:,T) = xa_f(:,T); %filtered and smoothed solution are same !
        for t = T-1:-1:1
            textprogressbar( (T-t)/T*100 );
            % noise index for t , t+1
            ni_tp1 = min(t+1, length(sQ) );
            ni_t = min(t, length(sR) );
            %inputs at time t-1 and t
            try inputs_tp1 = inputs(:, tp1); catch e,  inputs_tp1 = inputs;      end
            try inputs_t = inputs(:, t);    catch e,  inputs_t = inputs;    end
            if Nf > na % Havlicek's method - aug with noise also.
                % Xa{t+1} contains propagated cub pts X_{t+1|t}
                S_tp1t = Tria([Xa{t+1}] );                  %compute  S_SCK{x_t+1|t}
                Xa_tt = sA{t}*cub_pts(1:na,:)/sqrt(nCub);   %cubature pts sampled from P_{t|t}
                P_ttp1 = Xa_tt*Xa{t+1}';            %P_{t,t+1|t}
                Gt = (P_ttp1/S_tp1t')/S_tp1t; %P_{t,t+1|t}/P_{t+1|t}
                % store S_{t|N}
                sA{t} = Tria( [Xa_tt - Gt * Xa{t+1}, Gt * sA{t+1}] );
            else % Arasaratam's method
                Xa_tt = sA{t}*cub_pts/sqrt(nCub);   %cubature pts to replicate P_{t|t}
                %Xa{t+1} contains propagated cubature pts X*_{t+1|t}
                UU = Tria( ...
                    [Xa{t+1}, blkdiag(sQ{ni_tp1}, sV{ni_tp1}, sW{ni_tp1}); ...
                     Xa_tt,   zeros(na, na) ] );
                U11 = UU(1:na, 1:na);
                U21 = UU(na+1:end, 1:na);
                U22 = UU(na+1:end, na+1:end);
                Gt = U21/U11;
                % store S_{t|N}
                sA{t} = Tria( [U22, Gt * sA{t+1}] );
            end
            % -- apply smoothing correction --
            % x_{t+1|t}  was stored here
            x_tp1t = xa_s(:,t);
            % now store x_{t|N}
            xa_s(:,t) = xa_f(:,t) + Gt*( xa_s(:,t+1) - x_tp1t );
            % -- optional step --- compute a new y
            % the problem with this is
            %  a) it is an overfit
            %  b) doesn't always work (ill conditioned cov matrix)
            %  c) time consuming
            %  d) unnecessary !
            % generate new cubature points for y from S_{t|t-1} - which
            % incorporate state noise.
            if smooth_est
                cpts = sA{t}*cub_pts(1:na,:);   %eqn 13 of Arasaratam2011 called X_{k+1|k}
                if parallel_flag   % a parallel call to g
                    %display (['t = ', num2str(t)])
                    for m = 1:nCub
                        x_t(:,m)     = xa_s(x_ind,t) + cpts(x_ind,m);
                        u_t(:,m)     = xa_s(u_ind,t) + cpts(u_ind,m);
                        theta_t(:,m) = xa_s(theta_ind,t) + cpts(theta_ind,m);
                        if ntheta > 0
                            theta_t(:,m) = min( max( theta_t(:,m), lb_theta), ub_theta) ;  % constraints
                        end
                    end
                    Y_t = parallel_g(x_t, u_t, theta_t, inputs_t, t);
                else
                    for m = 1:nCub
                        x_t     = xa_s(x_ind,t) + cpts(x_ind,m);
                        u_t     = xa_s(u_ind,t) + cpts(u_ind,m);
                        theta_t = xa_s(theta_ind,t) + cpts(theta_ind,m);
                        if ntheta > 0
                            theta_t = min( max( theta_t, lb_theta), ub_theta) ;  % constraints
                        end
                        Y_t(:,m) = g(x_t, u_t, theta_t, inputs_t);
                    end
                end
                y_s(:,t) = mean(Y_t,2); %y_{t|t-1}
                resid_s(:,t) = y(:,t) - y_s(:,t);
                Y_t = (Y_t - repmat( y_s(:,t), 1, nCub))/sqrt(nCub); %mean shifted cubature points
                Syy = Tria(Y_t); % empirical noise precision
                % compute the log-likihood of this observation
                if ~isinf ( log(det(Syy)) )
                    r_t = Syy\resid_s(:,t);
                    loglik_bck = loglik_bck - 2*log(det(Syy)) - r_t'*r_t;
                end
            end
            callback(t, xa_s(x_ind,t), xa_s(u_ind,t), xa_s(theta_ind,t), inputs_t, y_s(:,t) );
        end
        if loglik_fwd > loglik_bck
            warning ('your smoothing has failed !');
        end
        textprogressbar(' Smoothing ... Done ');
    else % smooth_flg
       display(sprintf(' Skipping smoothing ', it));
    end
    display(sprintf(' F-B pass %d ... DONE, log_like_fwd = %g, log_like_max = %g', it, loglik_fwd, loglik_max));
    toc;
end % it
%
out.x_s = xa_s(x_ind,:); out.x_f = xa_f(x_ind,:);
out.u_s = xa_s(u_ind,:); out.u_f = xa_f(u_ind,:);
out.theta_s = xa_s(theta_ind,:); out.theta_f = xa_f(theta_ind,:);
out.resid_f = resid_f; out.resid_s = resid_s;
out.y_f = y_f; out.y_s = y_s;
out.loglik_f = loglik_max;
out.loglik_s = loglik_bck;
out.sA = sA;
out.sQ = sQ;
out.sR = sR;
out.sW = sW;
out.sV = sV;



function M = mysqrt( A)
    % replace sqrtm with this. for SPD matrices only !!
    [V, D] = eig(A); D(D<eps) = 0;
    M = V*sqrt(D);


