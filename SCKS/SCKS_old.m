function [model, Xout] = SCKS( model, data )
% Square-root Cubature Smoother solves the augmented system of Havlicek
%
%__________________________________________________________________________
% -- State Model --
% x_t = f(x_t, theta_t, u_t, t) + q_t ~ Q
% -- Measurement Model --
% y_t = g(x_t, theta_t, u_t, t) + r_t ~ R
% -- Parameter Model --
% theta_t = theta_{t-1} + w  ~ W
% -- Input Model
% u_t = u_t^obs + v_t ~ V
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
%__________________________________________________________________________
% Input fields 
% -- data --
% .Y - observations (p x T )
% .x - states (initial)
% .u - inputs (if known)
%
%  -- model --
% .f = @f(x,u, w)  (state, input, params)
% .g = @g(x,u, w)  (state, input, params)
% .xP = state error cov mtx  -- prior on states Px
% .uP = input error cov mtx 
% .wP = param error cov mtx 
%
% .pE = prior mean of params p 
% .pC = prior cov of params p
% .pP = prior noise covariances of model parametrs p
% .ip = indices of params to estimate
% .cb = parameter bounds ([lower, upper] )
%
% .Q = precision cmpnts on observation noise ?? --> should be renamed R
% .V = input noise precision (fixed) --> should remain Vinv.
% .W = precision of additive state noise -- should be renamed Qinv.
%
% .n = states dimension
% .l = ops dimension
% .m = ips dimension
%
% .Qf      = form of measarument noise cov estimate:
%                  'auto'(=default),'min','mean'
% .Alg.nN    = number of SCKF-SCKS algorithm iterations
% .Alg.Itol  = tolerance value for SCKF-SCKS convergence 
%
%
%   conditional moments of model-states - q(u)
%--------------------------------------------------------------------------
%   qU.x{1}  = Conditional expectation of hidden states (backward estimate)
%   qU.x{2}  = Conditional expectation of hidden states (forward estimate)
%   qU.v{1}  = Conditional expectation of input (backward estimate)
%   qU.v{2}  = Conditional expectation of input (forward estimate)
%   qU.r{1}  = Conditional prediction of response 
%   qU.z{1}  = Conditional prediction error 
%   qU.S{1}  = Conditional covariance: cov(x) (states - backward estimate)
%   qU.S{2}  = Conditional covariance: cov(x) (states - forward estimate)
%   qU.C{1}  = Conditional covariance: cov(u) (input - backward estimate)
%   qU.C{2}  = Conditional covariance: cov(u) (input - forward estimate)
%
% conditional moments of model-parameters - q(p)
%--------------------------------------------------------------------------
%   qP.P    = Conditional expectation
%   qP.C    = Conditional covariance
%
%      F    = negative log-likelihood
%__________________________________________________________________________


% run some checks x
%--------------------------------------------------------------------------
model specification
M  = model;
M.f  = fcnchk(M.f,'x','u','w');
M.g  = fcnchk(M.g,'x','u','w');

%% INITIALISATION:
% =========================================================================
% interpolate observation according to integration step
%--------------------------------------------------------------------------
y     = data.Y;            % observations
if size(y,1)>size(y,2)     % check the dimensions
    y = y';
end
T    = size(y,2);            % number of time points 
if M.l ~= size(y,1)
    warning('check failed');
    M.l = size(y,1);
end

% inital conditions
try 
    x= data.x
    if M.n ~= size(1,1)
        warning('check failed');
        M.n = size(x,1)
    end        
catch
    x = zeros( n, 1 );
end


u       = M.u;               % input %note this change here !
pE      = M.pE;             % model parameter prior
ip      = M.ip;            % parameter indices (must be specified or empty [] to avoid parameter estimation)


try cb  = M.cb;               catch, cb = []; end; % parameter constrains  
try tE  = M.tE.P{1};             catch, tE = []; end; % true paramters for display (if available)

% covariances (square-roots)
%--------------------------------------------------------------------------
sR      = cell(1,T);
[sR{:}] = deal(sparse(((M.V))));   % observation noise variance-- should be Q??
               
% process error covariances (square-roots)
%--------------------------------------------------------------------------
Sx      = sparse(real(spm_sqrtm(M.xP)));             %state error covariance matrix (prior)
sQ      = sparse(real(spm_sqrtm(inv(M.W))));         %variance of noise q added to state eqn.

if trace(M.uP)~=0  
    Su  = sparse(real(spm_sqrt(M.uP)));            % prior of the input variance u
    sV  = sparse(real(spm_sqrtm(inv(M.V))));       % variance of input noise v (innovations)
else
    Su  = [];
    sV  = [];
    u   = []; % no given input
end

if ~isempty(ip)  % need to estimate params
    % do svn on the parameters to decorrelate them.
    qp.u0   = spm_svd(M.pC,exp(-32));              % basis for parameters
    M.p  = size(qp.u0,2);                          % number of qp.p
    theta   = sparse(M.p,1);                       % initial deviates
    % compute cov-sqrt of parameters w
    Sw      = sparse(real(spm_sqrtm(qp.u0'*M.wP*qp.u0)));
    PB      = (qp.u0'*M.pC*qp.u0);              % parameter prior covariance in w-space
    sW      = cell(1,T-1);                      % propagated covs.
    % initialize with parameter noise variance prior
    [sW{:}] = deal(sparse(real(spm_sqrtm(qp.u0'*M.pP*qp.u0))));   
    dv      = diag(sW{1}); 
else
    Sw = [];
    sW = [];
end
  


ap1 = [1e2 1e8];
ap2 = [5e2 1e8];
% number of states, inputs and parameters:
%--------------------------------------------------------------------------
nx      = size(sQ,1);             % dim of states
nu      = size(sV,1);             % dim of inputs
nw      = size(sW{1},1);          % dim of paramters
no      = size(sR{1},1);          % dim of observations

% concatenate state vector and square-root error covariance:
%--------------------------------------------------------------------------
xc      = [x(:); u(:); theta(:)];  % initial state
xx      = zeros(nx+nu+nw,T);  % keep track of state estimates
xx(:,1) = xc;
Sc      = cell(1,T);
[Sc{:}] = deal(sparse(nx+nu+nw,nx+nu+nw));   % S^a{t} - aug state cov sqrt
Sc{1}   = blkdiag(Sx,Su,Sw);                 % initialize with priors
 
% get vector indices for components of concatenated state vector
xmask   = [ones(1,nx),ones(1,nu)*2,ones(1,nw)*3];
xind    = find(xmask==1);
uind    = find(xmask==2);
wind    = find(xmask==3);

%this is dcm specific where you have nx dims-states/node and N nodes. not
%irrelevant for me
xind2   = spm_vec((repmat([1:nx/N:length(xind)],nx/N,1) + repmat([0:1:nx/N-1]',1,N))');
xind1   = spm_vec((repmat([1:N:length(xind)],N,1) + repmat([0:1:N-1]',1,nx/N))');
clear xmask;

% Precalculate cubature points: 
%--------------------------------------------------------------------------
n          = nx + nu + nw;                % total state vector dimension
nPts       = 2*n;                         % number of cubature points
CubPtArray = sqrt(n)*[eye(n) -eye(n)];    % cubature points array

% some VB shit
iter  = 1;
RR    = repmat(diag(sR{1}),1,T-1);        %sR sqrt of Cov[r]

% augment paramter matrix by number of cubature points: -- not sure why --
pE0  = sparse(pE(:,ones(1,nPts)));
% kron a nPts x nPts with matrix parameter basis 
qp.u = kron(speye(nPts),qp.u0);

% % prepare matrix template for integration by Local linearization scheme:
% %--------------------------------------------------------------------------
% r      = [1:12];
% EXstep = (N*nPts)./r;
% EXstep = r(mod(EXstep,1)==0);
% EXstep = EXstep(end);
% EXPm   = sparse([ones(nx/N),2*ones(nx/N,1);zeros(1,nx/N+1)]);
% EXPm   = kron(speye(EXstep),EXPm);
% xt     = repmat([zeros(1,nx/N) 1]',EXstep,1);
% 
% OnesNpts = ones(1,nPts);
% XiR      = zeros(nx/N,nPts*N);
% dx       = zeros(nx,nPts);
% dx0      = zeros(nx,1);
% xPred    = zeros(n,nPts);
% yPred    = zeros(N,nPts);
% Xf       = cell(1,T-1);
% [Xf{:}]  = deal(sparse(nx+nu+nw,nPts));
% x1f      = zeros(nx+nu+nw,T-1);
% Jm       = zeros(((nx/N)^2)*N,nPts);
% % Initialize display:
%--------------------------------------------------------------------------
try M.nograph; catch, M.nograph = 1; end
if ~M.nograph
    f1 = spm_figure('Create','Graphics','SCKF-SCKS estimates');
    movegui(f1,'northeast'); drawnow;
    f2 = spm_figure('Create','Graphics','SCKF-SCKS estimates');
    movegui(f2,'northwest'); drawnow;
end



% ==================================================================
% Main iteration scheme:
% ==================================================================
% get maximum number of iterations and tolerance:
try  Itol   = M.E.Itol;   catch,  Itol   = 1e-3;      end
try  RUN    = M.E.nN;     catch,  RUN    = 60;        end

MLdiff0  = 1e-4;
mloglik0 = 0;
ML       = [];
VBrun    = RUN; 
EXEC     = 0;
t0  = tic;


UU = [];
ww = [];
qq = [];
lastwarn('');


% =========================================================================
% Iteration loop (until convergence)
% =========================================================================

dq2 = diag(sQ); % sqrt var[q]
sQ0 = sQ; 

RRR = [];

tt  =1;
for run = 1:RUN
     t1 = tic;
     mloglik    = 0;
     SSSb = [];
     SS = [];
     
    % ==================================================================
    %   Forward pass:
    % ==================================================================
    for t = 2:T,
       
                
        S = Sc{t-1};
        Xi            =  xc(:,OnesNpts) + S*CubPtArray;
        %------------------------------------------------------------------
        % PREDICTION STEP:
        %------------------------------------------------------------------
        xPred(uind,:) = Xi(uind,:);   % 
        xPred(wind,:) = Xi(wind,:);   % 
        
        % parameter constrain:
        if ~isempty(cb) && ~isempty(ip)
            xPred(wind,:) = min(cb(:,2*OnesNpts),xPred(wind,:)); % upper constrain
            xPred(wind,:) = max(cb(:,1*OnesNpts),xPred(wind,:)); % lower constrain
        end

        pE                = pE0 + (spm_unvec(qp.u*spm_vec(xPred(wind,:)),pE0));
   
        % propagation of cubature points through nonlinear function:
        %------------------------------------------------------------------

        XiR(:)        = Xi(xind1,:);
        f             = M(1).f(XiR,xPred(uind,:),pE);
 
        % integration by local-linearization scheme:
        %------------------------------------------------------------------
        dfdx           = spm_diff_all(M(1).f,XiR,xPred(uind,:),pE,1);
        [dx(:) dq]     = expmall2(dfdx,f,dt,xt,EXstep,EXPm,sQ0*sQ0',nx,Jm,xind1);
        xPred(xind,:)  = Xi(xind,:) + dx(xind2,:);
      
        % mean prediction:
        %------------------------------------------------------------------
        x1            = sum(xPred,2)/nPts;
        X0            = (xPred-x1(:,OnesNpts))/sqrt(nPts);
        Xf{t-1}       = X0;  % store for the backwards run (then no need to propaget through the nonlinear fcn)
        x1f(:,t-1)    = x1;  % store for the backwards run
        sQ0           = spm_sqrtm(dq(xind2,xind2));
        sQ(N+1:end,N+1:end) = sQ0(N+1:end,N+1:end); 
    
    
        S             = spm_qr([X0 blkdiag(sQ,sV,sW{t-1})]');
      
        Xi            = x1(:,OnesNpts) + S*CubPtArray;
        X             = (Xi-x1(:,OnesNpts))/sqrt(nPts);
        
        %------------------------------------------------------------------
        % UPDATE STEP:
        %------------------------------------------------------------------
        pE            = pE0 + spm_unvec(qp.u*spm_vec(Xi(wind,:)),pE0);
        XiR(:)        = Xi(xind1,:);
        % propagate cubature points through observation function:
        yPred(:)      = M(1).g(XiR,xPred(uind,:),pE);
        y1            = sum(yPred,2)/nPts;
        Y             = (yPred-y1(:,OnesNpts))/sqrt(nPts);

        resid         = y(:,t) - y1;         % innovations 
        RES(:,t)    = resid;
 
        for it = 1:iter,   
            % VB part - update of square-root measurement noise covarinace:
            if ~isempty(M(1).Q)
                Rtype = 'mean';
                switch(Rtype)
                    case('full')
                        sR{t-1} = spm_sqrtm(beta./alpha);
                        rr      = diag(diag(beta./alpha));
                    case('diag')
                        sR{t-1} = diag(sqrt(diag(beta./alpha)));
                        rr      = diag(diag(beta./alpha));
                    case('mean')
                        sR{t-1} = mean(sqrt(diag(beta./alpha)))*eye(no);
                        rr      = diag(diag(beta./alpha));
                end
            end
             
            SA  = spm_qr([Y sR{t-1};  X zeros(nx+nu+nw,N)]');
        
            
            Sy  = diag(diag(SA(1:no,1:no)));
            Sxy = SA(no+1:end,1:no);
            S   = SA(no+1:end,no+1:end);
            K   = Sxy/Sy;
           
            xc  = x1 + K*(resid);

            
            % VB part:
            if ~isempty(M(1).Q) 
                Xi0      = xc(:,OnesNpts) + S*CubPtArray;
                XiR(:)   = Xi0(xind1,:);
                yPred(:) = M(1).g(XiR,xPred(uind,:),pE); 
                D        = (y(:,t*OnesNpts)-yPred)/sqrt(nPts);
                beta     = beta0 + (D*D');
            end
            
        end
        
          
        if ~isempty(ip)
            d      = 1/ap1(1);
            subKG  = K(wind,:);
            try
                dv     = sqrt((1-d)*(dv.^2) + d*(diag(subKG*(subKG*(rr))')));
            catch
                dv     = sqrt((1-d)*(dv.^2) + d*(diag(subKG*(subKG*(resid*resid'))')));
            end
            sW{t}  = diag(dv);
            ap1(1) = min(ap1(1)+0.001,ap1(2));
        end
        
        d2     = 1/ap2(1);
        subKG  = K(xind,:);
        try
            dq2    = sqrt((1-d2)*(dq2.^2) + d2*(diag(subKG*(subKG*(rr))') ));
        catch
            dq2    = sqrt((1-d2)*(dq2.^2) + d2*(diag(subKG*(subKG*(resid*resid'))') ));
        end
        sQ0(N+1:end,N+1:end)   = diag((dq2(N+1:end)));
        ap2(1)                 = min(ap2(1)+0.001,ap2(2));
            
        Sc{t}    = S;
        xx(:,t)  = xc;

        
        % Maximum log-Likelihood
        %------------------------------------------------------------------
        mloglik  =  mloglik - log(2.*pi).*(no/2)- log(det(Sy*Sy'))/2 - resid'/(Sy*Sy')*resid/2; 
        %------------------------------------


        if ~isempty(M(1).Q)
            RR(:,t-1) = diag(sR{t-1});
        end
        
        % stop the loop if warning appears
        if ~isempty(lastwarn), error(lastwarn); 
          disp('ERROR')
          return;
        end
        
    end
 
    sW{1}  = sW{T};
    xxf = xx;
    Sf  = Sc;

   

    %----------------------------------------------------------------------
    % END of forward pass
    % ---------------------------------------------------------------------
    if run>4, mxw = mean(xx(wind,:),2); mxwi = mxw>(max(abs(mxw))*0.5);
        xx(wind(mxwi),end) = mxw(mxwi);
    end
