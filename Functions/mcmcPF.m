function [newpart, wnew, info] = mcmcPF(k, xpartminus, wminus, yobv, fDyn, fParm, hLike, hParm, gCon, gParm)

% xpartminus - particles at k-1
% wminus - weights at k-1
% yobv - observation, or measurement, e.g., y_k = h(x_k) + v_k
% fDyn - function handle for state propagation model, i.e., syste dynamics
% fDyn = @(x) f(x)
% hLike - function handle for the likelihood function 
% hLike = @(x,y) f(x,y,parm)
% gCon - funciton handle for the state constraints

% ====== Initialisation ===================================================
ct = tic; % overall calculation time;
N = length(wminus); % number of particles
n = 8; %length(fieldnames(xpartminus));

% xest = zeros(n,1);


% ========== propagation of state =========================================

% x_k+1 = f(x_k) + w_k 


% xpart = fDyn(xpartminus,fParm);

xpart = xpartminus;
%=================== update ===========================================

wupdate = hLike(xpart,yobv,hParm);

if isempty(gCon)
    wnew = wminus.*wupdate;
else
    wcon = gCon(xpart);
    wnew = wminus.*wupdate.*wcon;
end


wnew = wnew / sum(wnew);

%==========================================================================


ess = ESS(wnew);



% ===================== Resampling ========================================
% Resampling with MCMC step is based on Branko's book
% see Chapte 3 in Beyond the Kalman Filter

if ess < 0.5*N

    State = [xpart.x xpart.y xpart.z xpart.Q xpart.u xpart.phi xpart.ci xpart.cii]';

    % [State, wnew, ~] = resampling(State,wnew);

    avgState = sum(ones(n,1)*wnew'.*State,2);
    % covState = (State - avgState*ones(1,N))*diag(wnew)*(State - avgState*ones(1,N))';
    % D = cholcov(covState)';

    % assume the indepedence of differernt terms
    covPos = (State(1:3,:) - avgState(1:3)*ones(1,N))*diag(wnew)*(State(1:3,:) - avgState(1:3)*ones(1,N))';
    covQ = (State(4,:) - avgState(4)*ones(1,N))*diag(wnew)*(State(4,:) - avgState(4)*ones(1,N))';
    covWind = (State(5:6,:) - avgState(5:6)*ones(1,N))*diag(wnew)*(State(5:6,:) - avgState(5:6)*ones(1,N))';
    covDiff = (State(7:8,:) - avgState(7:8)*ones(1,N))*diag(wnew)*(State(7:8,:) - avgState(7:8)*ones(1,N))';
    

    Dpos = cholcov(covPos);
    Dq = cholcov(covQ);
    Dwind = cholcov(covWind);
    Ddiff = cholcov(covDiff);

    [wnew, index] = resamplingIndex(wnew);
    State = State(:,index);    
    wnew = wnew'; % convert to column vector


    A=(4/(n+2))^(1/(n+4)); %InstantaneousGaussian
    hopt=A*(N^(-1/(n+4)));
    
    idx = true(1,N);
    newState = State;

    for jj = 1:4
        % newState(:,idx) = State(:,idx) + hopt*D*randn(n,sum(idx));

        newState(1:3,idx) = State(1:3,idx) + hopt*Dpos*randn(3,sum(idx));

        newState(4,idx) = State(4,idx) + hopt*Dq*randn(1,sum(idx));

        newState(5:6,idx) = State(5:6,idx) + hopt*Dwind*randn(2,sum(idx));
        newState(7:8,idx) = State(7:8,idx) + hopt*Ddiff*randn(2,sum(idx));

        idx = gCon(newState)~=1;
        if sum(idx) == 0
            break;
        else 
            newState(:,idx) = State(:,idx);
        end  
    end

    newerr = newState-State;
    SIG = hopt^2*blkdiag(covPos, covQ, covWind, covDiff);
    % probnewerr= mvnpdf(newerr', zeros(1,8), SIG);
    % probx = mvnpdf(zeros(N,8), zeros(1,8), SIG);
    logratio = -0.5 * sum((newerr'/SIG)'.*newerr,1) + 0.5*sum((zeros(n,N)'/SIG)'.*zeros(n,N),1); 
    

    xupdate = hLike(State,yobv,hParm); 
    xnewupdate = hLike(newState,yobv,hParm); 

    % xprior = [xpartminus.x xpartminus.y xpartminus.z xpartminus.Q xpartminus.u xpartminus.phi xpartminus.ci xpartminus.cii]';
    % SIG = diag([fParm.x fParm.y fParm.z fParm.Q fParm.u fParm.phi fParm.ci fParm.cii].^2);
    % logratio = -0.5 * sum(((newState-xprior)'/SIG)'.*(State-xprior),1) + 0.5*sum(((State-xprior)'/SIG)'.*(State-xprior),1); 



    % alpha = xnewupdate./xupdate.*(probnewerr./probx);
    alpha = xnewupdate./xupdate.*exp(logratio)';


    mcrand = rand(N,1);
    accept = alpha>=mcrand;
    reject = alpha<mcrand;
    % sum(reject)
    newState(reject) = State(reject); 

    newpart.x = newState(1,:)';
    newpart.y = newState(2,:)';
    newpart.z = newState(3,:)';
    newpart.Q = newState(4,:)';
    newpart.u = newState(5,:)';
    newpart.phi = newState(6,:)';
    newpart.ci = newState(7,:)';
    newpart.cii = newState(8,:)';

else
    newpart = xpart;
end



% =========================================================================

time = toc(ct); % overall calculation time;

% ============= Outputs ==================================================

% xest = sum(newpart.*(ones(n,1)*wnew), 2);
xest = [];
info.ess = ess;
info.avgSampling = 0;
info.time = time;

% =======================================================================

end



function ess = ESS(w)

% sumw2 = w*w';
sumw2 = sum(w.^2);
ess = 1/sumw2;

end

