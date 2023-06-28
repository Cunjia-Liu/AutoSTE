function splus = fDyn(s, varargin)

% k - time steps
% s - state at k (source term)
% varargin - process noise 
% splus - state at k+1


optargin = size(varargin,2);
[~, N] = size(s.x);

if optargin == 1
    sigma = varargin{1};

    % source location
    splus.x = s.x + sigma.x*randn(N,1);
    splus.y = s.y + sigma.y*randn(N,1);
    splus.z = s.z + sigma.z*randn(N,1);

    % release rate 
    splus.Q = s.Q + sigma.Q*randn(N,1); % Q should larger than zero
    idx = find(splus.Q<0);
    while ~isempty(idx)
        splus.Q(idx,1) = s.Q(idx,1) + sigma.Q*randn(length(idx),1);
        idx = find(splus.Q<=0);
    end
    
    % Wind speed
    splus.u = s.u + sigma.u*randn(N,1);
    % Wind direction
    splus.phi = s.phi + sigma.phi*randn(N,1);
    
    % Diffusivity
    splus.ci = s.ci + sigma.ci*randn(N,1);
    idx = find(splus.ci<=0);
    while ~isempty(idx)
        splus.ci(idx,1) = s.ci(idx,1) + sigma.ci*randn(length(idx),1);
        idx = find(splus.ci<=0);
    end
    
    % Lifetime
    splus.cii = s.cii + sigma.cii*randn(N,1);
    idx = find(splus.cii<=0);
    while ~isempty(idx)
        splus.cii(idx,1) = s.cii(idx,1) + sigma.cii*randn(length(idx),1);
        idx = find(splus.cii<=0);
    end


else

    splus = s;
end





end