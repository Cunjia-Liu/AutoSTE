function likelihood = hLikePlume(xpart, yObv, pos, m)
% The likelihood function is based on the formulation in the follwoing
% paper, but without considering the first term in eq(12), i.e., no background sensor noise
% Hutchinson, M., Liu, C., & Chen, W. H. (2019). 
% Source term estimation of a hazardous airborne release using an unmanned aerial vehicle. 
% Journal of Field Robotics, 36(4), 797-817.

conc = plumeModel(xpart, pos);
 
% m.sigma - Standard deviation of sensor noise 
sigma0 = m.thresh; % or m.sig

sigmaN = m.sig_pct*conc+m.sig; 
sigma0 = sigmaN;

if yObv <= m.thresh
    likelihood = m.Pd*1/2*(1 + erf((m.thresh-conc)./sigma0/sqrt(2))) + (1-m.Pd); 
    % probably should consider the truncated normal distribution
else
    likelihood = 1./sigmaN/sqrt(2*pi).*exp(-(yObv-conc).^2./(2*sigmaN.^2));
end



end 
