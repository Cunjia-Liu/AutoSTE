function consTrue = gCon(theta)
% check if the constraint g(x) >= 0 is satisfied


if isstruct(theta)
    gVal = ones(length(theta.Q),4);
    gVal(:,1) = theta.Q >= 0;
    gVal(:,2) = theta.u >= 0;
    gVal(:,3) = theta.ci > 0;
    gVal(:,4) = theta.cii >0 ;
    
    consTrue = prod(gVal, 2);
else
    gVal = ones(size(theta,2),4);
    gVal(:,1) = theta(4,:)' >= 0;
    gVal(:,2) = theta(5,:)' >= 0;
    gVal(:,3) = theta(7,:)' > 0;
    gVal(:,4) = theta(8,:)' > 0 ;
    
    consTrue = prod(gVal, 2);

end