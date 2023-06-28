close all
clear all
% clc

addpath('Functions')




% Simulated source parameters
% true source
s.Q = 5; % Release rate per time unit
% source coodinates
s.x = 40;
s.y = 60;
s.z= 1;
s.u = 4; % wind speed
s.phi = 90 * pi/180; % wind direction 
s.ci = 1;
s.cii = 8;


% Create rectangular domain area and draw the plume 
xmin = 0;
xmax = 75;
ymin = 0;
ymax = 75;
zmin = 0;
zmax = 4;
domain = [xmin xmax ymin ymax zmin zmax]; % Size of search area

% Plot example dispersion from true source
fig1 = figure;
hold on;
hdl_plume = drawPlume(fig1, s, domain); 

% sensor model parameters
m.thresh = 5e-4; % sensor threshold
m.Pd = 0.7; % probability of detection
m.sig = 1e-4; % minimun sensor noise
m.sig_pct = 0.5; % the standard deviation is a percentage of the concentration level

% process noise parameters (not used)
sigma.x = 0.2;
sigma.y = 0.2;
sigma.z = 0.1;
sigma.Q = 0.2;
sigma.u = 0.2;
sigma.phi = 2*pi/180;
sigma.ci = 0.1;
sigma.cii = 0.5;


% Initialisation and parameters of the mobile sensor
StartingPosition = [2 2 4]; % Starting position [x,y,z]
moveDist = 2; % How far to move for one step


P_k = StartingPosition; % Current robot/sensor position
P_k_store = P_k;

pos.x_matrix = P_k(1);
pos.y_matrix = P_k(2);
pos.z_matrix = P_k(3);

D=[]; % store sensor readings



% initialise PF
N = 20000; % number of particles 

% Uniform prior for location
theta.x = xmin + (xmax-xmin) * rand(N,1);
theta.y = ymin + (ymax-ymin) * rand(N,1);
theta.z = 0+5*rand(N,1);

% Gamma prior for release rate Q
a = ones(N,1)*2;
b = ones(N,1)*5;
theta.Q = gamrnd(a,b);%200*rand(N,1);%

theta.u =s.u + randn(N,1)*2;%2+6*rand(N,1);%0.75+0.5*rand(N,1);0 + randn(N,1)*0.5;%
theta.phi = s.phi*0.9 + randn(N,1)*10.*pi/180;%(10 + 30*rand(N,1)).*pi/180;
theta.ci = s.ci+2*rand(N,1);%0.12+0.1*rand(N,1);
theta.cii = s.cii + 2*rand(N,1) - 2;%0.5+ 0.1*rand(N,1);

Wpnorm = ones(N,1)/N;


fig2 = figure;
hold on
preprocess(s,theta);


for i = 1:100
    
    % generate sensor data with added noise and miss-detection
    Dsim = sensorModel(s, pos, m);


    D(i)=Dsim;

    f_dyn = @fDyn;
    h_likeli = @(s, yObv, m) hLikePlume(s, yObv, pos, m);
    g_const = @gCon; 

    [theta, Wpnorm, info] = mcmcPF(i, theta, Wpnorm, Dsim, f_dyn, sigma, h_likeli, m, g_const,[]);

    
    figure(fig1)
    hold off
    drawPlume(fig1, s, domain);
    hold on
    S(i) = 5+ceil(D(i)*5e4); % size of the marker
    scatter3(theta.x,theta.y,theta.z,3,'g','filled');
    plot3(pos.x_matrix,pos.y_matrix,pos.z_matrix,'ro','MarkerFaceColor','r','MarkerSize',5);
    plot3(P_k_store(:,1),P_k_store(:,2),P_k_store(:,3),'r-');
    scatter3(P_k_store(:,1),P_k_store(:,2),P_k_store(:,3),S,'r','MarkerFaceColor','red');
    view(0,90)
    
    drawnow



    % define the action set
    ynew = [[0,moveDist,-moveDist,0] 2*[0,moveDist,-moveDist,0] 3*[0,moveDist,-moveDist,0]];
    xnew = [[moveDist,0,0,-moveDist] 2*[moveDist,0,0,-moveDist] 3*[moveDist,0,0,-moveDist]];
    znew = [0,0,0,0,0,0,0,0,0,0,0,0];
    
    % 
    Xneighbour = zeros(1,size(xnew,2));
    Yneighbour = zeros(1,size(ynew,2));
    Zneighbour = zeros(1,size(znew,2));

    
    
    Nz= 25; % down sample the source term particles (theta_i, i=1,...N) from N to Nz for generating the hypothetical measurements
    MM = 1; % the number of hypothetical measurements for each source term particle due to measurement noise
    
    % down sample the source term particles
    [~, indx_z]= resamplingIndex(Wpnorm,Nz);



    reward = zeros(1, size(xnew,2));
   
    for k = 1:size(xnew,2)

        Xneighbour(k) = pos.x_matrix+xnew(k);
        Yneighbour(k) = pos.y_matrix+ynew(k);
        Zneighbour(k) = pos.z_matrix+znew(k);
        
        if pos.x_matrix+xnew(k)<xmin || pos.x_matrix+xnew(k)>xmax || pos.y_matrix+ynew(k)<ymin || pos.y_matrix+ynew(k)>ymax || pos.z_matrix+znew(k)<zmin || pos.z_matrix+znew(k)>zmax
            reward(k)=NaN;
            continue
        end


        npos.x_matrix = pos.x_matrix+xnew(k);
        npos.y_matrix = pos.y_matrix+ynew(k);
        npos.z_matrix = pos.z_matrix+znew(k);

        infoGain=0;

        for jj = 1:Nz
            
            d.x = theta.x(indx_z(jj));
            d.y = theta.y(indx_z(jj));
            d.z = theta.z(indx_z(jj));
            d.Q = theta.Q(indx_z(jj));
            d.u = theta.u(indx_z(jj));
            d.phi = theta.phi(indx_z(jj));
            d.ci = theta.ci(indx_z(jj));
            d.cii = theta.cii(indx_z(jj));


            for jjj = 1:MM
                
                z = sensorModel(d, npos, m); % hypothetical measurements
                
                zUpate = hLikePlume(theta, z, npos, m);
                zWp = Wpnorm.*zUpate;
                zWpnorm = zWp./sum(zWp);
            
                WW = zWpnorm./Wpnorm;
                WW(WW<=0)=1;
                WW(isinf(WW))=1;
                WW(isnan(WW))=1;
                
                % Caculate the information gain 
                % comment/uncomment to choose one of those information gains
                % Note: here we used the sum rather than the averaged value

                % ---------------------------------------------------------
                % KLD
                infoGain = infoGain + (sum(zWpnorm.*log(WW)));
                %----------------------------------------------------------
                
                %----------------------------------------------------------
                % Entropy
                % infoGain = infoGain - (-sum(zWpnorm.*log2(zWpnorm+(zWpnorm==0))));
                %----------------------------------------------------------

                %----------------------------------------------------------
                % Dual control reward
                % 
                % [~, indx] = resamplingIndex(zWpnorm,round(N/5)); % downsample for quick calculation 
                % posPlus = [theta.x(indx), theta.y(indx)]';
                % posPlus_avg = mean(posPlus,2);
                % covPlus = cov(posPlus');
                % 
                % err_x = posPlus_avg(1)-npos.x_matrix;
                % err_y = posPlus_avg(2)-npos.y_matrix; 
                % 
                % infoGain = infoGain - ((err_x^2+err_y^2) + trace(covPlus));
                %----------------------------------------------------------
            end
        end

         
        reward(k) = infoGain;
        
    end


    [val,ind] = max(reward);
    
    pos.x_matrix = Xneighbour(ind);
    pos.y_matrix = Yneighbour(ind);
    pos.z_matrix = Zneighbour(ind);

    P_k = [pos.x_matrix pos.y_matrix pos.z_matrix];

    P_k_store = [P_k_store; P_k];
   


    % stop criteria
    Covar = cov(theta.x,theta.y);
    Spread = sqrt(trace(Covar));
    
    if Spread<5 % 3.5? 4?
        break
    end

end



[~, indx] = resamplingIndex(Wpnorm);

theta.x = theta.x(indx);
theta.y = theta.y(indx);
theta.z = theta.z(indx);
theta.Q = theta.Q(indx);
theta.u = theta.u(indx);
theta.phi = theta.phi(indx);
theta.ci = theta.ci(indx);
theta.cii = theta.cii(indx);

figure(fig2)
hold on
preprocess(s,theta);








