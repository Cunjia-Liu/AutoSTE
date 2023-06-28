function sensorData = sensorModel(s,pos,m)
% to genreated simulated sensor data based on the source term, the sensor
% position and sensor characteristics.

conc = plumeModel(s,pos);

% add noise 
datasize = size(conc);
error = m.sig_pct * conc .* randn(datasize); % add noise or fluctuations 

sensorData = conc + error;

% not detection if below the threshold
sensorData(sensorData<m.thresh) = 0;

% not detectin due to the mis-detection rate
sensorData(rand(datasize)<(1-m.Pd)) = 0;

end

