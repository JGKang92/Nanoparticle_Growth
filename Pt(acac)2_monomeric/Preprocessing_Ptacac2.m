%% Load data and set initial values
clear
clc
addpath(fullfile("..","functions"))
% load data containing time-dependent size and position(x, y)
nanoparticle_data = load('Monomeric_pixel.csv'); 
nanoparticle_data = nanoparticle_data(8:end,:);     % First 7 seconds does not have data
pixel_len = 0.53; % pixel length in nm
rho_s = 63.2825; % Pt monomer number density in atom/nm^-3
dt = 1; % Time interval in second
%% Convert Pixel data to radius and monomer number data
timedata = nanoparticle_data(:,1);
radidata = pixel_len.*sqrt(nanoparticle_data(:,2:3:end) ./ pi);
sigmas = (3 / (4 * pi * rho_s))^(1/3);
% Convert to monomer number: assuming the spherical geomerty
numberdata = round((radidata./sigmas).^3); 
radidata(radidata == 0) = nan;
numberdata(numberdata == 0) = nan;
[lent, ntrj] = size(radidata);
% Raw data information
meannumberdata = mean(numberdata, 2, 'omitnan');
varnumberdata = var(numberdata, 1, 2, 'omitnan');
maxnum = max(numberdata,[],2,'omitnan');
minnum = min(numberdata,[],2,'omitnan');

Preprocessing_in_situ(timedata, radidata, rho_s, dt,[], "expdata_Ptacac2.mat", ...
    trim = 200)


