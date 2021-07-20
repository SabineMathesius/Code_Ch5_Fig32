%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code for creating IPCC AR6 Figure 5.32 
%%% CDR pulse impact: CO2 concentration, land carbon change, ocean carbon change
%%%
%%% author: Sabine Mathesius
%%% institution: Simon Fraser University
%%% e-mail: sabine_mathesius@sfu.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input data availability
% Required CMIP6 model NetCDF files are available at https://esgf-node.llnl.gov/search/cmip6/
% Required UVic ESCM data available at https://zenodo.org/record/4641434 (doi: 10.5281/zenodo.4641434)

%% MIROC 
area_o = ncread('areacello_Ofx_MIROC-ES2L_esm-piControl_r1i1p1f2_gn.nc', 'areacello');
area_a = ncread('areacella_fx_MIROC-ES2L_esm-piControl_r1i1p1f2_gn.nc', 'areacella');

netAtmosLandCO2Flux_cdr_pulse = ncread('netAtmosLandCO2Flux_Emon_MIROC-ES2L_esm-pi-cdr-pulse_r1i1p1f2_gn_186001-205912.nc', 'netAtmosLandCO2Flux',[1 1 1],[128 64 1200]);
netAtmosLandCO2Flux_CO2pulse = ncread('netAtmosLandCO2Flux_Emon_MIROC-ES2L_esm-pi-CO2pulse_r1i1p1f2_gn_186001-205912.nc', 'netAtmosLandCO2Flux',[1 1 1],[128 64 1200]);
netAtmosLandCO2Flux_piControl = ncread('netAtmosLandCO2Flux_Emon_MIROC-ES2L_esm-piControl_r1i1p1f2_gn_185001-204912.nc', 'netAtmosLandCO2Flux',[1 1 121],[128 64 1200]);

nbp_cdr_pulse = ncread('nbp_Lmon_MIROC-ES2L_esm-pi-cdr-pulse_r1i1p1f2_gn_186001-205912.nc', 'nbp',[1 1 1],[128 64 1200]);
nbp_CO2pulse = ncread('nbp_Lmon_MIROC-ES2L_esm-pi-CO2pulse_r1i1p1f2_gn_186001-205912.nc', 'nbp',[1 1 1],[128 64 1200]);
nbp_piControl = ncread('nbp_Lmon_MIROC-ES2L_esm-piControl_r1i1p1f2_gn_185001-204912.nc', 'nbp',[1 1 121],[128 64 1200]);

a2o_flux_cdr_pulse = ncread('fgco2_Oyr_MIROC-ES2L_esm-pi-cdr-pulse_r1i1p1f2_gn_1860-2059.nc', 'fgco2',[1 1 1],[360 256 100]);
a2o_flux_CO2pulse = ncread('fgco2_Oyr_MIROC-ES2L_esm-pi-CO2pulse_r1i1p1f2_gn_1860-2059.nc', 'fgco2',[1 1 1],[360 256 100]);
a2o_flux_piControl = ncread('fgco2_Oyr_MIROC-ES2L_esm-piControl_r1i1p1f2_gn_1850-2049.nc', 'fgco2',[1 1 11],[360 256 100]);

co2_CO2pulse = ncread('co2_AERmon_MIROC-ES2L_esm-pi-CO2pulse_r1i1p1f2_gn_186001-205912.nc', 'co2',[1 1 1 1],[128 64 1 1200]);
co2_cdr_pulse = ncread('co2_AERmon_MIROC-ES2L_esm-pi-cdr-pulse_r1i1p1f2_gn_186001-205912.nc', 'co2',[1 1 1 1],[128 64 1 1200]);
co2_piControl = ncread('co2_AERmon_MIROC-ES2L_esm-piControl_r1i1p1f2_gn_185001-204912.nc', 'co2',[1 1 1 121],[128 64 1 1200]);

ppm_to_pgc = 2.124;

sftlf = ncread('sftlf_fx_MIROC-ES2L_esm-pi-CO2pulse_r1i1p1f2_gn (1).nc', 'sftlf');
year_in_sec = 60 .* 60 .* 24 .* 365;  %unit: s

%% netAtmosLandCO2Flux
time_max_2 = length(netAtmosLandCO2Flux_cdr_pulse(1,1,:));
netAtmosLandCO2Flux_cdr_pulse_ts = nan(time_max_2,1);
netAtmosLandCO2Flux_CO2pulse_ts  = nan(time_max_2,1);
netAtmosLandCO2Flux_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   netAtmosLandCO2Flux_cdr_pulse_ts(t,1) = nansum(nansum(netAtmosLandCO2Flux_cdr_pulse(:,:,t).*area_a.*sftlf./100)); 
   netAtmosLandCO2Flux_CO2pulse_ts(t,1)  = nansum(nansum(netAtmosLandCO2Flux_CO2pulse(:,:,t) .*area_a.*sftlf./100)); 
   netAtmosLandCO2Flux_piControl_ts(t,1) = nansum(nansum(netAtmosLandCO2Flux_piControl(:,:,t).*area_a.*sftlf./100));
end   

% calculate PgC/yr (yearly mean C flux multiplied by 365)
netAtmosLandCO2Flux_CO2pulse_ts_yrly  = nan(100,1);
netAtmosLandCO2Flux_cdr_pulse_ts_yrly = nan(100,1);
netAtmosLandCO2Flux_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  netAtmosLandCO2Flux_cdr_pulse_ts_yrly(t) = mean(netAtmosLandCO2Flux_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  netAtmosLandCO2Flux_CO2pulse_ts_yrly(t)  = mean(netAtmosLandCO2Flux_CO2pulse_ts(m:m2))  ./ 10^12 * year_in_sec;
  netAtmosLandCO2Flux_piControl_ts_yrly(t) = mean(netAtmosLandCO2Flux_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end
Net_Land_Flux_anomaly = netAtmosLandCO2Flux_cdr_pulse_ts_yrly-netAtmosLandCO2Flux_piControl_ts_yrly;
Net_Land_Flux_anomaly_CO2pulse = netAtmosLandCO2Flux_CO2pulse_ts_yrly-netAtmosLandCO2Flux_piControl_ts_yrly;

cum_Net_Land_Flux_anomaly = Net_Land_Flux_anomaly;
cum_Net_Land_Flux_anomaly_CO2pulse = Net_Land_Flux_anomaly_CO2pulse;
for i=2:length(Net_Land_Flux_anomaly)
cum_Net_Land_Flux_anomaly_CO2pulse(i) = Net_Land_Flux_anomaly_CO2pulse(i) + cum_Net_Land_Flux_anomaly_CO2pulse(i-1);
cum_Net_Land_Flux_anomaly(i) = Net_Land_Flux_anomaly(i) + cum_Net_Land_Flux_anomaly(i-1);
end

%% nbp 
time_max_2 = length(nbp_cdr_pulse(1,1,:));
nbp_CO2pulse_ts = nan(time_max_2,1);
nbp_cdr_pulse_ts = nan(time_max_2,1);
nbp_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   nbp_CO2pulse_ts(t,1) = nansum(nansum(nbp_CO2pulse(:,:,t).*area_a.*sftlf./100)); 
   nbp_cdr_pulse_ts(t,1) = nansum(nansum(nbp_cdr_pulse(:,:,t).*area_a.*sftlf./100)); 
   nbp_piControl_ts(t,1) = nansum(nansum(nbp_piControl(:,:,t).*area_a.*sftlf./100));
end   

% (yearly mean multiplied by 365)
nbp_CO2pulse_ts_yrly = nan(100,1);
nbp_cdr_pulse_ts_yrly = nan(100,1);
nbp_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  nbp_CO2pulse_ts_yrly(t) = mean(nbp_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_cdr_pulse_ts_yrly(t) = mean(nbp_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_piControl_ts_yrly(t) = mean(nbp_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

nbp_anomaly = nbp_cdr_pulse_ts_yrly-nbp_piControl_ts_yrly;
nbp_anomaly_CO2pulse = nbp_CO2pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_nbp_anomaly_CO2pulse = nbp_anomaly_CO2pulse;
cum_nbp_anomaly = nbp_anomaly;
for i=2:length(nbp_anomaly)
cum_nbp_anomaly(i) = nbp_anomaly(i) + cum_nbp_anomaly(i-1);
cum_nbp_anomaly_CO2pulse(i) = nbp_anomaly_CO2pulse(i) + cum_nbp_anomaly_CO2pulse(i-1);
end

% nbp and netCO2flux are identical (as expected)

%% Ocean C change
time_max_2 = length(a2o_flux_cdr_pulse(1,1,:));
a2o_flux_CO2pulse_ts = nan(time_max_2,1);
a2o_flux_cdr_pulse_ts = nan(time_max_2,1);
a2o_flux_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   a2o_flux_CO2pulse_ts(t,1) = nansum(nansum(a2o_flux_CO2pulse(:,:,t).*area_o)); 
   a2o_flux_cdr_pulse_ts(t,1) = nansum(nansum(a2o_flux_cdr_pulse(:,:,t).*area_o)); 
   a2o_flux_piControl_ts(t,1) = nansum(nansum(a2o_flux_piControl(:,:,t).*area_o));
end   
a2o_flux_CO2pulse_ts_PgC = a2o_flux_CO2pulse_ts ./ 10^12 * year_in_sec; 
a2o_flux_cdr_pulse_ts_PgC = a2o_flux_cdr_pulse_ts ./ 10^12 * year_in_sec; 
a2o_flux_piControl_ts_PgC = a2o_flux_piControl_ts ./ 10^12 * year_in_sec; 

cum_Net_Land_Flux_anomaly_MIROC = cum_Net_Land_Flux_anomaly;
cum_Net_Land_Flux_anomaly_CO2pulse_MIROC = cum_Net_Land_Flux_anomaly_CO2pulse;

Ocean_C_Flux_anomaly_cdr_pulse_MIROC = a2o_flux_cdr_pulse_ts_PgC(1:100)-a2o_flux_piControl_ts_PgC(1:100);
Ocean_C_Flux_anomaly_CO2pulse_MIROC = a2o_flux_CO2pulse_ts_PgC(1:100)-a2o_flux_piControl_ts_PgC(1:100);

cum_Ocean_C_Flux_anomaly_cdr_pulse_MIROC = Ocean_C_Flux_anomaly_cdr_pulse_MIROC;
cum_Ocean_C_Flux_anomaly_CO2pulse_MIROC = Ocean_C_Flux_anomaly_CO2pulse_MIROC;

for i=2:length(Ocean_C_Flux_anomaly_cdr_pulse_MIROC)
  cum_Ocean_C_Flux_anomaly_cdr_pulse_MIROC(i) = Ocean_C_Flux_anomaly_cdr_pulse_MIROC(i) + cum_Ocean_C_Flux_anomaly_cdr_pulse_MIROC(i-1);
  cum_Ocean_C_Flux_anomaly_CO2pulse_MIROC(i) = Ocean_C_Flux_anomaly_CO2pulse_MIROC(i) + cum_Ocean_C_Flux_anomaly_CO2pulse_MIROC(i-1);
end

Net_Land_Flux_anomaly_cdr_pulse_MIROC = Net_Land_Flux_anomaly;
Net_Land_Flux_anomaly_CO2pulse_MIROC = Net_Land_Flux_anomaly_CO2pulse;

%% atmospheric CO2 concentration to mass
%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
co2_control_z1_ts = nan(1200,1);
max_length=length(co2_piControl(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_piControl(i,j,1,t))) 
            sum1 = sum1 + (co2_piControl(i,j,1,t).*area_a(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a(i,j);
         end
     end
    end
   co2_control_z1_ts(t)=sum1./sum2; 
   end
  co2_control_z1_ts_MIROC = co2_control_z1_ts;
co2_control_z1_PgC = co2_control_z1_ts * ppm_to_pgc;

%%%%%%%%%% CO2 PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_CO2pulse(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_CO2pulse(i,j,1,t))) 
            sum1 = sum1 + (co2_CO2pulse(i,j,1,t).*area_a(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a(i,j);
         end
     end
    end
   co2_CO2pulse_z1_ts(t)=sum1./sum2; 
  end
co2_CO2pulse_z1_PgC = co2_CO2pulse_z1_ts * ppm_to_pgc;

%%%%%%%%%% CDR PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
co2_cdr_pulse_z1_ts = nan(1200,1);
max_length=length(co2_cdr_pulse(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_cdr_pulse(i,j,1,t))) 
            sum1 = sum1 + (co2_cdr_pulse(i,j,1,t).*area_a(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a(i,j);
         end
     end
    end
   co2_cdr_pulse_z1_ts(t)=sum1./sum2; 
   end
  co2_cdr_pulse_z1_ts_MIROC = co2_cdr_pulse_z1_ts;
co2_cdr_pulse_z1_PgC = co2_cdr_pulse_z1_ts * ppm_to_pgc;

%%% calculate annual average
co2_PgC_cdr_pulse_z1_yr = nan(100,1);
co2_PgC_CO2pulse_z1_yr = nan(100,1);
co2_PgC_piControl_z1_yr = nan(100,1);

for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  if (m2 <= length(co2_control_z1_PgC))
  co2_PgC_piControl_z1_yr(t) = mean(co2_control_z1_PgC(m:m2)) ;
  end
  if m2<=length(co2_CO2pulse_z1_PgC) 
  co2_PgC_CO2pulse_z1_yr(t) = mean(co2_CO2pulse_z1_PgC(m:m2));
  end
  if m2<=length(co2_cdr_pulse_z1_PgC) 
      co2_PgC_cdr_pulse_z1_yr(t) = mean(co2_cdr_pulse_z1_PgC(m:m2));
  end
end

co2_PgC_atm_anomaly_CO2pulse = co2_PgC_CO2pulse_z1_yr-co2_PgC_piControl_z1_yr;
co2_PgC_atm_anomaly_cdr_pulse = co2_PgC_cdr_pulse_z1_yr-co2_PgC_piControl_z1_yr;


%% cumulative emissions
% sum of -cum_Ocean_flux, -cum_Land_flux and atmospheric change is defined as
% cumulative emissions at given point in time.

% atmospheric CO2 in PgC, year 100
atm_CO2_PgC_CO2pulse_100yr_MIROC = co2_PgC_CO2pulse_z1_yr(100)-co2_PgC_piControl_z1_yr(100);
atm_CO2_PgC_cdr_pulse_100yr_MIROC = co2_PgC_cdr_pulse_z1_yr(100)-co2_PgC_piControl_z1_yr(100);

% cumulative emissions of CO2 in PgC, year 100
cum_emissions_cdr_pulse_100yr_MIROC = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_MIROC(100))...
                               +abs(cum_Net_Land_Flux_anomaly_MIROC(100))...
                               +abs(co2_PgC_cdr_pulse_z1_yr(100)-co2_PgC_piControl_z1_yr(100));
 
cum_emissions_CO2pulse_100yr_MIROC = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_MIROC(100))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_MIROC(100))...
                               +abs(co2_PgC_CO2pulse_z1_yr(100)-co2_PgC_piControl_z1_yr(100));
 
% airborne fraction, year 100
airborne_fraction_CO2pulse_MIROC  = atm_CO2_PgC_CO2pulse_100yr_MIROC/cum_emissions_CO2pulse_100yr_MIROC;
airborne_fraction_cdr_pulse_MIROC = atm_CO2_PgC_cdr_pulse_100yr_MIROC/cum_emissions_cdr_pulse_100yr_MIROC;

%%
% atmospheric CO2 in PgC, year 100
atm_CO2_PgC_CO2pulse_100yr_MIROC = co2_PgC_CO2pulse_z1_yr(1:100)-co2_PgC_piControl_z1_yr(1:100);
atm_CO2_PgC_cdr_pulse_100yr_MIROC = co2_PgC_cdr_pulse_z1_yr(1:100)-co2_PgC_piControl_z1_yr(1:100);

% cumulative emissions of CO2 in PgC, year 100
cum_emissions_cdr_pulse_100yr_MIROC = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_MIROC(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_MIROC(1:100))...
                               +abs(co2_PgC_cdr_pulse_z1_yr(1:100)-co2_PgC_piControl_z1_yr(1:100));
 
cum_emissions_CO2pulse_100yr_MIROC = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_MIROC(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_MIROC(1:100))...
                               +abs(co2_PgC_CO2pulse_z1_yr(1:100)-co2_PgC_piControl_z1_yr(1:100));
 
%% delete variables that are not necessary for following calculations/plots
clearvars -except *ACCESS *MIROC *CanESM* *NorESM *UKESM *UVic *CanOE *_OE *GFDL
 
%% CanESM5 
area_o_CanESM5 = ncread('areacello_Ofx_CanESM5_esm-pi-cdr-pulse_r1i1p2f1_gn.nc', 'areacello');
area_a_CanESM5 = ncread('areacella_fx_CanESM5_esm-piControl_r1i1p1f1_gn.nc', 'areacella');

netAtmosLandCO2Flux_cdr_pulse_CanESM5 = ncread('netAtmosLandCO2Flux_Emon_CanESM5_esm-pi-cdr-pulse_r1i1p2f1_gn_540101-560012.nc', 'netAtmosLandCO2Flux',[1 1 1],[128 64 1200]);
netAtmosLandCO2Flux_CO2pulse_CanESM5 = ncread('netAtmosLandCO2Flux_Emon_CanESM5_esm-pi-CO2pulse_r1i1p2f1_gn_540101-560012.nc', 'netAtmosLandCO2Flux',[1 1 1],[128 64 1200]);
netAtmosLandCO2Flux_piControl_CanESM5 = ncread('netAtmosLandCO2Flux_Emon_CanESM5_esm-piControl_r1i1p1f1_gn_540101-560012.nc', 'netAtmosLandCO2Flux',[1 1 1],[128 64 1200]);

% nbp_cdr_pulse = ncread('nbp_Lmon_CanESM5_esm-pi-cdr-pulse_r1i1p2f1_gn_540101-560012.nc', 'nbp',[1 1 1],[128 64 1200]);
% nbp_CO2pulse = ncread('nbp_Lmon_CanESM5_esm-pi-CO2pulse_r1i1p2f1_gn_540101-560012.nc', 'nbp',[1 1 1],[128 64 1200]);
% nbp_piControl = ncread('nbp_Lmon_CanESM5_esm-piControl_r1i1p1f1_gn_540101-560012.nc', 'nbp',[1 1 1],[128 64 1200]);

a2o_flux_cdr_pulse_CanESM5 = ncread('fgco2_Oyr_CanESM5_esm-pi-cdr-pulse_r1i1p2f1_gn_5401-5600.nc', 'fgco2',[1 1 1],[360 291 100]);
a2o_flux_CO2pulse_CanESM5 = ncread('fgco2_Oyr_CanESM5_esm-pi-CO2pulse_r1i1p2f1_gn_5401-5600.nc', 'fgco2',[1 1 1],[360 291 100]);
a2o_flux_piControl_CanESM5 = ncread('fgco2_Oyr_CanESM5_esm-piControl_r1i1p1f1_gn_5401-5600.nc', 'fgco2',[1 1 1],[360 291 100]);

co2_CO2pulse_CanESM5 = ncread('co2_Amon_CanESM5_esm-pi-CO2pulse_r1i1p2f1_gn_540101-560012.nc', 'co2',[1 1 1 1],[128 64 1 1200]);
co2_cdr_pulse_CanESM5 = ncread('co2_Amon_CanESM5_esm-pi-cdr-pulse_r1i1p2f1_gn_540101-560012.nc', 'co2',[1 1 1 1],[128 64 1 1200]);
co2_piControl_CanESM5 = ncread('co2_Amon_CanESM5_esm-piControl_r1i1p1f1_gn_540101-560012.nc', 'co2',[1 1 1 1],[128 64 1 1200]);

sftlf_CanESM5 = ncread('sftlf_fx_CanESM5_esm-pi-CO2pulse_r1i1p2f1_gn.nc', 'sftlf');

ppm_to_pgc = 2.124;

%% Ocean C change
year_in_sec = 60 .* 60 .* 24 .* 365;  %unit: s
time_max_2 = length(a2o_flux_cdr_pulse_CanESM5(1,1,:));
a2o_flux_CO2pulse_ts = nan(time_max_2,1);
a2o_flux_cdr_pulse_ts = nan(time_max_2,1);
a2o_flux_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   a2o_flux_CO2pulse_ts(t,1) = nansum(nansum(a2o_flux_CO2pulse_CanESM5(:,:,t).*area_o_CanESM5)); 
   a2o_flux_cdr_pulse_ts(t,1) = nansum(nansum(a2o_flux_cdr_pulse_CanESM5(:,:,t).*area_o_CanESM5)); 
   a2o_flux_piControl_ts(t,1) = nansum(nansum(a2o_flux_piControl_CanESM5(:,:,t).*area_o_CanESM5));
end   
a2o_flux_CO2pulse_ts_PgC = a2o_flux_CO2pulse_ts ./ 10^12 * year_in_sec; 
a2o_flux_cdr_pulse_ts_PgC = a2o_flux_cdr_pulse_ts ./ 10^12 * year_in_sec; 
a2o_flux_piControl_ts_PgC = a2o_flux_piControl_ts ./ 10^12 * year_in_sec; 

%% Land C change
time_max_2 = length(netAtmosLandCO2Flux_cdr_pulse_CanESM5(1,1,:));
netAtmosLandCO2Flux_CO2pulse_ts = nan(time_max_2,1);
netAtmosLandCO2Flux_cdr_pulse_ts = nan(time_max_2,1);
netAtmosLandCO2Flux_piControl_ts = nan(time_max_2,1);
for(t=1:time_max_2)
   netAtmosLandCO2Flux_CO2pulse_ts(t,1) = nansum(nansum(netAtmosLandCO2Flux_CO2pulse_CanESM5(:,:,t).*area_a_CanESM5.*sftlf_CanESM5./100)); 
   netAtmosLandCO2Flux_cdr_pulse_ts(t,1) = nansum(nansum(netAtmosLandCO2Flux_cdr_pulse_CanESM5(:,:,t).*area_a_CanESM5.*sftlf_CanESM5./100)); 
   netAtmosLandCO2Flux_piControl_ts(t,1) = nansum(nansum(netAtmosLandCO2Flux_piControl_CanESM5(:,:,t).*area_a_CanESM5.*sftlf_CanESM5./100));
end   

% (yearly mean multiplied by 365)
netAtmosLandCO2Flux_CO2pulse_ts_yrly = nan(100,1);
netAtmosLandCO2Flux_cdr_pulse_ts_yrly = nan(100,1);
netAtmosLandCO2Flux_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  netAtmosLandCO2Flux_CO2pulse_ts_yrly(t) = mean(netAtmosLandCO2Flux_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  netAtmosLandCO2Flux_cdr_pulse_ts_yrly(t) = mean(netAtmosLandCO2Flux_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  netAtmosLandCO2Flux_piControl_ts_yrly(t) = mean(netAtmosLandCO2Flux_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

Net_Land_Flux_anomaly = netAtmosLandCO2Flux_cdr_pulse_ts_yrly-netAtmosLandCO2Flux_piControl_ts_yrly;

cum_Net_Land_Flux_anomaly = Net_Land_Flux_anomaly;
for i=2:length(Net_Land_Flux_anomaly)
cum_Net_Land_Flux_anomaly(i) = Net_Land_Flux_anomaly(i) + cum_Net_Land_Flux_anomaly(i-1);
end

Net_Land_Flux_anomaly_CO2pulse = netAtmosLandCO2Flux_CO2pulse_ts_yrly-netAtmosLandCO2Flux_piControl_ts_yrly;

cum_Net_Land_Flux_anomaly_CO2pulse = Net_Land_Flux_anomaly_CO2pulse;
for i=2:length(Net_Land_Flux_anomaly_CO2pulse)
cum_Net_Land_Flux_anomaly_CO2pulse(i) = Net_Land_Flux_anomaly_CO2pulse(i) + cum_Net_Land_Flux_anomaly_CO2pulse(i-1);
end

%% UPTKAE/RELEASE/AIRBORNE FRACTIONS, 100 YEARS AFTER PULSE
cum_Net_Land_Flux_anomaly_CanESM = cum_Net_Land_Flux_anomaly;
cum_Net_Land_Flux_anomaly_CO2pulse_CanESM = cum_Net_Land_Flux_anomaly_CO2pulse;

Ocean_C_Flux_anomaly_cdr_pulse_CanESM = a2o_flux_cdr_pulse_ts_PgC(1:100)-a2o_flux_piControl_ts_PgC(1:100);
Ocean_C_Flux_anomaly_CO2pulse_CanESM = a2o_flux_CO2pulse_ts_PgC(1:100)-a2o_flux_piControl_ts_PgC(1:100);

cum_Ocean_C_Flux_anomaly_cdr_pulse_CanESM = Ocean_C_Flux_anomaly_cdr_pulse_CanESM;
cum_Ocean_C_Flux_anomaly_CO2pulse_CanESM = Ocean_C_Flux_anomaly_CO2pulse_CanESM;

for i=2:length(Ocean_C_Flux_anomaly_cdr_pulse_CanESM)
  cum_Ocean_C_Flux_anomaly_cdr_pulse_CanESM(i) = Ocean_C_Flux_anomaly_cdr_pulse_CanESM(i) + cum_Ocean_C_Flux_anomaly_cdr_pulse_CanESM(i-1);
  cum_Ocean_C_Flux_anomaly_CO2pulse_CanESM(i) = Ocean_C_Flux_anomaly_CO2pulse_CanESM(i) + cum_Ocean_C_Flux_anomaly_CO2pulse_CanESM(i-1);
end

Net_Land_Flux_anomaly_cdr_pulse_CanESM = Net_Land_Flux_anomaly;
Net_Land_Flux_anomaly_CO2pulse_CanESM = Net_Land_Flux_anomaly_CO2pulse;

%% nbp (confirmation that nbp gives identical result as variable netAtmosLandCO2Flux)
% time_max_2 = length(nbp_cdr_pulse(1,1,:))
% nbp_CO2pulse_ts = nan(time_max_2,1);
% nbp_cdr_pulse_ts = nan(time_max_2,1);
% nbp_piControl_ts = nan(time_max_2,1);
% 
% for(t=1:time_max_2)
%    nbp_CO2pulse_ts(t,1) = nansum(nansum(nbp_CO2pulse(:,:,t).*area_a.*sftlf./100)); 
%    nbp_cdr_pulse_ts(t,1) = nansum(nansum(nbp_cdr_pulse(:,:,t).*area_a.*sftlf./100)); 
%    nbp_piControl_ts(t,1) = nansum(nansum(nbp_piControl(:,:,t).*area_a.*sftlf./100));
% end   
% 
% % (yearly mean multiplied by 365)
% nbp_CO2pulse_ts_yrly = nan(100,1);
% nbp_cdr_pulse_ts_yrly = nan(100,1);
% nbp_piControl_ts_yrly = nan(100,1);
% for(t=1:100)
%     m = (t-1)*12+1;
%     m2 = m + 11;
%   nbp_CO2pulse_ts_yrly(t) = mean(nbp_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
%   nbp_cdr_pulse_ts_yrly(t) = mean(nbp_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
%   nbp_piControl_ts_yrly(t) = mean(nbp_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
% end
% 
% Net_Land_Flux2_anomaly = nbp_cdr_pulse_ts_yrly-nbp_piControl_ts_yrly;
% cLand_anomaly = cLand_cdr_pulse_ts_PgC - cLand_piControl_ts_PgC
% 
% cum_Net_Land_Flux2_anomaly = Net_Land_Flux2_anomaly
% for i=2:length(Net_Land_Flux2_anomaly)
% cum_Net_Land_Flux2_anomaly(i) = Net_Land_Flux2_anomaly(i) + cum_Net_Land_Flux2_anomaly(i-1)
% end
%  
% Net_Land_Flux2_anomaly_CO2pulse = nbp_CO2pulse_ts_yrly-nbp_piControl_ts_yrly;
% cLand_anomaly_CO2pulse = cLand_CO2pulse_ts_PgC - cLand_piControl_ts_PgC
% 
% cum_Net_Land_Flux2_anomaly_CO2pulse = Net_Land_Flux2_anomaly_CO2pulse
% for i=2:length(Net_Land_Flux2_anomaly_CO2pulse)
% cum_Net_Land_Flux2_anomaly_CO2pulse(i) = Net_Land_Flux2_anomaly_CO2pulse(i) + cum_Net_Land_Flux2_anomaly_CO2pulse(i-1)
% end

%% atmospheric CO2 concentration to mass
%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
co2_control_z1_ts_CanESM5 = nan(1200,1);
max_length=length(co2_piControl_CanESM5(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_piControl_CanESM5(i,j,1,t))) 
            sum1 = sum1 + (co2_piControl_CanESM5(i,j,1,t).*area_a_CanESM5(i,j)); %original unit: ppm
            sum2 = sum2 + area_a_CanESM5(i,j);
         end
     end
    end
   co2_control_z1_ts_CanESM5(t)=sum1./sum2; 
   end
  co2_control_z1_ts_CanESM = co2_control_z1_ts_CanESM5;
co2_control_z1_PgC_CanESM = co2_control_z1_ts_CanESM5 * ppm_to_pgc;


%%%%%%%%%% CO2 PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_CO2pulse_CanESM5(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_CO2pulse_CanESM5(i,j,1,t))) 
            sum1 = sum1 + (co2_CO2pulse_CanESM5(i,j,1,t).*area_a_CanESM5(i,j)); %original unit: mol/mol
            sum2 = sum2 + area_a_CanESM5(i,j);
         end
     end
    end
   co2_CO2pulse_z1_ts_CanESM5(t)=sum1./sum2; 
  end
co2_CO2pulse_z1_PgC_CanESM = co2_CO2pulse_z1_ts_CanESM5 * ppm_to_pgc;


%%%%%%%%%% CDR PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_cdr_pulse_CanESM5(1,1,1,:));
co2_cdr_pulse_z1_ts_CanESM5 = nan(1200,1);
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_cdr_pulse_CanESM5(i,j,1,t))) 
            sum1 = sum1 + (co2_cdr_pulse_CanESM5(i,j,1,t).*area_a_CanESM5(i,j)); %original unit: mol/mol
            sum2 = sum2 + area_a_CanESM5(i,j);
         end
     end
    end
   co2_cdr_pulse_z1_ts_CanESM5(t)=sum1./sum2; 
  end
  co2_cdr_pulse_z1_ts_CanESM = co2_cdr_pulse_z1_ts_CanESM5;
co2_cdr_pulse_z1_PgC_CanESM = co2_cdr_pulse_z1_ts_CanESM5 * ppm_to_pgc;

%%% calculate annual average
co2_PgC_cdr_pulse_z1_yr_CanESM5 = nan(100,1);
co2_PgC_CO2pulse_z1_yr_CanESM5 = nan(100,1);
co2_PgC_piControl_z1_yr_CanESM5 = nan(100,1);

for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  if (m2 <= length(co2_control_z1_PgC_CanESM))
  co2_PgC_piControl_z1_yr_CanESM5(t) = mean(co2_control_z1_PgC_CanESM(m:m2)) ;
  end
  if m2<=length(co2_CO2pulse_z1_PgC_CanESM) 
  co2_PgC_CO2pulse_z1_yr_CanESM5(t) = mean(co2_CO2pulse_z1_PgC_CanESM(m:m2));
  end
  if m2<=length(co2_cdr_pulse_z1_PgC_CanESM) 
      co2_PgC_cdr_pulse_z1_yr_CanESM5(t) = mean(co2_cdr_pulse_z1_PgC_CanESM(m:m2));
  end
end

%% cumulative emissions
% sum of -cum_Ocean_flux, -cum_Land_flux and atmospheric change is defined as
% cumulative emissions at given point in time.

% atmospheric CO2 in PgC, year 100
atm_CO2_PgC_CO2pulse_100yr_CanESM = co2_PgC_CO2pulse_z1_yr_CanESM5(1:100)-co2_PgC_piControl_z1_yr_CanESM5(1:100);
atm_CO2_PgC_cdr_pulse_100yr_CanESM = co2_PgC_cdr_pulse_z1_yr_CanESM5(1:100)-co2_PgC_piControl_z1_yr_CanESM5(1:100);

% cumulative emissions of CO2 in PgC, year 100
cum_emissions_cdr_pulse_100yr_CanESM = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_CanESM(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CanESM(1:100))...
                               +abs(co2_PgC_cdr_pulse_z1_yr_CanESM5(1:100)-co2_PgC_piControl_z1_yr_CanESM5(1:100));
 
cum_emissions_CO2pulse_100yr_CanESM = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_CanESM(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_CanESM(1:100))...
                               +abs(co2_PgC_CO2pulse_z1_yr_CanESM5(1:100)-co2_PgC_piControl_z1_yr_CanESM5(1:100));
 
% airborne fraction, year 100
airborne_fraction_CO2pulse_CanESM  = atm_CO2_PgC_CO2pulse_100yr_CanESM./cum_emissions_CO2pulse_100yr_CanESM;
airborne_fraction_cdr_pulse_CanESM = atm_CO2_PgC_cdr_pulse_100yr_CanESM./cum_emissions_cdr_pulse_100yr_CanESM;

%% delete variables that are not necessary for following calculations/plots
clearvars -except *ACCESS *MIROC *CanESM* *NorESM *UKESM *UVic *CanOE *_OE *GFDL
 
%% CanESM5-CanOE 
area_o_CanESM5_CanOE = ncread('areacello_Ofx_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn.nc', 'areacello');
area_a_CanESM5_CanOE = ncread('areacella_fx_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn.nc', 'areacella');
lat_CanESM5_CanOE = ncread('areacella_fx_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn.nc', 'lat');
lon_CanESM5_CanOE = ncread('areacella_fx_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn.nc', 'lon');

co2_cdr_pulse_CanESM5_CanOE = ncread('co2_Amon_CanESM5-CanOE_esm-pi-cdr-pulse_r1i1p2f1_gn_540101-560012.nc', 'co2',[1 1 1 1],[128 64 1 1200]); %
co2_CO2pulse_CanESM5_CanOE = ncread('co2_Amon_CanESM5-CanOE_esm-pi-CO2pulse_r1i1p2f1_gn_540101-560012.nc', 'co2',[1 1 1 1],[128 64 1 1200]);
co2_piControl_CanESM5_CanOE = ncread('co2_Amon_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn_600101-620012.nc', 'co2',[1 1 1 1],[128 64 1 1200]);
ppm_to_pgc = 2.124;

nbp_cdr_pulse_CanESM5_CanOE = ncread('nbp_Lmon_CanESM5-CanOE_esm-pi-cdr-pulse_r1i1p2f1_gn_540101-560012.nc', 'nbp',[1 1 1],[128 64 1200]);
nbp_CO2pulse_CanESM5_CanOE = ncread('nbp_Lmon_CanESM5-CanOE_esm-pi-CO2pulse_r1i1p2f1_gn_540101-560012.nc', 'nbp',[1 1 1],[128 64 1200]);
nbp_piControl_CanESM5_CanOE = ncread('nbp_Lmon_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn_600101-620012.nc', 'nbp',[1 1 1],[128 64 1200]);

a2o_flux_cdr_pulse_CanESM5_CanOE = ncread('fgco2_Omon_CanESM5-CanOE_esm-pi-cdr-pulse_r1i1p2f1_gn_540101-560012.nc', 'fgco2',[1 1 1],[360 291 1200]);
a2o_flux_CO2pulse_CanESM5_CanOE = ncread('fgco2_Omon_CanESM5-CanOE_esm-pi-CO2pulse_r1i1p2f1_gn_540101-560012.nc', 'fgco2',[1 1 1],[360 291 1200]);
a2o_flux_piControl_CanESM5_CanOE = ncread('fgco2_Omon_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn_600101-620012.nc', 'fgco2',[1 1 1],[360 291 1200]);

sftlf_CanESM5_CanOE = ncread('sftlf_fx_CanESM5-CanOE_esm-piControl_r1i1p2f1_gn.nc', 'sftlf');

%% Ocean C change
year_in_sec = 60 .* 60 .* 24 .* 365;  %unit: s
time_max_2 = length(a2o_flux_cdr_pulse_CanESM5_CanOE(1,1,:));
a2o_flux_CO2pulse_ts = nan(time_max_2,1);
a2o_flux_cdr_pulse_ts = nan(time_max_2,1);
a2o_flux_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   a2o_flux_CO2pulse_ts(t,1) = nansum(nansum(a2o_flux_CO2pulse_CanESM5_CanOE(:,:,t).*area_o_CanESM5_CanOE)); 
   a2o_flux_cdr_pulse_ts(t,1) = nansum(nansum(a2o_flux_cdr_pulse_CanESM5_CanOE(:,:,t).*area_o_CanESM5_CanOE)); 
   a2o_flux_piControl_ts(t,1) = nansum(nansum(a2o_flux_piControl_CanESM5_CanOE(:,:,t).*area_o_CanESM5_CanOE));
end   

% (yearly mean multiplied by 365)
a2o_flux_CO2pulse_ts_yrly = nan(100,1);
a2o_flux_cdr_pulse_ts_yrly = nan(100,1);
a2o_flux_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  a2o_flux_CO2pulse_ts_yrly(t) = mean(a2o_flux_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  a2o_flux_cdr_pulse_ts_yrly(t) = mean(a2o_flux_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  a2o_flux_piControl_ts_yrly(t) = mean(a2o_flux_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

a2o_flux2_anomaly = a2o_flux_cdr_pulse_ts_yrly-a2o_flux_piControl_ts_yrly;
cum_a2o_flux2_anomaly_CanESM5_CanOE = a2o_flux2_anomaly;
for i=2:length(a2o_flux2_anomaly)
cum_a2o_flux2_anomaly_CanESM5_CanOE(i) = a2o_flux2_anomaly(i) + cum_a2o_flux2_anomaly_CanESM5_CanOE(i-1);
end

a2o_flux2_anomaly_CO2pulse = a2o_flux_CO2pulse_ts_yrly-a2o_flux_piControl_ts_yrly;
cum_a2o_flux2_anomaly_CO2pulse_CanESM5_CanOE = a2o_flux2_anomaly_CO2pulse;
for i=2:length(a2o_flux2_anomaly_CO2pulse)
cum_a2o_flux2_anomaly_CO2pulse_CanESM5_CanOE(i) = a2o_flux2_anomaly_CO2pulse(i) + cum_a2o_flux2_anomaly_CO2pulse_CanESM5_CanOE(i-1);
end

%% nbp 
time_max_2 = length(nbp_cdr_pulse_CanESM5_CanOE(1,1,:));
nbp_CO2pulse_ts = nan(time_max_2,1);
nbp_cdr_pulse_ts = nan(time_max_2,1);
nbp_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   nbp_CO2pulse_ts(t,1) = nansum(nansum(nbp_CO2pulse_CanESM5_CanOE(:,:,t).*area_a_CanESM5_CanOE.*sftlf_CanESM5_CanOE./100)); 
   nbp_cdr_pulse_ts(t,1) = nansum(nansum(nbp_cdr_pulse_CanESM5_CanOE(:,:,t).*area_a_CanESM5_CanOE.*sftlf_CanESM5_CanOE./100)); 
   nbp_piControl_ts(t,1) = nansum(nansum(nbp_piControl_CanESM5_CanOE(:,:,t).*area_a_CanESM5_CanOE.*sftlf_CanESM5_CanOE./100));
end   

% (yearly mean multiplied by 365)
nbp_CO2pulse_ts_yrly = nan(100,1);
nbp_cdr_pulse_ts_yrly = nan(100,1);
nbp_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  nbp_CO2pulse_ts_yrly(t) = mean(nbp_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_cdr_pulse_ts_yrly(t) = mean(nbp_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_piControl_ts_yrly(t) = mean(nbp_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

Net_Land_Flux2_anomaly = nbp_cdr_pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_CanESM5_CanOE = Net_Land_Flux2_anomaly;
for i=2:length(Net_Land_Flux2_anomaly)
cum_Net_Land_Flux2_anomaly_CanESM5_CanOE(i) = Net_Land_Flux2_anomaly(i) + cum_Net_Land_Flux2_anomaly_CanESM5_CanOE(i-1);
end

Net_Land_Flux2_anomaly_CO2pulse = nbp_CO2pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_CO2pulse_CanESM5_CanOE = Net_Land_Flux2_anomaly_CO2pulse;
for i=2:length(Net_Land_Flux2_anomaly_CO2pulse)
cum_Net_Land_Flux2_anomaly_CO2pulse_CanESM5_CanOE(i) = Net_Land_Flux2_anomaly_CO2pulse(i) + cum_Net_Land_Flux2_anomaly_CO2pulse_CanESM5_CanOE(i-1);
end

%% atmospheric CO2 concentration to mass
%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
co2_control_z1_ts = nan(1200,1);
max_length=length(co2_piControl_CanESM5_CanOE(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_piControl_CanESM5_CanOE(i,j,1,t))) 
            sum1 = sum1 + (co2_piControl_CanESM5_CanOE(i,j,1,t).*area_a_CanESM5_CanOE(i,j)); %original unit: ppm
            sum2 = sum2 + area_a_CanESM5_CanOE(i,j);
         end
     end
    end
   co2_control_z1_ts(t)=sum1./sum2; 
   end
  co2_control_z1_ts_CanESM5_CanOE = co2_control_z1_ts;
co2_control_z1_PgC_CanESM5_CanOE = co2_control_z1_ts * ppm_to_pgc;

%%%%%%%%%% CO2 PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
co2_CO2pulse_z1_ts = nan(1200,1);
max_length=length(co2_CO2pulse_CanESM5_CanOE(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_CO2pulse_CanESM5_CanOE(i,j,1,t))) 
            sum1 = sum1 + (co2_CO2pulse_CanESM5_CanOE(i,j,1,t).*area_a_CanESM5_CanOE(i,j)); %original unit: mol/mol
            sum2 = sum2 + area_a_CanESM5_CanOE(i,j);
         end
     end
    end
   co2_CO2pulse_z1_ts(t)=sum1./sum2; 
  end
  co2_CO2pulse_z1_ts_CanESM5_CanOE = co2_CO2pulse_z1_ts;
co2_CO2pulse_z1_PgC_CanESM5_CanOE = co2_CO2pulse_z1_ts * ppm_to_pgc;

%%%%%%%%%% CDR PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_cdr_pulse_CanESM5_CanOE(1,1,1,:));
co2_cdr_pulse_z1_ts = nan(1200,1);
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:128)
     for(j=1:64)
         if(~isnan(co2_cdr_pulse_CanESM5_CanOE(i,j,1,t))) 
            sum1 = sum1 + (co2_cdr_pulse_CanESM5_CanOE(i,j,1,t).*area_a_CanESM5_CanOE(i,j)); %original unit: mol/mol
            sum2 = sum2 + area_a_CanESM5_CanOE(i,j);
         end
     end
    end
   co2_cdr_pulse_z1_ts(t)=sum1./sum2; 
  end
co2_cdr_pulse_z1_ts_CanESM5_CanOE = co2_cdr_pulse_z1_ts;
co2_cdr_pulse_z1_PgC_CanESM5_CanOE = co2_cdr_pulse_z1_ts * ppm_to_pgc;

%%% calculate annual average
co2_PgC_cdr_pulse_z1_yr_CanESM5_CanOE = nan(100,1);
co2_PgC_CO2pulse_z1_yr_CanESM5_CanOE = nan(100,1);
co2_PgC_piControl_z1_yr_CanESM5_CanOE = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  if (m2 <= length(co2_control_z1_PgC_CanESM5_CanOE))
  co2_PgC_piControl_z1_yr_CanESM5_CanOE(t) = mean(co2_control_z1_PgC_CanESM5_CanOE(m:m2));
  end
  if m2<=length(co2_CO2pulse_z1_PgC_CanESM5_CanOE) 
  co2_PgC_CO2pulse_z1_yr_CanESM5_CanOE(t) = mean(co2_CO2pulse_z1_PgC_CanESM5_CanOE(m:m2));
  end
  if m2<=length(co2_cdr_pulse_z1_PgC_CanESM5_CanOE) 
      co2_PgC_cdr_pulse_z1_yr_CanESM5_CanOE(t) = mean(co2_cdr_pulse_z1_PgC_CanESM5_CanOE(m:m2));
  end
end

%% UPTKAE/RELEASE/AIRBORNE FRACTIONS, TIMESERIES
cum_Net_Land_Flux_anomaly_CanESM5_CanOE = cum_Net_Land_Flux2_anomaly_CanESM5_CanOE;
cum_Net_Land_Flux_anomaly_CO2pulse_CanESM5_CanOE = cum_Net_Land_Flux2_anomaly_CO2pulse_CanESM5_CanOE;

Net_Land_Flux_anomaly_cdr_pulse_CanESM5_CanOE = Net_Land_Flux2_anomaly;
Net_Land_Flux_anomaly_CO2pulse_CanESM5_CanOE = Net_Land_Flux2_anomaly_CO2pulse;

Ocean_C_Flux_anomaly_cdr_pulse_CanESM5_CanOE = a2o_flux2_anomaly;
Ocean_C_Flux_anomaly_CO2pulse_CanESM5_CanOE = a2o_flux2_anomaly_CO2pulse;

cum_Ocean_C_Flux_anomaly_cdr_pulse_CanESM5_CanOE = cum_a2o_flux2_anomaly_CanESM5_CanOE;
cum_Ocean_C_Flux_anomaly_CO2pulse_CanESM5_CanOE = cum_a2o_flux2_anomaly_CO2pulse_CanESM5_CanOE;


%% calculate annual average
co2_PgC_cdr_pulse_z1_yr_CanESM5_CanOE = nan(100,1);
co2_PgC_CO2pulse_z1_yr_CanESM5_CanOE = nan(100,1);
co2_PgC_piControl_z1_yr_CanESM5_CanOE = nan(100,1);

for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  if (m2 <= length(co2_control_z1_PgC_CanESM5_CanOE))
  co2_PgC_piControl_z1_yr_CanESM5_CanOE(t) = mean(co2_control_z1_PgC_CanESM5_CanOE(m:m2)) ;
  end
  if m2<=length(co2_CO2pulse_z1_PgC_CanESM5_CanOE) 
  co2_PgC_CO2pulse_z1_yr_CanESM5_CanOE(t) = mean(co2_CO2pulse_z1_PgC_CanESM5_CanOE(m:m2));
  end
  if m2<=length(co2_cdr_pulse_z1_PgC_CanESM5_CanOE) 
      co2_PgC_cdr_pulse_z1_yr_CanESM5_CanOE(t) = mean(co2_cdr_pulse_z1_PgC_CanESM5_CanOE(m:m2));
  end
end

%% cumulative emissions
% sum of -cum_Ocean_flux, -cum_Land_flux and atmospheric change is defined as
% cumulative emissions at given point in time.

% atmospheric CO2 in PgC, year 100
atm_CO2_PgC_CO2pulse_100yr_CanESM5_CanOE = co2_PgC_CO2pulse_z1_yr_CanESM5_CanOE(1:100)-co2_PgC_piControl_z1_yr_CanESM5_CanOE(1:100);
atm_CO2_PgC_cdr_pulse_100yr_CanESM5_CanOE = co2_PgC_cdr_pulse_z1_yr_CanESM5_CanOE(1:100)-co2_PgC_piControl_z1_yr_CanESM5_CanOE(1:100);

% cumulative emissions of CO2 in PgC, year 100
cum_emissions_cdr_pulse_100yr_CanESM5_CanOE = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_CanESM5_CanOE(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CanESM5_CanOE(1:100))...
                               +abs(co2_PgC_cdr_pulse_z1_yr_CanESM5_CanOE(1:100)-co2_PgC_piControl_z1_yr_CanESM5_CanOE(1:100));
 
cum_emissions_CO2pulse_100yr_CanESM5_CanOE = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_CanESM5_CanOE(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_CanESM5_CanOE(1:100))...
                               +abs(co2_PgC_CO2pulse_z1_yr_CanESM5_CanOE(1:100)-co2_PgC_piControl_z1_yr_CanESM5_CanOE(1:100));
 
% airborne fraction, year 100
airborne_fraction_CO2pulse_CanESM5_CanOE  = atm_CO2_PgC_CO2pulse_100yr_CanESM5_CanOE./cum_emissions_CO2pulse_100yr_CanESM5_CanOE;
airborne_fraction_cdr_pulse_CanESM5_CanOE = atm_CO2_PgC_cdr_pulse_100yr_CanESM5_CanOE./cum_emissions_cdr_pulse_100yr_CanESM5_CanOE;

%% delete variables that are not necessary for following calculations/plots
clearvars -except *ACCESS *MIROC *CanESM* *NorESM *UKESM *UVic *CanOE *_OE *GFDL
 
%% NorESM 
area_o_NorESM = ncread('areacello_Ofx_NorESM2-LM_esm-piControl_r1i1p1f1_gn.nc', 'areacello');
area_a_NorESM = ncread('areacella_fx_NorESM2-LM_esm-piControl_r1i1p1f1_gn.nc', 'areacella');
lat_NorESM = ncread('areacella_fx_NorESM2-LM_esm-piControl_r1i1p1f1_gn.nc', 'lat');
lon_NorESM = ncread('areacella_fx_NorESM2-LM_esm-piControl_r1i1p1f1_gn.nc', 'lon');

co2_cdr_pulse_NorESM = squeeze(cat(4,ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_185001-185912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_186001-186912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_187001-187912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_188001-188912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_189001-189912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_190001-190912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_191001-191912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_192001-192912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_193001-193912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_194001-194912.nc', 'co2',[1 1 1 1],[144 96 1 120])));

co2_CO2pulse_NorESM = squeeze(cat(4,ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_185001-185912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_186001-186912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_187001-187912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_188001-188912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_189001-189912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_190001-190912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_191001-191912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_192001-192912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_193001-193912.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_194001-194912.nc', 'co2',[1 1 1 1],[144 96 1 120])));

co2_piControl_NorESM = squeeze(cat(4,ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_185101-186012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_186101-187012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_187101-188012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_188101-189012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_189101-190012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_190101-191012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_191101-192012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_192101-193012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_193101-194012.nc', 'co2',[1 1 1 1],[144 96 1 120]),... 
    ncread('co2_Amon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_194101-195012.nc', 'co2',[1 1 1 1],[144 96 1 120])));


nbp_cdr_pulse_NorESM = cat(3,ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_185001-185912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_186001-186912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_187001-187912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_188001-188912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_189001-189912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_190001-190912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_191001-191912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_192001-192912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_193001-193912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_194001-194912.nc', 'nbp'));

nbp_CO2pulse_NorESM = cat(3,ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_185001-185912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_186001-186912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_187001-187912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_188001-188912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_189001-189912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_190001-190912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_191001-191912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_192001-192912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_193001-193912.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_194001-194912.nc', 'nbp'));

nbp_piControl_NorESM = cat(3,ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_185101-186012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_186101-187012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_187101-188012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_188101-189012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_189101-190012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_190101-191012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_191101-192012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_192101-193012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_193101-194012.nc', 'nbp'),... 
    ncread('nbp_Lmon_NorESM2-LM_esm-piControl_r1i1p1f1_gn_194101-195012.nc', 'nbp'));

a2o_flux_cdr_pulse_NorESM = cat(3,ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1850-1859.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1860-1869.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1870-1879.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1880-1889.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1890-1899.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1900-1909.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1910-1919.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1920-1929.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1930-1939.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-cdr-pulse_r1i1p1f1_gn_1940-1949.nc', 'fgco2'));

a2o_flux_CO2pulse_NorESM = cat(3,ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1850-1859.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1860-1869.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1870-1879.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1880-1889.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1890-1899.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1900-1909.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1910-1919.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1920-1929.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1930-1939.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-pi-CO2pulse_r1i1p1f1_gn_1940-1949.nc', 'fgco2'));

a2o_flux_piControl_NorESM = cat(3,ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1851-1860.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1861-1870.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1871-1880.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1881-1890.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1891-1900.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1901-1910.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1911-1920.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1921-1930.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1931-1940.nc', 'fgco2'),... 
    ncread('fgco2_Oyr_NorESM2-LM_esm-piControl_r1i1p1f1_gn_1941-1950.nc', 'fgco2'));

sftlf_NorESM = ncread('sftlf_fx_NorESM2-LM_esm-piControl_r1i1p1f1_gn.nc', 'sftlf');

ppm_to_pgc = 2.124;

%% Ocean C change
year_in_sec = 60 .* 60 .* 24 .* 365;  %unit: s
time_max_2 = length(a2o_flux_cdr_pulse_NorESM(1,1,:));
a2o_flux_CO2pulse_ts = nan(time_max_2,1);
a2o_flux_cdr_pulse_ts = nan(time_max_2,1);
a2o_flux_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   a2o_flux_CO2pulse_ts(t,1) = nansum(nansum(a2o_flux_CO2pulse_NorESM(:,:,t).*area_o_NorESM)); 
   a2o_flux_cdr_pulse_ts(t,1) = nansum(nansum(a2o_flux_cdr_pulse_NorESM(:,:,t).*area_o_NorESM)); 
   a2o_flux_piControl_ts(t,1) = nansum(nansum(a2o_flux_piControl_NorESM(:,:,t).*area_o_NorESM));
end   

a2o_flux_CO2pulse_ts_yrly = a2o_flux_CO2pulse_ts ./ 10^12 * year_in_sec;
a2o_flux_cdr_pulse_ts_yrly = a2o_flux_cdr_pulse_ts ./ 10^12 * year_in_sec;
a2o_flux_piControl_ts_yrly = a2o_flux_piControl_ts ./ 10^12 * year_in_sec;

a2o_flux2_anomaly_NorESM = a2o_flux_cdr_pulse_ts_yrly-a2o_flux_piControl_ts_yrly;
cum_a2o_flux2_anomaly = a2o_flux2_anomaly_NorESM;
for i=2:length(a2o_flux2_anomaly_NorESM)
cum_a2o_flux2_anomaly(i) = a2o_flux2_anomaly_NorESM(i) + cum_a2o_flux2_anomaly(i-1);
end

a2o_flux2_anomaly_CO2pulse_NorESM = a2o_flux_CO2pulse_ts_yrly-a2o_flux_piControl_ts_yrly;
cum_a2o_flux2_anomaly_CO2pulse = a2o_flux2_anomaly_CO2pulse_NorESM;
for i=2:length(a2o_flux2_anomaly_CO2pulse_NorESM)
cum_a2o_flux2_anomaly_CO2pulse(i) = a2o_flux2_anomaly_CO2pulse_NorESM(i) + cum_a2o_flux2_anomaly_CO2pulse(i-1);
end

%% nbp 
time_max_2 = length(nbp_cdr_pulse_NorESM(1,1,:));
nbp_CO2pulse_ts = nan(time_max_2,1);
nbp_cdr_pulse_ts = nan(time_max_2,1);
nbp_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   nbp_CO2pulse_ts(t,1) = nansum(nansum(nbp_CO2pulse_NorESM(:,:,t).*area_a_NorESM.*sftlf_NorESM./100)); 
   nbp_cdr_pulse_ts(t,1) = nansum(nansum(nbp_cdr_pulse_NorESM(:,:,t).*area_a_NorESM.*sftlf_NorESM./100)); 
   nbp_piControl_ts(t,1) = nansum(nansum(nbp_piControl_NorESM(:,:,t).*area_a_NorESM.*sftlf_NorESM./100));
end   

% (yearly mean multiplied by 365)
nbp_CO2pulse_ts_yrly = nan(100,1);
nbp_cdr_pulse_ts_yrly = nan(100,1);
nbp_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  nbp_CO2pulse_ts_yrly(t) = mean(nbp_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_cdr_pulse_ts_yrly(t) = mean(nbp_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_piControl_ts_yrly(t) = mean(nbp_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

Net_Land_Flux2_anomaly = nbp_cdr_pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_NorESM = Net_Land_Flux2_anomaly;
for i=2:length(Net_Land_Flux2_anomaly)
cum_Net_Land_Flux2_anomaly_NorESM(i) = Net_Land_Flux2_anomaly(i) + cum_Net_Land_Flux2_anomaly_NorESM(i-1);
end

Net_Land_Flux2_anomaly_CO2pulse = nbp_CO2pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_CO2pulse_NorESM = Net_Land_Flux2_anomaly_CO2pulse;
for i=2:length(Net_Land_Flux2_anomaly_CO2pulse)
cum_Net_Land_Flux2_anomaly_CO2pulse_NorESM(i) = Net_Land_Flux2_anomaly_CO2pulse(i) + cum_Net_Land_Flux2_anomaly_CO2pulse_NorESM(i-1);
end

%% UPTKAE/RELEASE/AIRBORNE FRACTIONS, TIMESERIES
cum_Net_Land_Flux_anomaly_NorESM = cum_Net_Land_Flux2_anomaly_NorESM;
cum_Net_Land_Flux_anomaly_CO2pulse_NorESM = cum_Net_Land_Flux2_anomaly_CO2pulse_NorESM;

Net_Land_Flux_anomaly_cdr_pulse_NorESM = Net_Land_Flux2_anomaly;
Net_Land_Flux_anomaly_CO2pulse_NorESM = Net_Land_Flux2_anomaly_CO2pulse;

Ocean_C_Flux_anomaly_cdr_pulse_NorESM = a2o_flux2_anomaly_NorESM;
Ocean_C_Flux_anomaly_CO2pulse_NorESM = a2o_flux2_anomaly_CO2pulse_NorESM;

cum_Ocean_C_Flux_anomaly_cdr_pulse_NorESM = cum_a2o_flux2_anomaly;
cum_Ocean_C_Flux_anomaly_CO2pulse_NorESM = cum_a2o_flux2_anomaly_CO2pulse;

%% co2
% LOWEST layer global average co2
co2_cdr_pulse_z1_ts_NorESM = nan(1200,1);
max_length=length(co2_cdr_pulse_NorESM(1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:144)
     for(j=1:96)
         if(~isnan(co2_cdr_pulse_NorESM(i,j,t))) 
            sum1 = sum1 + (co2_cdr_pulse_NorESM(i,j,t).*area_a_NorESM(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_NorESM(i,j);
         end
     end
    end
   co2_cdr_pulse_z1_ts_NorESM(t)=sum1./sum2; 
   end
co2_cdr_pulse_z1_PgC_NorESM = co2_cdr_pulse_z1_ts_NorESM * ppm_to_pgc;

%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_piControl_NorESM(1,1,:));
co2_control_z1_ts_NorESM = nan(1200,1);
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:144)
     for(j=1:96)
         if(~isnan(co2_piControl_NorESM(i,j,t))) 
            sum1 = sum1 + (co2_piControl_NorESM(i,j,t).*area_a_NorESM(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_NorESM(i,j);
         end
     end
    end
   co2_control_z1_ts_NorESM(t)=sum1./sum2; 
   end
co2_control_z1_PgC_NorESM = co2_control_z1_ts_NorESM * ppm_to_pgc;

%%%%%%%%% CO2 PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_CO2pulse_NorESM(1,1,:));
co2_CO2pulse_NorESM_z1_ts = nan(1200,1);
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:144)
     for(j=1:96)
         if(~isnan(co2_CO2pulse_NorESM(i,j,t))) 
            sum1 = sum1 + (co2_CO2pulse_NorESM(i,j,t).*area_a_NorESM(i,j))*10^6;%*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_NorESM(i,j);
         end
     end
    end
   co2_CO2pulse_NorESM_z1_ts(t)=sum1./sum2; 
   end
co2_CO2pulse_z1_PgC_NorESM = co2_CO2pulse_NorESM_z1_ts * ppm_to_pgc;

%% calculate annual CO2 average
co2_PgC_cdr_pulse_z1_yr_NorESM = nan(100,1);
co2_PgC_CO2pulse_z1_yr_NorESM = nan(100,1);
co2_PgC_piControl_z1_yr_NorESM = nan(100,1);

for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  if (m2 <= length(co2_control_z1_PgC_NorESM))
  co2_PgC_piControl_z1_yr_NorESM(t) = mean(co2_control_z1_PgC_NorESM(m:m2)) ;
  end
  if m2<=length(co2_CO2pulse_z1_PgC_NorESM) 
  co2_PgC_CO2pulse_z1_yr_NorESM(t) = mean(co2_CO2pulse_z1_PgC_NorESM(m:m2));
  end
  if m2<=length(co2_cdr_pulse_z1_PgC_NorESM) 
      co2_PgC_cdr_pulse_z1_yr_NorESM(t) = mean(co2_cdr_pulse_z1_PgC_NorESM(m:m2));
  end
end

%% cumulative emissions
% sum of -cum_Ocean_flux, -cum_Land_flux and atmospheric change is defined as
% cumulative emissions at given point in time.

% atmospheric CO2 in PgC
atm_CO2_PgC_CO2pulse_100yr_NorESM = co2_PgC_CO2pulse_z1_yr_NorESM(1:100)-co2_PgC_piControl_z1_yr_NorESM(1:100);
atm_CO2_PgC_cdr_pulse_100yr_NorESM = co2_PgC_cdr_pulse_z1_yr_NorESM(1:100)-co2_PgC_piControl_z1_yr_NorESM(1:100);

% cumulative emissions of CO2 in PgC
cum_emissions_cdr_pulse_100yr_NorESM = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_NorESM(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_NorESM(1:100))...
                               +abs(co2_PgC_cdr_pulse_z1_yr_NorESM(1:100)-co2_PgC_piControl_z1_yr_NorESM(1:100));
 
cum_emissions_CO2pulse_100yr_NorESM = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_NorESM(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_NorESM(1:100))...
                               +abs(co2_PgC_CO2pulse_z1_yr_NorESM(1:100)-co2_PgC_piControl_z1_yr_NorESM(1:100));
 
% airborne fraction
airborne_fraction_CO2pulse_NorESM  = atm_CO2_PgC_CO2pulse_100yr_NorESM./cum_emissions_CO2pulse_100yr_NorESM;
airborne_fraction_cdr_pulse_NorESM = atm_CO2_PgC_cdr_pulse_100yr_NorESM./cum_emissions_cdr_pulse_100yr_NorESM;

%% delete variables that are not necessary for following calculations/plots
clearvars -except *ACCESS *MIROC *CanESM* *NorESM *UKESM *UVic *CanOE *_OE *GFDL

%% Asymmetry UVic
A_co2_CO2pulse_UVic = ncread('ts_intgrls_1xco2equil_100.nc', 'A_co2');
O_totcarb_CO2pulse_UVic = ncread('ts_intgrls_1xco2equil_100.nc', 'O_totcarb');
L_totcarb_CO2pulse_UVic = ncread('ts_intgrls_1xco2equil_100.nc', 'L_totcarb');
A_totcarb_CO2pulse_UVic = ncread('ts_intgrls_1xco2equil_100.nc', 'A_totcarb');

A_co2_cdr_UVic = ncread('ts_intgrls_1xco2equil_neg100.nc', 'A_co2');
O_totcarb_cdr_UVic = ncread('ts_intgrls_1xco2equil_neg100.nc', 'O_totcarb');
L_totcarb_cdr_UVic = ncread('ts_intgrls_1xco2equil_neg100.nc', 'L_totcarb');
A_totcarb_cdr_UVic = ncread('ts_intgrls_1xco2equil_neg100.nc', 'A_totcarb');

A_co2_control_UVic = ncread('ts_intgrls_1xco2equil_ctrl.nc', 'A_co2');
O_totcarb_control_UVic = ncread('ts_intgrls_1xco2equil_ctrl.nc', 'O_totcarb');
L_totcarb_control_UVic = ncread('ts_intgrls_1xco2equil_ctrl.nc', 'L_totcarb');
A_totcarb_control_UVic = ncread('ts_intgrls_1xco2equil_ctrl.nc', 'A_totcarb');

A_co2_spinup_UVic = ncread('ts_intgrls_1xco2spinup.nc', 'A_co2');

A_co2_controlled_UVic = A_co2_cdr_UVic - A_co2_control_UVic + A_co2_spinup_UVic(250);
O_totcarb_controlled_UVic = O_totcarb_cdr_UVic - O_totcarb_control_UVic;% + O_totcarb_spinup(250);
L_totcarb_controlled_UVic = L_totcarb_cdr_UVic - L_totcarb_control_UVic;% + L_totcarb_spinup(250);

%% a2o flux
a2o_flux_cdr_pulse_UVic = ncread('ts_intgrls_1xco2equil_neg100.nc', 'F_carba2o');
a2o_flux_piControl_UVic = ncread('ts_intgrls_1xco2equil_ctrl.nc', 'F_carba2o');
a2o_flux_cdr_controlled_UVic = a2o_flux_cdr_pulse_UVic(1:100) - a2o_flux_piControl_UVic(1:100);

a2o_cum_cdr_controlled_UVic = a2o_flux_cdr_controlled_UVic;
for i=2:100
a2o_cum_cdr_controlled_UVic(i) = a2o_flux_cdr_controlled_UVic(i) + a2o_cum_cdr_controlled_UVic(i-1);
end

%% a2l flux
a2l_flux_cdr_pulse_UVic = ncread('ts_intgrls_1xco2equil_neg100.nc', 'F_carba2l');
a2l_flux_piControl_UVic = ncread('ts_intgrls_1xco2equil_ctrl.nc', 'F_carba2l');
a2l_flux_cdr_controlled_UVic = a2l_flux_cdr_pulse_UVic(1:100) - a2l_flux_piControl_UVic(1:100);

a2l_cum_cdr_controlled_UVic = a2l_flux_cdr_controlled_UVic;
for i=2:100
a2l_cum_cdr_controlled_UVic(i) = a2l_flux_cdr_controlled_UVic(i) + a2l_cum_cdr_controlled_UVic(i-1);
end

%% cumulative emissions
% sum of -cum_Ocean_flux, -cum_Land_flux and atmospheric change is defined as
% cumulative emissions at given point in time.

% atmospheric CO2 in PgC, year 1:100
atm_CO2_PgC_CO2pulse_100yr_UVic = A_totcarb_CO2pulse_UVic(1:100)-A_totcarb_control_UVic(1:100);
atm_CO2_PgC_cdr_pulse_100yr_UVic = A_totcarb_cdr_UVic(1:100)-A_totcarb_control_UVic(1:100);

% cumulative emissions of CO2 in PgC, year 1:100
cum_emissions_cdr_pulse_100yr_UVic = abs(A_totcarb_cdr_UVic(1:100)-A_totcarb_control_UVic(1:100))...
                               +abs(L_totcarb_cdr_UVic(1:100)-L_totcarb_control_UVic(1:100))...
                               +abs(O_totcarb_cdr_UVic(1:100)-O_totcarb_control_UVic(1:100));
 
% cumulative emissions of CO2 in PgC, year 1:100
cum_emissions_CO2pulse_100yr_UVic = abs(A_totcarb_CO2pulse_UVic(1:100)-A_totcarb_control_UVic(1:100))...
                               +abs(L_totcarb_CO2pulse_UVic(1:100)-L_totcarb_control_UVic(1:100))...
                               +abs(O_totcarb_CO2pulse_UVic(1:100)-O_totcarb_control_UVic(1:100));
 
% airborne fraction, year 100
airborne_fraction_CO2pulse_UVic  = atm_CO2_PgC_CO2pulse_100yr_UVic./cum_emissions_CO2pulse_100yr_UVic;
airborne_fraction_cdr_pulse_UVic = atm_CO2_PgC_cdr_pulse_100yr_UVic./cum_emissions_cdr_pulse_100yr_UVic;

%% delete variables that are not necessary for following calculations/plots
clearvars -except *ACCESS *MIROC *CanESM* *NorESM *UKESM *UVic *CanOE *_OE *GFDL

%% UKESM1-0-LL 
%% same as in ACCESS: first 10 years are still control (even in pulse data), second dataset required
area_o_UKESM = ncread('areacello_Ofx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc', 'areacello');
area_a_UKESM = ncread('areacella_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc', 'areacella');

nbp_cdr_pulse_UKESM = cat(3,ncread('nbp_Lmon_UKESM1-0-LL_esm-pi-cdr-pulse_r1i1p1f2_gn_185001-194912.nc', 'nbp',[1 1 121],[192 144 1080]), ...
                      ncread('nbp_Lmon_UKESM1-0-LL_esm-pi-cdr-pulse_r1i1p1f2_gn_195001-198912.nc', 'nbp',[1 1 1],[192 144 120]));
nbp_CO2pulse_UKESM = cat(3,ncread('nbp_Lmon_UKESM1-0-LL_esm-pi-CO2pulse_r1i1p1f2_gn_185001-194912.nc', 'nbp',[1 1 121],[192 144 1080]), ...
                      ncread('nbp_Lmon_UKESM1-0-LL_esm-pi-CO2pulse_r1i1p1f2_gn_195001-196912.nc', 'nbp',[1 1 1],[192 144 120]));
nbp_piControl_UKESM = cat(3,ncread('nbp_Lmon_UKESM1-0-LL_esm-piControl_r1i1p1f2_gn_185001-194912.nc', 'nbp',[1 1 121],[192 144 1080]), ...
                      ncread('nbp_Lmon_UKESM1-0-LL_esm-piControl_r1i1p1f2_gn_195001-204912.nc', 'nbp',[1 1 1],[192 144 120]));

a2o_flux_cdr_pulse_UKESM = cat(3,ncread('fgco2_Omon_UKESM1-0-LL_esm-pi-cdr-pulse_r1i1p1f2_gn_185001-194912.nc', 'fgco2',[1 1 121],[360 330 1080]), ...
                      ncread('fgco2_Omon_UKESM1-0-LL_esm-pi-cdr-pulse_r1i1p1f2_gn_195001-198912.nc', 'fgco2',[1 1 1],[360 330 120]));
a2o_flux_CO2pulse_UKESM = cat(3,ncread('fgco2_Omon_UKESM1-0-LL_esm-pi-CO2pulse_r1i1p1f2_gn_185001-194912.nc', 'fgco2',[1 1 121],[360 330 1080]), ...
                      ncread('fgco2_Omon_UKESM1-0-LL_esm-pi-CO2pulse_r1i1p1f2_gn_195001-196912.nc', 'fgco2',[1 1 1],[360 330 120]));
a2o_flux_piControl_UKESM = cat(3,ncread('fgco2_Omon_UKESM1-0-LL_esm-piControl_r1i1p1f2_gn_185001-194912.nc', 'fgco2',[1 1 121],[360 330 1080]), ...
                      ncread('fgco2_Omon_UKESM1-0-LL_esm-piControl_r1i1p1f2_gn_195001-204912.nc', 'fgco2',[1 1 1],[360 330 120]));

sftlf_UKESM = ncread('sftlf_fx_UKESM1-0-LL_piControl_r1i1p1f2_gn.nc', 'sftlf');

co2mass_CO2pulse_UKESM =  cat(1,ncread('co2mass_Amon_UKESM1-0-LL_esm-pi-CO2pulse_r1i1p1f2_gm_185001-194912.nc', 'co2mass',[121],[1080]), ...
                          ncread('co2mass_Amon_UKESM1-0-LL_esm-pi-CO2pulse_r1i1p1f2_gm_195001-196912.nc', 'co2mass',[1],[120]) );
co2mass_cdr_pulse_UKESM =  cat(1,ncread('co2mass_Amon_UKESM1-0-LL_esm-pi-cdr-pulse_r1i1p1f2_gm_185001-194912.nc', 'co2mass',[121],[1080]), ...
                           ncread('co2mass_Amon_UKESM1-0-LL_esm-pi-cdr-pulse_r1i1p1f2_gm_195001-198912.nc', 'co2mass',[1],[120]) );
co2mass_piControl_UKESM =  cat(1,ncread('co2mass_Amon_UKESM1-0-LL_esm-piControl_r1i1p1f2_gm_185001-194912.nc', 'co2mass',[121],[1080]), ...
                           ncread('co2mass_Amon_UKESM1-0-LL_esm-piControl_r1i1p1f2_gm_195001-204912.nc', 'co2mass',[1],[120]) );

ppm_to_pgc = 2.124;

%% Ocean C change
year_in_sec = 60 .* 60 .* 24 .* 365;  %unit: s
time_max_2 = length(a2o_flux_cdr_pulse_UKESM(1,1,:));
a2o_flux_CO2pulse_ts = nan(time_max_2,1);
a2o_flux_cdr_pulse_ts = nan(time_max_2,1);
a2o_flux_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   a2o_flux_CO2pulse_ts(t,1) = nansum(nansum(a2o_flux_CO2pulse_UKESM(:,:,t).*area_o_UKESM)); 
   a2o_flux_cdr_pulse_ts(t,1) = nansum(nansum(a2o_flux_cdr_pulse_UKESM(:,:,t).*area_o_UKESM)); 
   a2o_flux_piControl_ts(t,1) = nansum(nansum(a2o_flux_piControl_UKESM(:,:,t).*area_o_UKESM));
end   

% (yearly mean multiplied by 365)
a2o_flux_CO2pulse_ts_yrly = nan(100,1);
a2o_flux_cdr_pulse_ts_yrly = nan(100,1);
a2o_flux_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  a2o_flux_CO2pulse_ts_yrly(t) = mean(a2o_flux_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  a2o_flux_cdr_pulse_ts_yrly(t) = mean(a2o_flux_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  a2o_flux_piControl_ts_yrly(t) = mean(a2o_flux_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

a2o_flux2_anomaly = a2o_flux_cdr_pulse_ts_yrly-a2o_flux_piControl_ts_yrly;
cum_a2o_flux2_anomaly_UKESM = a2o_flux2_anomaly;
for i=2:length(a2o_flux2_anomaly)
cum_a2o_flux2_anomaly_UKESM(i) = a2o_flux2_anomaly(i) + cum_a2o_flux2_anomaly_UKESM(i-1);
end

a2o_flux2_anomaly_CO2pulse = a2o_flux_CO2pulse_ts_yrly-a2o_flux_piControl_ts_yrly;
cum_a2o_flux2_anomaly_CO2pulse_UKESM = a2o_flux2_anomaly_CO2pulse;
for i=2:length(a2o_flux2_anomaly_CO2pulse)
cum_a2o_flux2_anomaly_CO2pulse_UKESM(i) = a2o_flux2_anomaly_CO2pulse(i) + cum_a2o_flux2_anomaly_CO2pulse_UKESM(i-1);
end

%% nbp 
time_max_2 = length(nbp_cdr_pulse_UKESM(1,1,:));
nbp_CO2pulse_ts = nan(time_max_2,1);
nbp_cdr_pulse_ts = nan(time_max_2,1);
nbp_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   nbp_CO2pulse_ts(t,1) = nansum(nansum(nbp_CO2pulse_UKESM(:,:,t).*area_a_UKESM.*sftlf_UKESM./100)); 
   nbp_cdr_pulse_ts(t,1) = nansum(nansum(nbp_cdr_pulse_UKESM(:,:,t).*area_a_UKESM.*sftlf_UKESM./100)); 
   nbp_piControl_ts(t,1) = nansum(nansum(nbp_piControl_UKESM(:,:,t).*area_a_UKESM.*sftlf_UKESM./100));
end   

% (yearly mean multiplied by 365)
nbp_CO2pulse_ts_yrly = nan(100,1);
nbp_cdr_pulse_ts_yrly = nan(100,1);
nbp_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  nbp_CO2pulse_ts_yrly(t) = mean(nbp_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_cdr_pulse_ts_yrly(t) = mean(nbp_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_piControl_ts_yrly(t) = mean(nbp_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

Net_Land_Flux2_anomaly = nbp_cdr_pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_UKESM = Net_Land_Flux2_anomaly;
for i=2:length(Net_Land_Flux2_anomaly)
cum_Net_Land_Flux2_anomaly_UKESM(i) = Net_Land_Flux2_anomaly(i) + cum_Net_Land_Flux2_anomaly_UKESM(i-1);
end

Net_Land_Flux2_anomaly_CO2pulse = nbp_CO2pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_CO2pulse_UKESM = Net_Land_Flux2_anomaly_CO2pulse;
for i=2:length(Net_Land_Flux2_anomaly_CO2pulse)
cum_Net_Land_Flux2_anomaly_CO2pulse_UKESM(i) = Net_Land_Flux2_anomaly_CO2pulse(i) + cum_Net_Land_Flux2_anomaly_CO2pulse_UKESM(i-1);
end

%% UPTKAE/RELEASE/AIRBORNE FRACTIONS
cum_Net_Land_Flux_anomaly_UKESM = cum_Net_Land_Flux2_anomaly_UKESM;
cum_Net_Land_Flux_anomaly_CO2pulse_UKESM = cum_Net_Land_Flux2_anomaly_CO2pulse_UKESM;

Net_Land_Flux_anomaly_cdr_pulse_UKESM = Net_Land_Flux2_anomaly;
Net_Land_Flux_anomaly_CO2pulse_UKESM = Net_Land_Flux2_anomaly_CO2pulse;

Ocean_C_Flux_anomaly_cdr_pulse_UKESM = a2o_flux2_anomaly;
Ocean_C_Flux_anomaly_CO2pulse_UKESM = a2o_flux2_anomaly_CO2pulse;

cum_Ocean_C_Flux_anomaly_cdr_pulse_UKESM = cum_a2o_flux2_anomaly_UKESM;
cum_Ocean_C_Flux_anomaly_CO2pulse_UKESM = cum_a2o_flux2_anomaly_CO2pulse_UKESM;

%% atmospheric CO2 mass
%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%
co2_control_z1_PgC_UKESM = co2mass_piControl_UKESM/10^12/(44/12);
co2_control_z1_ts_UKESM = co2_control_z1_PgC_UKESM/ppm_to_pgc;

%%%%%%%%%% CO2 PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
co2_CO2pulse_z1_PgC_UKESM = co2mass_CO2pulse_UKESM/10^12/(44/12);

%%%%%%%%%% CDR PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
co2_cdr_pulse_z1_PgC_UKESM = co2mass_cdr_pulse_UKESM/10^12/(44/12);
co2_cdr_pulse_z1_ts_UKESM = co2_cdr_pulse_z1_PgC_UKESM/ppm_to_pgc;

%%% calculate annual average
co2_PgC_cdr_pulse_z1_yr_UKESM = nan(100,1);
co2_PgC_CO2pulse_z1_yr_UKESM = nan(100,1);
co2_PgC_piControl_z1_yr_UKESM = nan(100,1);

for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  if (m2 <= length(co2_control_z1_PgC_UKESM))
  co2_PgC_piControl_z1_yr_UKESM(t) = mean(co2_control_z1_PgC_UKESM(m:m2)) ;
  end
  if m2<=length(co2_CO2pulse_z1_PgC_UKESM) 
  co2_PgC_CO2pulse_z1_yr_UKESM(t) = mean(co2_CO2pulse_z1_PgC_UKESM(m:m2));
  end
  if m2<=length(co2_cdr_pulse_z1_PgC_UKESM) 
      co2_PgC_cdr_pulse_z1_yr_UKESM(t) = mean(co2_cdr_pulse_z1_PgC_UKESM(m:m2));
  end
end

%% cumulative emissions
% sum of -cum_Ocean_flux, -cum_Land_flux and atmospheric change is defined as
% cumulative emissions at given point in time.

% atmospheric CO2 in PgC, year 100
atm_CO2_PgC_CO2pulse_100yr_UKESM = co2_PgC_CO2pulse_z1_yr_UKESM(1:100)-co2_PgC_piControl_z1_yr_UKESM(1:100);
atm_CO2_PgC_cdr_pulse_100yr_UKESM = co2_PgC_cdr_pulse_z1_yr_UKESM(1:100)-co2_PgC_piControl_z1_yr_UKESM(1:100);

% cumulative emissions of CO2 in PgC, year 100
cum_emissions_cdr_pulse_100yr_UKESM = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_UKESM(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_UKESM(1:100))...
                               +abs(co2_PgC_cdr_pulse_z1_yr_UKESM(1:100)-co2_PgC_piControl_z1_yr_UKESM(1:100));
 
cum_emissions_CO2pulse_100yr_UKESM = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_UKESM(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_UKESM(1:100))...
                               +abs(co2_PgC_CO2pulse_z1_yr_UKESM(1:100)-co2_PgC_piControl_z1_yr_UKESM(1:100));
 
% airborne fraction, year 100
airborne_fraction_CO2pulse_UKESM  = atm_CO2_PgC_CO2pulse_100yr_UKESM(1:100)./cum_emissions_CO2pulse_100yr_UKESM(1:100);
airborne_fraction_cdr_pulse_UKESM = atm_CO2_PgC_cdr_pulse_100yr_UKESM(1:100)./cum_emissions_cdr_pulse_100yr_UKESM(1:100);

%% delete variables that are not necessary for following calculations/plots
clearvars -except *ACCESS *MIROC *CanESM* *NorESM *UKESM *UVic *CanOE *_OE *GFDL

%% ACCESS-ESM1-5 
%% NOTE: Pulse starts after 10 years (month 121), data end 90(!) years after pulse. No further data available.
area_o_ACCESS = ncread('areacello_Ofx_ACCESS-ESM1-5_esm-pi-cdr-pulse_r1i1p1f1_gn.nc', 'areacello');
area_a_ACCESS = ncread('areacella_fx_ACCESS-ESM1-5_esm-pi-cdr-pulse_r1i1p1f1_gn.nc', 'areacella');

nbp_cdr_pulse_ACCESS = ncread('nbp_Lmon_ACCESS-ESM1-5_esm-pi-cdr-pulse_r1i1p1f1_gn_027101-037012.nc', 'nbp');
nbp_CO2pulse_ACCESS = ncread('nbp_Lmon_ACCESS-ESM1-5_esm-pi-CO2pulse_r1i1p1f1_gn_027101-037012.nc', 'nbp');
nbp_piControl_ACCESS = ncread('nbp_Lmon_ACCESS-ESM1-5_esm-piControl_r1i1p1f1_gn_027101-057012.nc', 'nbp',[1 1 1],[192 145 1200]);

a2o_flux_cdr_pulse_ACCESS = ncread('fgco2_Omon_ACCESS-ESM1-5_esm-pi-cdr-pulse_r1i1p1f1_gn_027101-037012.nc', 'fgco2');
a2o_flux_CO2pulse_ACCESS = ncread('fgco2_Omon_ACCESS-ESM1-5_esm-pi-CO2pulse_r1i1p1f1_gn_027101-037012.nc', 'fgco2');
a2o_flux_piControl_ACCESS = ncread('fgco2_Omon_ACCESS-ESM1-5_esm-piControl_r1i1p1f1_gn_027101-057012.nc', 'fgco2',[1 1 1],[360 300 1200]);

sftlf_ACCESS = ncread('sftlf_fx_ACCESS-ESM1-5_esm-pi-CO2pulse_r1i1p1f1_gn.nc', 'sftlf');

co2_cdr_pulse_ACCESS = squeeze(ncread('co2_Amon_ACCESS-ESM1-5_esm-pi-cdr-pulse_r1i1p1f1_gn_027101-037012.nc', 'co2',[1 1 1 1],[192 145 1 1200]));
co2_piControl_ACCESS = squeeze(cat(4,ncread('co2_Amon_ACCESS-ESM1-5_esm-piControl_r1i1p1f1_gn_027101-037012.nc', 'co2',[1 1 1 1],[192 145 1 984]),...
    ncread('co2_Amon_ACCESS-ESM1-5_esm-piControl_r1i1p1f1_gn_027101-037012.nc', 'co2',[1 1 1 984],[192 145 1 1]),...
    ncread('co2_Amon_ACCESS-ESM1-5_esm-piControl_r1i1p1f1_gn_027101-037012.nc', 'co2',[1 1 1 984],[192 145 1 1]),...
    ncread('co2_Amon_ACCESS-ESM1-5_esm-piControl_r1i1p1f1_gn_027101-037012.nc','co2',[1 1 1 987],[192 145 1 214])));
co2_CO2pulse_ACCESS = squeeze(cat(4,ncread('co2_Amon_ACCESS-ESM1-5_esm-pi-CO2pulse_r1i1p1f1_gn_027101-037012.nc', 'co2',[1 1 1 1],[192 145 1 1183]),...
    ncread('co2_Amon_ACCESS-ESM1-5_esm-pi-CO2pulse_r1i1p1f1_gn_027101-037012.nc', 'co2',[1 1 1 1183],[192 145 1 1]),...
    ncread('co2_Amon_ACCESS-ESM1-5_esm-pi-CO2pulse_r1i1p1f1_gn_027101-037012.nc', 'co2',[1 1 1 1183],[192 145 1 1]),...
    ncread('co2_Amon_ACCESS-ESM1-5_esm-pi-CO2pulse_r1i1p1f1_gn_027101-037012.nc', 'co2',[1 1 1 1186],[192 145 1 15])));

ppm_to_pgc = 2.124;

%% Ocean C change
year_in_sec = 60 .* 60 .* 24 .* 365;  %unit: s
time_max_2 = length(a2o_flux_cdr_pulse_ACCESS(1,1,:));
a2o_flux_CO2pulse_ts = nan(time_max_2,1);
a2o_flux_cdr_pulse_ts = nan(time_max_2,1);
a2o_flux_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   a2o_flux_CO2pulse_ts(t,1) = nansum(nansum(a2o_flux_CO2pulse_ACCESS(:,:,t).*area_o_ACCESS)); 
   a2o_flux_cdr_pulse_ts(t,1) = nansum(nansum(a2o_flux_cdr_pulse_ACCESS(:,:,t).*area_o_ACCESS)); 
   a2o_flux_piControl_ts(t,1) = nansum(nansum(a2o_flux_piControl_ACCESS(:,:,t).*area_o_ACCESS));
end   

% (yearly mean multiplied by 365)
a2o_flux_CO2pulse_ts_yrly = nan(100,1);
a2o_flux_cdr_pulse_ts_yrly = nan(100,1);
a2o_flux_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  a2o_flux_CO2pulse_ts_yrly(t) = mean(a2o_flux_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  a2o_flux_cdr_pulse_ts_yrly(t) = mean(a2o_flux_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  a2o_flux_piControl_ts_yrly(t) = mean(a2o_flux_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

a2o_flux2_anomaly = a2o_flux_cdr_pulse_ts_yrly-a2o_flux_piControl_ts_yrly;
cum_a2o_flux2_anomaly = a2o_flux2_anomaly;
for i=2:length(a2o_flux2_anomaly)
cum_a2o_flux2_anomaly(i) = a2o_flux2_anomaly(i) + cum_a2o_flux2_anomaly(i-1);
end

a2o_flux2_anomaly_CO2pulse = a2o_flux_CO2pulse_ts_yrly-a2o_flux_piControl_ts_yrly;
cum_a2o_flux2_anomaly_CO2pulse = a2o_flux2_anomaly_CO2pulse;
for i=2:length(a2o_flux2_anomaly_CO2pulse)
cum_a2o_flux2_anomaly_CO2pulse(i) = a2o_flux2_anomaly_CO2pulse(i) + cum_a2o_flux2_anomaly_CO2pulse(i-1);
end

%% nbp 
time_max_2 = length(nbp_cdr_pulse_ACCESS(1,1,:));
nbp_CO2pulse_ts = nan(time_max_2,1);
nbp_cdr_pulse_ts = nan(time_max_2,1);
nbp_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   nbp_CO2pulse_ts(t,1) = nansum(nansum(nbp_CO2pulse_ACCESS(:,:,t).*area_a_ACCESS.*sftlf_ACCESS./100)); 
   nbp_cdr_pulse_ts(t,1) = nansum(nansum(nbp_cdr_pulse_ACCESS(:,:,t).*area_a_ACCESS.*sftlf_ACCESS./100)); 
   nbp_piControl_ts(t,1) = nansum(nansum(nbp_piControl_ACCESS(:,:,t).*area_a_ACCESS.*sftlf_ACCESS./100));
end   

% (yearly mean multiplied by 365)
nbp_CO2pulse_ts_yrly = nan(100,1);
nbp_cdr_pulse_ts_yrly = nan(100,1);
nbp_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  nbp_CO2pulse_ts_yrly(t) = mean(nbp_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_cdr_pulse_ts_yrly(t) = mean(nbp_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_piControl_ts_yrly(t) = mean(nbp_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

Net_Land_Flux2_anomaly_ACCESS = nbp_cdr_pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_ACCESS = Net_Land_Flux2_anomaly_ACCESS;
for i=2:length(Net_Land_Flux2_anomaly_ACCESS)
cum_Net_Land_Flux2_anomaly_ACCESS(i) = Net_Land_Flux2_anomaly_ACCESS(i) + cum_Net_Land_Flux2_anomaly_ACCESS(i-1);
end

Net_Land_Flux2_anomaly_CO2pulse = nbp_CO2pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_CO2pulse_ACCESS = Net_Land_Flux2_anomaly_CO2pulse;
for i=2:length(Net_Land_Flux2_anomaly_CO2pulse)
cum_Net_Land_Flux2_anomaly_CO2pulse_ACCESS(i) = Net_Land_Flux2_anomaly_CO2pulse(i) + cum_Net_Land_Flux2_anomaly_CO2pulse_ACCESS(i-1);
end

%% UPTKAE/RELEASE/AIRBORNE FRACTIONS
cum_Net_Land_Flux_anomaly_ACCESS = cum_Net_Land_Flux2_anomaly_ACCESS;
cum_Net_Land_Flux_anomaly_CO2pulse_ACCESS = cum_Net_Land_Flux2_anomaly_CO2pulse_ACCESS;

Net_Land_Flux_anomaly_cdr_pulse_ACCESS = Net_Land_Flux2_anomaly_ACCESS;
Net_Land_Flux_anomaly_CO2pulse_ACCESS = Net_Land_Flux2_anomaly_CO2pulse;

Ocean_C_Flux_anomaly_cdr_pulse_ACCESS = a2o_flux2_anomaly;
Ocean_C_Flux_anomaly_CO2pulse_ACCESS = a2o_flux2_anomaly_CO2pulse;

cum_Ocean_C_Flux_anomaly_cdr_pulse_ACCESS = cum_a2o_flux2_anomaly;
cum_Ocean_C_Flux_anomaly_CO2pulse_ACCESS = cum_a2o_flux2_anomaly_CO2pulse;

%% co2
% LOWEST layer global average co2
co2_cdr_pulse_z1_ts = nan(1200,1);
max_length=length(co2_cdr_pulse_ACCESS(1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:192)
     for(j=1:145)
         if(~isnan(co2_cdr_pulse_ACCESS(i,j,t))) 
            sum1 = sum1 + (co2_cdr_pulse_ACCESS(i,j,t).*area_a_ACCESS(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_ACCESS(i,j);
         end
     end
    end
   co2_cdr_pulse_z1_ts(t)=sum1./sum2; 
   end
  co2_cdr_pulse_z1_ts_ACCESS = co2_cdr_pulse_z1_ts;
co2_cdr_pulse_z1_PgC_ACCESS = co2_cdr_pulse_z1_ts * ppm_to_pgc;

%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_piControl_ACCESS(1,1,:));
co2_control_z1_ts = nan(1200,1);
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:192)
     for(j=1:145)
         if(~isnan(co2_piControl_ACCESS(i,j,t))) 
            sum1 = sum1 + (co2_piControl_ACCESS(i,j,t).*area_a_ACCESS(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_ACCESS(i,j);
         end
     end
    end
   co2_control_z1_ts(t)=sum1./sum2; 
   end
  co2_control_z1_ts_ACCESS = co2_control_z1_ts;
co2_control_z1_PgC_ACCESS = co2_control_z1_ts * ppm_to_pgc;

%%%%%%%%%% CO2 PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_CO2pulse_ACCESS(1,1,:));
co2_CO2pulse_z1_ts_ACCESS = nan(1200,1);
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:192)
     for(j=1:145)
         if(~isnan(co2_CO2pulse_ACCESS(i,j,t))) 
            sum1 = sum1 + (co2_CO2pulse_ACCESS(i,j,t).*area_a_ACCESS(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_ACCESS(i,j);
         end
     end
    end
   co2_CO2pulse_z1_ts_ACCESS(t)=sum1./sum2; 
  end
co2_CO2pulse_z1_PgC_ACCESS = co2_CO2pulse_z1_ts_ACCESS * ppm_to_pgc;

%%% calculate annual average
co2_PgC_cdr_pulse_z1_yr_ACCESS = nan(100,1);
co2_PgC_CO2pulse_z1_yr_ACCESS = nan(100,1);
co2_PgC_piControl_z1_yr_ACCESS = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  if (m2 <= length(co2_control_z1_PgC_ACCESS))
  co2_PgC_piControl_z1_yr_ACCESS(t) = mean(co2_control_z1_PgC_ACCESS(m:m2)) ;
  end
  if m2<=length(co2_CO2pulse_z1_PgC_ACCESS) 
  co2_PgC_CO2pulse_z1_yr_ACCESS(t) = mean(co2_CO2pulse_z1_PgC_ACCESS(m:m2));
  end
  if m2<=length(co2_cdr_pulse_z1_PgC_ACCESS) 
      co2_PgC_cdr_pulse_z1_yr_ACCESS(t) = mean(co2_cdr_pulse_z1_PgC_ACCESS(m:m2));
  end
end

%% cumulative emissions
% sum of -cum_Ocean_flux, -cum_Land_flux and atmospheric change is defined as
% cumulative emissions at given point in time.

% atmospheric CO2 in PgC
atm_CO2_PgC_CO2pulse_100yr_ACCESS = co2_PgC_CO2pulse_z1_yr_ACCESS(1:100)-co2_PgC_piControl_z1_yr_ACCESS(1:100);
atm_CO2_PgC_cdr_pulse_100yr_ACCESS = co2_PgC_cdr_pulse_z1_yr_ACCESS(1:100)-co2_PgC_piControl_z1_yr_ACCESS(1:100);

% cumulative emissions of CO2 in PgC
cum_emissions_cdr_pulse_100yr_ACCESS = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_ACCESS(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_ACCESS(1:100))...
                               +abs(co2_PgC_cdr_pulse_z1_yr_ACCESS(1:100)-co2_PgC_piControl_z1_yr_ACCESS(1:100));
 
cum_emissions_CO2pulse_100yr_ACCESS = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_ACCESS(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_ACCESS(1:100))...
                               +abs(co2_PgC_CO2pulse_z1_yr_ACCESS(1:100)-co2_PgC_piControl_z1_yr_ACCESS(1:100));
 
% airborne fraction
airborne_fraction_CO2pulse_ACCESS  = atm_CO2_PgC_CO2pulse_100yr_ACCESS./cum_emissions_CO2pulse_100yr_ACCESS;
airborne_fraction_cdr_pulse_ACCESS = atm_CO2_PgC_cdr_pulse_100yr_ACCESS./cum_emissions_cdr_pulse_100yr_ACCESS;

%% delete variables that are not necessary for following calculations/plots
clearvars -except *ACCESS *MIROC *CanESM* *NorESM *UKESM *UVic *CanOE *_OE *GFDL

%% GFDL-ESM4
area_o_GFDL = ncread('areacello_Ofx_GFDL-ESM4_piControl_r1i1p1f1_gr.nc', 'areacello');
area_a_GFDL = ncread('areacella_fx_GFDL-ESM4_piControl_r1i1p1f1_gr1.nc', 'areacella');

nbp_cdr_pulse_GFDL = ncread('nbp_Lmon_GFDL-ESM4_esm-pi-cdr-pulse_r1i1p1f1_gr1_010101-020012.nc', 'nbp');
nbp_CO2pulse_GFDL = ncread('nbp_Lmon_GFDL-ESM4_esm-pi-CO2pulse_r1i1p1f1_gr1_010101-020012.nc', 'nbp');
nbp_piControl_GFDL = ncread('nbp_Lmon_GFDL-ESM4_esm-piControl_r1i1p1f1_gr1_010101-020012.nc', 'nbp');

a2o_flux_cdr_pulse_GFDL = ncread('fgco2_Oyr_GFDL-ESM4_esm-pi-cdr-pulse_r1i1p1f1_gr_0101-0200.nc', 'fgco2');
a2o_flux_CO2pulse_GFDL = ncread('fgco2_Oyr_GFDL-ESM4_esm-pi-CO2pulse_r1i1p1f1_gr_0101-0200.nc', 'fgco2');
a2o_flux_piControl_GFDL = ncread('fgco2_Oyr_GFDL-ESM4_esm-piControl_r1i1p1f1_gr_0101-0200.nc', 'fgco2');

sftlf_GFDL = ncread('sftlf_fx_GFDL-ESM4_historical_r1i1p1f1_gr1.nc', 'sftlf');

co2_cdr_pulse_GFDL = ncread('co2_Amon_GFDL-ESM4_esm-pi-cdr-pulse_r1i1p1f1_gr1_010101-020012.nc', 'co2',[1 1 1 1],[288 180 1 1200]);
co2_CO2pulse_GFDL  = ncread('co2_Amon_GFDL-ESM4_esm-pi-CO2pulse_r1i1p1f1_gr1_010101-020012.nc', 'co2',[1 1 1 1],[288 180 1 1200]);
co2_piControl_GFDL = ncread('co2_Amon_GFDL-ESM4_esm-piControl_r1i1p1f1_gr1_010101-020012.nc', 'co2',[1 1 1 1],[288 180 1 1176]);

ppm_to_pgc = 2.124;

%% Ocean C change
year_in_sec = 60 .* 60 .* 24 .* 365;  %unit: s
time_max_2 = length(a2o_flux_cdr_pulse_GFDL(1,1,:));
a2o_flux_CO2pulse_ts = nan(time_max_2,1);
a2o_flux_cdr_pulse_ts = nan(time_max_2,1);
a2o_flux_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   a2o_flux_CO2pulse_ts(t,1) = nansum(nansum(a2o_flux_CO2pulse_GFDL(:,:,t).*area_o_GFDL)); 
   a2o_flux_cdr_pulse_ts(t,1) = nansum(nansum(a2o_flux_cdr_pulse_GFDL(:,:,t).*area_o_GFDL)); 
   a2o_flux_piControl_ts(t,1) = nansum(nansum(a2o_flux_piControl_GFDL(:,:,t).*area_o_GFDL));
end   
a2o_flux_CO2pulse_ts_PgC = a2o_flux_CO2pulse_ts ./ 10^12 * year_in_sec; 
a2o_flux_cdr_pulse_ts_PgC = a2o_flux_cdr_pulse_ts ./ 10^12 * year_in_sec; 
a2o_flux_piControl_ts_PgC = a2o_flux_piControl_ts ./ 10^12 * year_in_sec; 

% TIMESERIES ANOMALY
C_change_ocean_100yr = sum(a2o_flux_cdr_pulse_ts_PgC(1:100)-a2o_flux_piControl_ts_PgC(1:100));
C_change_ocean_100yr_CO2pulse = sum(a2o_flux_CO2pulse_ts_PgC(1:100)-a2o_flux_piControl_ts_PgC(1:100));

%% nbp 
time_max_2 = length(nbp_cdr_pulse_GFDL(1,1,:));
nbp_CO2pulse_ts = nan(time_max_2,1);
nbp_cdr_pulse_ts = nan(time_max_2,1);
nbp_piControl_ts = nan(time_max_2,1);

for(t=1:time_max_2)
   nbp_CO2pulse_ts(t,1) = nansum(nansum(nbp_CO2pulse_GFDL(:,:,t).*area_a_GFDL.*sftlf_GFDL./100)); 
   nbp_cdr_pulse_ts(t,1) = nansum(nansum(nbp_cdr_pulse_GFDL(:,:,t).*area_a_GFDL.*sftlf_GFDL./100)); 
   nbp_piControl_ts(t,1) = nansum(nansum(nbp_piControl_GFDL(:,:,t).*area_a_GFDL.*sftlf_GFDL./100));
end   

% (yearly mean multiplied by 365)
nbp_CO2pulse_ts_yrly = nan(100,1);
nbp_cdr_pulse_ts_yrly = nan(100,1);
nbp_piControl_ts_yrly = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  nbp_CO2pulse_ts_yrly(t) = mean(nbp_CO2pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_cdr_pulse_ts_yrly(t) = mean(nbp_cdr_pulse_ts(m:m2)) ./ 10^12 * year_in_sec;
  nbp_piControl_ts_yrly(t) = mean(nbp_piControl_ts(m:m2)) ./ 10^12 * year_in_sec;
end

Net_Land_Flux2_anomaly = nbp_cdr_pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly = Net_Land_Flux2_anomaly;
for i=2:length(Net_Land_Flux2_anomaly)
cum_Net_Land_Flux2_anomaly(i) = Net_Land_Flux2_anomaly(i) + cum_Net_Land_Flux2_anomaly(i-1);
end

Net_Land_Flux2_anomaly_CO2pulse = nbp_CO2pulse_ts_yrly-nbp_piControl_ts_yrly;
cum_Net_Land_Flux2_anomaly_CO2pulse = Net_Land_Flux2_anomaly_CO2pulse;
for i=2:length(Net_Land_Flux2_anomaly_CO2pulse)
cum_Net_Land_Flux2_anomaly_CO2pulse(i) = Net_Land_Flux2_anomaly_CO2pulse(i) + cum_Net_Land_Flux2_anomaly_CO2pulse(i-1);
end

%% UPTKAE/RELEASE/AIRBORNE FRACTIONS
cum_Net_Land_Flux_anomaly_GFDL = cum_Net_Land_Flux2_anomaly;
cum_Net_Land_Flux_anomaly_CO2pulse_GFDL = cum_Net_Land_Flux2_anomaly_CO2pulse;

Ocean_C_Flux_anomaly_cdr_pulse_GFDL = a2o_flux_cdr_pulse_ts_PgC(1:100)-a2o_flux_piControl_ts_PgC(1:100);
Ocean_C_Flux_anomaly_CO2pulse_GFDL = a2o_flux_CO2pulse_ts_PgC(1:100)-a2o_flux_piControl_ts_PgC(1:100);

cum_Ocean_C_Flux_anomaly_cdr_pulse_GFDL = Ocean_C_Flux_anomaly_cdr_pulse_GFDL;
cum_Ocean_C_Flux_anomaly_CO2pulse_GFDL = Ocean_C_Flux_anomaly_CO2pulse_GFDL;

for i=2:length(Ocean_C_Flux_anomaly_cdr_pulse_GFDL)
  cum_Ocean_C_Flux_anomaly_cdr_pulse_GFDL(i) = Ocean_C_Flux_anomaly_cdr_pulse_GFDL(i) + cum_Ocean_C_Flux_anomaly_cdr_pulse_GFDL(i-1);
  cum_Ocean_C_Flux_anomaly_CO2pulse_GFDL(i) = Ocean_C_Flux_anomaly_CO2pulse_GFDL(i) + cum_Ocean_C_Flux_anomaly_CO2pulse_GFDL(i-1);
end

Net_Land_Flux_anomaly_cdr_pulse_GFDL = Net_Land_Flux2_anomaly;
Net_Land_Flux_anomaly_CO2pulse_GFDL = Net_Land_Flux2_anomaly_CO2pulse;

%% atmospheric CO2 concentration to mass
%%%%%%%%%% CONTROL %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_piControl_GFDL(1,1,1,:));
co2_control_z1_ts = nan(1200,1);
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:288)
     for(j=1:180)
         if(~isnan(co2_piControl_GFDL(i,j,1,t))) 
            sum1 = sum1 + (co2_piControl_GFDL(i,j,1,t).*area_a_GFDL(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_GFDL(i,j);
         end
     end
    end
   co2_control_z1_ts(t)=sum1./sum2; 
   end
  co2_control_z1_ts_GFDL = co2_control_z1_ts;
co2_control_z1_PgC_GFDL = co2_control_z1_ts * ppm_to_pgc;

%%%%%%%%%% CO2 PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
max_length=length(co2_CO2pulse_GFDL(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:288)
     for(j=1:180)
         if(~isnan(co2_CO2pulse_GFDL(i,j,1,t))) 
            sum1 = sum1 + (co2_CO2pulse_GFDL(i,j,1,t).*area_a_GFDL(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_GFDL(i,j);
         end
     end
    end
   co2_CO2pulse_z1_ts(t)=sum1./sum2; 
  end
co2_CO2pulse_z1_PgC_GFDL = co2_CO2pulse_z1_ts * ppm_to_pgc;

%%%%%%%%%% CDR PULSE %%%%%%%%%%%%%%%%%%
% LOWEST layer global average co2
co2_cdr_pulse_z1_ts = nan(1200,1);
max_length=length(co2_cdr_pulse_GFDL(1,1,1,:));
   for(t=1:max_length)
sum1 = 0;
sum2 = 0;
    for(i=1:288)
     for(j=1:180)
         if(~isnan(co2_cdr_pulse_GFDL(i,j,1,t))) 
            sum1 = sum1 + (co2_cdr_pulse_GFDL(i,j,1,t).*area_a_GFDL(i,j))*10^6; %original unit: mol/mol
            sum2 = sum2 + area_a_GFDL(i,j);
         end
     end
    end
   co2_cdr_pulse_z1_ts(t)=sum1./sum2; 
   end
  co2_cdr_pulse_z1_ts_GFDL = co2_cdr_pulse_z1_ts;
co2_cdr_pulse_z1_PgC_GFDL = co2_cdr_pulse_z1_ts * ppm_to_pgc;

%%% calculate annual average
co2_PgC_cdr_pulse_z1_yr_GFDL = nan(100,1);
co2_PgC_CO2pulse_z1_yr_GFDL = nan(100,1);
co2_PgC_piControl_z1_yr_GFDL = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  if (m2 <= length(co2_control_z1_PgC_GFDL))
  co2_PgC_piControl_z1_yr_GFDL(t) = mean(co2_control_z1_PgC_GFDL(m:m2));
  end
  if m2<=length(co2_CO2pulse_z1_PgC_GFDL) 
  co2_PgC_CO2pulse_z1_yr_GFDL(t) = mean(co2_CO2pulse_z1_PgC_GFDL(m:m2));
  end
  if m2<=length(co2_cdr_pulse_z1_PgC_GFDL) 
      co2_PgC_cdr_pulse_z1_yr_GFDL(t) = mean(co2_cdr_pulse_z1_PgC_GFDL(m:m2));
  end
end

co2_PgC_atm_anomaly_CO2pulse_GFDL = co2_PgC_CO2pulse_z1_yr_GFDL-co2_PgC_piControl_z1_yr_GFDL;
co2_PgC_atm_anomaly_cdr_pulse_GFDL = co2_PgC_cdr_pulse_z1_yr_GFDL-co2_PgC_piControl_z1_yr_GFDL;

%% cumulative emissions
% sum of -cum_Ocean_flux, -cum_Land_flux and atmospheric change is defined as
% cumulative emissions at given point in time.

% atmospheric CO2 in PgC
atm_CO2_PgC_CO2pulse_100yr_GFDL = co2_PgC_CO2pulse_z1_yr_GFDL(1:100)-co2_PgC_piControl_z1_yr_GFDL(1:100);
atm_CO2_PgC_cdr_pulse_100yr_GFDL = co2_PgC_cdr_pulse_z1_yr_GFDL(1:100)-co2_PgC_piControl_z1_yr_GFDL(1:100);

% cumulative emissions of CO2 in PgC
cum_emissions_cdr_pulse_100yr_GFDL = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_GFDL(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_GFDL(1:100))...
                               +abs(co2_PgC_cdr_pulse_z1_yr_GFDL(1:100)-co2_PgC_piControl_z1_yr_GFDL(1:100));
 
cum_emissions_CO2pulse_100yr_GFDL = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_GFDL(1:100))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_GFDL(1:100))...
                               +abs(co2_PgC_CO2pulse_z1_yr_GFDL(1:100)-co2_PgC_piControl_z1_yr_GFDL(1:100));

%% CO2 only until year 98 after pulse
% atmospheric CO2 in PgC
atm_CO2_PgC_CO2pulse_98yr_GFDL = co2_PgC_CO2pulse_z1_yr_GFDL(1:98)-co2_PgC_piControl_z1_yr_GFDL(1:98);
atm_CO2_PgC_cdr_pulse_98yr_GFDL = co2_PgC_cdr_pulse_z1_yr_GFDL(1:98)-co2_PgC_piControl_z1_yr_GFDL(1:98);

% cumulative emissions of CO2 in PgC
cum_emissions_cdr_pulse_98yr_GFDL = abs(cum_Ocean_C_Flux_anomaly_cdr_pulse_GFDL(1:98))...
                               +abs(cum_Net_Land_Flux_anomaly_GFDL(1:98))...
                               +abs(co2_PgC_cdr_pulse_z1_yr_GFDL(1:98)-co2_PgC_piControl_z1_yr_GFDL(1:98));
 
cum_emissions_CO2pulse_98yr_GFDL = abs(cum_Ocean_C_Flux_anomaly_CO2pulse_GFDL(1:98))...
                               +abs(cum_Net_Land_Flux_anomaly_CO2pulse_GFDL(1:98))...
                               +abs(co2_PgC_CO2pulse_z1_yr_GFDL(1:98)-co2_PgC_piControl_z1_yr_GFDL(1:98));

%% delete variables that are not necessary for following calculations/plots
clearvars -except *ACCESS *MIROC *CanESM* *NorESM *UKESM *UVic *CanOE *_OE *GFDL

%% FIGURE TIMESERIES
%% CO2 Concentration, Land C Change, Ocean C Change
co2_CanESM = co2_cdr_pulse_z1_ts_CanESM - co2_control_z1_ts_CanESM + co2_control_z1_ts_CanESM(1);
co2_MIROC = co2_cdr_pulse_z1_ts_MIROC - co2_control_z1_ts_MIROC + co2_control_z1_ts_MIROC(1);
co2_GFDL = co2_cdr_pulse_z1_ts_GFDL - co2_control_z1_ts_GFDL + co2_control_z1_ts_GFDL(1);
co2_UKESM = co2_cdr_pulse_z1_ts_UKESM - co2_control_z1_ts_UKESM + co2_control_z1_ts_UKESM(1);
co2_CanOE = co2_cdr_pulse_z1_ts_CanESM5_CanOE - co2_control_z1_ts_CanESM5_CanOE + co2_control_z1_ts_CanESM5_CanOE(1);
co2_NorESM = co2_cdr_pulse_z1_ts_NorESM - co2_control_z1_ts_NorESM + co2_control_z1_ts_NorESM(1);
co2_ACCESS = co2_cdr_pulse_z1_ts_ACCESS - co2_control_z1_ts_ACCESS + co2_control_z1_ts_ACCESS(1);

% calculate yearly mean
co2_MIROC_yr = nan(100,1);
co2_CanESM_yr = nan(100,1);
co2_GFDL_yr = nan(100,1);
co2_UKESM_yr = nan(100,1);
co2_NorESM_yr = nan(100,1);
co2_CanOE_yr = nan(100,1);
co2_ACCESS_yr = nan(100,1);
for(t=1:100)
    m = (t-1)*12+1;
    m2 = m + 11;
  co2_MIROC_yr(t) = mean(co2_MIROC(m:m2));
  co2_CanESM_yr(t) = mean(co2_CanESM(m:m2));
  co2_GFDL_yr(t) = mean(co2_GFDL(m:m2));   % NOTE: DATA END IN YEAR 98!
  co2_UKESM_yr(t) = mean(co2_UKESM(m:m2));
  co2_NorESM_yr(t) = mean(co2_NorESM(m:m2));
  co2_CanOE_yr(t) = mean(co2_CanOE(m:m2));
  co2_ACCESS_yr(t) = mean(co2_ACCESS(m:m2));
end

%% IPCC Figure
width_singlecolumn = 3.5;
width_doublecolumn = 7.08;
max_height = 10;

fig = figure;
set(gcf,'units','inch','position',[0,0,width_singlecolumn,max_height])
xwidth = 0.16;
xtitle = -0.2;

subplot(3,1,1)
hold on
text(0.67,0.52,'MIROC-ES2L','Units','normalized','FontSize',11, 'FontName', 'Arial','Color',[0.721568627 0.37254902 0.71372549])
text(0.67,0.42,'CanESM5','Units','normalized','FontSize',11, 'FontName', 'Arial','Color',[0.117647059 0.298039216 0.141176471])
text(0.67,0.32,'GFDL-ESM4','Units','normalized','FontSize',11, 'FontName', 'Arial','Color',[0.137254902 0.211764706 0.42745098])
text(0.67,0.22,'UKESM1-0-LL','Units','normalized','FontSize',11, 'FontName', 'Arial','Color',[0.643137255 0.235294118 0.439215686])
text(0.67,0.11,'UVic ESCM','Units','normalized','FontSize',11, 'FontName', 'Arial','Color',[0.9290, 0.6940, 0.1250])
text(0.15,0.11,'CanESM5-CanOE','Units','normalized','FontSize',11, 'FontName', 'Arial','Color',[0.117647059	0.298039216	0.141176471])
text(0.15,0.33,'NorESM2-LM','Units','normalized','FontSize',11, 'FontName', 'Arial','Color',[0.945098039	0.22745098	0.654901961])
text(0.15,0.22,'ACCESS-ESM1-5','Units','normalized','FontSize',11, 'FontName', 'Arial','Color',[0	0.690196078	0.31372549])

plot(0:100,[co2_control_z1_ts_MIROC(1);co2_MIROC_yr],'LineWidth',1,'Color', [0.721568627 0.37254902 0.71372549])
plot(0:100,[co2_control_z1_ts_CanESM(1);co2_CanESM_yr],'LineWidth',1,'Color', [0.117647059 0.298039216 0.141176471])
plot(0:100,[co2_control_z1_ts_CanESM5_CanOE(1);co2_CanOE_yr],'LineWidth',1,'Color', [0.117647059 0.298039216 0.141176471])
plot(0:100,[co2_control_z1_ts_GFDL(1);co2_GFDL_yr],'LineWidth',1,'Color', [0.137254902 0.211764706 0.42745098]) % NOTE: data only UNTIL YEAR 98! (in removal fraction in year 100, control of year 98 was used)
plot(0:100,[co2_control_z1_ts_UKESM(1);co2_UKESM_yr],'LineWidth',1,'Color', [0.643137255 0.235294118 0.439215686]) % NOTE: CO2 conc. calculated from co2mass
plot(0:100,[co2_control_z1_ts_NorESM(1);co2_NorESM_yr],'LineWidth',1,'Color', [0.945098039	0.22745098	0.654901961]) % NOTE: CO2 conc. calculated from co2mass
plot(0:100,[co2_control_z1_ts_ACCESS(1);co2_ACCESS_yr(11:100);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN],'LineWidth',1,'Color', [0	0.690196078	0.31372549]) % NOTE: CO2 conc. calculated from co2mass
plot(0:100,[A_co2_spinup_UVic(250);A_co2_controlled_UVic(1:100)],'LineWidth',1,'Color',[0.9290, 0.6940, 0.1250]) 
xlabel('')
grid off
ylim([230 290])
set(gcf,'color','w')
set(gca,'TickDir','out');
set(gca,'FontSize',11)
set(gca, 'FontName', 'Arial')
ylabel('(ppm)')
text(xtitle,1.12,'(a)  Atmospheric CO_2 concentration','Units','normalized','FontSize',11, 'FontName', 'Arial')

    pos = get(gca, 'Position');
    pos(1) = xwidth;
    pos(4) = 0.2;
    set(gca, 'Position', pos)

subplot(3,1,2)
hold on
plot(0:100,[0;cum_Net_Land_Flux_anomaly_MIROC],'LineWidth',1,'Color', [0.721568627 0.37254902 0.71372549])
plot(0:100,[0;cum_Net_Land_Flux_anomaly_CanESM],'LineWidth',1,'Color', [0.117647059 0.298039216 0.141176471])
plot(0:100,[0;cum_Net_Land_Flux_anomaly_CanESM5_CanOE],'LineWidth',1,'Color', [0.117647059 0.298039216 0.141176471])
plot(0:100,[0;cum_Net_Land_Flux_anomaly_GFDL],'LineWidth',1,'Color', [0.137254902 0.211764706 0.42745098]) % NOTE: Land-to-Ocean-C flux not accounted for
plot(0:100,[0;cum_Net_Land_Flux_anomaly_UKESM],'LineWidth',1,'Color', [0.643137255 0.235294118 0.439215686]) % NOTE: Land-to-Ocean-C flux not accounted for
plot(0:100,[0;cum_Net_Land_Flux_anomaly_NorESM],'LineWidth',1,'Color', [0.945098039	0.22745098	0.654901961]) % NOTE: Land-to-Ocean-C flux not accounted for
plot(0:100,[0;cum_Net_Land_Flux_anomaly_ACCESS(11:100);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN],'LineWidth',1,'Color', [0	0.690196078	0.31372549]) % NOTE: Land-to-Ocean-C flux not accounted for
plot(0:100,[0;a2l_cum_cdr_controlled_UVic(1:100)],'LineWidth',1,'Color', [0.9290, 0.6940, 0.1250]) 
xlabel('')
ylim([-80 0])
set(gca,'TickDir','out');
set(gca,'FontSize',11)
set(gca, 'FontName', 'Arial')
text(xtitle,1.12,'(b)  Land carbon change','Units','normalized','FontSize',11, 'FontName', 'Arial')
ylabel('(PgC)')

    pos = get(gca, 'Position');
    pos(1) = xwidth;
    pos(4) = 0.2;
    set(gca, 'Position', pos)

subplot(3,1,3)
hold on
plot(0:100,[0;cum_Ocean_C_Flux_anomaly_cdr_pulse_MIROC],'LineWidth',1,'Color', [0.721568627 0.37254902 0.71372549])
plot(0:100,[0;cum_Ocean_C_Flux_anomaly_cdr_pulse_CanESM],'LineWidth',1,'Color', [0.117647059 0.298039216 0.141176471])
plot(0:100,[0;cum_Ocean_C_Flux_anomaly_cdr_pulse_CanESM5_CanOE],'LineWidth',1,'Color', [0.117647059 0.298039216 0.141176471])
plot(0:100,[0;cum_Ocean_C_Flux_anomaly_cdr_pulse_GFDL],'LineWidth',1,'Color', [0.137254902 0.211764706 0.42745098])
plot(0:100,[0;cum_Ocean_C_Flux_anomaly_cdr_pulse_UKESM],'LineWidth',1,'Color', [0.643137255 0.235294118 0.439215686])
plot(0:100,[0;cum_Ocean_C_Flux_anomaly_cdr_pulse_NorESM],'LineWidth',1,'Color', [0.945098039	0.22745098	0.654901961])
plot(0:100,[0;cum_Ocean_C_Flux_anomaly_cdr_pulse_ACCESS(11:100);NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN],'LineWidth',1,'Color', [0	0.690196078	0.31372549])
plot(0:100,[0;a2o_cum_cdr_controlled_UVic(1:100)],'LineWidth',1,'Color', [0.9290, 0.6940, 0.1250])
xlabel('Years since pulse removal')
ylabel('Ocean C Change (PgC)')
set(gca,'TickDir','out');
ylim([-50 0])
set(gca,'FontSize',11)
set(gca, 'FontName', 'Arial')
text(xtitle,1.12,'(c)  Ocean carbon change','Units','normalized','FontSize',11, 'FontName', 'Arial')
ylabel('(PgC)')

    pos = get(gca, 'Position');
    pos(1) = xwidth;
    pos(4) = 0.2;
    set(gca, 'Position', pos)

    