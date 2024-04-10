function [time_c1,time_p2,time_l1,time_l2,truerange_sat,elevation_sat,azi_sat] = obstorange(obs,nav,elev_mask,rcv_pos)
% Calculate Total Electron Content (TEC)
% Inputs:
%        obs     = observation data
%        nav     = navigation data
%        satb    = Satellite bias
%        S_path  = results path
% Save files:
%       TEC      = Total Electron Content (STEC VTEC STEC_with_rcvbias STEC_with_bias STECp STECl )
%       ROTI     = Rate Of Change TEC Index
%       DCB      = Satellite and receiver Differential Code Bias (DCB)
%       prm      = parameters (elevation angle)
%       refpos   = reference position


% ======================================
% Ref position (Read PPP from PPPindex.txt)
f1 = 1575.42*10^6;          %   f1 = 1575.42 MHz (L1)
f2 = 1227.60*10^6;          %   f2 = 1227.60 MHz (L2)
c  = 299792458;             %   light speed = 299792458 m/s
lambda1 = c/f1;             %   wave length of f1
lambda2 = c/f2;             %   wave length of f2
center_E = [0 0 0];         %   Center of Earth
% 1. Prepare matrix
% NaN
Times                = nan(86400,32); % Time
time_elevation       = nan(86400,32); % elevation angle
azi_sat              = nan(86400,32); % azi angle
time_c1              = nan(86400,32);
time_p2              = nan(86400,32);
time_l1              = nan(86400,32);
time_l2              = nan(86400,32);
truerange_sat        = nan(86400,32);

% 2. TEC calculation
%== Satellte index
Sat_obs = unique(nav.index);
lla = ecef2lla(rcv_pos');
Rcv_lat  = lla(1);
Rcv_long = lla(2);
Rcv_h    = lla(3);
R = [-sind(Rcv_long)                cosd(Rcv_long)                            0;
    -sind(Rcv_lat)*cosd(Rcv_long) -sind(Rcv_lat)*sind(Rcv_long)  cosd(Rcv_lat);
    cosd(Rcv_lat)*cosd(Rcv_long)  cosd(Rcv_lat)*sind(Rcv_long)  sind(Rcv_lat)];
for i = 1 : length(Sat_obs) % GPS 1 - 32
    disp(['PRN# ... ' num2str(Sat_obs(i)) ' ...'])
    PRN = Sat_obs(i);
    Sat  = find(obs.index == PRN);
    if(~isempty(Sat))
        Time = obs.epoch(Sat)+((obs.date(4)*60*60)+(obs.date(5)*60)+obs.date(6));
        
        % 2.1 Read pseudorange from observation file
        % Pseudorange C/A code    :L1   (m)
        C1   = obs.data(Sat,ismember(obs.type,'C1'));
        % Pseudorange P code      :L2   (m)
        P2 = obs.data(Sat,ismember(obs.type,'P2'));
        % Carrier Phase in length :L1   (m)
        L1 = lambda1*obs.data(Sat,ismember(obs.type,'L1'));
        % Carrier Phase in length :L2   (m)
        L2 = lambda2*obs.data(Sat,ismember(obs.type,'L2'));
        
        % 2.2 Calculate elevation angle
        [satpos,~]  = satpos_xyz_sbias(Time,PRN,nav.eph,nav.index,C1);
        vector_s  = satpos-rcv_pos';
        truerange = sqrt(sum(vector_s'.^2));
        vector_r2 = rcv_pos'-center_E;
        vector_r  = repmat(vector_r2,length(vector_s),1);
        % elevation angle
        time_elevation(Time+1,PRN) = 90-acosd(dot(vector_s,vector_r,2)./(vecnorm(vector_s')'...
            .*vecnorm(vector_r')'));
        Xs = satpos;
        Xr = repmat(rcv_pos',length(vector_s),1);
        Rs = [Xs(:,1)-Xr(:,1) Xs(:,2)-Xr(:,2) Xs(:,3)-Xr(:,3)];
        RL = (R*Rs')';
        Xl = RL(:,1);
        Yl = RL(:,2);
        azi_sat(Time+1,PRN) = atan2(Xl,Yl)*180/pi;
        
        Times(Time+1,PRN) = Time+1;
        time_c1(Time+1,PRN)= C1;
        time_p2(Time+1,PRN)= P2;
        time_l1(Time+1,PRN)= L1;
        time_l2(Time+1,PRN)= L2;
        truerange_sat(Time+1,PRN)= truerange;
    end
end

% 2.4 elevation angle cutoff <15 degree
mask                 = time_elevation;
mask(mask<elev_mask) = NaN;
mask(~isnan(mask))   = 1;
% angle cut STEC
elevation_sat = mask.*time_elevation;
time_c1= mask.*time_c1;
time_p2= mask.*time_p2;
time_l1= mask.*time_l1;
time_l2= mask.*time_l2;
truerange_sat=mask.*truerange_sat;
end

