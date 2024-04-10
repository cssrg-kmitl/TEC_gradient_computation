clear
clc

%put the grad results here
load('E:\GBAS_SF_NICE_program\Results\GRAD_KMIT-STFD_2019_01_07.mat')
%put the grad results here

window_for_std = 100*2; %left plus right
sample_test = nan(window_for_std+1,32);
std_grad = nan(86400,32);
for ind_test = window_for_std/2+1:86399-window_for_std/2+1
    sample_test = output_PRM.ion_fix_mm_km(ind_test-window_for_std/2:ind_test+window_for_std/2,:);
    std_all_sat = nanstd(sample_test);
    std_grad(ind_test,:) = std_all_sat;
end
start_azi = -170;
stop_azi = 180;
start_ele = 90;
stop_ele = 5;
step_size = 5;
value_for_round = start_ele:-step_size:stop_ele;
% elevation = round(prm_2019_01_09.elevation,-1);
elevation = interp1(value_for_round,value_for_round,output_PRM.elevation,'nearest','extrap');
value_for_round = start_azi:step_size:stop_azi;
% azi = round(prm_2019_01_09.azi,-1);
azi = interp1(value_for_round,value_for_round,output_PRM.azi,'nearest','extrap');
elevation_gpu = gpuArray(elevation);
azi_gpu = gpuArray(azi);
tic
x_plot = [360+start_azi:step_size:(360-step_size) 0:step_size:stop_azi];
y_plot = start_ele:-(step_size):stop_ele;
std_grad_ele_azi = nan(length(y_plot),length(x_plot));
find_ind = 1;
for azi_std = start_azi:step_size:stop_azi
    for ele_std = start_ele:-(step_size):stop_ele
        std_grad_ele_azi(find_ind) = median(std_grad(elevation == ele_std & azi == azi_std),'omitnan');
        find_ind = find_ind + 1;
    end
end
disp('CPU time')
toc
% h = heatmap(x_plot,y_plot,std_MP1_ele_azi,'colormap',turbo);
% figure
h = heatmap(x_plot,y_plot,std_grad_ele_azi,'colormap',turbo);
xlabel('Azimuth angle (degree)')
ylabel('Elevation angle (degree)')
h.ColorLimits = [0 prctile(std_grad_ele_azi(:),95)];
%for skyplot
std_MP1_ele_azi_hortz = std_grad_ele_azi(:);
ind_notnan = ~isnan(std_MP1_ele_azi_hortz);
std_MP1_ele_azi_hortz = std_MP1_ele_azi_hortz(ind_notnan);
% max_std = max(std_MP1_ele_azi(:));
% min_std = min(std_MP1_ele_azi(:));
max_std = 1;
min_std = 0.01;
step_size = (max_std - min_std)/6;
satAz = (x_plot'*ones(1,length(y_plot)))';
satAz = satAz(:);
satAz = satAz(ind_notnan);
satEl = y_plot'*ones(1,length(x_plot));
satEl = satEl(:);
satEl = satEl(ind_notnan);
value_for_round = min_std:step_size:max_std;
vRounded = interp1(value_for_round,value_for_round,std_MP1_ele_azi_hortz,'nearest','extrap');
% skyplot(satAz,satEl,'GroupData',constellationGroup) 
constellationGroup = categorical(cellstr(num2str(vRounded)));
figure
newcolors = {'c','b','g','y','r','m','k'};
colororder(newcolors)
skyplot(satAz,satEl,'GroupData',constellationGroup)
legend(num2str(value_for_round'))