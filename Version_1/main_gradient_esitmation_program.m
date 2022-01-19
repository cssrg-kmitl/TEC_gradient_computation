%% Ionospheric delay gradient estimation (grad) from RINEX 2.11
% Estimation based on single and dual-frequency approachs (GPS)
% Original by Jirapoom budtho
% Version 1.00 
% (2021-06-24) - First version
% 
% 1. Main program is main_gradient_esitmation_program.m
% 2. Put the rinex files in \RINEX_files
% 3. The gradient results are in \Results
% 
% 4. We have laboratory website, you can visit
% - http://iono-gnss.kmitl.ac.th/
% =================================================
% CSSRG Laboratory
% School of Telecommunication Engineering
% Faculty of Engineering
% King Mongkut's Institute of Technology Ladkrabang
% Bangkok, Thailand
% =================================================
% Output : data 1 day 
% output_PRM.elevation    = elevation angle
% output_PRM.azi          = azimuth angle
% output_PRM.TImes        = UTC hour time
% output_PRM.ion_fix_mm_km= ionospheric delay gradient (mmn\km)

clear
clc
p_path = [pwd '\'];
file_path = [p_path 'RINEX_files\'];
BaseObsfile = [file_path 'KMIT0070.19o'];
RoverObsfile = [file_path 'STFD0070.19o'];
Navfile = [file_path 'brdc0070.19n'];
cd functions

%please obtain position from PPP or RTK approach
AER1_pos = [-1157218.063;6088967.72;1500160.56]; %rtkposition
AER2_pos = [-1156298.928;6090212.175;1495840.958]; %rtkposition
AERO_pos = [-1157218.154;6088968.441;1500160.563]; %rtkposition
STFD_pos = [-1146423.909;6089932.380;1504580.968]; %rtkposition
KMIT_pos = [-1158319.147;6087918.882;1503747.417]; %rtkposition

%set up the estimation parameters here
input_PRM.mode = 1; %mode = 1 is a single-frequency, mode = 2 is a dual-frequency
input_PRM.elevation = 30; %elevation is the elevation cut-off angle
input_PRM.base_pos = KMIT_pos;
input_PRM.rover_pos = STFD_pos;
input_PRM.p_path = p_path;
input_PRM.S_path = [p_path 'Results\'];   % Results path
input_PRM.DCB_path = [p_path 'DCB\'];   % DCB path

time_index = (0:86399)./86400*24;
[output_PRM,satava_s,doy,Year,month,date,name] = sfgrad_no_rm_sat(BaseObsfile,RoverObsfile,Navfile,input_PRM);
plot(time_index,output_PRM.ion_fix_mm_km,'DisplayName','ion_fix_mm_km_KMST')
xlabel('UTC Time (Hr)')
ylabel('Ionospheric delay gradient (mm/km)')
cd ..
%% Save file
filename = [input_PRM.S_path 'GRAD_' name.basename '-' name.rovername '_' Year '_' month '_' date];
save(filename,'output_PRM')