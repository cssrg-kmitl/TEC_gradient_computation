% clear
% clc
% load('before KALMAN 174')
Re = 6378137;
mf = sqrt(1-((Re*cosd(elevation_sat))/(Re+350000)).^2);
ion_fix_mf = ion_fix_mm_km.*mf';
% fig = plot(timeindex,ion_fix,'DisplayName','timeindex');
% title('DOY 174 2015')
% ylabel('ionospheric delay gradient (mm/km)')
% xlabel('UTC Time (hr)')
% figure
% fig = plot(timeindex,ion_fix_mf,'DisplayName','timeindex');
% title('DOY 174 2015')
% ylabel('ionospheric delay gradient (mm/km)')
% xlabel('UTC Time (hr)')
indd = find(ratio(:)>2);
ion_plo = nan(32,86400);
ion_plo(:,indd) = ion_fix_mf(:,indd);
set_left = satava_s - sat_remove;
% ratio_num = length(find(ratio(:)>2))/86400*100;
% fig = plot(timeindex,ion_plo,'DisplayName','timeindex');
% title(['DOY 002 2015 with ' num2str(elev_mask) ' degree'])
% ylabel('ionospheric delay gradient (mm/km)')
% xlabel('UTC Time (hr)')
% ylim([-100 100])
% dim = [.6 .6 .3 .3];
% str = ['Number of cycle slip: ' num2str(slip_count)];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% dim = [.15 .6 .3 .3];
% str = ['Fixed rate by the ratio test: ' num2str(ratio_num) '%'];
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% figure
% fig2 = plot(timeindex,satava_s,'DisplayName','timeindex');
% title(['DOY 043 2015 with ' num2str(elev_mask) ' degree'])
% ylabel('Number of satellite')
% xlabel('UTC Time (hr)')
% ylim([0 20])