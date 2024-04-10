% clear
% clc
% load MP_test.mat
MP1 = cod1 - ((f1^2+f2^2)/(f1^2-f2^2)).*phe1 + ((2*(f2^2))/(f1^2-f2^2)).*phe2;
MP1 = MP1 - ones(86400,1) * nanmean(MP1);
MP1_diff = nan(86400,32);
MP1_diff(2:86400,:) = MP1(2:86400,:) - MP1(1:86399,:);
threshold = 2;
% length(mod(find(~isnan(ind_slip(:,slip_check))),86400))

index_slip = find((max(MP1) - min(MP1))>threshold);
index_noslip = find((max(MP1) - min(MP1))<threshold);
MP1_zero = nan(86400,32);
MP1_zero(:,index_noslip) = MP1(:,index_noslip) - ones(86400,1)*nanmean(MP1(:,index_noslip));
ind_slip = nan(86400,32);
for k = 1:length(index_slip)
    %disp(index_slip(k))
    %count cycle slip
    mea = 0;
    sat = index_slip(k);
    ii = 1;
    start = 1; %bound start with no cycle slip
    mea_bound = [];
    count = 0;
    for i = 1:86400
        if(~isnan(MP1(i,sat)))
            mea = ((ii-1)/ii)*mea + (1/ii)*MP1(i,sat);
            if(abs(mea - MP1(i,sat))>0.8) %<-bound value
                count = count + 1;
                start = [start i];
                mea_bound = [mea_bound mea];
                mea = 0;
                ii = 1;
                mea = ((ii-1)/ii)*mea + (1/ii)*MP1(i,sat);
            end
            ii = ii+1;
        end
    end
    start = [start 86401];
    mea_bound = [mea_bound mea]; % last mean of no cycle slip bound
    %end count cycle slip
    %start repair cycle slip
    for i = 1:length(mea_bound)
        MP1_zero(start(i):start(i+1)-1,sat) = MP1(start(i):start(i+1)-1,sat) - mea_bound(i);
    end
    %end repair cycle slip
    start(start == 1 | start == 86401) = [];
    ind_slip(start,index_slip(k)) = ones(length(start),1);
end
mod(find(~isnan(ind_slip)),86400)
threshold = 0.05;
MP1 = cod1 - ((f1^2+f2^2)/(f1^2-f2^2)).*phe1 + ((2*(f2^2))/(f1^2-f2^2)).*phe2;
MP1_diff = nan(86400,32);
MP1_diff(2:86400,:) = MP1(2:86400,:) - MP1(1:86399,:);
index_slip = find((max(MP1_diff) - min(MP1_diff))>(threshold*2));
index_noslip = find((max(MP1_diff) - min(MP1_diff))<(threshold*2));
for k = 1:length(index_slip)
    %disp(index_slip(k))
    %count cycle slip
    slip_ind = find(abs(MP1_diff(:,index_slip(k)))>threshold);
    ind_slip(slip_ind,index_slip(k)) = ones(length(slip_ind),1);
end

window_for_std = 100*2; %left plus right
sample_test = nan(window_for_std+1,32);
std_MP1 = nan(86400,32);
for ind_test = window_for_std/2+1:86399-window_for_std/2+1
    sample_test = MP1_zero(ind_test-window_for_std/2:ind_test+window_for_std/2,:);
    std_all_sat = nanstd(sample_test);
    std_MP1(ind_test,:) = std_all_sat;
end
% plot(std_MP1,'DisplayName','std_MP1')

LMW = (1/(f1-f2)*(f1*phe1-f2*phe2))-(1/(f1+f2)*(f1*cod1+f2*cod2));
LMW_diff = nan(86400,32);
LMW_diff(2:86400,:) = LMW(2:86400,:) - LMW(1:86399,:);
windows_size = 15; %window size of std for cycle-slip detection
threshold = 14; %time of std
LMW_std = nan(86400,32);
std_buff = movstd(abs(LMW_diff(1:86400,:)),[windows_size-1 0],'omitnan');
LMW_std(51:86400,:) = std_buff(50:86399,:);
ind_slip = nan(86400,32);
ind_slip(abs(LMW_diff) >= LMW_std.*threshold) = ones(length(find(abs(LMW_diff) >= LMW_std.*threshold)),1);
count_notnan = movsum(~isnan(LMW),[windows_size-1 0]);
ind_slip(count_notnan<10) = nan(length(find(count_notnan<10)),1);
mod(find(~isnan(ind_slip)),86400)