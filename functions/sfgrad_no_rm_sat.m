function [output_PRM,satava_s,doy,Year,month,date,name] = sfgrad_no_rm_sat(Basefile,Roverfile,navfile,input_PRM)
%============ Constant ================================
f1 = 1575.42*10^6;          %   f1 = 1575.42 MHz (L1)
f2 = 1227.60*10^6;          %   f2 = 1227.60 MHz (L2)
c  = 299792458;             %   light speed = 299792458 m/s
lambda1 = c/f1;             %   wave length of f1
lambda2 = c/f2;             %   wave length of f2
elev_mask = input_PRM.elevation;
wait_reset_filter = 15;
time_index = (0:86399)./86400*24;

%==============================================================================
if(input_PRM.mode == 1) %single-frequency approach
    disp('Readrinex');
    % Read RINEX
    [obs,nav,~,Year] = readrinex211(Basefile,navfile);
    name.basename = obs.station;
    YYYY  = num2str(obs.date(1));
    [time_c1,time_p2,time_l1,time_l2,truerange_sat_a,elevation_sat,azi_sat]=obstorange(obs,nav,elev_mask,input_PRM.base_pos);
    disp('Readrinexfinish');
    phea = time_l1;
    La = (time_l1 + time_c1)/2;
    cod = time_c1 - truerange_sat_a;
    codd = time_l1 - truerange_sat_a;
    phe1 = time_l1;
    cod2 = time_p2;
    phe2 = time_l2;
    cod1 = time_c1;
    multipath_zero_adj % Start the multipath estimation
    ind_slip_use = ind_slip;
    [obs,nav,doy,Year] = readrinex211(Roverfile,navfile);
    name.rovername = obs.station;
    YYYY  = num2str(obs.date(1));
    month = num2str(obs.date(2),'%.2d');
    date  = num2str(obs.date(3),'%.2d');
    range_error_a = time_c1 - truerange_sat_a;
    [time_c1,time_p2,time_l1,time_l2,truerange_sat_b,elevation_sat,azi_sat]=obstorange(obs,nav,elev_mask,input_PRM.rover_pos);
    pheb = time_l1;
    Lb = (time_l1 + time_c1)/2;
    cod = (cod - (time_c1 - truerange_sat_b )); %code_SD
    phe_SD = ((phea - pheb)-(truerange_sat_a - truerange_sat_b))'; %carrier
    cod = -1.*(cod' - phe_SD)./lambda1; %code minus carrier
    L_SD = ((La - Lb)-(truerange_sat_a - truerange_sat_b))'; %m iono free combination
    phe1 = time_l1;
    cod2 = time_p2;
    phe2 = time_l2;
    cod1 = time_c1;
    range_error_b = time_c1 - truerange_sat_b;
    multipath_zero_adj % Start the multipath estimation
    ind_slip_use(ind_slip == 1) = ones(length(find(ind_slip == 1)),1);
    %Start the kalman algorithm
    y_SD = [phe_SD ;L_SD];
    X_SD = nan(66,86400);
    flac = 0; %check num sat from last epoch and now are same
    d_t = 1;%receiver clock difference
    alpha = 1/30;
    t_ion = ones(32,1);
    f_clk = [d_t d_t;0 d_t];
    sig_b = 15;
    sig_ion = 10^(-3)*1; %m/km
    sig_bd = 1;
    sig_amb = 10^(-4)*1;
    beta = 1/30;
    % sig_b = 1;
    % sig_ion = 10^(-3); %m/km
    % sig_bd = 0.0424;
    % sig_amb = 10^(-6)*6;
    % beta = 1/30;
    P_minus = eye(2);
    satb = [];
    ratio = nan(86400,1);
    ratio2 = nan(86400,1);
    integeram_DD_int = nan(32,86400);
    integeram_DD_float = nan(32,86400);
    satava_s = nan(86400,1); %sat in view
    sat_remove = nan(86400,1);
    ion_fix = nan(32,86400);
    satellite_pair = nan(86400,32);
    wait = 0;
    %start KALMAN algorithm
    % delete_sat = find(max(L_SD') - min(L_SD') > 10);
    % L_SD(delete_sat,:) = nan(length(delete_sat),86400);
    for i = 1:length(L_SD)
        L_buffer_notnanindex    = find(~isnan(L_SD(:,i)));
        L_buffer_nanindex       = find(isnan(L_SD(:,i)));
        satava                  = length(L_buffer_notnanindex);
        L_buffer                = L_SD(L_buffer_notnanindex);
        phe_buffer              = phe_SD(L_buffer_notnanindex);
        ll = length(L_buffer_notnanindex);
        if (i == 1)
            satava_past     = satava;
            P_minus         = zeros(2);
            x_minus         = zeros(2,1);
        end
        %disp(i); %number of previous and present sats not same%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H               = [ones(satava,1) zeros(satava,1) -1.*eye(satava,satava) lambda1.*eye(satava,satava);
            ones(satava,1) zeros(satava,1) zeros(satava,satava) (lambda1/2).*eye(satava,satava);];
        if(i>2)
            check_sat_consistence = ~isempty(L_buffer_notnanindex)&&isequal(L_buffer_notnanindex,satb);
            if((ratio(i-1)<2)&&(ratio(i-2)<2)&&(i~=1)&&(wait>wait_reset_filter)&&check_sat_consistence)   %%filter reset
                disp([num2str(i) ' filter reset'])
                ll = length(L_buffer_notnanindex);
                P_minus_buf             = zeros(2+2*satava); %prepare with zero
                x_minus_buf         = zeros(satava*2+2,1); %prepare with zero
                indice = find(L_buffer_notnanindex); %indexs of new sat in new vector
                P_minus_buf(2+indice,2+indice) = diag((sig_ion^2).*ones(length(indice),1));
                P_minus_buf(2+ll+indice,2+ll+indice) = diag(10.*ones(length(indice),1));
                P_minus = P_minus_buf;
                x_minus_buf(2+ll+indice) = round(cod(L_buffer_notnanindex(indice),i));
                %             x_minus_buf(2+ll+indice) = zeros(length(L_buffer_notnanindex(~ismember(L_buffer_notnanindex,satb))),1);
                x_minus             = x_minus_buf;
                wait = 0;
            end
        end
        if(~isempty(L_buffer_notnanindex))
            if ~isequal(L_buffer_notnanindex,satb) %% old = satb , new = L_buffer_notnanindex
                if(length(L_buffer_notnanindex)>length(satb))   %%sat add
                    ll = length(L_buffer_notnanindex);
                    P_minus_buf             = zeros(2+2*satava); %prepare with zero
                    indice = [1; 2; 2+find(ismember(L_buffer_notnanindex,satb)); 2+ll+find(ismember(L_buffer_notnanindex,satb))]; %indexs of old sat in new vector
                    P_minus_buf(indice,indice) = P_minus;
                    x_minus_buf         = zeros(satava*2+2,1); %prepare with zero
                    x_minus_buf(indice) = x_minus;
                    indice = find(~ismember(L_buffer_notnanindex,satb)); %indexs of new sat in new vector
                    P_minus_buf(2+indice,2+indice) = diag((sig_ion^2).*ones(length(indice),1));
                    P_minus_buf(2+ll+indice,2+ll+indice) = diag(10.*ones(length(indice),1));
                    P_minus = P_minus_buf;
                    x_minus_buf(2+ll+indice) = round(cod(L_buffer_notnanindex(~ismember(L_buffer_notnanindex,satb)),i));
                    %             x_minus_buf(2+ll+indice) = zeros(length(L_buffer_notnanindex(~ismember(L_buffer_notnanindex,satb))),1);
                    x_minus             = x_minus_buf;
                    disp([num2str(i) ' add : ' num2str(L_buffer_notnanindex(indice)')])
                elseif(length(L_buffer_notnanindex)<length(satb))%%sat lost
                    ll = length(satb);
                    indice = [1; 2; 2+find(ismember(satb,L_buffer_notnanindex)); 2+ll+find(ismember(satb,L_buffer_notnanindex))];
                    P_minus         = P_minus(indice,indice);
                    x_minus         = x_minus(indice);
                    disp([num2str(i) ' lost : ' num2str(satb(~ismember(satb,L_buffer_notnanindex))')])
                else %%sat lost and add
                    ll = length(satb);
                    indice = [1; 2; 2+find(ismember(satb,L_buffer_notnanindex)); 2+ll+find(ismember(satb,L_buffer_notnanindex))];
                    P_minus         = P_minus(indice,indice);
                    x_minus         = x_minus(indice);
                    ll = length(L_buffer_notnanindex);
                    P_minus_buf             = zeros(2+2*satava);
                    indice = [1; 2; 2+find(ismember(L_buffer_notnanindex,satb)); 2+ll+find(ismember(L_buffer_notnanindex,satb))];
                    P_minus_buf(indice,indice) = P_minus;
                    x_minus_buf         = zeros(satava*2+2,1);
                    x_minus_buf(indice) = x_minus;
                    x_minus             = x_minus_buf;
                    indice = find(~ismember(L_buffer_notnanindex,satb));
                    P_minus_buf(2+indice,2+indice) = diag((sig_ion^2).*ones(length(indice),1));
                    P_minus_buf(2+ll+indice,2+ll+indice) = diag(10.*ones(length(indice),1));
                    P_minus = P_minus_buf;
                end
                f_ion           = eye(satava)*diag(exp(-alpha.*t_ion(L_buffer_notnanindex)));
                f_amb           = eye(satava);
                F_k             = blkdiag(f_clk,f_ion,f_amb);
                R_k             = nancov(y_SD([L_buffer_notnanindex; 32+L_buffer_notnanindex],:)');
                Qk_clk          = [((sig_b^2)+(sig_bd^2)/3) ((sig_bd^2)/2); (sig_bd^2)/2 sig_bd^2];
                Qk_ion          = diag(((ones(satava,1).*sig_ion).^2)./(2*beta).*(1-exp((-2*beta).*t_ion(L_buffer_notnanindex))));
                Qk_amb          = diag((sig_amb.^2).*t_ion(L_buffer_notnanindex));
                Q_k             = blkdiag(Qk_clk,Qk_ion,Qk_amb);
                flac = 1;
            end
            if flac == 2
                f_ion           = eye(satava)*diag(exp(-alpha.*t_ion(L_buffer_notnanindex)));
                F_k             = blkdiag(f_clk,f_ion,f_amb);
                Qk_ion          = diag(((ones(satava,1).*sig_ion).^2)./(2*beta).*(1-exp((-2*beta).*t_ion(L_buffer_notnanindex))));
                Qk_amb          = diag((sig_amb.^2).*t_ion(L_buffer_notnanindex));
                Q_k             = blkdiag(Qk_clk,Qk_ion,Qk_amb);
            end
            %check cycle slip %%%%%%%%%%%%%%%%%%%%%%%%%
            %     if(max(abs(x_minus(2+satava+[1:satava])-cod(L_buffer_notnanindex,i)))>10)
            %         disp('Cycle slip detect')
            %         slip = find(abs(x_minus(2+satava+[1:satava])-cod(L_buffer_notnanindex,i))>10);
            %         x_minus(2+satava+slip) = cod(L_buffer_notnanindex(slip),i);
            %         P_minus         = blkdiag(P_minus(1:2,1:2),(sig_ion^2).*eye(satava),10.*eye(satava));
            %     end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %check cycle slip %%%%%%%%%%%%%%%%%%%%%%%%%
            if(any(~isnan(ind_slip_use(i,L_buffer_notnanindex))))
                disp(['Cycle slip detect at: ' num2str(i)])
                slip = find(ind_slip_use(i,L_buffer_notnanindex) == 1);
                x_minus(2+satava+slip) = round(cod(L_buffer_notnanindex(slip),i));
                %         x_minus(2+satava+slip) = zeros(length(L_buffer_notnanindex(slip)),1);
                P_minus_buf(2+slip,2+slip) = diag((sig_ion^2).*ones(length(slip),1));
                P_minus_buf(2+satava+slip,2+satava+slip) = diag(10.*ones(length(slip),1));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Kk                                                              = P_minus*H'/(H*P_minus*H'+R_k);
            X_SD([1; 2; 2+L_buffer_notnanindex; 34+L_buffer_notnanindex],i) = x_minus + Kk*(y_SD([L_buffer_notnanindex; 32+L_buffer_notnanindex],i)-(H*x_minus));
            P_plus                                                          = P_minus-(Kk*H*P_minus);
            x_minus                                                         = F_k*X_SD([1; 2; 2+L_buffer_notnanindex; 34+L_buffer_notnanindex],i);
            P_minus                                                         = (F_k*P_plus*F_k')+Q_k;
            
            t_ion(L_buffer_nanindex)            = t_ion(L_buffer_nanindex)+1; %plus delta t only sat not available
            t_ion(L_buffer_notnanindex)         = ones(satava,1);
            satb                                = L_buffer_notnanindex;
            flac                                = flac+1;
            satava_past                         = satava;
            
            %LAMBDA
            flac1 = 1; %check ratio test more than 2 or not
            integeram_buff = X_SD(34+L_buffer_notnanindex,i);
            D_met = [-1.*ones(satava-1,1) diag(ones(1,satava-1))];
            P_met = P_plus(3+satava:end,3+satava:end);
            EL = nan(1,32);
            EL(L_buffer_notnanindex) = elevation_sat(i,L_buffer_notnanindex);
            %     flac2 = 1;
            flac2 = 0;
            D_met = [-1.*ones(satava-1,1) diag(ones(1,satava-1))];
            integeram_DD = D_met*integeram_buff; %integer am SD to DD
            integeram_DD_float(1:satava-1,i) = integeram_DD;
            P_DD = D_met*P_met*D_met'; %cov int SD to DD
            %LAMBDA METHOD
            %         try
            [afixed,sqnorm] = LAMBDA(integeram_DD,P_DD,1,'ncands',2);
            ratio(i) = ((integeram_DD-afixed(:,2))'/P_DD*(integeram_DD-afixed(:,2)))/((integeram_DD-afixed(:,1))'/P_DD*(integeram_DD-afixed(:,1)));
            ratio2(i) = 1/(sqnorm(1)/sqnorm(2));
            if(ratio(i)<2)
                %                 disp(ratio(i))
                wait = wait+1;
            end
            if(ratio(i)>2)
                wait = 0;
            end
            %         catch
            %             afixed = round(integeram_DD);
            %             ratio(i) = 0;
            %             ratio2(i) = 0;
            %         end
            %LAMBDA METHOD
            if(ratio2(i)>=2) % in case of pass ratio test
                integeram_DD_int(1:satava-1,i) = afixed(:,1); %fixed solution
                flac1 = 0;
            end
            
            if(satava<4&&satava>0) %sat lower than 4 (float solution)
                if(flac2)
                    satava = satava - 1;
                    P_DD = D_met*P_met*D_met';
                    integeram_DD = D_met*integeram_buff;
                    flac2 = 1;
                end
                integeram_DD_float(1:satava-1,i) = integeram_DD;
                integeram_DD_int(1:satava-1,i) = integeram_DD; %float solution
                ion_SD = X_SD(3:34,i);
                ion_SD = ion_SD(L_buffer_notnanindex(2:end));
                %         ion_fix(L_buffer_notnanindex(2:end),i) = ion_SD;
            end
            satava_s(i) = satava_past;
            sat_remove(i) = satava_past - satava;
            if(satava>1&&ratio2(i)>=2) %fix solution
                ion_SD = X_SD(3:34,i);
                ion_SD = ion_SD(L_buffer_notnanindex(2:end));
                P_IN = P_plus(4+sat_remove(i):satava_past+2,satava_past+sat_remove(i)+4:2*satava_past+2);
                P_N = P_plus(satava_past+sat_remove(i)+4:2*satava_past+2,satava_past+sat_remove(i)+4:2*satava_past+2);
                ion_fix(L_buffer_notnanindex(2:end),i) = ion_SD - P_IN/P_N*(integeram_DD-integeram_DD_int(1:satava-1,i));
                L_buffer_notnanindex_ori    = find(~isnan(L_SD(:,i)));
                X_SD(2+L_buffer_notnanindex(2:end),i) = ion_SD - P_IN/P_N*(integeram_DD-integeram_DD_int(1:satava-1,i));
                x_minus(2+find(ismember(L_buffer_notnanindex_ori,L_buffer_notnanindex(2:end)))) = ion_SD - P_IN/P_N*(integeram_DD-integeram_DD_int(1:satava-1,i));
            end
            satellite_pair(i,L_buffer_notnanindex(2:end)) = ones(satava_s(i)-1,1)*L_buffer_notnanindex(1);
        end
    end
    rec_clk([1 2],:) = X_SD([1 2],:);
    % ion_SD = X_SD(3:34,:)./sqrt(nansum((posa-posb).^2))*1e6; %mm/km
    ion_SD = X_SD(3:34,:); %m
    ion_fix_mm_km = ion_fix./sqrt(nansum((input_PRM.base_pos-input_PRM.rover_pos).^2))*1e6; %mm/km
    integeram = X_SD(35:end,:);
    test_mapping_function
    max_gradient = max(abs(ion_plo(:)));
    output_PRM.elevation = elevation_sat;
    output_PRM.MP1_zero = MP1_zero;
    output_PRM.rec_clk = rec_clk';
    output_PRM.ion_SD = ion_SD';
    output_PRM.ion_fix_mm_km = ion_fix_mm_km';
    output_PRM.integeram = integeram';
    output_PRM.ratio = ratio;
end
if(input_PRM.mode == 2) %dual-frequency approach
    %Base station TEC
    r_o_name = Basefile; % observation file's name SISK0010.16
    r_n_name = navfile; % navigation file's name (if blank, program will downloaded from IGS)
    % r_n_name = ''; % navigation file's name
    
    % Setting#1
    % =========== Program's path ==========================
%     p_path = [pwd '\'];             % Program path
%     S_path = [p_path 'Results\'];   % Results path
%     DCB_path   = [p_path 'DCB\'];   % DCB path
    path(path,[input_PRM.p_path]);
    
    % Read RINEX
    [obs,nav,doy,Year] = readrinex211(r_o_name,r_n_name);
    name.basename = obs.station;
    year  = num2str(obs.date(1));
    month = num2str(obs.date(2),'%.2d');
    date  = num2str(obs.date(3),'%.2d');
    % download satellite bias (DCB)
    [satb.P1C1,satb.P1P2] = dlsat(obs,input_PRM.DCB_path);
    
    cd(input_PRM.p_path)
    %% 2. Calculate Total Electron Content(TEC)
    [~,output_PRM,~] = TECcalculation(obs,nav,satb,input_PRM.S_path);
    %Rover station TEC
    baseTEC = output_PRM.TEC_slant;
    r_o_name = Roverfile; % observation file's name SISK0010.16
    path(path,[input_PRM.p_path]);
    [obs,nav,doy,Year] = readrinex211(r_o_name,r_n_name);
    name.rovername = obs.station;
    year  = num2str(obs.date(1));
    month = num2str(obs.date(2),'%.2d');
    date  = num2str(obs.date(3),'%.2d');
    [satb.P1C1,satb.P1P2] = dlsat(obs,input_PRM.DCB_path);
    cd(input_PRM.p_path)
    [~,output_PRM,~] = TECcalculation(obs,nav,satb,input_PRM.S_path);
    roverTEC = output_PRM.TEC_slant;
    baseline = sqrt(nansum((input_PRM.base_pos-input_PRM.rover_pos).^2));
    output_PRM.ion_fix_mm_km = (baseTEC - roverTEC)/baseline*0.16*10^6;
    output_PRM.ion_fix_mm_km = output_PRM.ion_fix_mm_km - nanmedian(output_PRM.ion_fix_mm_km,1);
    satava_s = nan(1,1);
end
end