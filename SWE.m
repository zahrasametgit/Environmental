clc;
clear;

%%%% Variable: SWE %%%%

%magnitude frequency plot - CanESM2 - Historical 1991-2020
%SWE - Montreal, QC (45.45671,-73.85838)

%% Reading ERA5 data
%%% Historical 1991-2020
%input_e=ncinfo('g3_stlawrenceX_004deg_360x360_ERA5_SD_1991_2020_daymean.nc');
SD_H_era=ncread('g3_stlawrenceX_004deg_360x360_ERA5_SD_1991_2020_daymean.nc','SD'); %cm
lat_era=ncread('g3_stlawrenceX_004deg_360x360_ERA5_SD_1991_2020_daymean.nc','lat');
lon_era=ncread('g3_stlawrenceX_004deg_360x360_ERA5_SD_1991_2020_daymean.nc','lon');

%leap year
leap=false;
if mod(size(SD_H_era,3),365)
    leap=true;
end
yrs=1991:2020;
days_yr=365*ones(size(yrs));
if leap
    days_yr=days_yr+(mod(yrs,4)==0);
end
%days matrix
day_h=cumsum([1 days_yr]);

SD = SD_H_era; % filling SD
for yr=1:numel(yrs)
    yr_SD=SD(:,:,day_h(yr):(day_h(yr+1)-1)); % get data for the current year
    for day=(days_yr(yr)-334):(days_yr(yr)-151) % Spring (MAR-APR-MAY)
        SD_yr_h(:,:,yr)=max(yr_SD,[],3); % maximum daily SD
    end
end

%pwVw=psVs -> Dw=Ds*ps/pw
SWE_H_era = SD_yr_h*300/1000; % converting Snow depth to snow water equivalent (cm)

% put it in an array for the extreme value analysis
SWE = squeeze(SWE_H_era (:,:,:));


%% Extreme value analysis and plot - GEV distribution - method of L-moments
% ERA5 1991-2020

data_ex1=sort(SWE); % sort small to large

n=numel(data_ex1);

% Define GEV distribution function

% Calculate L-moments
b0=mean(data_ex1);
T=0;
for i = 2:n
    T=T+((i-1)*data_ex1(i));
end
b1=(1/(n*(n-1)))*T;
T2=0;
for i = 3:n
    T2=T2+((i-1)*(i-2)*data_ex1(i));
end
b2=(1/(n*(n-1)*(n-2)))*T2;

l1=b0;
l2=(2*b1)-b0;
l3=(6*b2)-(6*b1)+b0;
t3=l3/l2;

% Estimate the parameters
c=(2/(3+t3))-(log(2)/log(3));
k=(7.859*c)+(2.9554*c^2);
a=(l2*k)/((1-2^-k)*gamma(1+k));
xi=l1-(a*(1-gamma(1+k))/k);

ex_xi1=xi;
ex_alpha1=a;
ex_kappa1=k;

% ERA5 GEV parameters
xi_era = ex_xi1; % Location parameter for ERA5
alpha_era = ex_alpha1; % Scale parameter for ERA5
kappa_era = ex_kappa1; % Shape parameter for ERA5

ex_GEV1=@(F) ex_xi1+ex_alpha1/ex_kappa1*(1-(-log(F)).^ex_kappa1); %function for GEV distribution

figure (1);
F_data=(1:numel(data_ex1))/(numel(data_ex1)+1); % plotting position (Weibull with a=0 used here)

%plot(-log(-log(F_data)),data_ex1,'ob') % plot sample data
p01=plot(-log(-log(F_data)),data_ex1,'or','DisplayName', 'ERA5 Sample');
hold on;

F_dist=0.01:0.01:0.99; % range of F for distribution (up to 100-yr return period)
%plot(-log(-log(F_dist)),ex_GEV1(F_dist),'b') % plot fitted distribution
p02=plot(-log(-log(F_dist)),ex_GEV1(F_dist),'r','DisplayName', 'ERA5 Model');

%ylim([0 8]) % range of y-axis (if needed)
xlabel('Gumbel Reduced Variate') % label for x-axis
ylabel('Max Snow Water Equivalent (cm)') % label for y-axis

title('GEV Distribution - Method of L-Moments');

% mark return levels (dotted lines)
text(-log(-log(0.1)),0.5,'T (years)','HorizontalAlignment','center')
xline(-log(-log(0.5)),':k'); text(-log(-log(0.5)),0.5,'2','HorizontalAlignment','center') % 2y
xline(-log(-log(0.8)),':k'); text(-log(-log(0.8)),0.5,'5','HorizontalAlignment','center') % 5y
xline(-log(-log(0.9)),':k'); text(-log(-log(0.9)),0.5,'10','HorizontalAlignment','center') % 10y
xline(-log(-log(0.95)),':k'); text(-log(-log(0.95)),0.5,'20','HorizontalAlignment','center') % 20y
xline(-log(-log(0.98)),':k'); text(-log(-log(0.98)),0.5,'50','HorizontalAlignment','center') % 50y
xline(-log(-log(0.99)),':k'); text(-log(-log(0.99)),0.5,'100','HorizontalAlignment','center') % 100y


%% Reading CanESM2 data - Snow Depth
%%% Historical 1991-2020
input=ncinfo('g3_stlawrenceX_004deg_360x405_PC-CanESM2_histo_r1i1p1_SD_1991-2020_daymean.nc');
SD_H=ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_histo_r1i1p1_SD_1991-2020_daymean.nc','SD'); %cm
lat_can=ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_histo_r1i1p1_SD_1991-2020_daymean.nc','lat');
lon_can=ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_histo_r1i1p1_SD_1991-2020_daymean.nc','lon');

%leap year
leap=false;
if mod(size(SD_H,3),365)
    leap=true;
end
yrs=1991:2020;
days_yr=365*ones(size(yrs));
if leap
    days_yr=days_yr+(mod(yrs,4)==0);
end
%days matrix
day_h=cumsum([1 days_yr]);

SD = SD_H; % filling SD
for yr=1:numel(yrs)
    yr_SD=SD(:,:,day_h(yr):(day_h(yr+1)-1)); % get data for the current year
    for day=(days_yr(yr)-334):(days_yr(yr)-151) % Spring (MAR-APR-MAY)
        SD_yr_h(:,:,yr)=max(yr_SD,[],3); % maximum daily SD
    end
end

%pwVw=psVs -> Dw=Ds*ps/pw
SWE_H = SD_yr_h*300/1000; % converting Snow depth to snow water equivalent (cm)

% put it in an array for the extreme value analysis
SWE = squeeze(SWE_H (:,:,:));



%% Extreme value analysis and plot - GEV distribution - method of L-moments
% CanESM2 1991-2020

data_ex1=sort(SWE); % sort small to large

n=numel(data_ex1);

% Define GEV distribution function

% Calculate L-moments
b0=mean(data_ex1);
T=0;
for i = 2:n
    T=T+((i-1)*data_ex1(i));
end
b1=(1/(n*(n-1)))*T;
T2=0;
for i = 3:n
    T2=T2+((i-1)*(i-2)*data_ex1(i));
end
b2=(1/(n*(n-1)*(n-2)))*T2;

l1=b0;
l2=(2*b1)-b0;
l3=(6*b2)-(6*b1)+b0;
t3=l3/l2;

% Estimate the parameters
c=(2/(3+t3))-(log(2)/log(3));
k=(7.859*c)+(2.9554*c^2);
a=(l2*k)/((1-2^-k)*gamma(1+k));
xi=l1-(a*(1-gamma(1+k))/k);

ex_xi1=xi;
ex_alpha1=a;
ex_kappa1=k;

xi_can = ex_xi1; % Location parameter for CanESM2
alpha_can = ex_alpha1; % Scale parameter for CanESM2
kappa_can = ex_kappa1; % Shape parameter for CanESM2

ex_GEV1=@(F) ex_xi1+ex_alpha1/ex_kappa1*(1-(-log(F)).^ex_kappa1); %function for GEV distribution

figure (1);
F_data=(1:numel(data_ex1))/(numel(data_ex1)+1); % plotting position (Weibull with a=0 used here)

%plot(-log(-log(F_data)),data_ex1,'ob') % plot sample data
p11=plot(-log(-log(F_data)),data_ex1,'ob','DisplayName', 'CanESM2 Historical Sample');
hold on;

F_dist=0.01:0.01:0.99; % range of F for distribution (up to 100-yr return period)
%plot(-log(-log(F_dist)),ex_GEV1(F_dist),'b') % plot fitted distribution
p12=plot(-log(-log(F_dist)),ex_GEV1(F_dist),'b','DisplayName', 'CanESM2 Historical Model');

%ylim([0 8]) % range of y-axis (if needed)
xlabel('Gumbel Reduced Variate') % label for x-axis
ylabel('Max Snow Water Equivalent (cm)') % label for y-axis

title('GEV Distribution - Method of L-Moments');

% mark return levels (dotted lines)
text(-log(-log(0.1)),0.5,'T (years)','HorizontalAlignment','center')
xline(-log(-log(0.5)),':k'); text(-log(-log(0.5)),0.5,'2','HorizontalAlignment','center') % 2y
xline(-log(-log(0.8)),':k'); text(-log(-log(0.8)),0.5,'5','HorizontalAlignment','center') % 5y
xline(-log(-log(0.9)),':k'); text(-log(-log(0.9)),0.5,'10','HorizontalAlignment','center') % 10y
xline(-log(-log(0.95)),':k'); text(-log(-log(0.95)),0.5,'20','HorizontalAlignment','center') % 20y
xline(-log(-log(0.98)),':k'); text(-log(-log(0.98)),0.5,'50','HorizontalAlignment','center') % 50y
xline(-log(-log(0.99)),':k'); text(-log(-log(0.99)),0.5,'100','HorizontalAlignment','center') % 100y


% %%% certain hot spell durations return period
% y-axis value (y0) for which you want to find the return period and Gumbel variate
y0 = 14;
% Calculate the exceedance probability (F)
F_y0 = exp(-(1-(ex_kappa1*(y0-ex_xi1))/(ex_alpha1))^(1/ex_kappa1));
% Calculate the return period (T)
T_y0 = 1 / (1 - F_y0);
% Calculate the probability (%)
P_y0 = 1 / T_y0;
% Calculate the Gumbel variate (u)
u_y0 = -log(-log(F_y0));
% plot with the return period and Gumbel variate
figure(1);
hold on;
plot(u_y0, y0, 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Plot the point
text(u_y0, y0, {[' T = ', num2str(T_y0, '%.2f'), ' years'], [' P = ', num2str(P_y0, '%.2f')]}, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'blue');
yline(y0,':k'); % 2y
text(-2,y0,num2str(y0),'HorizontalAlignment','center') % y0 display

% % Given return period
% T0 = 10; 
% % Calculate the exceedance probability (F)
% F_T0 = 1 - 1/T0;
% % Calculate the corresponding y-value (y0) using the Gumbel distribution function
% y0 = ex_xi1 + (ex_alpha1/ex_kappa1) * (1-(-log(F_T0))^ex_kappa1);
% figure(1);
% hold on;
% plot(-log(-log(F_T0)), y0, 'kx', 'MarkerSize', 10, 'LineWidth', 2); % Plot the point
% text(-log(-log(F_T0)), y0, [' x = ', num2str(y0, '%.2f'), ' cm'], 'VerticalAlignment', 'bottom', ...
%     'HorizontalAlignment', 'right','Color', 'black');
% % T0 on the axis
% xline(-log(-log(F_T0)),':k');
% text(-log(-log(F_T0)),0.5,num2str(T0),'HorizontalAlignment','center') % y0 display

legend([p01, p02,p11,p12],'Location', 'northwest');


%% Performance Error: CanESM-ERA5

% Define return periods
T = [2, 5, 10, 20, 50, 100]; % years

% Calculate return levels for ERA5 and CanESM2
return_levels_era = xi_era + (alpha_era / kappa_era) * (1 - (-log(1 - 1 ./ T)).^kappa_era);
return_levels_can = xi_can + (alpha_can / kappa_can) * (1 - (-log(1 - 1 ./ T)).^kappa_can);

% Calculate performance errors
relative_error = abs(return_levels_can - return_levels_era) ./ return_levels_era * 100; % Relative error (%)
absolute_error = abs(return_levels_can - return_levels_era); % Absolute error

% Plot results
figure (5);
subplot(2,1,1);
plot(T, relative_error, '-o');
xlabel('Return Period (Years)');
ylabel('Relative Error (%)');
title('Relative Error Between CanESM2 and ERA5');
grid on;

subplot(2,1,2);
plot(T, absolute_error, '-o');
xlabel('Return Period (Years)');
ylabel('Absolute Error (cm)');
title('Absolute Error Between CanESM2 and ERA5');
grid on;


%% Choosing Spring to focus on

% Monthly SD plots for 6 years
%To choose the right time span for analysis
years_to_plot = [1992, 1998, 2004, 2010, 2016]; % 6 years for plotting

% Loop over the years you want to plot
for yr_idx = 1:length(years_to_plot)
    yr = years_to_plot(yr_idx);
    
    % start and end day indices for the year
    start_day = day_h(yr - min(yrs)) + 1;
    end_day = day_h(yr - min(yrs) + 1);
    
    %the snow depth data for the given year
    SD_year = SD_H_era(:,:,start_day:end_day);
    
    % Compute the monthly averages
    months_avg = zeros(12,1);
    for month = 1:12
        month_start = (month - 1) * 30 + 1; % Approximate start of month
        month_end = month * 30; % Approximate end of month
        
        % Clip if it's the last month and doesn't have exactly 30 days
        if month == 12
            month_end = size(SD_year,3); % Adjust the end for December
        end
        
        % Average over the days of the month
        months_avg(month) = mean(SD_year(:,:,month_start:month_end), 'all');
    end
    
    % Plotting
    figure (10);
    plot(1:12, months_avg, '-o', 'DisplayName', sprintf('Year %d', yr));
    hold on;
end

% Formatting the plot
xlabel('Month');
ylabel('Average Snow Depth (cm)');
title('Monthly Snow Depth Changes - ERA5');
legend('show');
grid on;
% x-axis ticks to be from 1 to 12 (representing each month)
xticks(1:12); 
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'}); 

% Monthly Temperature plots for 6 years
%To choose the right time span for analysis
years_to_plot = [1992, 1998, 2004, 2010, 2016]; % 6 years for plotting

%Temperature
t_H_era =ncread('g3_stlawrenceX_004deg_360x405_ERA5_TT2M_1991_to_2020_DAYMEAN.nc','TT2M'); %C

% Loop over the years you want to plot
for yr_idx = 1:length(years_to_plot)
    yr = years_to_plot(yr_idx);
    
    % start and end day indices for the year
    start_day = day_h(yr - min(yrs)) + 1;
    end_day = day_h(yr - min(yrs) + 1);
    
    %the temp data for the given year
    SD_year = t_H_era(:,:,start_day:end_day);
    
    % Compute the monthly averages
    months_avg = zeros(12,1);
    for month = 1:12
        month_start = (month - 1) * 30 + 1; % Approximate start of month
        month_end = month * 30; % Approximate end of month
        
        % Clip if it's the last month and doesn't have exactly 30 days
        if month == 12
            month_end = size(SD_year,3); % Adjust the end for December
        end
        
        % Average over the days of the month
        months_avg(month) = mean(SD_year(:,:,month_start:month_end), 'all');
    end
    
    % Plotting
    figure (20);
    plot(1:12, months_avg, '-o', 'DisplayName', sprintf('Year %d', yr));
    hold on;
end

% Formatting the plot
xlabel('Month');
ylabel('Average Temperature (C)');
title('Monthly Temperature Changes - ERA5');
legend('show');
grid on;
% x-axis ticks to be from 1 to 12 (representing each month)
xticks(1:12); 
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'});
yline(0, '-.k', 'LineWidth', 1);



%% Reading CanESM2 data - Snow Depth
%%% RCP4.5 - 2040-2069
input_445=ncinfo('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp45_r1i1p1_SD_2040-2069_daymean.nc');
SD_445=ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp45_r1i1p1_SD_2040-2069_daymean.nc','SD'); %cm

%leap year
leap=false;
if mod(size(SD_445,3),365)
    leap=true;
end
yrs=2040:2069;
days_yr=365*ones(size(yrs));
if leap
    days_yr=days_yr+(mod(yrs,4)==0);
end
%days matrix
day_h=cumsum([1 days_yr]);

%day_h = ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp45_r1i1p1_SD_2040-2069_cdaymean.nc','time');

SD = SD_445; % filling SD
for yr=1:numel(yrs)
    yr_SD = SD(:,:,day_h(yr):(day_h(yr+1)-1)); % get data for the current year
    for day=(days_yr(yr)-334):(days_yr(yr)-151) % Spring (MAR-APR-MAY)
        SD_yr_h(:,:,yr)=max(yr_SD,[],3); % maximum daily SD
    end
end

%pwVw=psVs -> Dw=Ds*ps/pw
SWE_445 = SD_yr_h*300/1000; % converting Snow depth to snow water equivalent (cm)

% put it in an array for the extreme value analysis
SWE = squeeze(SWE_445 (:,:,:));


%% Extreme value analysis and plot - GEV distribution - method of L-moments
% RCP4.5 - 2040-2069

data_ex1=sort(SWE); % sort small to large

n=numel(data_ex1);

% Define GEV distribution function

% Calculate L-moments
b0=mean(data_ex1);
T=0;
for i = 2:n
    T=T+((i-1)*data_ex1(i));
end
b1=(1/(n*(n-1)))*T;
T2=0;
for i = 3:n
    T2=T2+((i-1)*(i-2)*data_ex1(i));
end
b2=(1/(n*(n-1)*(n-2)))*T2;

l1=b0;
l2=(2*b1)-b0;
l3=(6*b2)-(6*b1)+b0;
t3=l3/l2;

% Estimate the parameters
c=(2/(3+t3))-(log(2)/log(3));
k=(7.859*c)+(2.9554*c^2);
a=(l2*k)/((1-2^-k)*gamma(1+k));
xi=l1-(a*(1-gamma(1+k))/k);

ex_xi1=xi;
ex_alpha1=a;
ex_kappa1=k;

ex_GEV1=@(F) ex_xi1+ex_alpha1/ex_kappa1*(1-(-log(F)).^ex_kappa1); %function for GEV distribution

figure (2);
F_data=(1:numel(data_ex1))/(numel(data_ex1)+1); % plotting position (Weibull with a=0 used here)

p21=plot(-log(-log(F_data)),data_ex1,'ob','DisplayName', 'CanESM2 2040-2069 RCP4.5 Sample');
hold on;

F_dist=0.01:0.01:0.99; % range of F for distribution (up to 100-yr return period)
p22=plot(-log(-log(F_dist)),ex_GEV1(F_dist),'b','DisplayName', 'CanESM2 2040-2069 RCP4.5 Model');

xlabel('Gumbel Reduced Variate') % label for x-axis
ylabel('Max Snow Water Equivalent (cm)') % label for y-axis

title('GEV Distribution - Method of L-Moments');

% mark return levels (dotted lines)
text(-log(-log(0.1)),0.5,'T (years)','HorizontalAlignment','center')
xline(-log(-log(0.5)),':k'); text(-log(-log(0.5)),0.5,'2','HorizontalAlignment','center') % 2y
xline(-log(-log(0.8)),':k'); text(-log(-log(0.8)),0.5,'5','HorizontalAlignment','center') % 5y
xline(-log(-log(0.9)),':k'); text(-log(-log(0.9)),0.5,'10','HorizontalAlignment','center') % 10y
xline(-log(-log(0.95)),':k'); text(-log(-log(0.95)),0.5,'20','HorizontalAlignment','center') % 20y
xline(-log(-log(0.98)),':k'); text(-log(-log(0.98)),0.5,'50','HorizontalAlignment','center') % 50y
xline(-log(-log(0.99)),':k'); text(-log(-log(0.99)),0.5,'100','HorizontalAlignment','center') % 100y

% %%% certain hot spell durations return period
% y-axis value (y0) for which you want to find the return period and Gumbel variate
y0 = 14;
% Calculate the exceedance probability (F)
F_y0 = exp(-(1-(ex_kappa1*(y0-ex_xi1))/(ex_alpha1))^(1/ex_kappa1));
% Calculate the return period (T)
T_y0 = 1 / (1 - F_y0);
% Calculate the probability (%)
P_y0 = 1 / T_y0;
% Calculate the Gumbel variate (u)
u_y0 = -log(-log(F_y0));
% plot with the return period and Gumbel variate
figure(2);
hold on;
plot(u_y0, y0, 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Plot the point
text(u_y0, y0, {[' T = ', num2str(T_y0, '%.2f'), ' years'], [' P = ', num2str(P_y0, '%.2f')]}, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'blue');
yline(y0,':k'); % 2y
text(-2,y0,num2str(y0),'HorizontalAlignment','center') % y0 display


%% Reading CanESM2 data - Snow Depth
%%% RCP4.5 - 2070-2099
input_745=ncinfo('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp45_r1i1p1_SD_2070-2099_daymean.nc');
SD_745=ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp45_r1i1p1_SD_2070-2099_daymean.nc','SD'); %cm

%leap year
leap=false;
if mod(size(SD_745,3),365)
    leap=true;
end
yrs=2070:2099;
days_yr=365*ones(size(yrs));
if leap
    days_yr=days_yr+(mod(yrs,4)==0);
end
%days matrix
day_h=cumsum([1 days_yr]);

SD = SD_745; % filling SD
for yr=1:numel(yrs)
    yr_SD=SD(:,:,day_h(yr):(day_h(yr+1)-1)); % get data for the current year
    for day=(days_yr(yr)-334):(days_yr(yr)-151) % Spring (MAR-APR-MAY)
        SD_yr_h(:,:,yr)=max(yr_SD,[],3); % maximum daily SD
    end
end

%pwVw=psVs -> Dw=Ds*ps/pw
SWE_745 = SD_yr_h*300/1000; % converting Snow depth to snow water equivalent (cm)

% put it in an array for the extreme value analysis
SWE = squeeze(SWE_745 (:,:,:));


%% Extreme value analysis and plot - GEV distribution - method of L-moments
% RCP4.5 - 2070-2099

data_ex1=sort(SWE); % sort small to large

n=numel(data_ex1);

% Define GEV distribution function

% Calculate L-moments
b0=mean(data_ex1);
T=0;
for i = 2:n
    T=T+((i-1)*data_ex1(i));
end
b1=(1/(n*(n-1)))*T;
T2=0;
for i = 3:n
    T2=T2+((i-1)*(i-2)*data_ex1(i));
end
b2=(1/(n*(n-1)*(n-2)))*T2;

l1=b0;
l2=(2*b1)-b0;
l3=(6*b2)-(6*b1)+b0;
t3=l3/l2;

% Estimate the parameters
c=(2/(3+t3))-(log(2)/log(3));
k=(7.859*c)+(2.9554*c^2);
a=(l2*k)/((1-2^-k)*gamma(1+k));
xi=l1-(a*(1-gamma(1+k))/k);

ex_xi1=xi;
ex_alpha1=a;
ex_kappa1=k;

ex_GEV1=@(F) ex_xi1+ex_alpha1/ex_kappa1*(1-(-log(F)).^ex_kappa1); %function for GEV distribution

figure (2);
F_data=(1:numel(data_ex1))/(numel(data_ex1)+1); % plotting position (Weibull with a=0 used here)

p31=plot(-log(-log(F_data)),data_ex1,'ok','DisplayName', 'CanESM2 2070-2099 RCP4.5 Sample');
hold on;

F_dist=0.01:0.01:0.99; % range of F for distribution (up to 100-yr return period)
p32=plot(-log(-log(F_dist)),ex_GEV1(F_dist),'k','DisplayName', 'CanESM2 2070-2099 RCP4.5 Model');

% xlabel('Gumbel Reduced Variate') % label for x-axis
% ylabel('Max Snow Water Equivalent (cm)') % label for y-axis
% 
% title('GEV Distribution - Method of L-Moments');

% mark return levels (dotted lines)
text(-log(-log(0.1)),0.5,'T (years)','HorizontalAlignment','center')
xline(-log(-log(0.5)),':k'); text(-log(-log(0.5)),0.5,'2','HorizontalAlignment','center') % 2y
xline(-log(-log(0.8)),':k'); text(-log(-log(0.8)),0.5,'5','HorizontalAlignment','center') % 5y
xline(-log(-log(0.9)),':k'); text(-log(-log(0.9)),0.5,'10','HorizontalAlignment','center') % 10y
xline(-log(-log(0.95)),':k'); text(-log(-log(0.95)),0.5,'20','HorizontalAlignment','center') % 20y
xline(-log(-log(0.98)),':k'); text(-log(-log(0.98)),0.5,'50','HorizontalAlignment','center') % 50y
xline(-log(-log(0.99)),':k'); text(-log(-log(0.99)),0.5,'100','HorizontalAlignment','center') % 100y

% %%% certain hot spell durations return period
% y-axis value (y0) for which you want to find the return period and Gumbel variate
y0 = 14;
% Calculate the exceedance probability (F)
F_y0 = exp(-(1-(ex_kappa1*(y0-ex_xi1))/(ex_alpha1))^(1/ex_kappa1));
% Calculate the return period (T)
T_y0 = 1 / (1 - F_y0);
% Calculate the probability (%)
P_y0 = 1 / T_y0;
% Calculate the Gumbel variate (u)
u_y0 = -log(-log(F_y0));
% plot with the return period and Gumbel variate
figure(2);
hold on;
plot(u_y0, y0, 'kx', 'MarkerSize', 10, 'LineWidth', 2); % Plot the point
text(u_y0, y0, {[' T = ', num2str(T_y0, '%.2f'), ' years'], [' P = ', num2str(P_y0, '%.2f')]}, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'black');
yline(y0,':k'); % 2y
text(-2,y0,num2str(y0),'HorizontalAlignment','center') % y0 display



%% Reading CanESM2 data - Snow Depth
%%% RCP8.5 - 2040-2069
%input_485=ncinfo('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp85_r1i1p1_SD_2040-2069_daymean.nc');
SD_485=ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp85_r1i1p1_SD_2040-2069_daymean.nc','SD'); %cm

%leap year
leap=false;
if mod(size(SD_485,3),365)
    leap=true;
end
yrs=2040:2069;
days_yr=365*ones(size(yrs));
if leap
    days_yr=days_yr+(mod(yrs,4)==0);
end
%days matrix
day_h=cumsum([1 days_yr]);

%day_h = ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp45_r1i1p1_SD_2040-2069_cdaymean.nc','time');

SD = SD_485; % filling SD
for yr=1:numel(yrs)
    yr_SD = SD(:,:,day_h(yr):(day_h(yr+1)-1)); % get data for the current year
    for day=(days_yr(yr)-334):(days_yr(yr)-151) % Spring (MAR-APR-MAY)
        SD_yr_h(:,:,yr)=max(yr_SD,[],3); % maximum daily SD
    end
end

%pwVw=psVs -> Dw=Ds*ps/pw
SWE_485 = SD_yr_h*300/1000; % converting Snow depth to snow water equivalent (cm)

% put it in an array for the extreme value analysis
SWE = squeeze(SWE_485 (:,:,:));


%% Extreme value analysis and plot - GEV distribution - method of L-moments
% RCP8.5 - 2040-2069

data_ex1=sort(SWE); % sort small to large

n=numel(data_ex1);

% Define GEV distribution function

% Calculate L-moments
b0=mean(data_ex1);
T=0;
for i = 2:n
    T=T+((i-1)*data_ex1(i));
end
b1=(1/(n*(n-1)))*T;
T2=0;
for i = 3:n
    T2=T2+((i-1)*(i-2)*data_ex1(i));
end
b2=(1/(n*(n-1)*(n-2)))*T2;

l1=b0;
l2=(2*b1)-b0;
l3=(6*b2)-(6*b1)+b0;
t3=l3/l2;

% Estimate the parameters
c=(2/(3+t3))-(log(2)/log(3));
k=(7.859*c)+(2.9554*c^2);
a=(l2*k)/((1-2^-k)*gamma(1+k));
xi=l1-(a*(1-gamma(1+k))/k);

ex_xi1=xi;
ex_alpha1=a;
ex_kappa1=k;

ex_GEV1=@(F) ex_xi1+ex_alpha1/ex_kappa1*(1-(-log(F)).^ex_kappa1); %function for GEV distribution

figure (2);
F_data=(1:numel(data_ex1))/(numel(data_ex1)+1); % plotting position (Weibull with a=0 used here)

p41=plot(-log(-log(F_data)),data_ex1,'or','DisplayName', 'CanESM2 2040-2069 RCP8.5 Sample');
hold on;

F_dist=0.01:0.01:0.99; % range of F for distribution (up to 100-yr return period)
p42=plot(-log(-log(F_dist)),ex_GEV1(F_dist),'r','DisplayName', 'CanESM2 2040-2069 RCP8.5 Model');

xlabel('Gumbel Reduced Variate') % label for x-axis
ylabel('Max Snow Water Equivalent (cm)') % label for y-axis

title('GEV Distribution - Method of L-Moments');

% mark return levels (dotted lines)
text(-log(-log(0.1)),0.5,'T (years)','HorizontalAlignment','center')
xline(-log(-log(0.5)),':k'); text(-log(-log(0.5)),0.5,'2','HorizontalAlignment','center') % 2y
xline(-log(-log(0.8)),':k'); text(-log(-log(0.8)),0.5,'5','HorizontalAlignment','center') % 5y
xline(-log(-log(0.9)),':k'); text(-log(-log(0.9)),0.5,'10','HorizontalAlignment','center') % 10y
xline(-log(-log(0.95)),':k'); text(-log(-log(0.95)),0.5,'20','HorizontalAlignment','center') % 20y
xline(-log(-log(0.98)),':k'); text(-log(-log(0.98)),0.5,'50','HorizontalAlignment','center') % 50y
xline(-log(-log(0.99)),':k'); text(-log(-log(0.99)),0.5,'100','HorizontalAlignment','center') % 100y

% %%% certain hot spell durations return period
% y-axis value (y0) for which you want to find the return period and Gumbel variate
y0 = 14;
% Calculate the exceedance probability (F)
F_y0 = exp(-(1-(ex_kappa1*(y0-ex_xi1))/(ex_alpha1))^(1/ex_kappa1));
% Calculate the return period (T)
T_y0 = 1 / (1 - F_y0);
% Calculate the probability (%)
P_y0 = 1 / T_y0;
% Calculate the Gumbel variate (u)
u_y0 = -log(-log(F_y0));
% plot with the return period and Gumbel variate
figure(2);
hold on;
plot(u_y0, y0, 'rx', 'MarkerSize', 10, 'LineWidth', 2); % Plot the point
text(u_y0, y0, {[' T = ', num2str(T_y0, '%.2f'), ' years'], [' P = ', num2str(P_y0, '%.2f')]}, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'red');
yline(y0,':k'); % 2y
text(-2,y0,num2str(y0),'HorizontalAlignment','center') % y0 display


%% Reading CanESM2 data - Snow Depth
%%% RCP8.5 - 2070-2099
input_785=ncinfo('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp85_r1i1p1_SD_2070-2099_daymean.nc');
SD_785=ncread('g3_stlawrenceX_004deg_360x405_PC-CanESM2_rcp85_r1i1p1_SD_2070-2099_daymean.nc','SD'); %cm

%leap year
leap=false;
if mod(size(SD_785,3),365)
    leap=true;
end
yrs=2070:2099;
days_yr=365*ones(size(yrs));
if leap
    days_yr=days_yr+(mod(yrs,4)==0);
end
%days matrix
day_h=cumsum([1 days_yr]);

SD = SD_785; % filling SD
for yr=1:numel(yrs)
    yr_SD=SD(:,:,day_h(yr):(day_h(yr+1)-1)); % get data for the current year
    for day=(days_yr(yr)-334):(days_yr(yr)-151) % Spring (MAR-APR-MAY)
        SD_yr_h(:,:,yr)=max(yr_SD,[],3); % maximum daily SD
    end
end

%pwVw=psVs -> Dw=Ds*ps/pw
SWE_785 = SD_yr_h*300/1000; % converting Snow depth to snow water equivalent (cm)

% put it in an array for the extreme value analysis
SWE = squeeze(SWE_785 (:,:,:));


%% Extreme value analysis and plot - GEV distribution - method of L-moments
% RCP8.5 - 2070-2099

data_ex1=sort(SWE); % sort small to large

n=numel(data_ex1);

% Define GEV distribution function

% Calculate L-moments
b0=mean(data_ex1);
T=0;
for i = 2:n
    T=T+((i-1)*data_ex1(i));
end
b1=(1/(n*(n-1)))*T;
T2=0;
for i = 3:n
    T2=T2+((i-1)*(i-2)*data_ex1(i));
end
b2=(1/(n*(n-1)*(n-2)))*T2;

l1=b0;
l2=(2*b1)-b0;
l3=(6*b2)-(6*b1)+b0;
t3=l3/l2;

% Estimate the parameters
c=(2/(3+t3))-(log(2)/log(3));
k=(7.859*c)+(2.9554*c^2);
a=(l2*k)/((1-2^-k)*gamma(1+k));
xi=l1-(a*(1-gamma(1+k))/k);

ex_xi1=xi;
ex_alpha1=a;
ex_kappa1=k;

ex_GEV1=@(F) ex_xi1+ex_alpha1/ex_kappa1*(1-(-log(F)).^ex_kappa1); %function for GEV distribution

figure (2);
F_data=(1:numel(data_ex1))/(numel(data_ex1)+1); % plotting position (Weibull with a=0 used here)

p51=plot(-log(-log(F_data)),data_ex1,'om','DisplayName', 'CanESM2 2070-2099 RCP8.5 Sample');
hold on;

F_dist=0.01:0.01:0.99; % range of F for distribution (up to 100-yr return period)
p52=plot(-log(-log(F_dist)),ex_GEV1(F_dist),'m','DisplayName', 'CanESM2 2070-2099 RCP8.5 Model');

% xlabel('Gumbel Reduced Variate') % label for x-axis
% ylabel('Max Snow Water Equivalent (cm)') % label for y-axis
% 
% title('GEV Distribution - Method of L-Moments');

% mark return levels (dotted lines)
text(-log(-log(0.1)),0.5,'T (years)','HorizontalAlignment','center')
xline(-log(-log(0.5)),':k'); text(-log(-log(0.5)),0.5,'2','HorizontalAlignment','center') % 2y
xline(-log(-log(0.8)),':k'); text(-log(-log(0.8)),0.5,'5','HorizontalAlignment','center') % 5y
xline(-log(-log(0.9)),':k'); text(-log(-log(0.9)),0.5,'10','HorizontalAlignment','center') % 10y
xline(-log(-log(0.95)),':k'); text(-log(-log(0.95)),0.5,'20','HorizontalAlignment','center') % 20y
xline(-log(-log(0.98)),':k'); text(-log(-log(0.98)),0.5,'50','HorizontalAlignment','center') % 50y
xline(-log(-log(0.99)),':k'); text(-log(-log(0.99)),0.5,'100','HorizontalAlignment','center') % 100y

% %%% certain hot spell durations return period
% y-axis value (y0) for which you want to find the return period and Gumbel variate
y0 = 14;
% Calculate the exceedance probability (F)
F_y0 = exp(-(1-(ex_kappa1*(y0-ex_xi1))/(ex_alpha1))^(1/ex_kappa1));
% Calculate the return period (T)
T_y0 = 1 / (1 - F_y0);
% Calculate the probability (%)
P_y0 = 1 / T_y0;
% Calculate the Gumbel variate (u)
u_y0 = -log(-log(F_y0));
% plot with the return period and Gumbel variate
figure(2);
hold on;
plot(u_y0, y0, 'mx', 'MarkerSize', 10, 'LineWidth', 2); % Plot the point
text(u_y0, y0, {[' T = ', num2str(T_y0, '%.2f'), ' years'], [' P = ', num2str(P_y0, '%.2f')]}, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', 'magenta');
yline(y0,':k'); % 2y
text(-2,y0,num2str(y0),'HorizontalAlignment','center') % y0 display

legend([p21, p22, p31, p32 ,p41, p42, p51, p52],'Location', 'northwest');