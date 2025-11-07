%% default setting

params_mouse=struct( ...
    'eta',{2.3e-4},... 
    'beta',{0.15}, ... 
    'epsilon',{0.16}, ...
    'X_c',{17}, ...
    'X_d',{10}, ... % artificial.. should be checked
    'kappa',{0.5});

params_yeast=struct( ...
    'T_eta',{24.31},... 
    'beta',{0}, ... 
    'epsilon',{223.651}, ...
    'X_c',{185.273}, ...
    'X_d',{150}, ... % artificial.. should be checked
    'kappa',{0});

%% yeast
power=1;
params_yeast.eta=params_yeast.T_eta^(-(power+1))*params_yeast.X_c;
[X_paths_yeast,t]=SR_simu_power(params_yeast,0.001,5000,100,0.005,power);     % for avoiding /0, X0!=0
[death_time,~,t]=SR_time(X_paths_yeast,params_yeast,t);
% Compertz law
[h,t]=hazard(death_time,1);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
set(gca, 'XScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("power=1: LOG hazard v.s. LOG time")

power=1.420;
params_yeast.eta=params_yeast.T_eta^(-(power+1))*params_yeast.X_c;
[X_paths_yeast,t]=SR_simu_power(params_yeast,0.001,5000,100,0.005,power);     % for avoiding /0, X0!=0
[death_time,~,t]=SR_time(X_paths_yeast,params_yeast,t);
% Compertz law
[h,t]=hazard(death_time,1);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
set(gca, 'XScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("power=1.420: LOG hazard v.s. LOG time")


power=1.425;
params_yeast.eta=params_yeast.T_eta^(-(power+1))*params_yeast.X_c;
[X_paths_yeast,t]=SR_simu_power(params_yeast,0.001,5000,100,0.005,power);     % for avoiding /0, X0!=0
[death_time,~,t]=SR_time(X_paths_yeast,params_yeast,t);
% Compertz law
[h,t]=hazard(death_time,1);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
set(gca, 'XScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("power=1.425: LOG hazard v.s. LOG time")
%%
power=2;
params_yeast.eta=params_yeast.T_eta^(-(power+1))*params_yeast.X_c;
[X_paths_yeast,t]=SR_simu_power(params_yeast,0.001,5000,100,0.005,power);     % for avoiding /0, X0!=0
[death_time,~,t]=SR_time(X_paths_yeast,params_yeast,t);
% Compertz law
[h,t]=hazard(death_time,1);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
set(gca, 'XScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("power=2: LOG hazard v.s. LOG time")

power=3;
params_yeast.eta=params_yeast.T_eta^(-(power+1))*params_yeast.X_c;
[X_paths_yeast,t]=SR_simu_power(params_yeast,0.001,5000,100,0.005,power);     % for avoiding /0, X0!=0
[death_time,~,t]=SR_time(X_paths_yeast,params_yeast,t);
% Compertz law
[h,t]=hazard(death_time,1);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
set(gca, 'XScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("power=3: LOG hazard v.s. LOG time")

power=4;
params_yeast.eta=params_yeast.T_eta^(-(power+1))*params_yeast.X_c;
[X_paths_yeast,t]=SR_simu_power(params_yeast,0.001,5000,50,0.001,power);     % for avoiding /0, X0!=0
[death_time,~,t]=SR_time(X_paths_yeast,params_yeast,t);
% Compertz law
[h,t]=hazard(death_time,1);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
set(gca, 'XScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("power=4: LOG hazard v.s. LOG time")


%%
power=0.5;
params_yeast.eta=params_yeast.T_eta^(-(power+1))*params_yeast.X_c;
[X_paths_yeast,t]=SR_simu_power(params_yeast,0.001,5000,300,0.01,power);     % for avoiding /0, X0!=0
[death_time,~,t]=SR_time(X_paths_yeast,params_yeast,t);
%% 
% Gompertz law
[h,t]=hazard(death_time,5);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("power=0.5: LOG hazard v.s. time")


%% mice
[X_paths_mice,t]=SR_simu_power(params_mouse,0.001,5000,1500,0.1,1);     % for avoiding /0, X0!=0
[death_time,~,t]=SR_time(X_paths_mice,params_mouse,t);
[f,x]=ecdf(death_time);
figure
stairs(x, 1-f, 'LineWidth', 0.75);
xlabel('Time');
ylabel('Survival probability');
ylim([0 1.05]);
grid on;
% Gompertz law
[h,t]=hazard(death_time,25);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("mice: LOG hazard v.s. time")

%% functions

function [X_paths,t]=SR_simu_power(param_struct,X0,M,EndTime,dt,power)
% param_struct must contain element "eta","beta","kappa","epsilon"
t=0:dt:EndTime;
N=EndTime/dt;
X_paths = X0.*ones(M, N+1);
rng(666)
for n = 1:N
    t_seq=t(n)*ones(M,1);
    mu = param_struct.eta.*(t_seq.^power) - param_struct.beta .* (X_paths(:, n) ./ (X_paths(:, n) + param_struct.kappa));
    dW = sqrt(dt) .* randn(M,1);
    X_paths(:,n+1) = X_paths(:,n) + mu*dt+ sqrt(2*param_struct.epsilon).*dW;
    X_paths(X_paths(:,n+1)<0,n+1)=0.001;
end

end


function [DeathTime,DiseaseTime,t]=SR_time(X,param_struct,t)
% param_struct must contain element "X_c","X_d"
n_sample=length(X(:,1));
n_step=length(X(1,:));

DeathTime=zeros(n_sample,1);
DiseaseTime=zeros(n_sample,1);
X_temp=X;
% for first-crossing
for m = 1:n_sample
    q=0;    % a switch for disease
    for n = 1:n_step
        if X_temp(m,n)> param_struct.X_d && q==0
            DiseaseTime(m)=t(n);
            q=1;
        end
        if X_temp(m,n)> param_struct.X_c
            DeathTime(m)=t(n);
            break
        end
    end
end
end

function [h,t]=hazard(times, binWidth)
    times = sort(times(:));
    n = length(times);
    edges = min(times):binWidth:max(times);
    deaths_per_bin = histcounts(times, edges);  % 每个时间区间的死亡人数
    survivors = n - [0, cumsum(deaths_per_bin(1:end-1))];
    hazard = deaths_per_bin ./ (survivors * binWidth);
    binCenters = edges(1:end-1) + binWidth/2;
    t=binCenters;
    h=hazard;
end


