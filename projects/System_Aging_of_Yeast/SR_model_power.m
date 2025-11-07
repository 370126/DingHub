%% default setting
params_ml_0=struct( ...
    'T_eta',{24.31},... 
    'beta',{0}, ... 
    'epsilon',{223.651}, ...
    'X_c',{185.273}, ...
    'X_d',{100}, ... % artificial.. should be checked
    'kappa',{0});

power=1;
params_ml_0.eta=params_ml_0.T_eta^(-(power+1))*params_ml_0.X_c;

[X_paths_yeast,t]=SR_simu_power(params_ml_0,0.001,5000,100,0.005,power);     % for avoiding /0, X0!=0
figure
plot(t, X_paths_yeast)
hold on
plot(t,params_ml_0.X_c.*ones(length(t),1),'LineWidth',1.5,'Color','r')
xlabel('T/day')
ylabel('X')
[death_time,~,t]=SR_time(X_paths_yeast,params_ml_0,t);
figure 
histogram(death_time,15)

[f,x]=ecdf(death_time);
figure
stairs(x, 1-f, 'LineWidth', 0.75);
xlabel('Time');
ylabel('Survival probability');
ylim([0 1.05]);
grid on;
% clear("disease_time","death_time")

%% Compertz law?
[h,t]=hazard(death_time,1);
figure
plot(t,h)
set(gca, 'YScale', 'log');  % y 轴对数
xlabel("time")
ylabel("hazard")
title("LOG hazard v.s. time")



%% parameter experiment
params_ml_0=struct( ...
    'T_eta',{24.31},... 
    'beta',{0}, ... 
    'epsilon',{223.651}, ...
    'X_c',{185.273}, ...
    'X_d',{150}, ... % artificial.. should be checked
    'kappa',{0});
X0=0.001;
M=10000;

power_list=[0.9,1.1,1.2,1.3,1.425];
% power_list=1.001;
power_result_list=cell(length(power_list),1);
for p=1:length(power_list)
power=power_list(p);
power_result_list{p}=struct();
params_ml_0.eta=params_ml_0.T_eta^(-(power+1))*params_ml_0.X_c;
fprintf("For power %.2f, eta = %.3f \n",power,params_ml_0.eta)
T_simu=round(300*(1/(power+0.5)),-1);
t_step=0.1;
fprintf("power %.2f start! \n",power)
% experiment setting
experi_para=["eta","epsilon","X_c"];
experi_span=[-1.5:0.3:1.5];  % instead of -0.75:0.15:0.75, for smoother curves

% document the result
% row:parameters; column: effect of change
steepness=zeros(length(experi_para),length(experi_span));
lifespan=steepness;lifespan_v=steepness;
sickspan=steepness;sickspan_v=steepness;
skewness=steepness;

for i=1:length(experi_para)
    para=experi_para(i);
    para_struct=params_ml_0;
    for j=1:length(experi_span)
        k=experi_span(j);
        % how to sample the parameter space ?
        para_experi=params_ml_0.(para) * exp(k);
        para_struct.(para)=para_experi;
            % T_simu=2*para_struct.beta/para_struct.eta;
        % run simulation
        [X_experi,t_experi]=SR_simu_power(para_struct,X0,M,T_simu,t_step,power);
        [DeathTime,DiseaseTime]=SR_time(X_experi,para_struct,t_experi);
        % documenting...
        % lifespan
        [lifespan_sdv,lifespan_mean]=std(DeathTime);
        % sickspan
        sickspan=DeathTime-DiseaseTime;
        [sickspan_sdv,sickspan_mean]=std(sickspan);
        % steepness
% "Steepness was defined by removing the 10% shortest lifespans and
%  computing the mean lifespan divided by the standard deviation of
%  lifespans"
        num_keep=0.9*length(DeathTime);
        temp=sort(DeathTime,"descend");
        lifespan_cutoff=temp(1:num_keep);
        steepness_temp=mean(lifespan_cutoff)/std(lifespan_cutoff);
        % NEW: skewness
        skewness_temp=skew(DeathTime);
        % save reult
        lifespan(i,j)=lifespan_mean;
        steepness(i,j)=steepness_temp;
        skewness(i,j)=skewness_temp;
        index_ref=find(experi_span==0);
        lifespan_rela=lifespan./lifespan(:,index_ref);
        steepness_rela=steepness./steepness(:,index_ref);
        skewness_rela=skewness./skewness(:,index_ref);
        fprintf('%s of %.2f done \n',para,k)
    end
end
power_result_list{p}.lifespan=lifespan_rela;
power_result_list{p}.steepness=steepness_rela;
power_result_list{p}.skewness=skewness_rela;
fprintf("power %.3f done! \n",power)

end

clear("DiseaseTime","DeathTime","X_experi","t_experi","para_experi")
save("myWorkspace.mat")
%% plot
load('myWorkspace.mat')

markers = {'+','^','d','v','>','<','p','h'};
% steepness v.s. lifespan
for i=1:length(power_result_list)
    result=power_result_list{i};
    figure
    for j=1:length(experi_para)
        x=result.lifespan(j,:);
        y=result.steepness(j,:);
        marker_idx = mod(j-1, length(markers)) + 1;
        plot(x, y, 'LineWidth', 1, ...
            'DisplayName', experi_para{j});
        hold on
        scatter(x, y, 36, ...
            'Marker', markers{marker_idx}, ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', 'k', ...
            'HandleVisibility', 'off');
    end
    hold off
    % axis([0.8,1.5,0.8,1.2])
    legend show
    xlabel('relative lifespan')
    ylabel('relative steepness')
    title(sprintf("power p=%.3f",power_list(i)))
end
%% plot
% skewness v.s. steepness
figure
for i=1:length(power_result_list)
    result=power_result_list{i};
    for j=1:length(experi_para)
        x=result.steepness(j,:);
        y=result.skewness(j,:);
        marker_idx = mod(j-1, length(markers)) + 1;
        plot(x, y, 'LineWidth', 1, ...
            'DisplayName', experi_para{j});
        hold on
        scatter(x, y, 36, ...
            'Marker', markers{marker_idx}, ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', 'k', ...
            'HandleVisibility', 'off');
    end
    hold on
    axis([0,3,-5,10])
    legend show
    xlabel('relative steepness')
    ylabel('relative skewness')
    title(sprintf("power p=%.3f",power_list(i)))
end

%% fitting
% eta_x=lifespan_rela(1,:);
% eta_y=steepness_rela(1,:);
% eta_z=skewness_rela(1,:);
% epsilon_x=lifespan_rela(2,:);
% epsilon_y=steepness_rela(2,:);
% epsilon_z=skewness_rela(2,:);
% xc_x=lifespan_rela(3,:);
% xc_y=steepness_rela(3,:);
% xc_z=skewness_rela(3,:);
% 
% all_x=lifespan_rela(:);
% all_y=steepness_rela(:);
% all_z=skewness_rela(:);
% 


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



function t=name_maxtrix(data,rownames,colnames)
row=cellstr(rownames);
col=cellstr(colnames);
t = array2table(data, 'VariableNames', col);
t.Properties.RowNames = row;
end


function skewness=skew(vec)
  n=length(vec);
  [s,mean_vec]=std(vec);
  vec_3=((vec-mean_vec)./s).^3;
  skewness=(n/((n-1)*(n-2)))*sum(vec_3);
end

function fieldName = num2name(k)
    str = sprintf('%.2f', k); 
    str = strrep(str, '.', '_');
    str = strrep(str, '-', 'm');
    fieldName = ['power_', str];
end
