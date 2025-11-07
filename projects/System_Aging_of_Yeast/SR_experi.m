%% default setting

% first dimension: default; second dimension: the power of the parameter
params_ref=struct( ...
    'eta',{12.5}, ... % 12.5
    'beta',{0}, ... % should be 0.11? beta=0 for much easier analytical solving
    'epsilon',{83}, ...
    'X_c',{176}, ...
    'X_d',{150}, ... % artificial.. should be checked
    'kappa',{0});

params_ml=struct( ...
    'eta',{0.314}, ... 
    'beta',{0.121}, ... 
    'epsilon',{223.651}, ...
    'X_c',{185.273}, ...
    'X_d',{150}, ... % artificial.. should be checked
    'kappa',{4.146*0.121});

params_ml_0=struct( ...
    'eta',{0.314}, ... 
    'beta',{0}, ... 
    'epsilon',{223.651}, ...
    'X_c',{185.273}, ...
    'X_d',{150}, ... % artificial.. should be checked
    'kappa',{4.146*0.121});

params_mo=struct( ...
    'eta',{0.0284}, ...
    'beta',{0.000572}, ... 
    'epsilon',{1.983}, ...
    'X_c',{185.273}, ...
    'X_d',{150}, ... % artificial.. should be checked
    'kappa',{87.402*0.000572});

params_mo_0=struct( ...
    'eta',{0.0284}, ...
    'beta',{0}, ... 
    'epsilon',{1.983}, ...
    'X_c',{185.273}, ...
    'X_d',{150}, ... % artificial.. should be checked
    'kappa',{87.402*0.000572});

params_strange=struct( ...
    'eta',{0.314*exp(-1.8)}, ...  % strange
    'beta',{0}, ... 
    'epsilon',{223.651}, ...
    'X_c',{185.273}, ...
    'X_d',{150}, ... % artificial.. should be checked
    'kappa',{4.146*0.121});



params_cells={params_ref,params_ml,params_ml_0,params_mo,params_mo_0};
names_cells={"reference","maximum likelihood","maximum likelihood(beta=0)",...
    "mode overall","mode overall(beta=0)"};

%%
% % 
% % for i=1:length(params_cells)
% % function [X_paths,t]=SR_simu(param_struct,X0,M,EndTime,dt)
% [X_paths_yeast,t]=SR_simu2(params_strange,0.001,10000,250,0.02);     % for avoiding /0, X0!=0
% figure
% plot(t, X_paths_yeast)
% hold on
% plot(t,params_ml.X_c.*ones(length(t),1),'LineWidth',1.5,'Color','r')
% xlabel('T/day')
% ylabel('X')
% % title(names_cells{i})
% % function [DeathTime,DiseaseTime,t]=SR_time(X,param_struct,t)
% [death_time,~,t]=SR_time(X_paths_yeast,params_strange,t);
% figure 
% histogram(death_time,15)
% % xlim([0,150])
% % title(names_cells{i})
% % fprintf("%s done\n",names_cells{i})
% clear("disease_time","death_time")
% % end
%% parameter experiment
T_simu=300;
X0=0.001;
M=10000;
t_step=0.01;
% experiment setting
experi_para=["eta","epsilon","X_c"];
experi_span=[-1.8:0.1:1.5];  % instead of -0.75:0.15:0.75, for smoother curves
result=struct();

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
        [X_experi,t_experi]=SR_simu2(para_struct,X0,M,T_simu,t_step);
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
        lifespan_v(i,j)=lifespan_sdv;
        sickspan(i,j)=sickspan_mean;
        sickspan_v(i,j)=sickspan_sdv;
        steepness(i,j)=steepness_temp;
        skewness(i,j)=skewness_temp;
        % save result
        % result.(para).(num2name(k))={X_experi,t_experi,DeathTime,DiseaseTime}; 
        % (1):trajectory (2):time (3):death time (4)disease time
        fprintf('%s of %.2f done \n',para,k)
    end
end

clear("DiseaseTime","DeathTime","X_experi","t_experi","para_experi")

%% calculate properties for each trial
% 
% % row:parameters; column: effect of change
% steepness=zeros(length(experi_para),length(experi_span));
% lifespan=steepness;lifespan_v=steepness;
% sickspan=steepness;sickspan_v=steepness;
% % NEW
% skewness=steepness;
% 
% for i=1:length(experi_para)
%     para=experi_para(i);
%     for j=1:length(experi_span)
%         k=experi_span(j);
%             % X_temp=result.result.(para).(num2name(k))(1);   % #traj X time
%         DeathTime_temp=result.(para).(num2name(k)){3};
%         DiseaseTime_temp=result.(para).(num2name(k)){4};
%         sickspan_temp=DeathTime_temp-DiseaseTime_temp;
%         % lifespan
%         [lifespan_sdv,lifespan_mean]=std(DeathTime_temp);
%         % sickspan
%         [sickspan_sdv,sickspan_mean]=std(sickspan_temp);
%         % steepness
% % "Steepness was defined by removing the 10% shortest lifespans and
% %  computing the mean lifespan divided by the standard deviation of
% %  lifespans"
%         num_keep=0.9*length(DeathTime_temp);
%         temp=sort(DeathTime_temp,"descend");
%         lifespan_cutoff=temp(1:num_keep);
%         steepness_temp=mean(lifespan_cutoff)/std(lifespan_cutoff);
%         % NEW: skewness
%         skewness_temp=skew(DeathTime_temp);
%         % save reult
%         lifespan(i,j)=lifespan_mean;
%         lifespan_v(i,j)=lifespan_sdv;
%         sickspan(i,j)=sickspan_mean;
%         sickspan_v(i,j)=sickspan_sdv;
%         steepness(i,j)=steepness_temp;
%         skewness(i,j)=skewness_temp;
%     end
%     fprintf('%s documented \n',para)
% end
% 
% clear("DiseaseTime_temp","sickspan_temp",'DeathTime_temp', ...
%     "lifespan_cutoff",'st1_sdv','st2_mean','steepness_temp', ...
%     'lifespan_mean','lifespan_sdv','sickspan_mean','sickspan_sdv','skewness_temp')

% save("mat_data.mat")


%% plot
% load("mat_data.mat")

index_ref=find(experi_span==0);
lifespan_rela=lifespan./lifespan(:,index_ref);
sickspan_rela=sickspan./sickspan(:,index_ref);
steepness_rela=steepness./steepness(:,index_ref);
skewness_rela=skewness./skewness(:,index_ref);

figure
for i=1:length(experi_para)
    x_plot=lifespan_rela(i,:);
    y_plot=steepness_rela(i,:);
    plot(x_plot,y_plot,'LineWidth',2)
    hold on
end
hold off
title('Relative Steepness v.s Relative Lifespan')
legend(experi_para)
xlabel('lifespan_rela')
ylabel('steepness_rela')
axis on
grid on
% axis([0.95,1.6,0.95,1.6])


figure
for i=1:length(experi_para)
    x_plot=skewness_rela(i,:);
    y_plot=steepness_rela(i,:);
    plot(x_plot,y_plot,'LineWidth',2)
    hold on
    scatter(x_plot,y_plot,'+')
    hold on
end
hold off
title('Relative Skewness v.s Relative Steepness')
xlabel('skewness_rela')
ylabel('steepness_rela')
legend(repelem(experi_para,2))
axis on
grid on
% axis([0.95,1.6,0.95,1.6])


figure
for i=1:length(experi_para)
    x_plot=lifespan_rela(i,:);
    y_plot=steepness_rela(i,:);
    z_plot=skewness_rela(i,:);
    plot3(x_plot,y_plot,z_plot,'LineWidth',2)
    hold on
end
hold off
title('Relative Steepness, Relative Lifespan, Relative Skewness')
xlabel('lifespan_rela')
ylabel('steepness_rela')
zlabel('skewness_rela')
legend(experi_para)
grid on
axis on



% figure
% for i=1:length(experi_para)
%     x_plot=lifespan_rela(i,:);
%     y_plot=sickspan_rela(i,:);
%     plot(x_plot,y_plot,'LineWidth',2)
%     hold on
% end
% hold off
% title('Relative Sickspan v.s Relative Lifespan')
% legend(experi_para)
% axis on
% axis([1,1.2,0.7,1.2])

clear("x_plot",'y_plot')


%% save into .xlsx file
lifespan_t=name_maxtrix(lifespan,experi_para,string(experi_span));
lifespan_v_t=name_maxtrix(lifespan_v,experi_para,string(experi_span));
%sickspan_t=name_maxtrix(sickspan,experi_para,string(experi_span));
%sickspan_v_t=name_maxtrix(sickspan_v,experi_para,string(experi_span));
steepness_t=name_maxtrix(steepness,experi_para,string(experi_span));
skewness_t=name_maxtrix(skewness,experi_para,string(experi_span));

delete('./simu_result.xlsx')
writetable(lifespan_t, 'simu_result.xlsx', 'Sheet', 'lifespan_absolute','WriteRowNames', true, 'WriteVariableNames', true);
writetable(lifespan_v_t, 'simu_result.xlsx', 'Sheet', 'lifespan_std','WriteRowNames', true, 'WriteVariableNames', true);
%writetable(sickspan_t, 'simu_result.xlsx', 'Sheet', 'sickspan_absolute','WriteRowNames', true, 'WriteVariableNames', true);
%writetable(sickspan_v_t, 'simu_result.xlsx', 'Sheet', 'sickspan_std','WriteRowNames', true, 'WriteVariableNames', true);
writetable(steepness_t, 'simu_result.xlsx', 'Sheet', 'steepness','WriteRowNames', true, 'WriteVariableNames', true);
writetable(skewness_t,'simu_result.xlsx', 'Sheet', 'skewness','WriteRowNames', true, 'WriteVariableNames', true)

%% semi-analytical solve

% no way!
%% fitting
eta_x=lifespan_rela(1,:);
eta_y=steepness_rela(1,:);
eta_z=skewness_rela(1,:);
epsilon_x=lifespan_rela(2,:);
epsilon_y=steepness_rela(2,:);
epsilon_z=skewness_rela(2,:);
xc_x=lifespan_rela(3,:);
xc_y=steepness_rela(3,:);
xc_z=skewness_rela(3,:);

all_x=lifespan_rela(:);
all_y=steepness_rela(:);
all_z=skewness_rela(:);





%%

% function [X_paths,t]=SR_simu(param_struct,X0,M,EndTime,dt)
% % param_struct must contain element "eta","beta","kappa","epsilon"
% t=0:dt:EndTime;
% N=EndTime/dt;
% X_paths = X0.*ones(M, N+1);
% rng(666)
% for m = 1:M
%     for n = 1:N
%         mu = param_struct.eta * t(n) - param_struct.beta * (X_paths(m, n) / (X_paths(m, n) + param_struct.kappa));
%         dW = sqrt(dt) * randn;
%         X_paths(m,n+1) = X_paths(m,n) + mu*dt+ sqrt(2*param_struct.epsilon) * dW;
%         if X_paths(m,n+1)<=0
%             X_paths(m,n+1)=0.001;  % for avoiding /0
%         end
%     end
% end
% end

function [X_paths,t]=SR_simu2(param_struct,X0,M,EndTime,dt)
% param_struct must contain element "eta","beta","kappa","epsilon"
t=0:dt:EndTime;
N=EndTime/dt;
X_paths = X0.*ones(M, N+1);
rng(666)
for n = 1:N
    t_seq=t(n)*ones(M,1);
    mu = param_struct.eta .* t_seq- param_struct.beta .* (X_paths(:, n) ./ (X_paths(:, n) + param_struct.kappa));
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


function fieldName = num2name(k)
    str = sprintf('%.2f', k); 
    str = strrep(str, '.', '_');
    str = strrep(str, '-', 'm');
    fieldName = ['exp_', str];
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