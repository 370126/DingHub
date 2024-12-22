%%% 系统生物学导论作业9
%%% 221505023 张牧原

%% 
% X + Y →c1→ 2Y
% 
% 2Y →c2→ Z

%% SSA

%开始计时
tic

% parameters
c1_X=5;
c2=0.00125;
X0=100;Z0=100;
Y0=12000;
    % S/P/A(1):reaction 1  S/P/A(2):reaction 2
S=[1,1,0;
    0,2,0];
P=[0,2,0;
    0,0,1];

figure
% for repeat=1:10
C=[X0,Y0,Z0]; %C(1):X  C(2):Y  C(3):Z
ts=0;

step_end=140000;  % 设置步长

for step=1:step_end
    A=[c1_X*C(end,2);
        c2*1/2*C(end,2)*(C(end,2)-1)];  
    a_0=sum(A);

    p1=rand(1);  % time step
    dt=(1/a_0)*log(1/p1);
    p2=rand(1);  % lucky bar
    for i=1:length(A)
        bar=p2*a_0;
        if i==1
            if bar<sum(A(1:i))
                r_luck=i;
                break
            end
        else
            if bar<sum(A(1:i)) && bar>sum(A(1:i-1))
                r_luck=i;
                break
            end
        end
    end
    ts(end+1)=ts(end)+dt;
    % update Components
    C(end+1,:)=C(end,:)-S(r_luck,:)+P(r_luck,:);
end
% 结束计时
running_time_SSA1=toc
running_time_per_time_SSA1=running_time_SSA1/ts(end)
running_time_per_step_SSA1=running_time_SSA1/step_end
plot(ts,C(:,2))
hold on
% end
axis([0,5,0,12000])
xlabel('Time')
ylabel('number of Y moleculars')
% average_running_time=running_time/10


tic
% parameters
Y0=40;
    % S/P/A(1):reaction 1  S/P/A(2):reaction 2
% figure
% for repeat=1:10
C=[X0,Y0,Z0]; %C(1):X  C(2):Y  C(3):Z
ts=0;
step_end=100000;  % 设置步长

for step=1:step_end
    A=[c1_X*C(end,2);
        c2*1/2*C(end,2)*(C(end,2)-1)];  
    a_0=sum(A);

    p1=rand(1);  % time step
    dt=(1/a_0)*log(1/p1);
    p2=rand(1);  % lucky bar
    for i=1:length(A)
        bar=p2*a_0;
        if i==1
            if bar<sum(A(1:i))
                r_luck=i;
                break
            end
        else
            if bar<sum(A(1:i)) && bar>sum(A(1:i-1))
                r_luck=i;
                break
            end
        end
    end
    ts(end+1)=ts(end)+dt;
    % update Components
    C(end+1,:)=C(end,:)-S(r_luck,:)+P(r_luck,:);
end
% 结束计时
running_time_SAA2=toc
running_time_per_time_SAA2=running_time_SAA2/ts(end)
running_time_per_step_SAA2=running_time_SAA2/step_end
plot(ts,C(:,2))
hold on
% end
axis([0,5,0,12000])
xlabel('Time')
ylabel('number of Y moleculars')
% average_running_time=running_time/10

%%
%% tau-leap

%开始计时
tic

% parameters
c1_X=5;
c2=0.00125;
X0=100;Z0=100;
Y0=12000;
    % S/P/A(1):reaction 1  S/P/A(2):reaction 2
S=[1,1,0;
    0,2,0];
P=[0,2,0;
    0,0,1];

C=[X0,Y0,Z0]; %C(1):X  C(2):Y  C(3):Z
ts=0;
step_end=10000*4.5;  % 设置步长
K=[0;        % number of each reaction that occur
    0];

figure
for step=1:step_end
    dt=0.0001;
    A=[c1_X*C(end,2);
            c2*1/2*C(end,2)*(C(end,2)-1)];  
    a_0=sum(A);
    Lamda=A.*dt;
    
    for i=1:length(A)
        lmd=Lamda(i);
        p_r=rand(1);
        p_sum=0;
        for n=0:1:1000
            p_n=lmd^n*exp(-lmd)/(factorial(n));
            p_sum=p_sum+p_n;
            if p_sum>p_r
                K(i)=n;
                break
            end
        end
    end
    ts(end+1)=ts(end)+dt;
    % update Components
    C_change=0*S(1,:);
    for k=1:length(K)
        Sk=S(k,:)*K(k);
        Pk=P(k,:)*K(k);
        C_change=C_change-Sk+Pk;
    end
    C(end+1,:)=C(end,:)+C_change;

end
running_time_tau1=toc
running_time_per_time_tau1=running_time_tau1/ts(end)
running_time_per_step_tau1=running_time_tau1/step_end

plot(ts,C(:,2))
hold on



tic
% parameters
Y0=40;
    % S/P/A(1):reaction 1  S/P/A(2):reaction 2
% figure
% for repeat=1:10
C=[X0,Y0,Z0]; %C(1):X  C(2):Y  C(3):Z
ts=0;
K=[0;        % number of each reaction that occur
    0];
for step=1:step_end
    dt=0.0001;
    A=[c1_X*C(end,2);
            c2*1/2*C(end,2)*(C(end,2)-1)];  
    a_0=sum(A);
    Lamda=A.*dt;
    
    for i=1:length(A)
        lmd=Lamda(i);
        p_r=rand(1);
        p_sum=0;
        for n=0:1:1000
            p_n=lmd^n*exp(-lmd)/(factorial(n));
            p_sum=p_sum+p_n;
            if p_sum>p_r
                K(i)=n;
                break
            end
        end
    end
    ts(end+1)=ts(end)+dt;
    % update Components
    C_change=0*S(1,:);
    for k=1:length(K)
        Sk=S(k,:)*K(k);
        Pk=P(k,:)*K(k);
        C_change=C_change-Sk+Pk;
    end
    C(end+1,:)=C(end,:)+C_change;

end
running_time_tau2=toc
running_time_per_time_tau2=running_time_tau2/ts(end)
running_time_per_step_tau2=running_time_tau2/step_end
plot(ts,C(:,2))
axis([0,5,0,12000])
xlabel('Time')
ylabel('number of Y moleculars')
hold on