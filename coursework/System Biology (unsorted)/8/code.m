%%% 系统生物学作业8
%%% 221505023 张牧原
%%% Noise

% 绘制噪声信号的图形
noise=stimu(20,60,120);
figure
plot([0:0.1:120],noise)
axis([0,120,-0.5,3])

%%
%% One slow loop
%parameters
t_s1=20;
t_s2=60;
t_end=120;
k_out_on=2;
k_out_off=0.3;
k_out_min=0.001;
k_min=0.01;
n=3;
ec50=0.35;
tao=0.008;
options=odeset('MaxStep',0.1);
dydt=@(t,y)[k_out_on*y(2)*(1-y(1))-k_out_off*y(1)+k_out_min;tao*(noise(ceil(round(t,1)*10))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min)];
[ts1,y_good1]=ode45(dydt,[1,120],[0.06562,0.01],options);
% for i=1:length(ts)
%     k=ceil(10*(round(ts(i),1)));
%     ts_test(i)=noise(k);
% end
% figure
% plot(ts,ts_test)

figure
subplot(3,1,1)
plot(ts1,y_good1(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.065,0.077])
title('OUT')
subplot(3,1,2)
plot(ts1,y_good1(:,2),'LineWidth',1,'Color','black')
title('A')
subplot(3,1,3)
plot([0:0.1:120],noise,'LineWidth',1,'Color','black')
title('noise')
axis([0,120,-0.5,4])

%%

%% One fast loop
%parameters
t_s1=20;
t_s2=60;
t_end=120;
k_out_on=4;
k_out_off=0.3;
k_out_min=0.001;
k_min=0.01;
n=3;
ec50=0.35;
tao=5;
options=odeset('MaxStep',0.1);
dydt=@(t,y)[k_out_on*y(2)*(1-y(1))-k_out_off*y(1)+k_out_min;tao*(noise(ceil(round(t,1)*10))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min)];
[ts3,y_good3]=ode45(dydt,[1,120],[0.1202,0.01],options);
% for i=1:length(ts)
%     k=ceil(10*(round(ts(i),1)));
%     ts_test(i)=noise(k);
% end
% figure
% plot(ts,ts_test)

figure
subplot(3,1,1)
plot(ts3,y_good3(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.1,1])
title('OUT')
subplot(3,1,2)
plot(ts3,y_good3(:,2),'LineWidth',1,'Color','black')
title('A')
subplot(3,1,3)
plot([0:0.1:120],noise,'LineWidth',1,'Color','black')
title('noise')
axis([0,120,-0.5,4])

%%
%% Two fast loop
%parameters
t_s1=20;
t_s2=60;
t_end=120;
k_out_on=4;
k_out_off=0.3;
k_out_min=0.001;
k_min=0.01;
n=3;
ec50=0.35;
tao_A=5;
tao_B=5;
options=odeset('MaxStep',0.1);
dydt=@(t,y)[k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(round(t,1)*10))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(round(t,1)*10))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts4,y_good4]=ode45(dydt,[1,120],[0.21296,0.01,0.01],options);
figure
subplot(4,1,1)
plot(ts4,y_good4(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.15,1])
title('OUT')
subplot(4,1,2)
plot(ts4,y_good4(:,2),'LineWidth',1,'Color','black')
title('A')
subplot(4,1,3)
plot(ts4,y_good4(:,3),'LineWidth',1,'Color','black')
title('B')
subplot(4,1,4)
plot([0:0.1:120],noise,'LineWidth',1,'Color','black')
title('noise')
axis([0,120,-0.5,4])
%%
%% Two slow loop
%parameters
t_s1=20;
t_s2=60;
t_end=120;
k_out_on=2;
k_out_off=0.3;
k_out_min=0.001;
k_min=0.01;
n=3;
ec50=0.35;
tao_A=0.008;
tao_B=0.008;
options=odeset('MaxStep',0.1);
dydt=@(t,y)[k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(round(t,1)*10))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(round(t,1)*10))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts2,y_good2]=ode45(dydt,[1,120],[0.1171,0.0097,0.0097],options);
% for i=1:length(ts)
%     k=ceil(10*(round(ts(i),1)));
%     ts_test(i)=noise(k);
% end
% figure
% plot(ts,ts_test)

figure
subplot(4,1,1)
plot(ts2,y_good2(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.1,0.47])
title('OUT')
subplot(4,1,2)
plot(ts2,y_good2(:,2),'LineWidth',1,'Color','black')
title('A')
subplot(4,1,3)
plot(ts2,y_good2(:,3),'LineWidth',1,'Color','black')
title('B')
subplot(4,1,4)
plot([0:0.1:120],noise,'LineWidth',1,'Color','black')
title('noise')
axis([0,120,-0.5,4])
%%
%% One fast and one slow loop
%parameters
t_s1=20;
t_s2=60;
t_end=120;
k_out_on=2;
k_out_off=0.3;
k_out_min=0.001;
k_min=0.01;
n=3;
ec50=0.35;
tao_A=5;
tao_B=0.008;
options=odeset('MaxStep',0.1);
dydt=@(t,y)[k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(round(t,1)*10))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(round(t,1)*10))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts5,y_good5]=ode45(dydt,[1,120],[0.1171,0.0097,0.0097],options);

figure
subplot(4,1,1)
plot(ts5,y_good5(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.08,0.9])
title('OUT')
subplot(4,1,2)
plot(ts5,y_good5(:,2),'LineWidth',1,'Color','black')
title('A')
subplot(4,1,3)
plot(ts5,y_good5(:,3),'LineWidth',1,'Color','black')
title('B')
subplot(4,1,4)
plot([0:0.1:120],noise,'LineWidth',1,'Color','black')
title('noise')
axis([0,120,-0.5,4])
%%
%% RESULT
figure
subplot(6,1,1)
plot(ts1,y_good1(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.065,0.077])
title('One slow loop')
axis off
subplot(6,1,2)
plot(ts2,y_good2(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.1,0.47])
title('Two slow loop')
axis off
subplot(6,1,3)
plot(ts3,y_good3(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.1,1])
title('One fast loop')
axis off
subplot(6,1,4)
plot(ts4,y_good4(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.15,1])
title('Two fast loop')
axis off
subplot(6,1,5)
plot(ts5,y_good5(:,1),'LineWidth',1,'Color','black')
axis([0,120,0.08,0.9])
title('Two loops,dual time')
axis off
subplot(6,1,6)
plot([0:0.1:120],noise,'LineWidth',1,'Color','black')
axis([0,120,-0.5,3.75])
axis off
title('stimulus')
xlabel('Time')
%%
function output=stimu(t_s1,t_s2,t_end)
    %parameters
    % t_s1=20;
    % t_s2=60;
    % t_end=120;
    k=1;
        %noise
    mu0=0.31;
    mu1=1;
    % dev=1;
        %probability
    p1=0.01;
    p2=0.70;
        %persist
    t0=0.75;
    t1=0.75;
 
    % tspan=length([0:t_end]);
    output=0*[0:0.1:t_end];
    %静息状态下，有~p1的概率产生均值为mu0,sd为1的噪声
    t=0.1;
    while t<t_s1
        
        if rand(1)<p1
            temp=mu0+1*randn;
            if temp>0
                k=ceil(t*10);
            for i=t:0.1:t+t0
                if (t+t0)<=t_s1
                output(k)=temp;
                k=k+1;
                end
            end
            t=t+t0;
            end
        end
        t=t+0.1;
    end
    
    t=t_s1;
    %刺激状态下，有~p2的概率产生均值为mu1,sd为0.75的噪声
    while t<=t_s2
        temp=mu1+0.75*randn;
        if rand(1)<p2 && temp>=0       
            k=floor(t*10);
                for i=t:0.1:t+t1
                    output(k)=temp;
                    k=k+1;
                end
            t=t+t1;
        else
            for i=t:0.1:t+t1
                output(k)=0;
                k=k+1;
            end
            t=t+t1;
        end
        t=t+0.1;
    end

    t=t_s2;
    while t_s2<=t && t<=t_end 
        
        if rand(1)<p1
            temp=mu0+1*randn;
            if temp>0
                k=ceil(t*10);
            for i=t:0.1:t+t0
                if (t+t0)<=t_end
                output(k)=temp;
                k=k+1;
                end
            end
            t=t+t0;
            end
        end
        t=t+0.1;
    end

end