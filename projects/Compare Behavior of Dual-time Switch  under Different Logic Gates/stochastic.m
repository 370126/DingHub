%%% essay
%%% 221505023 张牧原
%%% stochastic model

% 绘制噪声信号的图形
noise=stimu(200,600,1200);


%% Two fast loop
%parameters
% t_s1=200;
% t_s2=600;
% t_end=1200;
k_out_on=2;
k_out_off=0.1;
k_out_min=0.01;
k_min=0.01;
n=3;
ec50=0.35;
tao_A=0.5;
tao_B=0.5;
dydt=@(t,y)[k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(t))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(t))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts4o,y_good4]=ode45(dydt,[1,1200],[0.3571,0.01,0.01]);
dydt=@(t,y)[k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(t))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(t))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts4,y_good4_new]=ode45(dydt,[1,1200],[26/251,0.01,0.01]);
%%
%% Two slow loop
tao_A=0.008;
tao_B=0.008;
dydt=@(t,y)[k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(t))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(t))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts2o,y_good2]=ode45(dydt,[1,1200],[0.3571,0.0097,0.0097]);
dydt=@(t,y)[k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(t))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(t))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts2,y_good2_new]=ode45(dydt,[1,1200],[26/251,0.0097,0.0097]);
%%
%% One fast and one slow loop
tao_A=0.5;
tao_B=0.008;
dydt=@(t,y)[k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(t))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(t))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts5o,y_good5]=ode45(dydt,[1,1200],[0.3571,0.0097,0.0097]);
dydt=@(t,y)[k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;tao_A*(noise(ceil(t))*y(1)^n*(1-y(2))/(y(1)^n+ec50^n)-y(2)+k_min);tao_B*(noise(ceil(t))*y(1)^n*(1-y(3))/(y(1)^n+ec50^n)-y(3)+k_min)];
[ts5,y_good5_new]=ode45(dydt,[1,1200],[26/251,0.0097,0.0097]);

%%
%% RESULT
figure
subplot(4,1,1)
plot(ts4,y_good4_new(:,1),'LineWidth',1,'Color','black')
title('TWO fast loop')
axis off
subplot(4,1,2)
plot(ts2,y_good2_new(:,1),'LineWidth',1,'Color','black')
title('Two SLOW loop')
axis off
subplot(4,1,3)
plot(ts5,y_good5_new(:,1),'LineWidth',1,'Color','black')
title('Two loops,dual time')
axis off
subplot(4,1,4)
plot([0:1200],noise,'LineWidth',1,'Color','black')
axis off
title('stimulus')
xlabel('Time')

figure
title('OUT A+B')
subplot(4,1,1)
plot(ts4o,y_good4(:,1),'LineWidth',1,'Color','black')
title('TWO fast loop')
axis off
subplot(4,1,2)
plot(ts2o,y_good2(:,1),'LineWidth',1,'Color','black')
title('Two SLOW loop')
axis off
subplot(4,1,3)
plot(ts5o,y_good5(:,1),'LineWidth',1,'Color','black')
title('Two loops,dual time')
axis off
subplot(4,1,4)
plot([0:1200],noise,'LineWidth',1,'Color','black')
axis off
title('stimulus')
xlabel('Time')
%%
function output=stimu(t_s1,t_s2,t_end)
    %parameters
    % t_s1=200;
    % t_s2=600;
    % t_end=1200;
    k=1;
        %noise
    mu0=0.31;
    mu1=1;
    % dev=1;
        %probability
    p1=0.01;
    p2=0.70;
        %persist
    t0=5;
    t1=5;
     output=0*[0:1:t_end];
    %静息状态下，有p1的概率产生dev为0.75的噪声
    t=1;
    while t<t_s1
         if rand(1)<p1
            temp=mu0+1*randn;
            if temp>0
                k=ceil(t);
            for i=t:1:t+t0
                if (t+t0)<=t_s1
                output(k)=temp;
                k=k+1;
                end
            end
            t=t+t0;
            end
        end
        t=t+1;
    end
    t=t_s1;
    %刺激状态下，有p2的概率产生均值为mu1的噪声
    while t<=t_s2
        temp=mu1+0.75*randn;
        if rand(1)<p2 && temp>=0       
            k=floor(t);
                for i=t:1:t+t1
                    output(k)=temp;
                    k=k+1;
                end
            t=t+t1;
        else
            for i=t:1:t+t1
                output(k)=0;
                k=k+1;
            end
            t=t+t1;
        end
        t=t+1;
    end

    t=t_s2;
    while t_s2<=t && t<=t_end 
        
        if rand(1)<p1
            temp=mu0+1*randn;
            if temp>0
                k=ceil(t);
            for i=t:1:t+t0
                if (t+t0)<=t_end
                output(k)=temp;
                k=k+1;
                end
            end
            t=t+t0;
            end
        end
        t=t+1;
    end

end