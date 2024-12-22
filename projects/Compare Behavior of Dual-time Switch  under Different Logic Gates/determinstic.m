%%% essay 
%%% 221505023 уедат╜
%% change A+B to A*B, compare the results under two conditions
%% Determinstic

[ts1,y_slow]=ode45(@two_slow,[0,1200],[26/251,0.01,0.01]);
[ts1o,y_slow_old]=ode45(@two_slow_old,[0,1200],[0.3571,0.01,0.01]);
[ts2,y_fast]=ode45(@two_fast,[0,1200],[26/251,0.01,0.01]);
[ts2o,y_fast_old]=ode45(@two_fast_old,[0,1200],[0.3571,0.01,0.01]);
[ts3,y_faslow]=ode45(@slow_fast,[0,1200],[26/251,0.01,0.01]);
[ts3o,y_faslow_old]=ode45(@slow_fast_old,[0,1200],[0.3571,0.01,0.01]);

%%
sti=zeros(120,1);
sti([1:200])=0;
sti([600:1200])=0;
sti([200:600])=1;
%%
figure
subplot(4,1,2)
plot(ts1o,y_slow_old(:,1),'LineWidth',1,'Color','black')
title('TWO slow loop')
axis off
subplot(4,1,1)
plot(ts2o,y_fast_old(:,1),'LineWidth',1,'Color','black')
title('Two fast loop')
axis off
subplot(4,1,3)
plot(ts3o,y_faslow_old(:,1),'LineWidth',1,'Color','black')
title('Two loops,dual time')
axis off
subplot(4,1,4)
plot([1:1200],sti,'LineWidth',1,'Color','black')
axis([0,1200,-0.5,3.7])
title('Stimulus')
axis off
xlabel('time')

figure
subplot(4,1,2)
plot(ts1,y_slow(:,1),'LineWidth',1,'Color','black')
title('TWO slow loop')
axis off
subplot(4,1,1)
plot(ts2,y_fast(:,1),'LineWidth',1,'Color','black')
title('Two fast loop')
axis off
subplot(4,1,3)
plot(ts3,y_faslow(:,1),'LineWidth',1,'Color','black')
title('Two loops,dual time')
axis off
subplot(4,1,4)
plot([1:1200],sti,'LineWidth',1,'Color','black')
axis([0,1200,-0.5,3.7])
title('Stimulus')
axis off
xlabel('time')
%%
function dydt=two_slow(t,y)
dydt=zeros(3,1);

t_s1=200;
t_s2=600;
k_out_on=2;
k_out_off=0.1;
k_out_min=0.01;
k_min=0.01;
n=3;
ec50=0.35;
tau=0.008;
if (t<t_s1 || t>t_s2)
    S=0;
    dydt(1)=k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau;
else
    S=1;
    dydt(1)=k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau;
end
end

function dydt=two_slow_old(t,y)
dydt=zeros(3,1);
t_s1=200;
t_s2=600;
k_out_on=2;
k_out_off=0.1;
k_out_min=0.01;
k_min=0.01;
n=3;
ec50=0.35;
tau=0.008;
if (t<t_s1 || t>t_s2)
    S=0;
    dydt(1)=k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau;
else
    S=1;
    dydt(1)=k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau;
end
end


function dydt=two_fast(t,y)
dydt=zeros(3,1);
t_s1=200;
t_s2=600;
k_out_on=2;
k_out_off=0.1;
k_out_min=0.01;
k_min=0.01;
n=3;
ec50=0.35;
tau=0.5;
if (t<t_s1 || t>t_s2)
    S=0;
    dydt(1)=k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau;
else
    S=1;
    dydt(1)=k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau;
end
end

function dydt=two_fast_old(t,y)
dydt=zeros(3,1);
t_s1=200;
t_s2=600;
k_out_on=2;
k_out_off=0.1;
k_out_min=0.01;
k_min=0.01;
n=3;
ec50=0.35;
tau=0.5;
if (t<t_s1 || t>t_s2)
    S=0;
    dydt(1)=k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau;
else
    S=1;
    dydt(1)=k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau;
end
end



function dydt=slow_fast(t,y)
dydt=zeros(3,1);
t_s1=200;
t_s2=600;
k_out_on=2;
k_out_off=0.1;
k_out_min=0.01;
k_min=0.01;
n=3;
ec50=0.35;
tau_A=0.5;
tau_B=0.008;
if (t<t_s1 || t>t_s2)
    S=0;
    dydt(1)=k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau_A;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau_B;
else
    S=1;
    dydt(1)=k_out_on*(2*y(2)*y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau_A;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau_B;
end
end


function dydt=slow_fast_old(t,y)
dydt=zeros(3,1);
t_s1=200;
t_s2=600;
k_out_on=2;
k_out_off=0.1;
k_out_min=0.01;
k_min=0.01;
n=3;
ec50=0.35;
tau_A=0.5;
tau_B=0.008;
if (t<t_s1 || t>t_s2)
    S=0;
    dydt(1)=k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau_A;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau_B;
else
    S=1;
    dydt(1)=k_out_on*(y(2)+y(3))*(1-y(1))-k_out_off*y(1)+k_out_min;
    dydt(2)=(S*(1-y(2))*y(1)^n/(y(1)^n+ec50^n)-y(2)+k_min)*tau_A;
    dydt(3)=(S*(1-y(3))*y(1)^n/(y(1)^n+ec50^n)-y(3)+k_min)*tau_B;
end
end