%%% Assignment 4
%% Ultrasensitivity in bistability
%% 221505023 张牧原

load ('bifurcation.mat')

%% parameters
k1=1;
k_1_high=1;
k_1_low=0.3;
k2_high=2;
k2_low=0.8;
K=0.01;
n=5;

x_tot=1;  %
input=1;  %

%% functions for dephosphorylation
d_m = @(xp) k_1_high.*xp;
d_u = @(xp) k_1_low.*(xp./(xp+K));

xxp=linspace(0,1,500)';
%%
%% simulations

% B
% fsolve初值
x_solve=[0,0.5];
k=1;
%x_Bcross=zeros(length(x_solve)*7);
v_d=d_m(xxp);
figure
subplot(1,2,1)
plot(xxp,v_d,'red')
hold on
for i=0:0.2:1.2
    input=i;
    p_m_high = @(xp) (k1*input+k2_high.*xp).*(x_tot-xp);
    v_p=p_m_high(xxp);
    plot(xxp,v_p,'blue')
    hold on
    % 找交点
    fun= @ (xp) p_m_high (xp)-d_m (xp);
    for j=1:length(x_solve)
        [x0,~,exitflag]=fsolve(fun,x_solve(j));
        if exitflag>0 && x0>=0
            plot(x0,d_m(x0),'.','MarkerSize',10,'Color','black')
            hold on
            x_Bcross(k)=x0;
            k=k+1;
        end
    end
end
axis([0,1,0,1.3])
xlabel('Xp/X\_tot')
ylabel('phosphorylation or dephosphorylation rate')
set(gca,'PlotBoxAspectRatio',[1 1.4 1]);
title('Linear positive feedback (strong)')

subplot(1,2,2)
plot(x_B,y_B,'LineWidth',1)
hold on
for i=1:length(x_Bcross)
    input_Bcross=(x_Bcross(i)-2*x_Bcross(i)^2)/(x_Bcross(i)-1);
    plot(input_Bcross,x_Bcross(i),'.','MarkerSize',10,'Color','black')
    hold on
end
axis([0,5,0,1])
xlabel('Input')
ylabel('steady state xp/x\_tot')
set(gca,'PlotBoxAspectRatio',[1 1.4 1]);
hold off
%%
% syms xpx inputx
% eqn= (k1*inputx+k2_high.*xpx).*(x_tot-xpx)-k_1_high*xpx==0
% S=solve(eqn,inputx,'Real',true)
%%

% C
k=1;
% x_Ccross=zeros(7,1);
v_d=d_m(xxp);
figure
subplot(1,2,1)
plot(xxp,v_d,'red')
hold on
for i=0:0.2:1.2
    input=i;
    p_m_low = @(xp) (k1*input+k2_low.*xp).*(x_tot-xp);
    v_p=p_m_low(xxp);
    plot(xxp,v_p,'blue')
    hold on
    % 找交点
    fun= @ (xp) p_m_low (xp)-d_m (xp);
    x0=fsolve(fun,[0.5]);
    plot(x0,d_m(x0),'.','MarkerSize',10,'Color','black')
    hold on
    x_Ccross(k)=x0;
    k=k+1;
end
axis([0,1,0,1.3])
xlabel('Xp/X\_tot')
ylabel('phosphorylation or dephosphorylation rate')
set(gca,'PlotBoxAspectRatio',[1 1.4 1]);
title('Linear positive feedback (weak)')
x_Ccross;


subplot(1,2,2)
plot(x_C,y_C,'LineWidth',1)
hold on
for i=1:length(x_Ccross)
    input_Ccross=(-x_Ccross(i)-4*x_Ccross(i)^2)/(5*x_Ccross(i)-5);
    plot(input_Ccross,x_Ccross(i),'.','MarkerSize',10,'Color','black')
    hold on
end
axis([0,5,0,1])
xlabel('Input')
ylabel('steady state xp/x\_tot')
set(gca,'PlotBoxAspectRatio',[1 1.4 1]);
hold off
%%
 % syms xpx inputx
 % eqn= (k1*inputx+k2_low.*xpx).*(x_tot-xpx)-k_1_high*xpx==0
 % S=solve(eqn,inputx,'Real',true)
%%
% D
% fsolve初值
x_solve=[0,0.42,0.5,0.7,1];
k=1;
m=1;
%x_Dcross=zeros(length(x_solve)*7);
v_d=d_m(xxp);
figure
subplot(1,2,1)
plot(xxp,v_d,'red')
hold on
for i=0:0.2:1.2
    input=i;
    p_u = @(xp) (k1*input+k2_high.*(xp.^n./(0.5.^n+xp.^n))).*(x_tot-xp);
    v_p=p_u(xxp);
    plot(xxp,v_p,'blue')
    hold on
    % 找交点
    fun= @ (xp) p_u(xp)-d_m (xp);
    for j=1:length(x_solve)
    [x0,~,exitflag]=fsolve(fun,x_solve(j));
        if exitflag>0 && x0>=0
            if (x0>0.3 && x0<0.52)
                plot(x0,d_m(x0),'ro','MarkerSize',3,'Color','black')
                hold on
                x_Dcross2(m)=x0;
                m=m+1;
            else
                plot(x0,d_m(x0),'.','MarkerSize',10,'Color','black')
                hold on
                x_Dcross(k)=x0;
                k=k+1;
            end
        end
    end
end
axis([0,1,0,1.3])
xlabel('Xp/X\_tot')
ylabel('phosphorylation or dephosphorylation rate')
set(gca,'PlotBoxAspectRatio',[1 1.4 1]);
title('Ultrasensitive feedback')

subplot(1,2,2)
plot(x_D,y_D,'LineWidth',1)
hold on
for i=1:length(x_Dcross)
    input_Dcross=(x_Dcross(i)+2*x_Dcross(i)^5*(x_Dcross(i)-1)/(x_Dcross(i)^5+1/32))/(1-x_Dcross(i));
    plot(input_Dcross,x_Dcross(i),'.','MarkerSize',10,'Color','black')
    hold on
end
for i=1:length(x_Dcross2)
    input_Dcross2=(x_Dcross2(i)+2*x_Dcross2(i)^5*(x_Dcross2(i)-1)/(x_Dcross2(i)^5+1/32))/(1-x_Dcross2(i));
    plot(input_Dcross2,x_Dcross2(i),'ro','MarkerSize',4,'Color','black')
    hold on
end
axis([0,1,0,1])
xlabel('Input')
ylabel('steady state xp/x\_tot')
set(gca,'PlotBoxAspectRatio',[1 1.4 1]);
hold off
%%
 % syms xpx inputx
 % eqn= (k1*inputx+k2_high*(xpx^n/(0.5^n+xpx^n)))*(x_tot-xpx)-k_1_high*xpx==0
 % S=solve(eqn,inputx,'Real',true)
%%
% E
v_d=d_u(xxp);
figure
subplot(1,2,1)
plot(xxp,v_d,'red')
hold on
% fsolve初值
x_solve=[0,0.42,0.5,0.7,1];
k=1;
%x_Ecross=zeros(length(x_solve)*7);
for i=0:0.2:1.2
    input=i;
    p_m_low = @(xp) (k1*input+k2_low.*xp).*(x_tot-xp);
    v_p=p_m_low(xxp);
    plot(xxp,v_p,'blue')
    hold on
    % 找交点
    fun= @ (xp) p_m_low(xp)-d_u (xp);
    for j=1:length(x_solve)
    [x0,~,exitflag]=fsolve(fun,x_solve(j));
        if exitflag>0 && x0>=0
             if (x0>0.3 && x0<0.6)
                 plot(x0,d_u(x0),'ro','MarkerSize',3,'Color','black')
                 hold on
             else
                plot(x0,d_u(x0),'.','MarkerSize',10,'Color','black')
                hold on
             end
             x_Ecross(k)=x0;
             k=k+1;
        end
    end
end
axis([0,1,0,1.3])
xlabel('Xp/X\_tot')
ylabel('phosphorylation or dephosphorylation rate')
set(gca,'PlotBoxAspectRatio',[1 1.4 1]);
title('Saturated dephosphorylation')

subplot(1,2,2)
plot(x_E,y_E,'LineWidth',1)
hold on
for i=1:length(x_Ecross)
    input_Ecross=(3*x_Ecross(i)/(10*x_Ecross(i)+1/10)+4*x_Ecross(i)*(x_Ecross(i)-1)/5)/(1-x_Ecross(i));
    if x_Ecross(i)>0.3 && x_Ecross(i)<0.6
    plot(input_Ecross,x_Ecross(i),'ro','MarkerSize',4,'Color','black')
    else
    plot(input_Ecross,x_Ecross(i),'.','MarkerSize',10,'Color','black')
    end
    hold on
end
axis([0,1,0,1])
xlabel('Input')
ylabel('steady state xp/x\_tot')
set(gca,'PlotBoxAspectRatio',[1 1.4 1]);
hold off
%%
 % syms xpx inputx
 % eqn= (k1*inputx+k2_low*xpx)*(x_tot-xpx)-k_1_low*xpx/(K+xpx)==0
 % S=solve(eqn,inputx,'Real',true)