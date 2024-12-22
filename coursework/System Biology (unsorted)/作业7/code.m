%%% 系统生物学作业7
%% 221505023 张牧原
%% INTERLINKING POSITIVE AND NEGATIVE FEEDBACK

load ('data.mat')
% Stimu function
syms t
S(t)= piecewise((t>=500)&(t<600),5,1);
ts=[0:1:1500];
S_y=subs(S(t),ts);

figure
% a
subplot(3,2,1)
plot(bifur1_1x,bifur1_1y,'LineWidth',1.5,'Color','Black')
hold on
plot(bifur1_2x,bifur1_2y,'LineWidth',1.5,'Color','Black')
axis([0,5,0,4])
ylabel('[x]')
title('(a)')

% b
subplot(3,2,2)
yyaxis left
plot(ts1,x1,'LineWidth',1.5,'Color','Black')
set(gca,'ycolor','k')
ylim([0,4])
hold on
yyaxis right
plot(ts,S_y,'--','LineWidth',1.5,'Color','Black')
set(gca,'ycolor','k')
ylabel('S')
ylim([0,10])
title('(b)')

% c
subplot(3,2,3)
plot(bifur2_x,bifur2_y,'LineWidth',1.5,'Color','Black')
hold on
axis([0,5,0,2])
ylabel('[x]')
title('(c)')

% d
subplot(3,2,4)
yyaxis left
plot(ts2,x2,'LineWidth',1.5,'Color','Black')
set(gca,'ycolor','k')
ylim([0,4])
hold on
yyaxis right
plot(ts,S_y,'--','LineWidth',1.5,'Color','Black')
set(gca,'ycolor','k')
ylabel('S')
ylim([0,10])
title('(d)')

% e
subplot(3,2,5)
plot(bifur3_1x(1:4:end),bifur3_1y(1:4:end),'ro','MarkerSize',1,'Color','Black')
hold on
plot(bifur3_2x(1:4:end),bifur3_2y(1:4:end),'ro','MarkerSize',1,'Color','Black')
hold on
plot(bifur3x,bifur3y,'LineWidth',1.5,'Color','Black')
hold on
axis([0,5,0,2])
ylabel('[x]')
xlabel('S')
title('(e)')

% f 
S(t)= piecewise((t>=500)&(t<600),3.1,3.1);
ts=[0:1:1500];
S_y=subs(S(t),ts);

subplot(3,2,6)
yyaxis left
plot(ts3,x3,'LineWidth',1.5,'Color','Black')
set(gca,'ycolor','k')
ylim([0,2])
hold on
yyaxis right
plot(ts,S_y,'--','LineWidth',1.5,'Color','Black')
set(gca,'ycolor','k')
ylabel('S')
ylim([0,10])
xlabel('time')
title('(f)')