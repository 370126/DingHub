%% 系统生物学作业5
%% Ocsillation
%% 张牧原 221505023

load 'raw.mat'
% parameters
k1=0.05;
kd_x=0.05;kd_y=0.05;
p=4;
S=1;
ks_y=1;
Kd=1;
Km=0.1;Ki=2;k2=1;Et=1;

% functions
syms x y
f_fx=@(x,y) k1*S*Kd^p/(Kd^p+y^p)-kd_x*x==0;
f_fy=@(x,y) ks_y*x-kd_y*y-k2*Et*y/(Km+y+Ki*y^2)==0;
S_x=solve(f_fx,x)
S_y=solve(f_fy,x)
%%
%% figure b
prot=[0:0.01:4];
figure
plot(prot,subs(S_x,prot),'LineWidth',1.75,'color',[103/255,188/255,117/255])
hold on
plot(prot,subs(S_y,prot),'LineWidth',1.75,'color','red')
hold on
plot(phase_x(1:15000),phase_y(1:15000),'--','LineWidth',1.75,'Color',[0.7,0.7,0.7])
legend('mRNA','Protein');
axis([0,4,0,1.1])
set(gca,'color',[247/255,247/255,241/255])
xlabel('Protein')
ylabel('mRNA')
title('b','FontSize',12)
%%
%% figure c
figure
yyaxis left
plot(tsx_x,tsx_y,'LineWidth',1.75,'color',[103/255,188/255,117/255])
ylim([0,1])
set(gca,'ycolor','k')
hold on
yyaxis right
plot(tsy_x,tsy_y,'LineWidth',1.75,'color','red')
ylim([0,4])
set(gca,'ycolor','k')
xlim([0,150])
legend('mRNA','Protein')
set(gca,'color',[247/255,247/255,241/255])
xlabel('Time')
title('c','FontSize',12)
%%
%% figure d  
figure
plot(bir_x,bir_y,'LineWidth',1.5,'Color','red')
hold on
plot(limimax_x(1:3:end),limimax_y(1:3:end),'ro','MarkerSize',3.75,'Color',[161/255,177/255,184/255],'MarkerFaceColor',[214/255,222/255,226/255])
hold on
plot(limimin_x(1:3:end),limimin_y(1:3:end),'ro','MarkerSize',3.75,'Color',[161/255,177/255,184/255],'MarkerFaceColor',[214/255,222/255,226/255])
axis([0,4,0,6])
legend('Protein')
xlabel('Signal')
ylabel('Protein')
set(gca,'color',[247/255,247/255,241/255])
title('d','FontSize',12)
%%
%% figure e
figure
plot(e1_x(837:end),e1_y(837:end),'LineWidth',2,'Color',[206/255,81/255,151/255])
hold on
plot(e2_x(565:end),e2_y(565:end),'LineWidth',2,'Color',[206/255,81/255,151/255])
hold on
fill([e1_x(end:-1:837);e2_x(565:end)],[e1_y(end:-1:837);e2_y(565:end)],[224/255,157/255,194/255],'EdgeColor','none')
% [224/255,157/255,194/255]
axis([0,6,0,10])
xlabel('K_d·K_i')
ylabel('p')
set(gca,'color',[247/255,247/255,241/255])
title('e','FontSize',12)
%%
%% figure f: one possible solution
figure
plot(f_x,f_y,'LineWidth',2,'Color',[206/255,81/255,151/255])
hold on
fill([f_x;0.0642],[f_y;0.2],[224/255,157/255,194/255],'EdgeColor','none')
axis([0.0642,5,0,0.2])
xlabel('S')
ylabel('Turnover of mRNA')
title('f','FontSize',12)
set(gca,'color',[247/255,247/255,241/255])
%%
%% figure f: another possible solution
figure
plot([f2_x;f_x],[f2_y;f_y],'LineWidth',2,'Color',[206/255,81/255,151/255])
hold on
fill([f2_x;f_x],[f2_y;f_y],[224/255,157/255,194/255],'EdgeColor','none')
axis([0,3,0,0.3])
set(gca,'color',[247/255,247/255,241/255])
title('f','Fontsize',12)