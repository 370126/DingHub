%%% 系统生物学作业10
%%% 221505023 张牧原
%%% noise-induced stabilization of an unstable state

% strange parameters
    % 单位：/s
alpha_a=0.00875;
beta_a=7.5;
k_a=0.2;
lamda_a=1e-4;
alpha_r=0.025;
beta_r=2.5;
lamda_r=1e-4;
delta=4.1e-8;
n=2;
p=5;
Gamma=2.5e4;
K_a=k_a*Gamma;
% k_r=0.14;
% K_r=k_r*Gamma;

syms A R
f=alpha_a+beta_a*A^n/(K_a^n+A^n)-delta*A*R-lamda_a*A;
% g=alpha_r+beta_r*A^p/(K_r^p+A^p)-lamda_r*R;
% 零线
nuc_A=solve(f==0,R);
% nuc_R=solve(g==0,R);
% x_R=[2:30000];
x_A=[2:30000];
y_A=double(subs(nuc_A,[2:30000]));
% y_R=double(subs(nuc_R,[2:30000]));
%%
% different k_r
figure
loglog(x_A,y_A,'Color',[0.1294, 0.6078, 0.9916],'LineWidth',2)
hold on
k_r_set=[0.1,0.12,0.15,0.18,0.21,0.24];
for i=1:length(k_r_set)
k_r=k_r_set(i)
K_r=k_r*Gamma;
g=alpha_r+beta_r*A^p/(K_r^p+A^p)-lamda_r*R;
% 零线
nuc_R=solve(g==0,R);
x_R=[1000:30000];
y_R=double(subs(nuc_R,[1000:30000]));
% 找交点
h=matlabFunction(nuc_A-nuc_R);
cross=fsolve(h,[5000]);
A_f=cross    % focus
R_f=double(subs(nuc_R,cross))
loglog(x_R,y_R,'Color',[0.4863, 0.7216, 0.3490],'LineWidth',2)
hold on
    if k_r<0.1559
        plot(A_f,R_f,'d','MarkerEdgeColor','red','MarkerFaceColor','white','lineWidth',1,'MarkerSize',6)
        hold on
    else
        plot(A_f,R_f,'.','Color','red','MarkerSize',20)
        hold on
    end
end
k_r=0.1559;
K_r=k_r*Gamma;
g=alpha_r+beta_r*A^p/(K_r^p+A^p)-lamda_r*R;
% 零线
nuc_R=solve(g==0,R);
x_R=[1000:30000];
y_R=double(subs(nuc_R,[1000:30000]));
loglog(x_R,y_R,'Color','yellow','LineWidth',2)
hold on
loglog(x_R,y_R,'--','Color','black','LineWidth',1.5)

axis([800,25000,5000,30000])
xlabel('A (molec #)')
ylabel('R (molec #)')
%%
k_r=0.15;
K_r=k_r*Gamma;
g=alpha_r+beta_r*A^p/(K_r^p+A^p)-lamda_r*R;
% 零线
nuc_R=solve(g==0,R);
x_R=[2:20000];
y_R=double(subs(nuc_R,[2:20000]));

% 找交点
h=matlabFunction(nuc_A-nuc_R);
cross=fsolve(h,[100,400,3000]);
A_n=cross(1);    % node
R_n=double(subs(nuc_R,cross(2)));
A_s=cross(2);    % saddle
R_s=double(subs(nuc_R,cross(2)));
A_f=cross(3);    % focus
R_f=double(subs(nuc_R,cross(3)));

% time series
options=odeset('MaxStep',5);
[ts,y_ts]=ode45(@strange,[0,45*3600],[A_s+100,R_s-50],options);
yy_A=y_ts(:,1);
yy_R=y_ts(:,2);
figure
yyaxis left
plot(ts,yy_A,'Color',[0.1294, 0.6078, 0.9916],'LineWidth',2)
ylim([0,14e3])
ylabel('A (molec #)')
set(gca,'YColor',[0.1294, 0.6078, 0.9916])
hold on
yyaxis right
plot(ts,yy_R,'Color',[0.4863, 0.7216, 0.3490],'LineWidth',2)
ylim([0,2e4])
ylabel('R (molec #)')
set(gca,'YColor',[0.4863, 0.7216, 0.3490])
xlim([0,ts(end)])
xlabel('Time (s^-1)')

% nuclines & phase diagram
figure
loglog(x_A,y_A,'Color',[0.1294, 0.6078, 0.9916],'LineWidth',2)
hold on
loglog(x_R,y_R,'Color',[0.4863, 0.7216, 0.3490],'LineWidth',2)
hold on
plot(yy_A,yy_R,'Color','black','LineWidth',1)
hold on
plot(A_n,R_n,'.','Color','red','MarkerSize',20)
hold on
plot(A_s,R_s,'ro','MarkerEdgeColor','red','MarkerFaceColor','white','lineWidth',1)
hold on
plot(A_f,R_f,'d','MarkerEdgeColor','red','MarkerFaceColor','white','lineWidth',1,'MarkerSize',6)
axis([2,20000,30,50000])
xlabel('A (molec #)')
ylabel('R (molec #)')

%%
%% Discrete Stochastic Simulations

% 
% % SSA
% % strange parameters
%     % 单位：/s
% alpha_a=0.00875;
% beta_a=7.5;
% k_a=0.2;
% lamda_a=1e-4;
% alpha_r=0.025;
% beta_r=2.5;
% lamda_r=1e-4;
% delta=4.1e-8;
% n=2;
% p=5;
% Gamma=2.5e4;
% K_a=k_a*Gamma;
% k_r=0.15;
% K_r=k_r*Gamma;
% 
% A0=A_s+100;
% R0=R_s-50;
% ts=0;
% step_end=305000;
% 
% % reations
% S=[0,0;
%    0,0;
%    1,1;
%    1,0;
%    0,0;
%    0,0;
%    0,1];
% P=[1,0;
%    1,0;
%    0,1;
%    0,0;
%    0,1;
%    0,1;
%    0,0];
% 
% C=[A0,R0];
% 
% % SSA
% for step=1:step_end
% 
% A=[ alpha_a;
%     beta_a*C(end,1)^n/(K_a^n+C(end,1)^n);
%     delta*C(end,1)*C(end,2);
%     lamda_a*C(end,1);
%     alpha_r;
%     beta_r*C(end,1)^p/(K_r^p+C(end,1)^p);
%     lamda_r*C(end,2)  ];
% a_0=sum(A);
% p1=rand(1);  % time step
% dt=(1/a_0)*log(1/p1);
% p2=rand(1);  % lucky bar
% for i=1:length(A)
%     bar=p2*a_0;
%     if i==1
%         if bar<sum(A(1:i))
%             r_luck=i;
%             break
%         end
%     else
%         if bar<sum(A(1:i)) && bar>sum(A(1:i-1))
%             r_luck=i;
%             break
%         end
%     end
% end
% ts(end+1)=ts(end)+dt;
% % update Components
% C(end+1,:)=C(end,:)-S(r_luck,:)+P(r_luck,:);
% end
% figure
% yyaxis left
% plot(ts,C(:,1))
% hold on
% yyaxis right
% plot(ts,C(:,2))
% 
% 
% % nuclines & phase diagram
% figure
% loglog(x_A,y_A,'Color',[0.1294, 0.6078, 0.9916],'LineWidth',2)
% hold on
% loglog(x_R,y_R,'Color',[0.4863, 0.7216, 0.3490],'LineWidth',2)
% hold on
% plot(C(:,1),C(:,2),'Color','black','LineWidth',0.5)
% % hold on
% % plot(A_n,R_n,'.','Color','red','MarkerSize',20)
% % hold on
% % plot(A_s,R_s,'ro','MarkerEdgeColor','red','MarkerFaceColor','white','lineWidth',1)
% % hold on
% % plot(A_f,R_f,'d','MarkerEdgeColor','red','MarkerFaceColor','white','lineWidth',1,'MarkerSize',6)
% axis([2,20000,30,50000])
% xlabel('A (molec #)')
% ylabel('R (molec #)')
%%
% tau-leap

alpha_a=0.00875;
beta_a=7.5;
k_a=0.2;
lamda_a=1e-4;
alpha_r=0.025;
beta_r=2.5;
lamda_r=1e-4;
delta=4.1e-8;
n=2;
p=5;
Gamma=2.5e4;
K_a=k_a*Gamma;
k_r=0.15;
K_r=k_r*Gamma;

A0=A_s+100;
R0=R_s-50;
ts=0;
% reations
S=[0,0;
   0,0;
   1,1;
   1,0;
   0,0;
   0,0;
   0,1];
P=[1,0;
   1,0;
   0,1;
   0,0;
   0,1;
   0,1;
   0,0];


K=[ 0;        % number of each reaction that occur
    0;
    0;
    0;
    0;
    0;
    0];

step_end=60*3600/5;

% REPEAT! 

for repeat=1:1:10
C=[A0,R0];
ts=0;
for step=1:step_end
    dt=5;
    A=[ alpha_a;
    beta_a*C(end,1)^n/(K_a^n+C(end,1)^n);
    delta*C(end,1)*C(end,2);
    lamda_a*C(end,1);
    alpha_r;
    beta_r*C(end,1)^p/(K_r^p+C(end,1)^p);
    lamda_r*C(end,2)  ];

    Lamda=A.*dt;

    for i=1:length(A)
        lmd=Lamda(i);
        p_r=rand(1);
        p_sum=0;
        for n1=0:1:100
            p_n=lmd^n1*exp(-lmd)/(factorial(n1));
            p_sum=p_sum+p_n;
            if p_sum>p_r
                K(i)=n1;
                break
            end
        end
    end
    ts(end+1)=ts(end)+dt;
    % update Components
    C_change = -sum(S .* K) + sum(P .* K);
    C(end+1, :) = C(end, :) + C_change;
end

figure
loglog(x_A,y_A,'Color',[0.1294, 0.6078, 0.9916],'LineWidth',2)
hold on
loglog(x_R,y_R,'Color',[0.4863, 0.7216, 0.3490],'LineWidth',2)
hold on
plot(C(:,1),C(:,2),'Color','black','LineWidth',0.5)
axis([2,20000,30,50000])
xlabel('A (molec #)')
ylabel('R (molec #)')

 figure
 yyaxis left
 plot(ts,C(:,1),'Color',[0.1294, 0.6078, 0.9916],'LineWidth',2)
 ylim([0,14e3])
ylabel('A (molec #)')
 hold on
 yyaxis right
 plot(ts,C(:,2),'Color',[0.4863, 0.7216, 0.3490],'LineWidth',2)
 ylim([0,2e4])
ylabel('R (molec #)')

end





%%
function dydt=strange(t,y)
    dydt=zeros(2,1);
    % parameters
    alpha_a=0.00875;
    beta_a=7.5;
    k_a=0.2;
    lamda_a=1e-4;
    alpha_r=0.025;
    beta_r=2.5;
    lamda_r=1e-4;
    delta=4.1e-8;
    n=2;
    p=5;
    k_r=0.14;
    Gamma=2.5e4;
    K_a=k_a*Gamma;
    K_r=k_r*Gamma;
    dydt(1)=alpha_a+beta_a*y(1)^n/(K_a^n+y(1)^n)-delta*y(1)*y(2)-lamda_a*y(1);
    dydt(2)=alpha_r+beta_r*y(1)^p/(K_r^p+y(1)^p)-lamda_r*y(2);
end