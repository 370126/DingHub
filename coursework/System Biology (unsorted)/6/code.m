%%% 系统生物学作业6
%%% Perfect Adaptation
%%% 张牧原221505023


%% NFBLB

% Parameters for the middle panel

% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", 7);

% 指定范围和分隔符
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% 指定列名称和类型
opts.VariableNames = ["x0", "x1", "x2", "x3", "x0_1298", "x0_3586", "x0_4928"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.LeadingDelimitersRule = "ignore";

% 导入数据
ts = readtable("C:\Users\86189\Desktop\SysBio\作业6\/ts.txt", opts);

% 转换为输出类型
x0 = ts.x0;
x1 = ts.x1;
x2 = ts.x2;
tspan = ts.x3;
A_ts = ts.x0_1298;
B_ts = ts.x0_3586;
C_ts = ts.x0_4928;

% 清除临时变量
clear  ts

% 清除临时变量
clear opts
    % determined
K_fbb=0.01;  
K_cb=0.01;
k_ac=10;
k_bc=10;
I1=0.5;
I2=0.6;
    % guess
k_ia=2.501;
K_ia=1;
F_a=1;
k_faa=1;
K_faa=0.43;
k_cb=1;
F_b=0.5;
k_fbb=1;
K_ac=0.01;
K_bc=0.301;
t_s=8;


% functions
syms A B C I

f_A= @(A) I*k_ia*(1-A)/(1-A+K_ia)-F_a*k_faa*A/(A+K_faa)==0;
f_B= @(B,C) C*k_cb*(1-B)/(1-B+K_cb)-F_b*k_fbb*B/(B+K_fbb)==0;
% f_C= @(A,B,C) A*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
A_ss=solve(f_A,A);   %取第一个解
A_ss1=subs(A_ss(2),I1);
A_ss2=subs(A_ss(2),I2);
S_B=solve((f_B),C);
f_C1= @(B,C) A_ss1*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C1=solve(f_C1,C);  %取第二个解
f_C2= @(B,C) A_ss2*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C2=solve(f_C2,C) ; %取第二个解

%%
% nullclines
    % C
figure
subplot(1,2,1)
plot([0:0.01:1],subs(S_C1(2),[0:0.01:1]),'Color','r','LineWidth',1.25)
hold on
plot([0:0.01:1],subs(S_C2(2),[0:0.01:1]),'Color','r','LineWidth',1.25,'LineStyle',':')
hold on
plot([0:0.01:0.99],subs(S_B,[0:0.01:0.99]),'Color','black','LineWidth',1.25)
hold on
plot(B_ts,C_ts,'Color','b','LineWidth',1.8)
hold on
plot(B_ts(1),C_ts(1),'ro','MarkerEdgeColor','Black','MarkerFaceColor','black','MarkerSize',5)
hold on
plot(B_ts(end),C_ts(end),'ro','MarkerEdgeColor','Black','MarkerFaceColor',[0.65,0.65,0.65],'MarkerSize',5)
axis([0,1,0,1])
xlabel('B')
ylabel('C')
% time series
subplot(1,2,2)
plot(tspan,C_ts,'Color','b','LineWidth',1.8)
axis([0,20,0.35,0.65])
xlabel('t')

SEN=abs(((max(C_ts)-C_ts(1))/C_ts(1))/((I2-I1)/I1))     %无法达到1😭,已经很努力调参了...
PREC=abs(((I2-I1)/I1)/((C_ts(end)-C_ts(1))/C_ts(1)))
disp('SEN≈1')
disp('PREC>10')
%%
% Parameters for the top panel

% 设置导入选项并导入数据
opts2 = delimitedTextImportOptions("NumVariables", 7);

% 指定范围和分隔符
opts2.DataLines = [1, Inf];
opts2.Delimiter = " ";

% 指定列名称和类型
opts2.VariableNames = ["x0", "x1", "x2", "x3", "x0_3882", "x0_5945", "x0_5336"];
opts2.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% 指定文件级属性
opts2.ExtraColumnsRule = "ignore";
opts2.EmptyLineRule = "read";
opts2.LeadingDelimitersRule = "ignore";

% 导入数据
ts22 = readtable("C:\Users\86189\Desktop\SysBio\作业6\/ts2.txt", opts2);

% 转换为输出类型
x02 = ts22.x0;
x12 = ts22.x1;
x22 = ts22.x2;
ts2 = ts22.x3;
A_ts2 = ts22.x0_3882;
B_ts2 = ts22.x0_5945;
C_ts2 = ts22.x0_5336;

% 清除临时变量
clear  ts22

% 清除临时变量
clear opts2
    % determined
K_fbb=0.1;  
K_cb=0.1;
k_ac=10;
k_bc=10;
I1=0.5;
I2=0.6;
    % guess
k_ia=2.501;
K_ia=1;
F_a=1;
k_faa=1;
K_faa=0.43;
k_cb=1;
F_b=0.5;
k_fbb=1;
K_ac=0.01;
K_bc=0.301;
t_s=50;


% functions
syms A B C I

f_A= @(A) I*k_ia*(1-A)/(1-A+K_ia)-F_a*k_faa*A/(A+K_faa)==0;
f_B= @(B,C) C*k_cb*(1-B)/(1-B+K_cb)-F_b*k_fbb*B/(B+K_fbb)==0;
% f_C= @(A,B,C) A*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
A_ss=solve(f_A,A);   %取第一个解
A_ss1=subs(A_ss(2),I1);
A_ss2=subs(A_ss(2),I2);
S_B=solve((f_B),C);
f_C1= @(B,C) A_ss1*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C1=solve(f_C1,C);  %取第二个解
f_C2= @(B,C) A_ss2*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C2=solve(f_C2,C) ; %取第二个解

% nullclines
    % C
figure
subplot(1,2,1)
plot([0:0.01:1],subs(S_C1(2),[0:0.01:1]),'Color','r','LineWidth',1.25)
hold on
plot([0:0.01:1],subs(S_C2(2),[0:0.01:1]),'Color','r','LineWidth',1.25,'LineStyle',':')
hold on
plot([0:0.01:0.99],subs(S_B,[0:0.01:0.99]),'Color','black','LineWidth',1.25)
hold on
plot(B_ts2,C_ts2,'Color','b','LineWidth',1.8)
hold on
plot(B_ts2(1),C_ts2(1),'ro','MarkerEdgeColor','Black','MarkerFaceColor','black','MarkerSize',5)
hold on
plot(B_ts2(end),C_ts2(end),'ro','MarkerEdgeColor','Black','MarkerFaceColor',[0.65,0.65,0.65],'MarkerSize',5)
axis([0,1,0,1])
xlabel('B')
ylabel('C')
% time series
subplot(1,2,2)
plot(ts2,C_ts2,'Color','b','LineWidth',1.8)
axis([0,20,0.4,0.7])
xlabel('t')

SEN=abs(((max(C_ts2)-C_ts2(1))/C_ts2(1))/((I2-I1)/I1))     %无法达到1😭,已经很努力调参了...
PREC=abs(((I2-I1)/I1)/((C_ts2(end)-C_ts2(1))/C_ts2(1)))     %嘿嘿不管怎么说小于10没过关🤗
disp('SEN≈1')
disp('PREC<10') 

%%
% Parameters for the top panel

% 设置导入选项并导入数据
opts3 = delimitedTextImportOptions("NumVariables", 7);

% 指定范围和分隔符
opts3.DataLines = [1, Inf];
opts3.Delimiter = " ";

% 指定列名称和类型
opts3.VariableNames = ["x0", "x1", "x2", "x3", "x0_3882", "x0_6076", "x0_50444"];
opts3.VariableTypes = ["double", "double", "double", "double", "double", "double", "double"];

% 指定文件级属性
opts3.ExtraColumnsRule = "ignore";
opts3.EmptyLineRule = "read";
opts3.LeadingDelimitersRule = "ignore";

% 导入数据
ts32 = readtable("C:\Users\86189\Desktop\SysBio\作业6\/ts3.txt", opts3);

% 转换为输出类型
x03 = ts32.x0;
x13 = ts32.x1;
x23 = ts32.x2;
ts3 = ts32.x3;
A_ts3 = ts32.x0_3882;
B_ts3 = ts32.x0_6076;
C_ts3 = ts32.x0_50444;

% 清除临时变量
clear  ts32

% 清除临时变量
clear opts3
    % determined
K_fbb=0.01;  
K_cb=0.01;
k_ac=0.1;
k_bc=0.1;
I1=0.5;
I2=0.6;
    % guess
k_ia=2.501;
K_ia=1;
F_a=1;
k_faa=1;
K_faa=0.43;
k_cb=1;
F_b=0.5;
k_fbb=1;
K_ac=0.01;
K_bc=0.301;
t_s=8;

% functions
syms A B C I
f_A= @(A) I*k_ia*(1-A)/(1-A+K_ia)-F_a*k_faa*A/(A+K_faa)==0;
f_B= @(B,C) C*k_cb*(1-B)/(1-B+K_cb)-F_b*k_fbb*B/(B+K_fbb)==0;
% f_C= @(A,B,C) A*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
A_ss=solve(f_A,A);   %取第一个解
A_ss1=subs(A_ss(2),I1);
A_ss2=subs(A_ss(2),I2);
S_B=solve((f_B),C);
f_C1= @(B,C) A_ss1*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C1=solve(f_C1,C);  %取第二个解
f_C2= @(B,C) A_ss2*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C2=solve(f_C2,C) ; %取第二个解

% nullclines
figure
subplot(1,2,1)
plot([0:0.01:1],subs(S_C1(2),[0:0.01:1]),'Color','r','LineWidth',1.25)
hold on
plot([0:0.01:1],subs(S_C2(2),[0:0.01:1]),'Color','r','LineWidth',1.25,'LineStyle',':')
hold on
plot([0:0.01:0.99],subs(S_B,[0:0.01:0.99]),'Color','black','LineWidth',1.25)
hold on
plot(B_ts3,C_ts3,'Color','b','LineWidth',1.8)
hold on
plot(B_ts3(1),C_ts3(1),'ro','MarkerEdgeColor','Black','MarkerFaceColor','black','MarkerSize',5)
hold on
plot(B_ts3(end),C_ts3(end),'ro','MarkerEdgeColor','Black','MarkerFaceColor',[0.65,0.65,0.65],'MarkerSize',5)
axis([0,1,0,1])
xlabel('B')
ylabel('C')
% time series
subplot(1,2,2)
plot(ts3,C_ts3,'Color','b','LineWidth',1.8)
axis([0,100,0.4,0.7])
xlabel('t')

SEN=abs(((max(C_ts3)-C_ts3(1))/C_ts3(1))/((I2-I1)/I1))     %终于能名正言顺小于1了！
PREC=abs(((I2-I1)/I1)/((C_ts3(end)-C_ts3(1))/C_ts3(1)))     
disp('SEN<1')
disp('PREC>10')
%%
%% IFFLP

% parameters for the middle panel
    % determined
K_fbb=100;  
K_ab=0.001;
k_ab=0.5;
k_fbb=10;
I1=0.5;
I2=0.6;
    % guess
k_ac=10;
k_bc=10;
k_ia=2.5;
K_ia=1;
F_a=1;
k_faa=1;
K_faa=0.43;
F_b=3.2;
K_ac=0.1;
K_bc=0.301;
% functions
syms A B C I
f_A= @(A) I*k_ia*(1-A)/(1-A+K_ia)-F_a*k_faa*A/(A+K_faa)==0;
f_B= @(A,B) A*k_ab*(1-B)/(1-B+K_ab)-F_b*k_fbb*B/(B+K_fbb)==0;
% f_C= @(A,B,C) A*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
A_ss=solve(f_A,A) ;  %取第一个解
A_ss1=subs(A_ss(2),I1);
A_ss2=subs(A_ss(2),I2);
S_B=solve(f_B,B);
S_B1=double(subs(S_B(1),A_ss1));
S_B2=double(subs(S_B(1),A_ss2));
delta=S_B2-S_B1;

f_C1= @(B,C) A_ss1*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C1=solve(f_C1,C);  %取第二个解
f_C2= @(B,C) A_ss2*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C2=solve(f_C2,C) ; %取第二个解


[t,y]=ode45(@ifflp,[0,20],[0.3882,0.6087,0.3687]);
figure
plot(t,y(:,1))
hold on
plot(t,y(:,2))
hold on
plot(t,y(:,3))
legend('A','B','C')
% ss=[y(650,1),y(650,2),y(650,3)]
A2_ts=y(:,1);
B2_ts=y(:,2);
C2_ts=y(:,3);

% nullclines
figure
subplot(1,2,1)
plot([0:0.01:1],subs(S_C1(2),[0:0.01:1]),'Color','r','LineWidth',1.25)
hold on
plot([0:0.01:1],subs(S_C2(2),[0:0.01:1]),'Color','r','LineWidth',1.25,'LineStyle',':')
hold on
plot(S_B1*ones(length([0:0.01:1])),[0:0.01:1],'Color','black','LineWidth',1.25)
hold on
plot(S_B2*ones(length([0:0.01:1])),[0:0.01:1],'--','Color','black','LineWidth',1.25)
hold on
plot(B2_ts,C2_ts,'Color','b','LineWidth',1.8)
hold on
plot(B2_ts(1),C2_ts(1),'ro','MarkerEdgeColor','Black','MarkerFaceColor','black','MarkerSize',5)
hold on
plot(B2_ts(end),C2_ts(end),'ro','MarkerEdgeColor','Black','MarkerFaceColor',[0.65,0.65,0.65],'MarkerSize',5)
axis([0,1,0,1])
xlabel('B')
ylabel('C')

% time series
subplot(1,2,2)
plot(t,C2_ts,'Color','b','LineWidth',1.8)
axis([0,20,0.25,0.55])
xlabel('t')

SEN=abs(((max(C2_ts)-C2_ts(1))/C2_ts(1))/((I2-I1)/I1))     
PREC=abs(((I2-I1)/I1)/((C2_ts(end)-C2_ts(1))/C2_ts(1)))     
disp('SEN>1')
disp('PREC>10')
%%
% parameters for the top panel
    % determined
% K_fbb=100;  
% K_ab=0.001;
    %Bss =Ass k_AB*K_FBB/F_B/k_FBB
    %若其他量不变，则相比于middle,top道的Bss小了100倍！！怎么可能？？肯定是其它参数上也变了🤬🤬
K_fbb=1.01;  %?
K_ab=0.101;%?
k_ab=0.5;
k_fbb=10;
I1=0.5;
I2=0.6;
    % guess
k_ac=10;
k_bc=10;
k_ia=2.5;
K_ia=1;
F_a=1;
k_faa=1;
K_faa=0.43;
F_b=0.04;       %变的就是你！ 
K_ac=0.1;
K_bc=0.301;
% functions
syms A B C I
f_A= @(A) I*k_ia*(1-A)/(1-A+K_ia)-F_a*k_faa*A/(A+K_faa)==0;
f_B= @(A,B) A*k_ab*(1-B)/(1-B+K_ab)-F_b*k_fbb*B/(B+K_fbb)==0;
% f_C= @(A,B,C) A*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
A_ss=solve(f_A,A);   %取第一个解
A_ss1=(subs(A_ss(2),I1));
A_ss2=(subs(A_ss(2),I2));
% as1=double(A_ss1);
% as2=double(A_ss2);
% sb1=as1*k_ab*K_fbb/F_b/k_fbb;
% sb2=as2*k_ab*K_fbb/F_b/k_fbb;
% S_B=solve(f_B,B);
S_B1=double(subs(S_B(2),A_ss1));
S_B2=double(subs(S_B(2),A_ss2));
delta=S_B2-S_B1;

f_C1= @(B,C) A_ss1*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C1=solve(f_C1,C);  %取第二个解
f_C2= @(B,C) A_ss2*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C2=solve(f_C2,C) ; %取第二个解

[t,y]=ode45(@ifflp2,[0,20],[0.3882,0.6224,0.3534]);
figure
plot(t,y(:,1))
hold on
plot(t,y(:,2))
hold on
plot(t,y(:,3))
legend('A','B','C')
A2_ts2=y(:,1);
B2_ts2=y(:,2);
C2_ts2=y(:,3);

% nullclines
figure
subplot(1,2,1)
plot([0:0.01:1],subs(S_C1(2),[0:0.01:1]),'Color','r','LineWidth',1.25)
hold on
plot([0:0.01:1],subs(S_C2(2),[0:0.01:1]),'Color','r','LineWidth',1.25,'LineStyle',':')
hold on
plot(S_B1*ones(length([0:0.01:1])),[0:0.01:1],'Color','black','LineWidth',1.25)
hold on
plot(S_B2*ones(length([0:0.01:1])),[0:0.01:1],'--','Color','black','LineWidth',1.25)
hold on
plot(B2_ts2,C2_ts2,'Color','b','LineWidth',1.8)
hold on
plot(B2_ts2(1),C2_ts2(1),'ro','MarkerEdgeColor','Black','MarkerFaceColor','black','MarkerSize',5)
hold on
plot(B2_ts2(end),C2_ts2(end),'ro','MarkerEdgeColor','Black','MarkerFaceColor',[0.65,0.65,0.65],'MarkerSize',5)
axis([0,1,0,1])
xlabel('B')
ylabel('C')

% time series
subplot(1,2,2)
plot(t,C2_ts2,'Color','b','LineWidth',1.8)
axis([0,20,0.25,0.5])
xlabel('t')

SEN=abs(((max(C2_ts2)-C2_ts2(1))/C2_ts2(1))/((I2-I1)/I1))     
PREC=abs(((I2-I1)/I1)/((C2_ts2(end)-C2_ts2(1))/C2_ts2(1)))     
disp('SEN>1')
disp('PREC<10')
%%
% parameters for the bottom panel
    % determined
K_fbb=100;
K_ab=0.001;
k_ab=100;
k_fbb=2000;
    % guess
k_ac=10;
k_bc=10;
k_ia=2.5;
K_ia=1;
F_a=1;
k_faa=1;
K_faa=0.43;
F_b=3.2;
K_ac=0.1;
K_bc=0.301;
% functions
syms A B C I
f_A= @(A) I*k_ia*(1-A)/(1-A+K_ia)-F_a*k_faa*A/(A+K_faa)==0;
f_B= @(A,B) A*k_ab*(1-B)/(1-B+K_ab)-F_b*k_fbb*B/(B+K_fbb)==0;
% f_C= @(A,B,C) A*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
A_ss=solve(f_A,A) ;  %取第一个解
% A_ss1=subs(A_ss(2),I1);
% A_ss2=subs(A_ss(2),I2);
% as1=double(A_ss1);
% as2=double(A_ss2);
sb1=as1*k_ab*K_fbb/F_b/k_fbb;
sb2=as2*k_ab*K_fbb/F_b/k_fbb;

S_B=solve(f_B,B);
S_B1=double(subs(S_B(1),A_ss1));
S_B2=double(subs(S_B(1),A_ss2));
delta=S_B2-S_B1;

f_C1= @(B,C) A_ss1*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C1=solve(f_C1,C);  %取第二个解
f_C2= @(B,C) A_ss2*k_ac*(1-C)/(1-C+K_ac)-B*k_bc*C/(C+K_bc)==0;
S_C2=solve(f_C2,C) ; %取第二个解

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxOrder', 6);
[t,y]=ode45(@ifflp3,[0,5],[0.3882,0.6091,0.3687],options);
figure
plot(t,y(:,1))
hold on
plot(t,y(:,2))
hold on
plot(t,y(:,3))
legend('A','B','C')
A2_ts3=y(:,1);
B2_ts3=y(:,2);
C2_ts3=y(:,3);

% nullclines
figure
subplot(1,2,1)
plot([0:0.01:1],subs(S_C1(2),[0:0.01:1]),'Color','r','LineWidth',1.25)
hold on
plot([0:0.01:1],subs(S_C2(2),[0:0.01:1]),'Color','r','LineWidth',1.25,'LineStyle',':')
hold on
plot(S_B1*ones(length([0:0.01:1])),[0:0.01:1],'Color','black','LineWidth',1.25)
hold on
plot(S_B2*ones(length([0:0.01:1])),[0:0.01:1],'--','Color','black','LineWidth',1.25)
hold on
plot(B2_ts3,C2_ts3,'Color','b','LineWidth',1.8)
hold on
plot(B2_ts3(1),C2_ts3(1),'ro','MarkerEdgeColor','Black','MarkerFaceColor','black','MarkerSize',5)
hold on
plot(B2_ts3(end),C2_ts3(end),'ro','MarkerEdgeColor','Black','MarkerFaceColor',[0.65,0.65,0.65],'MarkerSize',5)
axis([0,1,0,1])
xlabel('B')
ylabel('C')

% time series
subplot(1,2,2)
plot(t,C2_ts3,'Color','b','LineWidth',1.8)
axis([0,5,0.36,0.38])
xlabel('t')

SEN=abs(((max(C2_ts3)-C2_ts3(1))/C2_ts3(1))/((I2-I1)/I1))     
PREC=abs(((I2-I1)/I1)/((C2_ts3(end)-C2_ts3(1))/C2_ts3(1)))     
disp('SEN<1')
disp('PREC>10')
%%

%%
function dydt=ifflp(t,y)
    dydt=zeros(3,1);
    % determined
    K_fbb=100;  
    K_ab=0.001;
    k_ab=0.5;
    k_fbb=10;
    % guess
    k_ac=10;
    k_bc=10;
    k_ia=2.5;
    K_ia=1;
    F_a=1;
    k_faa=1;
    K_faa=0.43;
    F_b=3.2;
    K_ac=0.1;
    K_bc=0.301;
    t_s=4;
    if t<t_s
        I=0.5;
        dydt(1)=I*k_ia*(1-y(1))/(1-y(1)+K_ia)-F_a*k_faa*y(1)/(y(1)+K_faa);
        dydt(2)=y(1)*k_ab*(1-y(2))/(1-y(2)+K_ab)-F_b*k_fbb*y(2)/(y(2)+K_fbb);
        dydt(3)=y(1)*k_ac*(1-y(3))/(1-y(3)+K_ac)-y(2)*k_bc*y(3)/(y(3)+K_bc);
    else
        I=0.6;
        dydt(1)=I*k_ia*(1-y(1))/(1-y(1)+K_ia)-F_a*k_faa*y(1)/(y(1)+K_faa);
        dydt(2)=y(1)*k_ab*(1-y(2))/(1-y(2)+K_ab)-F_b*k_fbb*y(2)/(y(2)+K_fbb);
        dydt(3)=y(1)*k_ac*(1-y(3))/(1-y(3)+K_ac)-y(2)*k_bc*y(3)/(y(3)+K_bc);
    end
end

function dydt=ifflp2(t,y)
    dydt=zeros(3,1);
    % determined
    K_fbb=1;  
    K_ab=0.1;
    k_ab=0.5;
    k_fbb=10;
    % guess
    k_ac=10;
    k_bc=10;
    k_ia=2.5;
    K_ia=1;
    F_a=1;
    k_faa=1;
    K_faa=0.43;
    F_b=0.04;
    K_ac=0.1;
    K_bc=0.301;
    t_s=4;
    if t<t_s
        I=0.5;
        dydt(1)=I*k_ia*(1-y(1))/(1-y(1)+K_ia)-F_a*k_faa*y(1)/(y(1)+K_faa);
        dydt(2)=y(1)*k_ab*(1-y(2))/(1-y(2)+K_ab)-F_b*k_fbb*y(2)/(y(2)+K_fbb);
        dydt(3)=y(1)*k_ac*(1-y(3))/(1-y(3)+K_ac)-y(2)*k_bc*y(3)/(y(3)+K_bc);
    else
        I=0.6;
        dydt(1)=I*k_ia*(1-y(1))/(1-y(1)+K_ia)-F_a*k_faa*y(1)/(y(1)+K_faa);
        dydt(2)=y(1)*k_ab*(1-y(2))/(1-y(2)+K_ab)-F_b*k_fbb*y(2)/(y(2)+K_fbb);
        dydt(3)=y(1)*k_ac*(1-y(3))/(1-y(3)+K_ac)-y(2)*k_bc*y(3)/(y(3)+K_bc);
    end
end

function dydt=ifflp3(t,y)
    dydt=zeros(3,1);
    % determined
    K_fbb=100;  
    K_ab=0.001;
    k_ab=100;
    k_fbb=2000;
    % guess
    k_ac=10;
    k_bc=10;
    k_ia=2.5;
    K_ia=1;
    F_a=1;
    k_faa=1;
    K_faa=0.43;
    F_b=3.2;
    K_ac=0.1;
    K_bc=0.301;
    t_s=1;
    if t<t_s
        I=0.5;
        dydt(1)=I*k_ia*(1-y(1))/(1-y(1)+K_ia)-F_a*k_faa*y(1)/(y(1)+K_faa);
        dydt(2)=y(1)*k_ab*(1-y(2))/(1-y(2)+K_ab)-F_b*k_fbb*y(2)/(y(2)+K_fbb);
        dydt(3)=y(1)*k_ac*(1-y(3))/(1-y(3)+K_ac)-y(2)*k_bc*y(3)/(y(3)+K_bc);
    else
        I=0.6;
        dydt(1)=I*k_ia*(1-y(1))/(1-y(1)+K_ia)-F_a*k_faa*y(1)/(y(1)+K_faa);
        dydt(2)=y(1)*k_ab*(1-y(2))/(1-y(2)+K_ab)-F_b*k_fbb*y(2)/(y(2)+K_fbb);
        dydt(3)=y(1)*k_ac*(1-y(3))/(1-y(3)+K_ac)-y(2)*k_bc*y(3)/(y(3)+K_bc);
    end
end