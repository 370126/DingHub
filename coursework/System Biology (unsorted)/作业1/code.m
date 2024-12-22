%%% 系统生物学导论 作业1 
%%% 221505023 张牧原
%%

%  Fast
tao_A_fast=0.5;
tao_B_fast=0.5;
%  Slow
tao_A_slow=0.008;
tao_B_slow=0.008;

% Parameters
k_out_on=2;
k_out_off=0.3;
k_out_min=0.001;
ec50=0.35;
k_min=0.01;
n=3;            %Hill factor
t_s1=20;        %刺激开始时间
t_s2=60;        %刺激结束时间
t_b=0;          %记录开始时间
t_f=120;        %记录结束时间
h=0.1;          %求解ODE步长

% steady state
A_slow_steady =0.0097;
OUT_slow_steady=0.1171;
B_slow_steady=0.0097;
A_fast_steady=0.0097;
OUT_fast_steady=0.1171;
B_fast_steady=0.0097;

A_1slow_steady =0.0100;
OUT_1slow_steady=0.0656;
A_1fast_steady=0.0100;
OUT_1fast_steady=0.0656;


% Stimulus function
syms t
syms stimu(t)
stimu(t)=piecewise((t>=t_s1)&(t<t_s2),1,0)
x=[t_b:0.1:t_f];
y=stimu(x);
plot(x,y)
axis([t_b,t_f,-0.5,1.5])
title("Stimulus")

% ODEs
% 1 loop
syms OUT A
syms f_out(OUT,A,t)
f_out(OUT,A,t)=k_out_on*A*(1-OUT)-k_out_off*OUT+k_out_min
syms f_A_fast(OUT,A,t)
f_A_fast(OUT,A,t)=(stimu(t)*OUT^n/(OUT^n+ec50^n)*(1-A)-A+k_min)*tao_A_fast
syms f_A_slow(OUT,A,t)
f_A_slow(OUT,A,t)=(stimu(t)*OUT^n/(OUT^n+ec50^n)*(1-A)-A+k_min)*tao_A_slow
% 2 loops
syms OUT A B
syms g_out(OUT,A,B,t)
g_out(OUT,A,B,t)=k_out_on*(A+B)*(1-OUT)-k_out_off*OUT+k_out_min
syms g_A_fast(OUT,A,B,t)
g_A_fast(OUT,A,B,t)=((stimu(t)*OUT^n/(OUT^n+ec50^n))*(1-A)-A+k_min)*tao_A_fast
syms g_B_fast(OUT,A,B,t)
g_B_fast(OUT,A,B,t)=((stimu(t)*OUT^n/(OUT^n+ec50^n))*(1-B)-B+k_min)*tao_B_fast
syms g_A_slow(OUT,A,B,t)
g_A_slow(OUT,A,B,t)=((stimu(t)*OUT^n/(OUT^n+ec50^n))*(1-A)-A+k_min)*tao_A_slow
syms g_B_slow(OUT,A,B,t)
g_B_slow(OUT,A,B,t)=((stimu(t)*OUT^n/(OUT^n+ec50^n))*(1-B)-B+k_min)*tao_B_slow

%%
%% 1 FAST
% Initial State
A_1fast=[t_b:h:t_f]*0;
OUT_1fast=[t_b:h:t_f]*0;
A0 = A_1fast_steady;         %初值
OUT0=OUT_1fast_steady;

A_1fast(1)=A0;
OUT_1fast(1)=OUT0;


i=1;
% SOLVE ODE
for t=t_b+h:h:t_f
    %OUT
    square_OUT=f_out(OUT_1fast(i),A_1fast(i),t)*h;
    OUT_1fast(i+1)=OUT_1fast(i)+square_OUT;
    %A
    square_A=f_A_fast(OUT_1fast(i),A_1fast(i),t)*h;
    A_1fast(i+1)=A_1fast(i)+square_A;
    i=i+1;
end

t_line=[t_b:h:t_f];
plot(t_line,OUT_1fast,'Color','Black','LineWidth',1.25)
axis([0,100,0,1])
title("OUT fast loop")
plot(t_line,A_1fast)
title("A fast loop")

%%

%% 1 SLOW
% Initial State
A_1slow=[t_b:h:t_f]*0;
OUT_1slow=[t_b:h:t_f]*0;
A0 = A_1slow_steady;         %初值
OUT0=OUT_1slow_steady;

A_1slow(1)=A0;
OUT_1slow(1)=OUT0;

i=1;
% SOLVE ODE
for t=t_b+h:h:t_f
    %OUT
    square_OUT=f_out(OUT_1slow(i),A_1slow(i),t)*h;
    OUT_1slow(i+1)=OUT_1slow(i)+square_OUT;
    %A
    square_A=f_A_slow(OUT_1slow(i),A_1slow(i),t)*h;
    A_1slow(i+1)=A_1slow(i)+square_A;
    i=i+1;
end

t_line=[t_b:h:t_f];
plot(t_line,OUT_1slow)
% axis([0,100,0,0.1])
title("OUT slow loop")
plot(t_line,A_1slow)
title("A slow loop")

%%
%% 2 FAST

% Initial State
A_fast=[t_b:h:t_f]*0;
OUT_fast=[t_b:h:t_f]*0;
B_fast=[t_b:h:t_f]*0;

A0 = A_fast_steady;         %初值
OUT0=OUT_fast_steady;
B0=B_fast_steady;

A_fast(1)=A0;
OUT_fast(1)=OUT0;
B_fast(1)=B0;


i=1;
% SOLVE ODE
for t=t_b+h:h:t_f
    %OUT
    square_OUT=g_out(OUT_fast(i),A_fast(i),B_fast(i),t)*h;
    OUT_fast(i+1)=OUT_fast(i)+square_OUT;
    %A
    square_A=g_A_fast(OUT_fast(i),A_fast(i),B_fast(i),t)*h;
    A_fast(i+1)=A_fast(i)+square_A;
    %B
    square_B=g_B_fast(OUT_fast(i),A_fast(i),B_fast(i),t)*h;
    B_fast(i+1)=B_fast(i)+square_B;
    
    i=i+1;
end

t_line=[t_b:h:t_f];
plot(t_line,OUT_fast)
title("OUT 2 fast loop")
plot(t_line,A_fast)
title("A 2 fast loop")
plot(t_line,B_fast)
title("B 2 fast loop")

%%

%% 2 SLOW

% Initial State
A_slow=[t_b:h:t_f]*0;
OUT_slow=[t_b:h:t_f]*0;
B_slow=[t_b:h:t_f]*0;

A_slow(1)=A0;
OUT_slow(1)=OUT0;
B_slow(1)=B0;

i=1;
% SOLVE ODE
for t=t_b+h:h:t_f
    %OUT
    square_OUT=g_out(OUT_slow(i),A_slow(i),B_slow(i),t)*h;
    OUT_slow(i+1)=OUT_slow(i)+square_OUT;
    %A
    square_A=g_A_slow(OUT_slow(i),A_slow(i),B_slow(i),t)*h;
    A_slow(i+1)=A_slow(i)+square_A;
    %B
    square_B=g_B_slow(OUT_slow(i),A_slow(i),B_slow(i),t)*h;
    B_slow(i+1)=B_slow(i)+square_B;
    
    i=i+1;
end

t_line=[t_b:h:t_f];
plot(t_line,OUT_slow)
title("OUT 2 slow loop")
plot(t_line,A_slow)
title("A 2 slow loop")
plot(t_line,B_slow)
title("B 2 slow loop")
%%
%% A SLOW B FAST

% Initial State
A_faslow=[t_b:h:t_f]*0;
OUT_faslow=[t_b:h:t_f]*0;
B_faslow=[t_b:h:t_f]*0;

A0 = A_slow_steady;         %初值
OUT0=OUT_fast_steady;
B0=B_fast_steady;

A_faslow(1)=A0;
OUT_faslow(1)=OUT0;
B_faslow(1)=B0;


i=1;
% SOLVE ODE
for t=t_b+h:h:t_f
    %OUT
    square_OUT=g_out(OUT_faslow(i),A_faslow(i),B_faslow(i),t)*h;
    OUT_faslow(i+1)=OUT_faslow(i)+square_OUT;
    %A
    square_A=g_A_slow(OUT_faslow(i),A_faslow(i),B_faslow(i),t)*h;
    A_faslow(i+1)=A_faslow(i)+square_A;
    %B
    square_B=g_B_fast(OUT_faslow(i),A_faslow(i),B_faslow(i),t)*h;
    B_faslow(i+1)=B_faslow(i)+square_B;
    
    i=i+1;
end

t_line=[t_b:h:t_f];
figure
plot(t_line,OUT_faslow)
title("OUT 1 fast 1 slow loop")
figure
plot(t_line,A_faslow)
title("A 1 fast 1 slow loop")
figure
plot(t_line,B_faslow)
title("B 1 fast 1 slow loop")

%%
% PLOT

subplot(6,1,6)
plot(x,y,'Color','Black','LineWidth',1.25)
axis([t_b,t_f,-0.5,3])
axis off
title("Stimulus")

subplot(6,1,4)
plot(t_line,OUT_fast,'Color','Black','LineWidth',1.25)
axis off
title("Two fast loops")

subplot(6,1,2)
plot(t_line,OUT_slow,'Color','Black','LineWidth',1.25)
%axis([0,120,0.1,0.55])
axis off
title("Two slow loops")

subplot(6,1,5)
plot(t_line,OUT_faslow,'Color','Black','LineWidth',1.25)
axis off
title("Two loops, dual time")


subplot(6,1,3)
plot(t_line,OUT_1fast,'Color','Black','LineWidth',1.25)
axis off
title("OUT fast loop")

subplot(6,1,1)
plot(t_line,OUT_1slow,'Color','Black','LineWidth',1.25)
%axis([0,120,0.0645,0.083])
axis off
title("OUT slow loop")