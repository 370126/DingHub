% OR gate
figure
x=0:0.01:1;
y=2.*x;
plot(x,y,'LineWidth',1.5,'Color','black')
hold on
plot(x,x,'--','LineWidth',1.5,'Color','black')
hold on
x=1:0.01:2;
y=2.*(2-x);
plot(x,y,'LineWidth',1.5,'Color','black')
hold on
plot(x,2-x,'--','LineWidth',1.5,'Color','black')
hold on
yy=0:0.05:2.1;
xx=ones(length(yy),1);
plot(xx,yy,'.','LineWidth',1,'Color','b')
legend('Output','Input')
title('OR gate')

% AND gate
figure
x=0:0.01:1;
y=2.*x.^2;
plot(x,y,'LineWidth',1.5,'Color','black')
hold on
plot(x,x,'--','LineWidth',1.5,'Color','black')
hold on
x=1:0.01:2;
y=2.*(2-x).^2;
plot(x,y,'LineWidth',1.5,'Color','black')
hold on
plot(x,2-x,'--','LineWidth',1.5,'Color','black')
hold on
yy=0:0.05:2.1;
xx=ones(length(yy),1);
plot(xx,yy,'.','LineWidth',1,'Color','b')
legend('Output','Input')
title('AND gate')