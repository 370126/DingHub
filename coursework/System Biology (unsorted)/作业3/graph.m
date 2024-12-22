figure
plot(x1,y1./1,'LineWidth',1.25)
hold on
plot(x125,y125./1.25,'LineWidth',1.25)
hold on
plot(x15,y15./1.5,'LineWidth',1.25)
hold on
plot(x175,y175./1.75,'LineWidth',1.25)
hold on
plot(x2,y2./2,'LineWidth',1.25)
hold on
plot(x225,y225./2.25,'LineWidth',1.25)
hold on
plot(x25,y25./2.5,'LineWidth',1.25)
axis([0,5,0,0.6])

xlabel("Kinase\_total concentration")
ylabel("XPP/X_total")

legend('x\_total=1','1.25','1.5','1.75','2','2.25','2.5')