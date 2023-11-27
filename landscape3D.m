%% ‘ÿ»Î ˝æ›
load("landscape_data.mat")

%% stage
xn = linspace(0,30,90);
yn = linspace(0,20,60);
[Xn,Yn] = meshgrid(xn,yn);
Zn = griddata(LUADdegsursomnew,LUADdegsursomnewS1,LUADdegsursomnewS2,Xn,Yn,'v4');
% linear cubic natural nearest v4
figure
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName4==0);
scatter3(VarName1(m),VarName2(m),VarName3(m),'k','filled','DisplayName','Normal');
hold on
m=find(VarName4==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'y','filled','DisplayName','lung I');
m=find(VarName4==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','lung II');
m=find(VarName4==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','lung III');
m=find(VarName4==4);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','lung IV');
legend
%%
subplot(4,3,1)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName4==0);
scatter3(VarName1(m),VarName2(m),VarName3(m),'k','filled','DisplayName','Normal');
hold on
m=find(VarName4==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'y','filled','DisplayName','lung I');
m=find(VarName4==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','lung II');
m=find(VarName4==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','lung III');
m=find(VarName4==4);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','lung IV');
legend

subplot(4,3,2)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName4==0);
scatter3(VarName1(m),VarName2(m),VarName3(m),'k','filled','DisplayName','Normal');

subplot(4,3,3)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName4==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'y','filled','DisplayName','lung I');

subplot(4,3,4)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on 
m=find(VarName4==0);
scatter(VarName1(m),VarName2(m),'k','filled','DisplayName','Normal');
hold on
m=find(VarName4==1);
scatter(VarName1(m),VarName2(m),'y','filled','DisplayName','lung I');
m=find(VarName4==2);
scatter(VarName1(m),VarName2(m),'g','filled','DisplayName','lung II');
m=find(VarName4==3);
scatter(VarName1(m),VarName2(m),'b','filled','DisplayName','lung III');
m=find(VarName4==4);
scatter(VarName1(m),VarName2(m),'r','filled','DisplayName','lung IV');
legend

subplot(4,3,5)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on 
m=find(VarName4==0);
scatter(VarName1(m),VarName2(m),'k','filled','DisplayName','Normal');

subplot(4,3,6)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on 
m=find(VarName4==1);
scatter(VarName1(m),VarName2(m),'y','filled','DisplayName','lung I');

subplot(4,3,7)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName4==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','lung II');

subplot(4,3,8)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName4==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','lung III');

subplot(4,3,9)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName4==4);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','lung IV');

subplot(4,3,10)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on 
m=find(VarName4==2);
scatter(VarName1(m),VarName2(m),'g','filled','DisplayName','lung II');


subplot(4,3,11)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on 
m=find(VarName4==3);
scatter(VarName1(m),VarName2(m),'b','filled','DisplayName','lung III');

subplot(4,3,12)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on 
m=find(VarName4==4);
scatter(VarName1(m),VarName2(m),'r','filled','DisplayName','lung IV');


%%
figure
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on 
m=find(VarName4==0);
scatter(VarName1(m),VarName2(m),'k','filled','DisplayName','Normal');
hold on
m=find(VarName4==1);
scatter(VarName1(m),VarName2(m),'y','filled','DisplayName','lung I');
m=find(VarName4==2);
scatter(VarName1(m),VarName2(m),'g','filled','DisplayName','lung II');
m=find(VarName4==3);
scatter(VarName1(m),VarName2(m),'b','filled','DisplayName','lung III');
m=find(VarName4==4);
scatter(VarName1(m),VarName2(m),'r','filled','DisplayName','lung IV');
legend

figure
mesh(Xn,Yn,Zn)
hold on
plot3(LUADdegsursomnew,LUADdegsursomnewS1,LUADdegsursomnewS2,'r+','MarkerSize',3)

%% subtype
figure
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName5==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','Proximal inflammatory');
m=find(VarName5==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','Proximal proliferative');
m=find(VarName5==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','Terminal respiratory unit');
legend

figure
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on
m=find(VarName5==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','Proximal inflammatory');
m=find(VarName5==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','Proximal proliferative');
m=find(VarName5==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','Terminal respiratory unit');
legend

%% subtype_all
subplot(4,2,1)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName5==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','Proximal inflammatory');
m=find(VarName5==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','Proximal proliferative');
m=find(VarName5==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','Terminal respiratory unit');
legend

subplot(4,2,2)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName5==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','Terminal respiratory unit');

subplot(4,2,3)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on
m=find(VarName5==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','Proximal inflammatory');
m=find(VarName5==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','Proximal proliferative');
m=find(VarName5==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','Terminal respiratory unit');
legend

subplot(4,2,4)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on
m=find(VarName5==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','Terminal respiratory unit');

subplot(4,2,5)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName5==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','Proximal proliferative');

subplot(4,2,6)
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName5==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','Proximal inflammatory');


subplot(4,2,7)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on
m=find(VarName5==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','Proximal proliferative');

subplot(4,2,8)
contour(Xn,Yn,Zn,'DisplayName','Contour')
hold on
m=find(VarName5==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','Proximal inflammatory');

%% smoke years
figure
surfc(Xn,Yn,Zn,'DisplayName','Surface')
hold on 
m=find(VarName6==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','Non-smoker');
m=find(VarName6==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','0-20 years');
m=find(VarName6==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','20-40 years');
m=find(VarName6==4);
scatter3(VarName1(m),VarName2(m),VarName3(m),'k','filled','DisplayName','40-60 years');
legend

figure
contourf(LUADdegsursomnew,LUADdegsursomnewS1,LUADdegsursomnewS2,'DisplayName','Contour')
hold on 
% m=find(VarName5==0);
% scatter3(VarName1(m),VarName2(m),VarName3(m),'k','filled','DisplayName','NA');
hold on
m=find(VarName6==1);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','Non-smoker');
m=find(VarName6==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','0-20 years');
m=find(VarName6==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','20-40 years');
m=find(VarName6==4);
scatter3(VarName1(m),VarName2(m),VarName3(m),'k','filled','DisplayName','40-60 years');
legend


figure
contourf(LUADdegsursomnew,LUADdegsursomnewS1,LUADdegsursomnewS2,'DisplayName','Contour')
hold on 
% m=find(VarName5==0);
% scatter3(VarName1(m),VarName2(m),VarName3(m),'k','filled','DisplayName','NA');
hold on
m=find(VarName7==2);
scatter3(VarName1(m),VarName2(m),VarName3(m),'g','filled','DisplayName','0-40 pack');
m=find(VarName7==3);
scatter3(VarName1(m),VarName2(m),VarName3(m),'b','filled','DisplayName','40-80 pack');
m=find(VarName7==4);
scatter3(VarName1(m),VarName2(m),VarName3(m),'r','filled','DisplayName','80-120 pack');
m=find(VarName7==5);
scatter3(VarName1(m),VarName2(m),VarName3(m),'k','filled','DisplayName','120- pack');
legend

