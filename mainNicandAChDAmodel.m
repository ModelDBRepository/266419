%mex Both_NicandAChDAmodel.cpp
%
TT=11*60000; dt=0.02; N=TT/dt; % total time; step; numver of steps
ggaba=repmat(2.5,1,N); gbarnmda1=repmat(4,1,N); gampa=repmat(0,1,N);
gl=0.18/6; gbarh=0.05;

% nicotinic input like in Tolu paper
t=0:1:10*60*10^3/dt; %samples (msec/dt)
tdecay=3*60*10^3/dt;
trise=0.5*60*10^3/dt;
Nicmax=0.5645; %uM
Nic1=(tdecay/trise)^(1/((tdecay/trise)-1));
Nic=Nicmax*Nic1*(exp(-t/tdecay)-exp(-t/trise));
nicinput1=[repmat(0,1,60*10^3/dt),Nic];
nn=10; % number of trials
for j=1:nn %(Ach and Glu inputs change for each j)
    j
TT=11*60000;

CinppoisAch=bimodal_input(12000,TT); % create bursty ACh input (bimodal)
achinputda=  2*CinppoisAch(1:TT/dt); % ACh input to the DA neuron
achinputgaba=10*CinppoisAch(1:TT/dt); % ACh input to the GABA neuron
%frach=length(find(achinputgaba>0))/(TT*50)*10^3
CinppoisAch=[];
% Glu input produced by a population of poisson diftributed spike trains 
lambda=0.05; TT1=670000;
[CinppoisGlu,st]=Glu_population_fun(lambda,TT1);
CinppoisGlu=CinppoisGlu(1:TT/dt);

% a4b2 no nAChR
gbarachgaba=0; gbarachda=0; 
[Vmno] = Both_NicandAChDAmodel(TT,CinppoisGlu,ggaba,gbarnmda1,gampa,achinputda,achinputgaba,nicinput1,gbarachgaba,gbarachda,gl,gbarh);
stdano{j}=find(diff(Vmno>-40)>0); % spiketimes
   
% a4b2 nAChR on DA neurons only
gbarachgaba=0; gbarachda=5; 
[VmDA] = Both_NicandAChDAmodel(TT,CinppoisGlu,ggaba,gbarnmda1,gampa,achinputda,achinputgaba,nicinput1,gbarachgaba,gbarachda,gl,gbarh);
stdaDA{j}=find(diff(VmDA>-40)>0); % spiketimes
% a4b2 nAChR on GABA neurons only
gbarachgaba=4; gbarachda=0;
[VmGABA] = Both_NicandAChDAmodel(TT,CinppoisGlu,ggaba,gbarnmda1,gampa,achinputda,achinputgaba,nicinput1,gbarachgaba,gbarachda,gl,gbarh);
stdaGABA{j}=find(diff(VmGABA>-40)>0); % spiketimes
% a4b2 nAChR on both DA and GABA neurons
gbarachgaba=1.5; gbarachda=10;  % gbarachgaba=0.15; gbarachda=10; 
[Vmboth] = Both_NicandAChDAmodel(TT,CinppoisGlu,ggaba,gbarnmda1,gampa,achinputda,achinputgaba,nicinput1,gbarachgaba,gbarachda,gl,gbarh);
stdaboth{j}=find(diff(Vmboth>-40)>0); % spiketimes
% "wild" type
% nicotine increases the frequency of the Glu input
lambdaglu1=0.05;
CinppoisGlu1=Glu_population_fun(lambdaglu1,70000);
CinppoisGlu1=CinppoisGlu1(1:60000/dt);
lambdaglu2=0.054;
[CinppoisGlu2,st]=Glu_population_fun(lambdaglu2,610000);
CinppoisGlu2=CinppoisGlu2(1:600000/dt);
CinppoisGlu=[CinppoisGlu1',CinppoisGlu2'];
gbarachgaba=1.5; gbarachda=10;
[Vmwild,allgababoth,achcurrent,nmdasig] = Both_NicandAChDAmodel(TT,CinppoisGlu,ggaba,gbarnmda1,gampa,achinputda,achinputgaba,nicinput1,gbarachgaba,gbarachda,gl,gbarh);
stdawild{j}=find(diff(Vmwild>-40)>0); % spiketimes
end

stda=[];
stda=[stdano; stdaGABA; stdaDA; stdaboth; stdawild]; % spike times of DA neuron for all the cases
% calculate firing rate and %SWB at baseline
meanFRno=[]; meanFRGABA=[]; meanFRDA=[]; meanFRboth=[]; meanFRwild=[];
meanSWBno=[]; meanSWBGABA=[]; meanSWBDA=[]; meanSWBboth=[]; meanSWBwild=[];
Tbaseline=1*50*10^3;
%
for j=1:nn
    stdanobaseline=stda{1,j}(stda{1,j}<Tbaseline/dt);
    stdagababaseline=stda{2,j}(stda{2,j}<Tbaseline/dt);
    stdadabaseline=stda{3,j}(stda{3,j}<Tbaseline/dt);
    stdabothbaseline=stda{4,j}(stda{4,j}<Tbaseline/dt);
    stdawildbaseline=stda{5,j}(stda{5,j}<Tbaseline/dt);
    
    [meanFRno(j),meanSWBno(j)]=FR_SWB(stdanobaseline,Tbaseline)
    [meanFRGABA(j),meanSWBGABA(j)]=FR_SWB(stdagababaseline,Tbaseline)
    [meanFRDA(j),meanSWBDA(j)]=FR_SWB(stdadabaseline,Tbaseline)
    [meanFRboth(j),meanSWBboth(j)]=FR_SWB(stdabothbaseline,Tbaseline)
    [meanFRwild(j),meanSWBwild(j)]=FR_SWB(stdawildbaseline,Tbaseline)
end
%
meanFR=[meanFRno; meanFRGABA; meanFRDA; meanFRboth; meanFRwild];
% FR and %SWB
FRnoseg=[]; FRDAseg=[]; FRGABAseg=[]; FRbothseg=[]; FRwildseg=[];
SWBnoseg=[]; SWBDAseg=[]; SWBGABAseg=[]; SWBbothseg=[];  SWBwildseg=[];
n=TT/15000;
for j=1:nn
    for jj=1:n % calculate firing rate and %SWB in overlapping 30 s bins
        [FRnoseg(j,jj),SWBnoseg(j,jj)]=...
            FR_SWB(stda{1,j}(stda{1,j}>(30*10^3/dt*(jj-1)-15*10^3/dt*(jj-1))...
            & stda{1,j}<(30*10^3/dt*jj-15*10^3/dt*(jj-1))),30*10^3);
        
        [FRDAseg(j,jj),SWBDAseg(j,jj)]=...
            FR_SWB(stda{3,j}(stda{3,j}>(30*10^3/dt*(jj-1)-15*10^3/dt*(jj-1))...
            & stda{3,j}<(30*10^3/dt*jj-15*10^3/dt*(jj-1))),30*10^3);
        
        [FRGABAseg(j,jj),SWBGABAseg(j,jj)]=...
            FR_SWB(stda{2,j}(stda{2,j}>(30*10^3/dt*(jj-1)-15*10^3/dt*(jj-1))...
            & stda{2,j}<(30*10^3/dt*jj-15*10^3/dt*(jj-1))),30*10^3);
        
        [FRbothseg(j,jj),SWBbothseg(j,jj)]=...
            FR_SWB(stda{4,j}(stda{4,j}>(30*10^3/dt*(jj-1)-15*10^3/dt*(jj-1))...
            & stda{4,j}<(30*10^3/dt*jj-15*10^3/dt*(jj-1))),30*10^3);
        
        [FRwildseg(j,jj),SWBwildseg(j,jj)]=...
            FR_SWB(stda{5,j}(stda{5,j}>(30*10^3/dt*(jj-1)-15*10^3/dt*(jj-1))...
            & stda{5,j}<(30*10^3/dt*jj-15*10^3/dt*(jj-1))),30*10^3);
    end
end

% plot baseline firing rate and %SWB
meanfrda=[mean(meanFRno); mean(meanFRDA); mean(meanFRGABA); mean(meanFRboth); mean(meanFRwild)];
meanswbda=[mean(meanSWBno); mean(meanSWBDA); mean(meanSWBGABA); mean(meanSWBboth); mean(meanSWBwild)];
stdfrda=[std(meanFRno); std(meanFRDA); std(meanFRGABA); std(meanFRboth); std(meanFRwild)];
stdswbda=[std(meanSWBno); std(meanSWBDA); std(meanSWBGABA); std(meanSWBboth); std(meanSWBwild)];

% percent change in FR after Nic from baseline
perchangeFRDA=[]; perchangeFRno=[]; perchangeFRGABA=[]; perchangeFRboth=[]; perchangeFRwild=[];
perchangeSWBDA=[]; perchangeSWBno=[]; perchangeSWBGABA=[]; perchangeSWBboth=[]; perchangeSWBwild=[];
%
for j=1:nn
perchangeFRDA(j,:)=(FRDAseg(j,:)-meanfrda(2))/meanfrda(2)*100;
perchangeFRno(j,:)=(FRnoseg(j,:)-meanfrda(1))/meanfrda(1)*100;
perchangeFRGABA(j,:)=(FRGABAseg(j,:)-meanfrda(3))/meanfrda(3)*100;
perchangeFRboth(j,:)=(FRbothseg(j,:)-meanfrda(4))/meanfrda(4)*100;
perchangeFRwild(j,:)=(FRwildseg(j,:)-meanfrda(5))/meanfrda(5)*100;
end

% plot FR and %SWB at baseline
figure(1); clf
subplot(1,2,1)
for i=1:length(meanfrda)
    if i==1
        colorcode = 'r';
        bar(i, meanfrda(i), colorcode); hold on
        errorbar(i,meanfrda(i),stdfrda(i)/sqrt(nn));
        %errorbar(i,meanfrda(i),stdfrda(i));
    elseif i==2
        colorcode = 'g';
        bar(i, meanfrda(i), colorcode); hold on
        errorbar(i,meanfrda(i),stdfrda(i)/sqrt(nn));
        %errorbar(i,meanfrda(i),stdfrda(i));
    elseif i==3
        colorcode = 'b';
        bar(i, meanfrda(i), colorcode); hold on
    errorbar(i,meanfrda(i),stdfrda(i)/sqrt(nn));
    %errorbar(i,meanfrda(i),stdfrda(i));
    elseif i==4
        colorcode = [.5 0 .5];
        bar(i, meanfrda(i), 'facecolor',colorcode);
        hold on
    errorbar(i,meanfrda(i),stdfrda(i)/sqrt(nn));
   % errorbar(i,meanfrda(i),stdfrda(i));
    else
    colorcode = 'k';
    bar(i, meanfrda(i), 'facecolor',colorcode); hold on
    errorbar(i,meanfrda(i),stdfrda(i)/sqrt(nn));
    %errorbar(i,meanfrda(i),stdfrda(i));
    end
end
set(gca,'Xtick',1:5)
set(gca,'XTickLabel',{'neither','DA', 'GABA', 'both','wild'})
ylabel('Firing rate, Hz')
title('FR')

subplot(1,2,2)
for i=1:length(meanswbda)
    if i==1
        colorcode = 'r';
        bar(i, 100*meanswbda(i), colorcode); hold on
        %errorbar(i,100*meanswbda(i),100*stdswbda(i)/sqrt(nn));
        errorbar(i,100*meanswbda(i),100*stdswbda(i));
    elseif i==2
        colorcode = 'g';
        bar(i, 100*meanswbda(i), colorcode); hold on
        %errorbar(i,100*meanswbda(i),100*stdswbda(i)/sqrt(nn));
        errorbar(i,100*meanswbda(i),100*stdswbda(i));
    elseif i==3
        colorcode = 'b';
        bar(i, 100*meanswbda(i), colorcode); hold on
        %errorbar(i,100*meanswbda(i),100*stdswbda(i)/sqrt(nn));
        errorbar(i,100*meanswbda(i),100*stdswbda(i));
    elseif i==4
        colorcode = [.5 0 .5];
        bar(i, 100*meanswbda(i), 'facecolor',colorcode); hold on
        %errorbar(i,100*meanswbda(i),100*stdswbda(i)/sqrt(nn));
        errorbar(i,100*meanswbda(i),100*stdswbda(i));
    else
        colorcode = 'k';
        bar(i, 100*meanswbda(i), 'facecolor',colorcode); hold on
        %errorbar(i,100*meanswbda(i),100*stdswbda(i)/(sqrt(nn))); 
        errorbar(i,100*meanswbda(i),100*stdswbda(i));
    end
end
set(gca,'Xtick',1:5)
%ylim([0 20])
set(gca,'XTickLabel',{'neiher','DA', 'GABA', 'both','wild'})
ylabel('%SWB')
title('%SWB')


% plot firing rate and %SWB of DA neuron after exposure to micotine as percent change from baseline
 g=0.07; l=0.07;
figure(2); clf;
n1=n-1;
x=1:15:length(perchangeFRDA(1,:))*15;
subtightplot(5,2,1,g,l,l)
%plot(x,(perchangeFRDA),'g','linewidth',1)
plot(x,mean(perchangeFRDA),'g','linewidth',1)
hold on, plot(x,mean(perchangeFRDA)+std(perchangeFRDA),'g','linewidth',1)
hold on, plot(x,mean(perchangeFRDA)-std(perchangeFRDA),'g','linewidth',1)
h=line([1 n1*15],[0 0]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[-20 200]); set(h,'linestyle','--','color','black')
title(' FR percent change nAChR on DA'); ylabel ('%')
set(gca, 'Fontsize', 14)
ylim([-40 110]); xlim([0 500])
%
subtightplot(5,2,3,g,l,l)
%plot(x,(perchangeFRno(:,:)),'r','linewidth',2)
plot(x,mean(perchangeFRno),'r','linewidth',1)
hold on, plot(x,mean(perchangeFRno)+std(perchangeFRno),'r','linewidth',1)
hold on, plot(x,mean(perchangeFRno)-std(perchangeFRno),'r','linewidth',1)
h=line([1 n1*15],[0 0]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[-30 200]); set(h,'linestyle','--','color','black')
title(' FR percent change no nAChR'); ylabel ('%')
set(gca, 'Fontsize', 14)
ylim([-40 110]); xlim([0 500])

subtightplot(5,2,5,g,l,l)
%plot(x,(perchangeFRGABA(:,:)),'b','linewidth',2)
plot(x,mean(perchangeFRGABA),'b','linewidth',1)
hold on, plot(x,mean(perchangeFRGABA)+std(perchangeFRGABA),'b','linewidth',1)
hold on, plot(x,mean(perchangeFRGABA)-std(perchangeFRGABA),'b','linewidth',1)
h=line([1 n1*15],[0 0]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[-100 20]); set(h,'linestyle','--','color','black')
title(' FR percent change nAChR on GABA'); ylabel ('%')
set(gca, 'Fontsize', 14)
ylim([-110 40]); xlim([0 500])

subtightplot(5,2,7,g,l,l)
%plot(x,(perchangeFRboth(:,:)),'m','linewidth',2)
plot(x,mean(perchangeFRboth),'m','linewidth',1)
hold on, plot(x,mean(perchangeFRboth)+std(perchangeFRboth),'m','linewidth',1)
hold on, plot(x,mean(perchangeFRboth)-std(perchangeFRboth),'m','linewidth',1)
h=line([1 n1*15],[0 0]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[-20 40]); set(h,'linestyle','--','color','black')
title(' FR percent change nAChR on both'); ylabel ('%')
set(gca, 'Fontsize', 14)
ylim([-40 110]); xlim([0 500])

subtightplot(5,2,9,g,l,l)
%plot(x,(perchangeFRwild(:,:)),'k','linewidth',2)
plot(x,mean(perchangeFRwild),'k','linewidth',1)
hold on, plot(x,mean(perchangeFRwild)+std(perchangeFRwild),'k','linewidth',1)
hold on, plot(x,mean(perchangeFRwild)-std(perchangeFRwild),'k','linewidth',1)
h=line([1 n1*15],[0 0]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[-30 200]); set(h,'linestyle','--','color','black')
title(' FR percent change "wild type"'); ylabel ('%'); xlabel('time, sec')
set(gca, 'Fontsize', 14)
ylim([-40 110]); xlim([0 500])

% plot %SWB
limy=35; limy1=-2;
n1=n-1;
x=1:15:length(SWBDAseg(1,:))*15;
subtightplot(5,2,2,g,l,l)
%plot(x,(SWBDAseg(:,:)*100),'g','linewidth',2)
plot(x,mean(SWBDAseg*100),'g','linewidth',1)
hold on, plot(x,mean(SWBDAseg*100)+std(SWBDAseg*100),'g','linewidth',1)
hold on, plot(x,mean(SWBDAseg*100)-std(SWBDAseg*100),'g','linewidth',1)
h=line([1 n1*15],[meanswbda(2)*100 meanswbda(2)*100]); set(h,'linestyle','--','color','black')
h=line([4*15 4*15],[2 12])
set(h,'linestyle','--','color','black')
title(' %SWB nAChR no DA only'); ylabel ('%SWB')
set(gca, 'Fontsize', 14)
ylim([-20 74]); xlim([0 500])
%
subtightplot(5,2,4,g,l,l)
%plot(x,(SWBnoseg(:,:))*100,'r','linewidth',2)
plot(x,mean(SWBnoseg*100),'r','linewidth',1)
hold on, plot(x,mean(SWBnoseg*100)+std(SWBnoseg*100),'r','linewidth',1)
hold on, plot(x,mean(SWBnoseg*100)-std(SWBnoseg*100),'r','linewidth',1)
h=line([1 n1*15],[meanswbda(1)*100 meanswbda(1)*100]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[0 10]); set(h,'linestyle','--','color','black')
title(' %SWB no nAChR'); ylabel ('%SWB')
set(gca, 'Fontsize', 14)
ylim([-20 74]); xlim([0 500])

subtightplot(5,2,6,g,l,l)
%plot(x,(SWBGABAseg(:,:))*100,'b','linewidth',2)
plot(x,mean(SWBGABAseg*100),'b','linewidth',1)
hold on, plot(x,mean(SWBGABAseg*100)+std(SWBGABAseg*100),'b','linewidth',1)
hold on, plot(x,mean(SWBGABAseg*100)-std(SWBGABAseg*100),'b','linewidth',1)
h=line([1 n1*15],[meanswbda(3)*100 meanswbda(3)*100]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[0 30]); set(h,'linestyle','--','color','black')
title(' %SWB no nAChR'); ylabel ('%SWB')
set(gca, 'Fontsize', 14)
ylim([-20 74]); xlim([0 500])

subtightplot(5,2,8,g,l,l)
%plot(x,(SWBbothseg(:,:))*100,'m','linewidth',1)
plot(x,mean(SWBbothseg*100),'m','linewidth',1)
hold on, plot(x,mean(SWBbothseg*100)+std(SWBbothseg*100),'m','linewidth',1)
hold on, plot(x,mean(SWBbothseg*100)-std(SWBbothseg*100),'m','linewidth',1)
h=line([1 n1*15],[meanswbda(4)*100 meanswbda(4)*100]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[0 30]); set(h,'linestyle','--','color','black')
xlim([1 n1*15])
title(' %SWB no nAChR'); ylabel ('%SWB')
set(gca, 'Fontsize', 14)
ylim([-20 74]); xlim([0 500])

subtightplot(5,2,10,g,l,l)
%plot(x,(SWBwildseg(:,:))*100,'k','linewidth',1)
plot(x,mean(SWBwildseg*100),'k','linewidth',1)
hold on, plot(x,mean(SWBwildseg*100)+std(SWBwildseg*100),'k','linewidth',1)
hold on, plot(x,mean(SWBwildseg*100)-std(SWBwildseg*100),'k','linewidth',1)
h=line([1 n1*15],[meanswbda(4)*100 meanswbda(4)*100]); set(h,'linestyle','--','color','black')
h=line([3*15 3*15],[0 50]); set(h,'linestyle','--','color','black')
title(' %SWB no nAChR'); ylabel ('%SWB'); xlabel('time, sec')
set(gca, 'Fontsize', 14)
ylim([-20 74]); xlim([0 500])

%% uncomment to plot DA neuron voltage traces during baseline
% xlim1=10; xlim2=15;
% x=dt/10^3:dt/10^3*10:TT/10^3;
% g=0.07; l=0.07; clf
% subtightplot(5,1,1,g,l,l),plot(x,Vmno(1:10:end),'r','linewidth',1)
% title('KO'); xlim([xlim1 xlim2])
% subtightplot(5,1,2,g,l,l),plot(x,VmDA(1:10:end),'g','linewidth',1)
% title('nAChR on DA neuron'); xlim([xlim1 xlim2])
% subtightplot(5,1,3,g,l,l),plot(x,VmGABA(1:10:end),'b','linewidth',1)
% title('nAChR on GABA neurons'); xlim([xlim1 xlim2])
% subtightplot(5,1,4,g,l,l),plot(x,Vmboth(1:10:end),'m','linewidth',1)
% title('nAChR on both DA and GABA neuron'); xlim([xlim1 xlim2])
% subtightplot(5,1,5,g,l,l),plot(x,Vmwild(1:10:end),'k','linewidth',1)
% title('"Wild type"'); xlim([xlim1 xlim2])
% xlabel('Time, s'); ylabel('Voltage, mV')

% uncomment to plot DA neuron voltage traces during nicotine exposure
% xlim1=40; xlim2=100;
% x=dt/10^3:dt/10^3*10:TT/10^3;
% g=0.07; l=0.07; clf
% subtightplot(5,1,1,g,l,l),plot(x,Vmno(1:10:end),'r','linewidth',0.1)
% title('KO'); xlim([xlim1 xlim2]);
% subtightplot(5,1,2,g,l,l),plot(x,VmDA(1:10:end),'g','linewidth',0.1)
% title('nAChR on DA neuron'); xlim([xlim1 xlim2])
% subtightplot(5,1,3,g,l,l),plot(x,VmGABA(1:10:end),'b','linewidth',0.1)
% title('nAChR on GABA neurons'); xlim([xlim1 xlim2])
% subtightplot(5,1,4,g,l,l),plot(x,Vmboth(1:10:end),'m','linewidth',0.1)
% title('nAChR on both DA and GABA neuron'); xlim([xlim1 xlim2])
% subtightplot(5,1,5,g,l,l),plot(x,Vmwild(1:10:end),'k','linewidth',0.1)
% title('"Wild type"'); xlim([xlim1 xlim2])
% xlabel('Time, s'); ylabel('Voltage, mV')

