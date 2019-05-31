%Runs 4-population model with and without ACH
clear all
mas=[5]; %max strength of ach


depols=[0 1 1 1].*mas;

T=300; dt=.5; T=T/dt;
n=2.2; % 2.2 lin/nonlin. 1=lin, >1=n for nonlin
maxopt=0; %0 for nonlin
k=.01; %.01
Ne=1; Npv=1; Nvip=1; Nsom=1; %one rate unit per cell type
Ntot=Ne+Npv+Nvip+Nsom; % # of cell groups
taus=[20 10 10 10]; %time constants, e i



W=[.017, -.956, -.045, -.512;
    .8535, -.99, -.09, -.307;
    2.104, -.184, 0, -.734;
    1.285, 0, -.14, 0];



NeI=[1]; NpvI=[2]; NvipI=[3];  NsomI=[4];

for tri=1:2
    %tri=1 ach, tri=2 no ach.
    depols=depols.*(tri==1);
    SpontRates=[12.8 24 8 8.9];
    
    InpRates=[93 74 0 0];
    
    BinpE=ones(1,T).*SpontRates(1); BinpPV=ones(1, T).*SpontRates(2); BinpVIP=ones(1, T).*SpontRates(3);  BinpSOM=ones(1, T).*SpontRates(4); %no reason for this to be over time....
    
    InpE=zeros(Ne,T); InpPV=zeros(Npv,T); InpVIP=zeros(Nvip,T); InpSOM=zeros(Nsom,T);
    endT=300; stT=200;
    InpE(:,stT:endT)=InpE(:,stT:endT)+InpRates(1); InpPV(:,stT:endT)=InpPV(:,stT:endT)+InpRates(2); InpVIP(:,stT:endT)=InpVIP(:,stT:endT)+InpRates(3); InpSOM(:,stT:endT)=InpSOM(:,stT:endT)+InpRates(4);
    
    Re(:,1)=zeros(Ne,1); Rpv(:,1)=zeros(Npv,1); Rvip(:,1)=zeros(Nvip,1); Rsom(:,1)=zeros(Nsom,1);
    
    for t=2:T %run simulation
        dRe=(-Re(:,t-1)+k.*(max((InpE(:,t-1)+depols(1)+BinpE(1,t-1)+ W(NeI,:)*[Re(:,t-1); Rpv(:,t-1); Rvip(:,t-1); Rsom(:,t-1)]),maxopt)).^n)./taus(1);
        dRpv=(-Rpv(:,t-1)+k.*(max((InpPV(:,t-1)+depols(2)+BinpPV(1,t-1)+ W(NpvI,:)*[Re(:,t-1); Rpv(:,t-1); Rvip(:,t-1); Rsom(:,t-1)]),maxopt)).^n)./taus(2);
        dRvip=(-Rvip(:,t-1)+k.*(max((InpVIP(:,t-1)+depols(3)+BinpVIP(1,t-1)+ W(NvipI,:)*[Re(:,t-1); Rpv(:,t-1); Rvip(:,t-1); Rsom(:,t-1)]),maxopt)).^n)./taus(3);
        dRsom=(-Rsom(:,t-1)+k.*(max((InpSOM(:,t-1)+depols(4)+BinpSOM(1,t-1)+ W(NsomI,:)*[Re(:,t-1); Rpv(:,t-1); Rvip(:,t-1); Rsom(:,t-1)]),maxopt)).^n)./taus(4);
        
        
        
        
        Re(:,t)= Re(:,t-1)+dRe*dt;
        Rpv(:,t)= Rpv(:,t-1)+dRpv*dt;
        Rvip(:,t)= Rvip(:,t-1)+dRvip*dt;
        Rsom(:,t)= Rsom(:,t-1)+dRsom*dt;
        
    end
    allOT(tri,:,:)=[Re; Rpv; Rvip; Rsom];
    
    %spontaneous rates:
    rezes(tri,1:4)=[mean(Re(:,stT-100:stT-1),2), mean(Rpv(:,stT-100:stT-1),2), mean(Rvip(:,stT-100:stT-1),2), mean(Rsom(:,stT-100:stT-1),2)];
    %transient average:
    rezesL(tri,1:4)=[mean(Re(:,stT:stT+100),2), mean(Rpv(:,stT:stT+100),2), mean(Rvip(:,stT:stT+100),2), mean(Rsom(:,stT:stT+100),2)];
end




OTdifs=[allOT(1,1,:)-allOT(2,1,:); allOT(1,2,:)-allOT(2,2,:); allOT(1,3,:)-allOT(2,3,:); allOT(1,4,:)-allOT(2,4,:) ];
OTpeaks=max(allOT,[],3);
OTpdifs=[OTpeaks(1,1)-OTpeaks(2,1); OTpeaks(1,2)-OTpeaks(2,2); OTpeaks(1,3)-OTpeaks(2,3); OTpeaks(1,4)-OTpeaks(2,4) ];


mr=squeeze(rezes);
mrL=squeeze(rezesL);
ssdif=OTdifs(:,1,endT);


figure; subplot(1,2,1); hold on; plot(squeeze(allOT(2,:,:))'); xlabel('Time (ms)'); ylabel('Firing Rate')
ax=gca;
set(ax,'XTick',[0 200 400 600 800])
set(ax,'XTickLabel',{'0', '100', '200', '300', '400'}); title('No ACH'); legend({'E', 'P', 'V', 'S'})

subplot(1,2,2); hold on; plot(squeeze(allOT(1,:,:))'); xlabel('Time (ms)'); ylabel('Firing Rate')
ax=gca;
set(ax,'XTick',[0 200 400 600 800])
set(ax,'XTickLabel',{'0', '100', '200', '300', '400'}); title('ACH')

figure
subplot(2,3,1)
bar(squeeze(OTpeaks(2,:)))
title('No ACH Evoked Rates (Peak)')
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})
subplot(2,3,2)
bar(squeeze(mrL(2,:)))
title('No ACH Evoked Rates (50ms)')
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})
subplot(2,3,3)
bar(squeeze(allOT(2,:,endT)))
title('No ACH Evoked Rates (SS)')
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})

subplot(2,3,4)
bar(squeeze(OTpdifs))
title('Differences at Peak')
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})
subplot(2,3,5)
ax=gca;
bar([mrL(1,1)-mrL(2,1)   mrL(1,2)-mrL(2,2)    mrL(1,3)-mrL(2,3)  mrL(1,4)-mrL(2,4) ])
title('Differences over 50ms')
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})
subplot(2,3,6)
bar(squeeze(ssdif))
title('Differences at SS')
ax=gca;
set(ax,'XTick',[1 2 3 4])
set(ax,'XTickLabel',{'E','P','V','S'})


