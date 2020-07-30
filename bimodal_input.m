function [inpX,spiketimeX]=bimodal_input(nsp,T)
%%
% gives input with bimodal distribution in samples 1/dt to fit in the modal
% nsp is a number of spikes in the input
% T is a total time
nb=nsp;
for i=1:10
q = [.7 .3];
m = [30 300];
s = [7 7];
distrib = struct('mu', m, 'sigma', s, 'weight', q);
%nb = 750;
%nb=2000;
X = randn(nb,length(q)).*repmat(s,nb,1)+repmat(m,nb,1);
rsel    = rand(nb,1);
idx1    = (repmat(rsel,1,length(q))>repmat(cumsum(q),nb,1));
idx2    = (repmat(rsel,1,length(q))<repmat(cumsum(q),nb,1));
idx1(:,2:end) = idx1(:,1:end-1).*idx2(:,2:end);
idx1(:,1) = idx2(:,1);
X = sum(X.*idx1,2);
[y,x]=ksdensity(X);

spiketimeX=cumsum(X);
spiketimeX={spiketimeX};
end
%%
[rasterX]=histc(spiketimeX{1},0:1:spiketimeX{1}(end))';
%% firing rate
meanfrX=nb/spiketimeX{1}(end);
%% T sec
rasterX=rasterX(1,1:T)';
%% saveinp to read in C
dt=0.02;
inputX=zeros(length(rasterX),1/dt);
for i=1:length(rasterX)
    for j=1:1/dt
   inputX(i,j)= rasterX(i,1);
    end
end;

inputX=inputX';
inpX=inputX(:);
inpX=inpX';
