%% Load data
load('data_cn_project_iii_a17.mat')

%% Question 1
% Auto correlation
tic
R = @(tou)(dot(Stimulus(51:end-50),Stimulus(51+tou:end-50+tou))/(length(Stimulus)-100) );

Tou = -50:1:50;
Auto = ones(size(Tou));

parfor i=1:101
    Auto(i) = R(Tou(i));
end
toc
figure;
plot(Tou,Auto);
%disp(Auto)

%% Question 2
% PSTH - Trail averaged spike rate
tic
rates = zeros(4,15000);

parfor i=1:4
    v=zeros(1,15000);
    for j=1:50
        for k=All_Spike_Times{i,j}         
            in = int64(k*1000);            
            if in<=15000
                v(in) = v(in) + 1;
            else
                break;
            end
        end
    end
    rates(i,:) = v/50;
end

rates = rates*100;

figure;

subplot(2,2,1);
plot(1:15000,rates(1,:));
title('Neuron 1')

subplot(2,2,2);
plot(1:15000,rates(2,:));
title('Neuron 2')

subplot(2,2,3);
plot(1:15000,rates(3,:));
title('Neuron 3')

subplot(2,2,4);
plot(1:15000,rates(4,:));
title('Neuron 4')
toc
%% Question 3
tic
valueset = {zeros(4,50,2001),zeros(4,50,1001),zeros(4,50,401),zeros(4,50,201),zeros(4,50,101),zeros(4,50,41)};
keyset = [10,20,50,100,200,500];
rateStore = containers.Map(keyset,valueset);

valueset2 = {zeros(4,2001),zeros(4,1001),zeros(4,401),zeros(4,201),zeros(4,101),zeros(4,41)};
meanStore = containers.Map(keyset,valueset2);
varStore = containers.Map(keyset,valueset2);

for i=keyset
    v1 = zeros(4,50,(20000/i)+1);
    parfor j=1:4
        v2 = zeros(1,50,(20000/i)+1);
        for k=1:50
            v3 = zeros(1,1,(20000/i)+1); 
            for l=All_Spike_Times{j,k}         
                in = int64(l*1000/(i))+1;                            
                v3(1,1,in) = v3(1,1,in) + 1;
            end 
            v2(1,k,:) = v3(1,1,:);
        end
        v1(j,:,:) = v2(1,:,:); 
    end
    rateStore(i) = v1;
end

for i = keyset
    meanStore(i) = mean(rateStore(i),2);
    varStore(i) = var(rateStore(i),0,2);
end

for j=1:4
    figure;
    index=0;
    for i = keyset
        index=index+1;
        m = meanStore(i);
        v = varStore(i);
        m1 = m(j,1,:);
        v1 = v(j,1,:);
        subplot(2,3,index);
        plot(m1(:),v1(:),'.')
    end
end
toc
%% Question 4 PART B
tic
costs = [10^-9,10^-6,0.0001,0.001,0.01,0.1,1,10,100,1000];
qvalues = length(costs);
s = RandStream('mlfg6331_64');

processedspikes = cell(4,50,200);

parfor i=1:4
    tempCell = cell(50,200);
    for j=1:50
        for k=All_Spike_Times{i,j}
            ind = ceil(k*1000/100);
            spiketrain = tempCell{j,ind};
            tempCell{j,ind} = [spiketrain k];
        end        
    end
    processedspikes(i,:,:) = tempCell;
end

save( 'samplespikes.mat', 'processedspikes');
toc

%%
tic
confusionMatrixTotal = zeros(100,4,8,8,qvalues); 

parfor i=1:100
    processedspi = load('samplespikes.mat');
    processedspikes = processedspi.processedspikes;
    
    y = datasample(s,1:200,8,'Replace',false);  
    
    confusionMatrix = zeros(4,8,8,qvalues);
    comparitions = ones(400,400,qvalues)*Inf;
    
    for j=1:4
        perNeuron = processedspikes(j,:,:);
        for k=1:400
                in1 = ceil(k/50);
                ST1 = perNeuron{1,in1,y(in1)};
                for m=1:400 % compare with which                     
                        if ~(k==m) 
                            if comparitions(k,m,:) == Inf
                                in2 = ceil(m/50);
                                value = spkd_qpara(ST1,perNeuron{1,in2,y(in2)},costs);                
                                comparitions(k,m,:) = value;
                                comparitions(m,k,:) = value;
                            end
                            %comparitions((m-1)*50+n,(k-1)*50+l,:) = value;                                                 
                        end
                end
                [argvalue, index] = min(reshape(comparitions(k,:,:),[400,qvalues]),[],1);
                index=ceil(index/50);                
                %confusionMatrix(j,in1,index(1,1,:),[1,2,3,4,5,6]) = confusionMatrix(j,in1,index(1,1,:),[1,2,3,4,5,6]) + 1;
                for o=1:length(index)
                    confusionMatrix(j,in1,index(o),o) = confusionMatrix(j,in1,index(o),o) + 1;
                end
                       
        end
    end
    confusionMatrixTotal(i,:,:,:,:) = confusionMatrix/50;
end
toc

%% Calculate mutual information
tic

MIData = cell(100,4);
parfor i=1:100
    for j=1:4
        temp = [];
        for k=1:qvalues
            P12 = reshape(confusionMatrixTotal(i,j,:,:,k),[8,8]);
            P1 = reshape(sum(P12,1)/8,[8,1]) ;
            P2 = reshape(sum(P12,2)/8,[1,8]) ;    
%             if i==1 && j ==1 && k==1
%                 disp(P1)
%                 disp(P2')
%                 disp(P12)
%             end
            Com = P12.*log((P12+10^-6)./(P1*P2+10^-6));
            
            MI = sum(Com);            
            MI = sum(MI);
            
            temp = [temp MI(1,1)];
        end
        MIData{i,j} = temp;
    end
end
toc
%%

for j=1:4
    figure;
    hold on
    ylim([0 100])
    xlim([-10 1000])


    for i=1:100
        %disp(squeeze(MIData{i,1}))
        plot(1./costs,squeeze(MIData{i,1}),'--')
    end

end
%%
function d=spkd_qpara(tli,tlj,costs)
    %
    % d=spkd(tli,tlj,costs) calculates the "spike time" distance
    % (Victor & Purpura 1996) for a single costs
    %
    % tli: vector of spike times for first spike train
    % tlj: vector of spike times for second spike train
    % costs: costs per unit time to move a spike
    %
    % Copyright (c) 1999 by Daniel Reich and Jonathan Victor.
    % Translated to Matlab by Daniel Reich from FORTRAN code by Jonathan Victor.

    nspi=length(tli);
    nspj=length(tlj);
    if costs==0
        d=abs(nspi-nspj);
        return
    elseif costs==Inf
        d=nspi+nspj;
        return
    end
    scr=zeros(nspi+1,nspj+1);
    scr(:,1)=(0:nspi)';
    scr(1,:)=(0:nspj);
    scr=repmat(shiftdim(scr,-1),[length(costs),1,1]);
    if nspi && nspj
        for i=2:nspi+1
            for j=2:nspj+1           
                scr(:,i,j)=min(cat(3,scr(:,i-1,j)+1,scr(:,i,j-1)+1,scr(:,i-1,j-1)+costs'*abs(tli(i-1)-tlj(j-1))),[],3);
            end
        end
    end
    d=scr(:,nspi+1,nspj+1);
end

%% 

function [clast,rmax,c]=spkdallq_recur(sa,sb)
% function [clast,rmax,c]=spkdallqk_recur(sa,sb) does the recursion for the DP
% algorithm for spike time distance.  It uses the method of spkdallq, but
% reformats and produces arguments to be parallel to spkdallqk_recur.
%
% sa: vector of spike times for first spike train
% sb: vector of spike times for second spike train
% c: array of critical lengthts of linkages. size(c)=[1+length(sa) 1+length(sb) 1+rmax]
% rmax: maximum number of linkages (=min(length(sa),length(sb))
% clast: vector of critical lengths for full spike trains, necessary for calculating all
%    distances, which is done in spkdallq_dist, via spkdallq_final.
%
%  Copyright (c) by Jonathan Victor.
%
%    See also SPKD, SPKDALLQ, SPKDALLQK_RECUR.
%
tli=sa; %reformat the input arguments
tlj=sb;
nspi=length(tli);
nspj=length(tlj);
nij=min(nspi,nspj);
lc=repmat(NaN,[nspi+1 nspj+1 nij]); % assume no length of movement
if (nij>0)
%
%     best linkage length is a recursion based on linkages for shorter trains
%
	for i=2:nspi+1
   	for j=2:nspj+1
         td=abs(tli(i-1)-tlj(j-1));
         li=squeeze(lc(i-1,j,:));
         lj=squeeze(lc(i,j-1,:));
         lij=td+[0;squeeze(lc(i-1,j-1,1:nij-1))];
         lc(i,j,:)=min([li,lj,lij],[],2);
   	end
	end
end
rmax=nij;
c=cat(3,zeros(nspi+1,nspj+1),lc);
clast=squeeze(c(end,end,:));
return
end

function dists=spkdallq_final(qlist,sa,sb,clast,rmax)
% dists=spkdallq_final(qlist,sa,sb,clast,rmax) does the 
% final portion of parallel DP spike tima algorithm, extended to simultaneous
%  calculation for all values of q.
%
%    dists: a column vector of the distances
%    qlist: a column vector of values of q. qklist(:,1): q-values
%    sa, sb: spike times on the two spike trains
%    clast, rmax: calculated by spkdallq_recur
%
% See spkdallqk.doc.
%
%  Copyright (c) by Jonathan Victor.
%
%    See also SPKDALLQ_RECUR, SPKDALLQ_DIST, SPKD, SPKDALLQ.
%
%
na=length(sa);
nb=length(sb);
nq=length(qlist);
qlist=reshape(qlist,[nq 1]);
%
%make column vector of na+nb-2r, indicating costs of deletions and insertions
%
nanbrs=(na+nb)*ones(1+rmax,1)-2*[0:rmax]';
%
% find the best strategy (all r (matched links) and all s (mismatched links)
%
clast=reshape(clast,[1+rmax 1]);
for iq=1:nq
   posscosts=qlist(iq,1)*clast+nanbrs;
   dists(iq,1)=min(min(posscosts));
end
return
end

function [dists,clast,rmax]=spkdallq_dist(qlist,sa,sb)
% function [dists,clast,rmax]=spkdallq_dist(qlist,sa,sb) does the 
% recursion and final portion of the parallel DP single unit algorithm, extended to simultaneous
%  calculation for all values of q.
%
%    qlist: a column vector of values of q. qklist(:,1): q-values.
%    sa, sb: spike times on the two spike trains
%    clast, rmax: calculated by spkdallq_recur.
%
% See spkdallqk.doc.
%
%  Copyright (c) by Jonathan Victor.
%
%    See also SPKDALLQ_RECUR, SPKDALLQ_DIST.
%
%do the recursion
[clast,rmax,c]=spkdallq_recur(sa,sb);
%do the final stage
dists=spkdallq_final(qlist,sa,sb,clast,rmax);
return
end