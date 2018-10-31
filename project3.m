%% Load data
load('data_cn_project_iii_a17.mat')

%% Question 1
% Auto correlation
R = @(tou)(dot(Stimulus(51:end-50),Stimulus(51+tou:end-50+tou))/(length(Stimulus)-100) );

Tou = -50:1:50;
Auto = ones(size(Tou));
tic
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
costs = [0,0.001,0.01,1,10,100];
s = RandStream('mlfg6331_64');

processedspikes = cell(4,50,201);

parfor i=1:4
    tempCell = cell(50,201);
    for j=1:50
        for k=All_Spike_Times{i,j}
            index = int64(k*1000/100)+1;
            spiketrain = tempCell{j,index};
            tempCell{j,index} = [spiketrain k];
        end        
    end
    processedspikes(i,:,:) = tempCell;
end

save( 'samplespikes.mat', 'processedspikes');
toc
%%
disp(processedspikes{1,1,1}(1));
processedspi = load('samplespikes.mat');
%%
tic
confusionMatrixTotal = zeros(100,4,6,8,8); 

processedspikes = processedspi.processedspikes;
%parfor i=1:100
    y = datasample(s,1:200,8,'Replace',false);  
    confusionMatrix = zeros(4,6,8,8);
    comparitions = zeros(400,400,6);
    for j=1:4
        perNeuron = processedspikes(j,1:end,1:end);
        for k=1:8
            for l=1:50
                ST1 = perNeuron{1,l,y(k)};
                for m=1:8 % compare with which 
                    for n=1:50
                        if ~(k==m && l==n)
                            value = spkd_qpara(ST1,perNeuron{1,n,y(m)},costs);                
                            comparitions((k-1)*50+l,(m-1)*50+n,:) = value;
                            %comparitions((m-1)*50+n,(k-1)*50+l,:) = value; 
                        end
                        if (k==m && l==n)
                            comparitions((k-1)*50+l,(m-1)*50+n,:) = 100000;
                        end
                    end
                end
                [argvalue, index] = min(comparitions((k-1)*50+l,:,:),[],2);
                index=ceil(index/50);
                for o=1:length(index)
                    confusionMatrix(j,o,k,index(1,1,o)) = confusionMatrix(j,o,k,index(1,1,o)) + 1;
                end
            end           
        end
    end
    confusionMatrixTotal(1,:,:,:,:) = confusionMatrix/50;
%end
toc
%%

arr = confusionMatrixTotal(1,1,1,1);
disp(arr)

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