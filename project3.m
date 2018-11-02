%% Load data
load('data_cn_project_iii_a17.mat')

%% Question 1
% Auto correlation
tic
R = @(tou)(dot(Stimulus(51:end-50),Stimulus(51+tou:end-50+tou))/(length(Stimulus)-100) );

Tou = -50:1:50;
Auto = ones(size(Tou));

for i=1:101
    Auto(i) = R(Tou(i));
end
toc
figure;

xlabel('time')
ylabel('Auto(t)')

plot(Tou,Auto);
%disp(Auto)

%% Question 2
% PSTH - Trail averaged spike rate
tic
PSTH = zeros(4,20000);


for i=1:4
    v=zeros(1,20000);
    for j=1:50
        for k=All_Spike_Times{i,j}         
            in = int64(k*1000);            
            %if in<=15000
                v(in) = v(in) + 1;
            %else
            %    break;
            %end
        end
    end
    PSTH(i,:) = v/50;
end

PSTH = PSTH*100;

figure;
xlabel('time')
ylabel('10*spikes/sec')

subplot(2,2,1);
plot(1:20000,PSTH(1,:));
title('Neuron 1')

subplot(2,2,2);
plot(1:20000,PSTH(2,:));
title('Neuron 2')

subplot(2,2,3);
plot(1:20000,PSTH(3,:));
title('Neuron 3')

subplot(2,2,4);
plot(1:20000,PSTH(4,:));
title('Neuron 4')
toc
%% Question 3
tic
valueset = {zeros(4,50,2000),zeros(4,50,1000),zeros(4,50,400),zeros(4,50,200),zeros(4,50,100),zeros(4,50,40)};
keyset = [10,20,50,100,200,500];
rateStore = containers.Map(keyset,valueset);

valueset2 = {zeros(4,2000),zeros(4,1000),zeros(4,400),zeros(4,200),zeros(4,100),zeros(4,40)};
meanStore = containers.Map(keyset,valueset2);
varStore = containers.Map(keyset,valueset2);

for i=keyset
    v1 = zeros(4,50,ceil(20000/i));
    for j=1:4
        v2 = zeros(1,50,ceil(20000/i));
        for k=1:50
            v3 = zeros(1,1,ceil(20000/i)); 
            for l=All_Spike_Times{j,k}         
                in = ceil(l*1000/(i));                            
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
    title(strcat('neuron',num2str(j)))
    for i = keyset
        index=index+1;
        m = meanStore(i);
        v = varStore(i);
        m1 = m(j,1,:);
        v1 = v(j,1,:);
        
        subplot(2,3,index);
        plot(m1(:),v1(:),'.')
        title(strcat('binsize ',num2str(i)))
        xlabel('mean')
        ylabel('variance')
    end
end
toc
%%
%  spiketrain = All_Spike_Times{1,1};
%         indices = find(spiketrain*1000>100);
%         filtertrain = spiketrain(indices)*1000;
%         disp(Stimulus(int64(filtertrain)))
%% Question 4
tic
STA = zeros(4,100);

for i=1:4
    for j=1:50
        spiketrain = All_Spike_Times{i,j};
        indices = find(spiketrain*1000>100);
        filtertrain = spiketrain(indices)*1000;
        for k=1:100
            STA(i,k) = STA(i,k) + sum(Stimulus(int64(filtertrain-k)))/length(filtertrain);
        end
        
    end
end
STA = STA/50;


figure;
for i=1:4
    subplot(2,2,i);
    
    plot(1:100,STA(i,:))
    title(strcat('neuron',num2str(i)))
    xlabel('time')
    ylabel('STA')
end
toc
%%
figure
subplot(2,1,1)
y = conv(Stimulus,STA(1,:),'same') ;
plot(Stimulus(1:200));
subplot(2,1,2)
plot(y(1:200));

figure
subplot(2,1,1)
y = conv(Stimulus,STA(2,:),'same') ;
plot(Stimulus(1:200));
subplot(2,1,2)
plot(y(1:200));

figure
subplot(2,1,1)
y = conv(Stimulus,STA(3,:),'same') ;
plot(Stimulus(1:200));
subplot(2,1,2)
plot(y(1:200));

figure
subplot(2,1,1)
y = conv(Stimulus,STA(4,:),'same') ;
plot(Stimulus(1:200));
subplot(2,1,2)
plot(y(1:200));


%% Question 5
tic

params = [];
for i=1:4
   
    binsize=100;
    y = conv(Stimulus,STA(i,:),'same');    
    [N,edges,bins] = histcounts(y,binsize);
    %disp(bins)
    
    %disp(PSTH(i,bins==4))
    available_data = [];
    bin_means = [];
    for n = 1:binsize
        
       if ~isempty(PSTH(i,bins==n))
        
        bin_means = [bin_means mean(PSTH(i,bins==n))];
        available_data = [available_data edges(n)];
      end
      %disp(bin_means(n))
    end
    %Xedges = edges(1,2:end);
    subplot(2,2,i);
    params = [params sigm_fit(available_data,bin_means)];
    hold on;
    plot(available_data,bin_means,'.');
end
toc
%% Question 6




%% PART B
tic
costs = [0,0.001,0.01,0.1,1,10,100];
qvalues = length(costs);
s = RandStream('mlfg6331_64');

processedspikes = cell(4,50,200);

for i=1:4
%    tempCell = cell(1,50,200);
    for j=1:50
        temparray = [];
        memory = 1;
        for k=All_Spike_Times{i,j}
            ind = ceil(k*1000/100);
            %spiketrain = tempCell{1,j,ind};
            %if isempty(tempCell{1,j,ind})
            %    tempCell{1,j,ind} = [k*1000];   
            %else            
            %    tempCell{1,j,ind} = [spiketrain k*1000];
            %end
            if ind == memory
                temparray = [temparray k];            
            else
                processedspikes{i,j,memory} = temparray;
                memory = ind;
                temparray = [];
            end
           
        end  
        processedspikes{i,j,memory} = temparray;
    end
    %processedspikes(i,:,:) = tempCell;
end

save( 'samplespikes.mat', 'processedspikes');
toc
%%
tic

confusionMatrixTotal = zeros(100,4,8,8,qvalues); 

parfor i=1:100
   processedspi = load('samplespikes.mat');
   processedspikes = processedspi.processedspikes;
    
    y = sort(datasample(s,1:200,8,'Replace',false));  
%     disp(y)
    confusionMatrix = zeros(4,8,8,qvalues);
    
    for j=1:4
        comparitions = ones(400,400,qvalues)*Inf;
        for k=1:400
                in1 = ceil(k/50);
                ST1 = processedspikes{j,k-(ceil(k/50)-1)*50,y(in1)};               
                for m=1:400 % compare with which                                    
                        if ~(k==m) 
                            if comparitions(k,m,:) == Inf
                                in2 = ceil(m/50);
                                value = spkd_qpara(ST1,processedspikes{j,m-(ceil(m/50)-1)*50,y(in2)},costs); 
%                                 if j==1 && k==52 && m==2
%                                     disp(ST1)
%                                     disp(processedspikes{j,m-(ceil(m/50)-1)*50,y(in2)})
%                                     disp(value)
%                                 end
                                comparitions(k,m,:) = value;
                                comparitions(m,k,:) = value;
%                                 if j==1 && k==2 && m==52
%                                     disp(comparitions(m,k,:));
%                                     disp(comparitions(k,m,:));
%                                 end
                            end
                            %comparitions((m-1)*50+n,(k-1)*50+l,:) = value;                                                 
                        end
                end
                [argvalue, index] = min(reshape(comparitions(k,:,:),[400, qvalues]),[],1);
%                 if k==51 && j==1
%                     disp(comparitions(k,:,3));
%                     disp(index(3))
%                 end
                index=ceil(index/50); 
                %if k == 157 && j==1
                %    disp(index)
                %end
                %confusionMatrix(j,in1,index(1,1,:),[1,2,3,4,5,6]) = confusionMatrix(j,in1,index(1,1,:),[1,2,3,4,5,6]) + 1;
                for q=1:qvalues
                    confusionMatrix(j,in1,index(q),q) = confusionMatrix(j,in1,index(q),q) + 1;                     
                end                       
        end
        %if j == 1
        %    for o=1:qvalues
        %        disp(reshape(confusionMatrix(j,:,:,o),[8,8]))
        %    end
        %end
    end
    confusionMatrixTotal(i,:,:,:,:) = confusionMatrix/50;
end

toc

%%
% 
% P12 = reshape(confusionMatrixTotal(1,1,:,:,3),[8,8]);
% P1 = reshape(sum(P12,1)/8,[8,1]) ;
% P2 = reshape(sum(P12,2)/8,[1,8]) ;
% disp(P1)
% disp(P2')
% disp(P12)
% 
% P12 = reshape(confusionMatrixTotal(1,2,:,:,3),[8,8]);
% P1 = reshape(sum(P12,1)/8,[8,1]) ;
% P2 = reshape(sum(P12,2)/8,[1,8]) ;
% disp(P1)
% disp(P2')
% disp(P12)
% 
% P12 = reshape(confusionMatrixTotal(1,3,:,:,3),[8,8]);
% P1 = reshape(sum(P12,1)/8,[8,1]) ;
% P2 = reshape(sum(P12,2)/8,[1,8]) ;
% disp(P1)
% disp(P2')
% disp(P12)

%% Calculate mutual information
tic
confusionMatrixNew = confusionMatrixTotal/8; 
MIData = zeros(100,4,qvalues);
parfor i=1:100
    for j=1:4
        temp = [];
        for k=1:qvalues
            P12 = reshape(confusionMatrixNew(i,j,:,:,k),[8,8]);
            P1 = reshape(sum(P12,1)/8,[8,1]) ;
            P2 = reshape(sum(P12,2)/8,[1,8]) ;    
%             if i==1 && j ==1 && k==1
%                 disp(P1)
%                 disp(P2')
%                 disp(P12)
%             end
            MI=0;
            
            for l=1:8
                for m=1:8
                    if ~ (P12(l,m) == 0)
                        MI = MI + P12(l,m)*log((P12(l,m))/(P1(m)*P2(l)));
                    end
                end
            end
            %Com(isnan(Com)) = 0;
            %disp(Com)
            
            %MI = sum(Com);            
            %MI = sum(MI);
%             disp(MI)
            temp = [temp MI];
        end
        MIData(i,j,:) = temp;
    end
end
toc
%%

%%

meanData = reshape(mean(MIData,1),[4,qvalues]);
stdData = reshape(std(MIData,1),[4,qvalues])/sqrt(100);
disp(meanData)


Ubound = meanData + 1.645*stdData;
Lbound = meanData - 1.645*stdData;

costsX = costs;
costsX(1) = 10^-9;
figure;
for j=1:4
    subplot(2,2,j);
    hold on;
    ylim([4.2 5.2])
    xlim([-0.1 1])
    
    xlabel('1/q')
    ylabel('Mutual Information')
    title('1/q vs Mutual Information')
    %xticks(sort(1./costsX))
    %for i=1:100
        %disp(squeeze(MIData{i,1}))

    plot((1./costsX),squeeze(meanData(j,:)))
    plot((1./costsX),squeeze(Ubound(j,:)))
    plot((1./costsX),squeeze(Lbound(j,:)))
    %legend('90% interval upper bound','mean','90% interval lower bound')
    %end

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

function [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
% Optimization of parameters of the sigmoid function
%
% Syntax:
%       [param]=sigm_fit(x,y)       
%
%       that is the same that
%       [param]=sigm_fit(x,y,[],[],[])     % no fixed_params, automatic initial_params
%
%       [param]=sigm_fit(x,y,fixed_params)        % automatic initial_params
%       [param]=sigm_fit(x,y,[],initial_params)   % use it when the estimation is poor
%       [param]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
%
% param = [min, max, x50, slope]
%
% if fixed_params=[NaN, NaN , NaN , NaN]        % or fixed_params=[]
% optimization of "min", "max", "x50" and "slope" (default)
%
% if fixed_params=[0, 1 , NaN , NaN]
% optimization of x50 and slope of a sigmoid of ranging from 0 to 1
%
%
% Additional information in the second output, STAT
% [param,stat]=sigm_fit(x,y,fixed_params,initial_params,plot_flag)
%
%
% Example:
% %% generate data vectors (x and y)
% fsigm = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)))
% param=[0 1 5 1];  % "min", "max", "x50", "slope"
% x=0:0.1:10;
% y=fsigm(param,x) + 0.1*randn(size(x));
%
% %% standard parameter estimation
% [estimated_params]=sigm_fit(x,y)
%
% %% parameter estimation with forced 0.5 fixed min
% [estimated_params]=sigm_fit(x,y,[0.5 NaN NaN NaN])
%
% %% parameter estimation without plotting
% [estimated_params]=sigm_fit(x,y,[],[],0)
%
%
% Doubts, bugs: rpavao@gmail.com
% Downloaded from http://www.mathworks.com/matlabcentral/fileexchange/42641-sigmoid-logistic-curve-fit

% warning off

x=x(:);
y=y(:);

if nargin<=1 %fail
    fprintf('');
    help sigm_fit
    return
end

automatic_initial_params=[quantile(y,0.05) quantile(y,0.95) NaN 1];
if sum(y==quantile(y,0.5))==0
    temp=x(y==quantile(y(2:end),0.5));    
else
    temp=x(y==quantile(y,0.5));
end
automatic_initial_params(3)=temp(1);

if nargin==2 %simplest valid input
    fixed_params=NaN(1,4);
    initial_params=automatic_initial_params;
    plot_flag=1;    
end
if nargin==3
    initial_params=automatic_initial_params;
    plot_flag=1;    
end
if nargin==4
    plot_flag=1;    
end

if exist('fixed_params','var')
    if isempty(fixed_params)
        fixed_params=NaN(1,4);
    end
end
if exist('initial_params','var')
    if isempty(initial_params)
        initial_params=automatic_initial_params;
    end
end
if exist('plot_flag','var')
    if isempty(plot_flag)
        plot_flag=1;
    end
end

%p(1)=min; p(2)=max-min; p(3)=x50; p(4)=slope como em Y=Bottom + (Top-Bottom)/(1+10^((LogEC50-X)*HillSlope))
%f = @(p,x) p(1) + (p(2)-p(1)) ./ (1 + 10.^((p(3)-x)*p(4)));

f_str='f = @(param,xval)';
free_param_count=0;
bool_vec=NaN(1,4);
for i=1:4;
    if isnan(fixed_params(i))
        free_param_count=free_param_count+1;
        f_str=[f_str ' param(' num2str(free_param_count) ')'];
        bool_vec(i)=1;
    else
        f_str=[f_str ' ' num2str(fixed_params(i))];
        bool_vec(i)=0;
    end
    if i==1; f_str=[f_str ' + (']; end
    if i==2;
        if isnan(fixed_params(1))            
            f_str=[f_str '-param(1) )./ (   1 + 10.^( (']; 
        else
            f_str=[f_str '-' num2str(fixed_params(1)) ')./ (1 + 10.^((']; 
        end
    end    
    if i==3; f_str=[f_str ' - xval ) *']; end
    if i==4; f_str=[f_str ' )   );']; end
end

eval(f_str)

[BETA,RESID,J,COVB,MSE] = nlinfit(x,y,f,initial_params(bool_vec==1));
stat.param=BETA';

% confidence interval of the parameters
stat.paramCI = nlparci(BETA,RESID,'Jacobian',J);

% confidence interval of the estimation
[stat.ypred,delta] = nlpredci(f,x,BETA,RESID,'Covar',COVB);
stat.ypredlowerCI = stat.ypred - delta;
stat.ypredupperCI = stat.ypred + delta;

% plot(x,y,'ko') % observed data
% hold on
% plot(x,ypred,'k','LineWidth',2)
% plot(x,[lower,upper],'r--','LineWidth',1.5)

free_param_count=0;
for i=1:4;
    if isnan(fixed_params(i))
        free_param_count=free_param_count+1;
        param(i)=BETA(free_param_count);
    else
        param(i)=fixed_params(i);
    end    
end
    
if plot_flag==1 
    x_vector=min(x):(max(x)-min(x))/100:max(x);
    plot(x,y,'k.',x_vector,f(param(isnan(fixed_params)),x_vector),'r-')
    xlim([min(x) max(x)])
end
end