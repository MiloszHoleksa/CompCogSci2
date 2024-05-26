clear; close all

% Analysis of VE and VAE biases from an experiment based on Triplets of Trials.
% Triplets were either 1) AV reset trial 2) AV trial (n-1) 3) AV trial (n)
% or   1) AV reset trial, 2) AV trial (n-1) 3) A trial (for VAE)

% we analyze the history dependence of VE for Triplets 1, the VE itself also for Triplets 2, and the VAE for Triplets 2

load('data_11.mat')

%  condMAt = [1 vispos, 2 soundpos, 3 DeltaAV(V-A), 4 modality (AV = 0, A = 1, V = 2, Reset = 4), 5 trialtype id with 1-16 integ and 17-32 recal]
%  6 resp  x,7 resp y,8 resp RT, 9resp whichButton, 10 % trileID, 11 SubId

D=1;
for s=1:length(Dataset{D})
  data = Dataset{D}{s};
  %------------------------------------------------------------
  % ventriloquist bias in AV trials
  % (2nd AV trial). Cond =[1:16]
  %------------------------------------------------------------
  
  biaseall=[];
  for c=1:16 % each stimulus configuration
    Seq_onset = find( (data(:,5)==c).*(data(:,4)==4));
    % take 2nd AV test trial
    biases=[];
    for l=1:length(Seq_onset) % each triplet
      % response bias in last trial
      resp_bias = data(Seq_onset(l)+2,6)-data(Seq_onset(l)+2,2); % resp - soundpos
      % deltaVA
      dva = data(Seq_onset(l)+2,3);
      % deltaVAprev
      dvaprev = data(Seq_onset(l)+1,3);
      tmp = [dva,dvaprev,resp_bias,data(Seq_onset(l)+2,6)];
      biaseall = cat(1,biaseall,tmp);
    end
  end
  
  % estimate slopes against DVA and DVAprev
  model = biaseall(:,[1:2]);
  model(:,3) = 1;
  [beta,Stats,stats_each] = ck_stat_glm(model,biaseall(:,3));
  VEbias{D}(s,:) = beta(1:2);
  
  %------------------------------------------------------------
  % recalibration bias (A trial)
  % Cond =[1:16]+16
  %------------------------------------------------------------
  biaseall=[];
  % compute mean response for each sound location as correction for vAE
  mresp=[];
  apos = unique(data(:,2)); apos = apos(~isnan(apos));
  for a=1:length(apos)
    j = find( (data(:,2)==apos(a)).*(data(:,5)>16).*(data(:,4)~=4));
    mresp(a) = mean(data(j,6));
  end
  
  for c=1:16 % each stimulus configuration
    Seq_onset = find( (data(:,5)==(c+16)).*(data(:,4)==4));
    % take 2nd AV test trial
    biases=[];
    for l=1:length(Seq_onset) % each triplet
      % response bias
      ai = find(apos==data(Seq_onset(l)+2,2));
      resp_bias = data(Seq_onset(l)+2,6)-mresp(ai);
      % deltaVAprev
      dvaprev = data(Seq_onset(l)+1,3);
      % response bias in AV trial prior to A trial
      resp_biasAV = data(Seq_onset(l)+1,6)-data(Seq_onset(l)+1,2); % resp - soundpos
      
      
      tmp = [dvaprev,resp_bias,resp_biasAV];
      biaseall = cat(1,biaseall,tmp);
    end
  end
  
  % estimate slopes against each predictor
  model = biaseall(:,[1]);
  model(:,2) = 1;
  [beta,Stats,stats_each] = ck_stat_glm(model,biaseall(:,2));
  VAEbias{D}(s,:) = beta(1);
  
  % VE in VAE trial estimate slopes against each predictor
  model = biaseall(:,[1]);
  model(:,2) = 1;
  [beta,Stats,stats_each] = ck_stat_glm(model,biaseall(:,3));
  VEbias2{D}(s,:) = beta(1);
  
  
end


%% ----------------------------------------------------
% display slopes

VElabel ={'a','b'};
VAElabel ={''};

figure(2);clf;
h  = subplot(1,2,1); hold on
line([0.5 4.5],[0 0],'LineWidth',2,'color',[0.5 0.5 0.5]);
tmp = cat(1,VEbias{1});
ckboxplotcompact(tmp(:,[2,1]),[1:2],1);
set(gca,'XTick',[1:2],'XTickLabel',VElabel,'YTick',[-0.2:0.2:1])
 ylabel('partial beta')
axis([0.5 2.5 -0.25 1])
xlabel('trial');
text(1.5,1.05,'VE','FontSize',14,'FontWeight','Bold')

% compute bootstrap p-values 
for t=1:2
  [Lo,Up,p]=ck_stat_bootCI_percentile(tmp(:,t),'median',0.01,8000,0);
  CIs_VE(t,:) =[median(tmp(:,t)) Lo Up];
end

fprintf('VE partial betas: \n');
for t=1:2
  fprintf('%1.2f [%1.2f, %1.2f]  \n',CIs_VE(t,1),CIs_VE(t,2),CIs_VE(t,3));
end
h = subplot(1,2,2); hold on
line([0.5 2],[0 0],'LineWidth',2,'color',[0.5 0.5 0.5]);
tmp = cat(1,VAEbias{1});

ckboxplotcompact(tmp,[1],1);
set(gca,'XTick',[],'YTick',[-0.2:0.2:1])
ylabel('beta')
axis([0.2 2.5 -0.25 1])
xlabel('');
text(1.2,1.05,'VAE','FontSize',14,'FontWeight','Bold')
set(gca,'Position',[0.58    0.15    0.2 0.7673])


[Lo,Up,p]=ck_stat_bootCI_percentile(tmp,'median',0.01,8000,0);
CIs_VAE =[median(tmp)  Lo Up];


fprintf('VAE partial beta: \n');
fprintf('%1.2f [%1.2f, %1.2f]  \n',CIs_VAE(1),CIs_VAE(2),CIs_VAE(3));

% VE effect in AV-A sequence
x = cat(1,VEbias2{1});
[Lo,Up,p]=ck_stat_bootCI_percentile(x,'median',0.01,8000,0);
CIs_VE2 =[median(x)  Lo Up];

fprintf('VE in AV-A sequence: \n');
fprintf('%1.2f [%1.2f, %1.2f]  \n',CIs_VE2(1),CIs_VE2(2),CIs_VE2(3));



