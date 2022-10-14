%ESTIMATE_WAKE_PROB Estimates wake probability from behavioral and physiological data
%
%   This code implements the Baysian framework, empirical wake probability model, and particle filter
%   described in Prerau et. al 2014, PLOS Computational Biology. Using behavioral (binary responses)
%   and physiological (EMG and EEG) data, we can estimate the instantanous probability that a subject
%   is awake, which is equivalent in this model to the instantanous probability of response.
%
%   Usage:
%   estimate_wake_prob() [GENERATES EXAMPLE DATA AND RUNS DEMO]
%   parameter_estimates=estimate_wake_prob(Fs, data, num_particles, ploton, prog_bar)
%
%   Input:
%       Fs: sampling frequency (in Hz)
%       data: 5xT matrix of simutaneously observed data, with rows:
%              1. Behavioral response (1=correct, 0=incorrect, NaN=missing at that time point)
%              2. Squeeze amplitude (in mV)
%              3. Alpha power (in dB)
%              4. Delta power (in dB)
%              5. Theta power (in dB)
%           All rows are sampled at Fs, with missing data represented by NaN.
%           A missing observation type is represented by a complete row of NaN values.
%       num_particles: Number of particles to use (Default: 5000)
%       ploton: 1 = Plot output graph, 0 = No output plot (Default: 1);
%       progbar: 1 = Display progress bar, 0 = No progress bar (Default: 1);
%
%   Output:
%       estimates: A structure with the 2.5, 50, and 97.5 percentiles of Pr(Wake), observation, and state estimates
%
%   Example:
%       %Runs demo and saves example data to workspace
%       parameter_estimates=estimate_wake_prob();
%
%   From the paper:
%   "Tracking the Sleep Onset Process: An Empirical Model of Behavioral and Physiological Dynamics"
%   Prerau, MJ Hartnack, KE, Obregón-Henao, G, Sampson, A, Merlino, M,
%   Gannon, K, Bianchi, MT, Ellenbogen, JM, Purdon, PL
%   PLOS Computational Biology, 2014
%
%   Copyright 2014 The General Hospital Corporation, authored by Michael J. Prerau, Ph.D.
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%   Last modified 8/19/2014
%********************************************************************

function estimates=estimate_wake_prob(Fs, data, num_particles, ploton, prog_bar)
%Run example
if nargin==0
    estimates=run_example();
    return;
end

%Set up default inputs
if nargin<2
    error('Must input sampling frequency and data');
end

if size(data,1)~=5
    warning('Data must be a 5xT matrix with rows: ');
    warning('   1. Behavioral response');
    warning('   2. Squeeze amplitude');
    warning('   3. Alpha power');
    warning('   4. Delta power');
    warning('   5. Theta power');
    warning('Missing data is indicated as NaN.');
    error('Invalid data size');
end

if nargin<3
    num_particles=5000;
end

if nargin<4
    ploton=1;
end

if nargin<5
    prog_bar=1;
end

disp(['Estimating Pr(Wake) using ' num2str(num_particles) ' particles...']);

%Number of time points
N=length(data(1,:));

%Number of parameters to estimate
num_params=24;

%Create <time>x<variable>x<particle> matrix
particles=zeros(N+1,num_params,num_particles);

%-----------------------------------------------
%Set Priors and state variances
%-----------------------------------------------
%STATE PRIORS
%X0
particles(1,1:2,:)=rand(2,num_particles)*4;
particles(1,3,:)=rand(1,num_particles)*-4;

%Variance multiplier
mult_fact=Fs/4;

%Sig2x: State variance
particles(1,4,:)=rand(1,num_particles)*.05*mult_fact;
particles(1,5:6,:)=rand(2,num_particles)*.05*mult_fact;

%EMG PRIORS
%Compute the max and min motor power
data(2,:)=log(data(2,:));
motor_rest=prctile(data(2,data(1,:)==0),2.5);
motor_smin=prctile(data(2,data(1,:)==1),2.5);

%EEG PRIORS
%Compute the max and min alpha power
alphabounds=prctile(data(3,:),[2.5 97.5]);
alpha_min=alphabounds(1);
alpha_max=alphabounds(2);

%Compute the max and min delta power
delta_bounds=prctile(data(4,:),[2.5 97.5]);
delta_min=delta_bounds(1);
delta_max=delta_bounds(2);

%Compute the max and min delta power
theta_bounds=prctile(data(5,:),[2.5 97.5]);
theta_min=theta_bounds(1);
theta_max=theta_bounds(2);

%motor_min
particles(1,7,:)=randn(1,num_particles)*.5+motor_rest;
%motor_max
particles(1,8,:)=randn(1,num_particles)*.5+motor_smin;

%alpha_min
particles(1,9,:)=randn(1,num_particles)*.5+alpha_min;
%alpha_max
particles(1,10,:)=randn(1,num_particles)*.5+alpha_max;

%delta_min
particles(1,11,:)=randn(1,num_particles)*.5+delta_min;
%delta_max
particles(1,12,:)=randn(1,num_particles)*.5+delta_max;

%theta_min
particles(1,13,:)=randn(1,num_particles)*.5+theta_min;
%theta_max
particles(1,14,:)=randn(1,num_particles)*.5+theta_max;

%Sig2m/a/d/t: Observation variance
particles(1,15,:)=rand(1,num_particles)/25;
particles(1,16:18,:)=rand(3,num_particles)*10;

%Scale factors
particles(1,19:22,:)=rand(4,num_particles)*5;

%EMG SQUEEZE PRIORS
%Mu1
particles(1,23,:)=rand(1,num_particles)*1;
%Sig2mu: Observation variance
particles(1,24,:)=rand(1,num_particles)/10;

%Start the progress bar
if prog_bar
    wh=waitbar(0,'Estimating wake probability...');
end

%Keep state variance under control
gamma=.9999;

%Iterate through all time
for t=2:N+1
    %Amount the state/observation variances can change
    sig2v=.01;
    
    %-----------------------------------------------
    %Update using the one step prediction equation
    %-----------------------------------------------
    %Sig2x: State variance
    particles(t,4,:)=abs(gamma*particles(t-1,4,:)+randn(1,1,num_particles)*sig2v/5);
    particles(t,5:6,:)=abs(gamma*particles(t-1,5:6,:)+randn(1,2,num_particles)*sig2v);
    %Xt=Gamma*Xt-1 + Ex
    particles(t,1:3,:)=gamma*particles(t-1,1:3,:)+randn(1,3,num_particles).*particles(t,4:6,:);
    
    %P_min
    particles(t,[7 9 11 13],:)=particles(t-1,[7 9 11 13],:)+randn(1,4,num_particles)*sig2v;
    %P_max
    particles(t,[8 10 12 14],:)=particles(t-1,[8 10 12 14],:)+randn(1,4,num_particles)*sig2v;
    
    %Sig2a/d/t: Observation variance
    particles(t,15,:)=max(gamma*particles(t-1,15,:)+randn(1,1,num_particles).*sig2v/50 ,1e-100);
    particles(t,16:18,:)=max(gamma*particles(t-1,16:18,:)+randn(1,3,num_particles).*sig2v*10,1e-100);
    
    %P_scale
    particles(t,19:22,:)=abs(particles(t-1,19:22,:)+randn(1,4,num_particles).*sig2v);
    
    %Mu1
    particles(t,23,:)=particles(t-1,23,:)+randn(1,1,num_particles)*sig2v/2.5;
    %Sig2mu: Observation variance
    particles(t,24,:)=max(gamma.*particles(t-1,24,:)+randn(1,1,num_particles).*sig2v/5,1e-100);
    
    %----------------------------------------------------------
    %Use the observation model to compute the estimated state
    %----------------------------------------------------------
    %Get all the particles at current time t
    particles_t=squeeze(particles(t,:,:))';
    
    %Extract the variables from the particles
    x_motor=particles_t(:,1);
    x_alpha=particles_t(:,2);
    x_deltatheta=particles_t(:,3);
    
    motor_rest=particles_t(:,7);
    motor_smin=particles_t(:,8);
    
    alpha_min=particles_t(:,9);
    alpha_max=particles_t(:,10);
    
    delta_min=particles_t(:,11);
    delta_max=particles_t(:,12);
    
    theta_min=particles_t(:,13);
    theta_max=particles_t(:,14);
    
    sig2_alpha=particles_t(:,16);
    sig2_delta=particles_t(:,17);
    sig2_theta=particles_t(:,18);
    
    motor_scale=particles_t(:,19);
    alpha_scale=particles_t(:,20);
    delta_scale=particles_t(:,21);
    theta_scale=particles_t(:,22);
    
    mu1=particles_t(:,23);
    sig2_motor=particles_t(:,24);
    
    %Estimated motor power
    mu0=motor_rest+(motor_smin-motor_rest).*exp(motor_scale.*x_motor)./(1+exp(motor_scale.*x_motor));
    motor_hat=mu0+mu1.*x_motor;
    
    %Estimated alpha power
    alpha_hat=alpha_min+(alpha_max-alpha_min).*exp(alpha_scale.*x_alpha)./(1+exp(alpha_scale.*x_alpha));
    %Estimated delta power
    delta_hat=delta_min+(delta_max-delta_min).*exp(delta_scale.*x_deltatheta)./(1+exp(delta_scale.*x_deltatheta));
    %Estimated theta value
    theta_hat=theta_min+(theta_max-theta_min).*exp(theta_scale.*x_deltatheta)./(1+exp(theta_scale.*x_deltatheta));
    
    %Estimated binomial probability
    x_wake=(x_motor+x_alpha-x_deltatheta)/3;
    %Compute sleep probability
    p_wake=exp(x_wake)./(1+exp(x_wake));
    
    %-----------------------------------------------
    %Compute the likelihood/weights
    %-----------------------------------------------
    loglikelihood=0;
    if ~isnan(data(1,t-1))
        bin_observation=repmat(data(1,t-1),num_particles,1);
        loglikelihood=loglikelihood+bin_observation.*log(p_wake)+(1-bin_observation).*log(1-p_wake);
    end
    
    if ~isnan(data(2,t-1))
        cont_observation=repmat(data(2,t-1),num_particles,1);
        loglikelihood=loglikelihood-((cont_observation-motor_hat).^2./(2*sig2_motor));
    end
    
    if ~isnan(data(3,t-1))
        cont_observation=repmat(data(3,t-1),num_particles,1);
        loglikelihood=loglikelihood-((cont_observation-alpha_hat).^2./(2*sig2_alpha));
    end
    
    if ~isnan(data(4,t-1))
        cont_observation=repmat(data(4,t-1),num_particles,1);
        loglikelihood=loglikelihood-((cont_observation-delta_hat).^2./(2*sig2_delta));
    end
    
    if ~isnan(data(5,t-1))
        cont_observation=repmat(data(5,t-1),num_particles,1);
        loglikelihood=loglikelihood-((cont_observation-theta_hat).^2./(2*sig2_theta));
    end
    
    %Resample if there is actual data
    if any(~isnan(data(:,t-1)))
        %Compute the weights
        pweights=sum(loglikelihood,2);
        
        %-----------------------------------------------
        %Resample Particles
        %-----------------------------------------------
        %Get distribution of weights
        weights=exp(pweights-max(pweights));
        weights(isnan(weights))=0;
        
        %Sample particles given the distribution of the weights
        [~,ind]=randsampleind(squeeze(particles(t,1,:)),num_particles,weights);
        particles(t,:,:)=squeeze(particles(t,:,ind));
    end
    
    %Show progress if wanted
    if prog_bar && ~mod(round((t-1)/N*100),5)
        waitbar((t-1)/N, wh, ['Estimating Pr(Wake): ' num2str(round(100*(t-1)/N)) '% Complete']);
    end
end
close(wh);

%----------------------------
%Compute the parameter estimates
%----------------------------
%Remove the prior
particles=particles(2:end,:,:);

%Extract the variables from the particles
x_motor=squeeze(particles(:,1,:))';
x_alpha=squeeze(particles(:,2,:))';
x_deltatheta=squeeze(particles(:,3,:))';

motor_rest=squeeze(particles(:,7,:))';
motor_smin=squeeze(particles(:,8,:))';

alpha_min=squeeze(particles(:,9,:))';
alpha_max=squeeze(particles(:,10,:))';

delta_min=squeeze(particles(:,11,:))';
delta_max=squeeze(particles(:,12,:))';

theta_min=squeeze(particles(:,13,:))';
theta_max=squeeze(particles(:,14,:))';

motor_scale=squeeze(particles(:,19,:))';
alpha_scale=squeeze(particles(:,20,:))';
delta_scale=squeeze(particles(:,21,:))';
theta_scale=squeeze(particles(:,22,:))';

mu1=squeeze(particles(:,23,:))';

%Estimated motor power
mu0=motor_rest+(motor_smin-motor_rest).*exp(motor_scale.*x_motor)./(1+exp(motor_scale.*x_motor));

%Estimated squeeze value
motor_hat=prctile(exp(mu0+mu1.*x_motor),[2.5, 50 97.5]);
alpha_hat=prctile(alpha_min+(alpha_max-alpha_min).*exp(alpha_scale.*x_alpha)./(1+exp(alpha_scale.*x_alpha)),[2.5, 50 97.5]);
delta_hat=prctile(delta_min+(delta_max-delta_min).*exp(delta_scale.*x_deltatheta)./(1+exp(delta_scale.*x_deltatheta)),[2.5, 50 97.5]);
theta_hat=prctile(theta_min+(theta_max-theta_min).*exp(theta_scale.*x_deltatheta)./(1+exp(theta_scale.*x_deltatheta)),[2.5, 50 97.5]);

%Compute the wake state
x_wake=(x_motor+x_alpha-x_deltatheta)/3;
p_wake=prctile(exp(x_wake)./(1+exp(x_wake)),[2.5, 50 97.5]);
x_wake_hat=prctile(x_wake,[2.5, 50 97.5]);
x_motor_hat=prctile(x_motor,[2.5, 50 97.5]);
x_alpha_hat=prctile(x_alpha,[2.5, 50 97.5]);
x_deltatheta_hat=prctile(x_deltatheta,[2.5, 50 97.5]);

%Prepare output
if nargout>0
    estimates.pr_wake=p_wake;
    estimates.observations.motor=motor_hat;
    estimates.observations.delta=delta_hat;
    estimates.observations.theta=theta_hat;
    estimates.observations.alpha=alpha_hat;
    estimates.states.x_wake=x_wake_hat;
    estimates.states.x_motor=x_motor_hat;
    estimates.states.x_alpha=x_alpha_hat;
    estimates.states.x_deltatheta=x_deltatheta_hat;
end

if ploton
    times=(1:size(data,2))/Fs/60;
    
    fig_h = figure('units', 'inches','papertype','usletter','paperorientation','portrait',...
        'units','normalized','position',[0 0 1 1],'color','w');
    % Create axes
    ax(1)=axes('Parent',fig_h,'Position',[0.04 0.8268888 0.9333335 0.1231112]);
    
    % Create axes
    ax(2)=axes('Parent',fig_h,'Position',[0.04 0.6326666 0.9333335 0.1231112]);
    
    % Create axes
    ax(3)=axes('Parent',fig_h,'Position',[0.04 0.4384444 0.9333335 0.1231112]);
    
    % Create axes
    ax(4)=axes('Parent',fig_h,...
        'Position',[0.04 0.31104033970276 0.9333335 0.05629306029724]);
    
    % Create axes
    ax(5)=axes('Parent',fig_h,'Position',[0.04 0.05 0.9333335 0.206900212314225]);
    
    linkaxes(ax,'x');
    set(get(gcf,'children'),'units','normalized');
    
    %Motor Data
    subplot(ax(1))
    hold on;
    shadebounds(times, motor_hat(2,:), motor_hat(3,:), motor_hat(1,:),'r',[1 .6 .6],'none');
    plot(times, exp(data(2,:)),'k.','markersize',10)
    axis tight
    title('EMG Squeeze Amplitude','fontname','helvetica');
    ylabel('Amplitude (mV)');
    
    %alpha axis
    subplot(ax(2))
    hold on
    shadebounds(times, alpha_hat(2,:), alpha_hat(3,:), alpha_hat(1,:),'b','c','none');
    plot(times, data(3,:),'k-')
    axis tight
    title('Alpha Power','fontname','helvetica');
    ylabel('Power (dB)');
    
    %Delta Theta axis
    subplot(ax(3))
    hold on
    shadebounds(times, theta_hat(2,:), theta_hat(3,:), theta_hat(1,:),'m',[1 .6 1 ],'none');
    plot(times, data(5,:),'k-')
    
    shadebounds(times, delta_hat(2,:), delta_hat(3,:), delta_hat(1,:),[0 .6 0],[.6 1 .6],'none');
    plot(times, data(4,:),'k-')
    
    axis tight
    title('Theta and Delta Power','fontname','helvetica');
    ylabel('Power (dB)');
    
    yl=[get(ax(2),'ylim') get(ax(3),'ylim')];
    set(ax(2:3),'ylim',[min(yl) max(yl)]);
    
    subplot(ax(4))
    binplot(times(~isnan(data(1,:))),logical(data(1,~isnan(data(1,:)))));
    set(gca,'ytick',[-.5 .5],'yticklabel',{'Incorr.','Correct'});
    axis tight
    title('Behavioral Responses','fontname','helvetica');
    
    subplot(ax(5))
    hold on
    shadebounds(times, p_wake(2,:), p_wake(3,:), p_wake(1,:),'b',[.6 .6 1],'none');
    ylabel('Pr(Wake), Pr(Response)');
    axis tight
    ylim([0 1]);
    title('Wake Probability Curve','fontname','helvetica');
    
    xlabel('Time (min)');
end

%Plot binary behavioral responses
function binplot(varargin)
times=varargin{1};
N=logical(varargin{2});

oinds=N;
xinds=~N;

hold on
stem(times(oinds),ones(1,sum(oinds)),'color', [0 .5 0],'marker','none')
stem(times(xinds),-ones(1,sum(xinds)),'color', [1 0 0],'marker','none')

ylim([-1.1 1.1]);
set(gca,'ytick',[-.5 .5],'yticklabel',{'Incorrect','Correct'});

%Plots shaded confidence bounds
function handle = shadebounds(x, mid, hi, lo, cmid, cbounds, edgecolor)

%Set defaults
if nargin==4
    cmid=[0 0 0];
    cbounds=[.95 .95 .95];
    edgecolor=[.9 .9 .9];
elseif nargin==5
    cbounds=[.95 .95 .95];
    edgecolor=[.9 .9 .9];
elseif nargin==6
    edgecolor=[.9 .9 .9];
end

%Make sure all vectors are pointing in the same direction
x=x(:)';
mid=mid(:)';
lo=lo(:)';
hi=hi(:)';

%Plot curve and bounds
hold on
handle = fill([x fliplr(x)],[lo fliplr(hi)],cbounds,'edgecolor',edgecolor);
plot(x,mid,'linewidth',2,'color',cmid);

%Run the demo
function estimates=run_example()
clc;

disp('Generating Data...');
%Simulate data
Fs=4;
T=Fs*60*10;
x_motor=zeros(1,T);
x_alpha=zeros(1,T);
x_deltatheta=zeros(1,T);
x_motor(1)=4;
x_alpha(1)=4;
x_deltatheta(1)=-4;

%Generate state processes
for i=2:T
    if i>T/3
        shift=.03;
    else
        shift=.0005;
    end
    
    x_motor(i)=x_motor(i-1)+randn*.05-shift;
    x_alpha(i)=x_alpha(i-1)+randn*.05-shift;
    x_deltatheta(i)=x_deltatheta(i-1)+randn*.05+shift;
end

x_wake=mean([x_motor;x_alpha;-x_deltatheta]);

pr_response=exp(x_wake)./(1+exp(x_wake));

response_inds=1:Fs:T;
bin=nan(1,T);
bin(response_inds)=rand(1,length(response_inds))<pr_response(response_inds);

motor=nan(1,T);
motor_rest=log(1); motor_smin=log(10); motor_scale=.0003; mu1=.4;
mu0=motor_rest+(motor_smin-motor_rest).*exp(motor_scale.*x_motor)./(1+exp(motor_scale.*x_motor))+randn(1,T)*.2;
motor(response_inds)=exp(mu0(response_inds)+mu1.*x_motor(response_inds));

alpha_min=-5; alpha_max=5;alpha_scale=5;
delta_min=-4; delta_max=5;delta_scale=5;
theta_min=-5; theta_max=2;theta_scale=5;

alpha=smooth(alpha_min+(alpha_max-alpha_min).*exp(alpha_scale.*x_alpha)./(1+exp(alpha_scale.*x_alpha))+randn(1,T)*5,10)';
delta=smooth(delta_min+(delta_max-delta_min).*exp(delta_scale.*x_deltatheta)./(1+exp(delta_scale.*x_deltatheta))+randn(1,T)*5,10)';
theta=smooth(theta_min+(theta_max-theta_min).*exp(theta_scale.*x_deltatheta)./(1+exp(theta_scale.*x_deltatheta))+randn(1,T)*5,10)';

%Create data matrix
data=[bin; motor; alpha; delta; theta];
disp(' ');
disp('> parameter_estimates=estimate_wake_prob(Fs, data, 1000);');
%Estimate data
estimates=estimate_wake_prob(Fs, data, 1000);

disp(' ');
assignin('base','example_data',data);
assignin('base','Fs',Fs);
assignin('base','parameter_estimates',estimates);
disp('Example data and parameter estimates saved to workspace.');


%RANDSAMPLEIND Take weighted samples with replacement and return indices.
% adapted from randomsample() by The MathWorks, Inc.
function [samples,  ind]= randsampleind(data, num_particles, weights)
emp_pdf = weights(:)' / sum(weights);

bin_edges = min([0 cumsum(emp_pdf)],1);
bin_edges(end) = 1;

[~, ind] = histc(rand(num_particles,1),bin_edges);

samples = data(ind);

