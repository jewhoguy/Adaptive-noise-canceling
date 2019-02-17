clear all; close all;
% %%
%
[c,fs] = audioread('clear_speech.wav'); %clear speech and take the sampling frequency
playClear = audioplayer(c,fs); %Save for listening
v = audioread('noise_source.wav'); %noise
s1 = audioread('speech_and_noise_through_room_1.wav'); %room 1
s2 = audioread('speech_and_noise_through_room_2.wav'); %room 2
playS1 = audioplayer(s1,fs);
playS2 = audioplayer(s2,fs);
%
v_s = audioread('structured_noise_source.wav'); %structured noise
s1_s = audioread('speech_and_structured_noise_through_room_1.wav'); %room 1 structured noise
playS1_s = audioplayer(s1_s,fs);
s2_s = audioread('speech_and_structured_noise_through_room_2.wav'); %room 2 structured noise
playS2_s = audioplayer(s2_s,fs);

% After code has run, use
% play(playNLMS,fs)
% play(playRLS,fs)
% play(playLS,fs)

% % % % %

%% % % % % % % % % % %
% % % % NLMS % % % %
% % % % % % % % % % %
tic
% [c,fs] = audioread('clear_speech.wav'); %clear speech and take the sampling frequency
% v = audioread('noise_source.wav'); %noise
% s2 = audioread('speech_and_noise_through_room_2.wav'); %room 2
h = waitbar(0,'NLMS');
% Initialization
mu = [0.05:0.05:1.95]; %µ - step-size/adaptation constant
M = 200; %number of filter taps
a = 1; %normalization variable for NLMS
AvgSE_NLMS = zeros(length(mu),1);

%for every µ
for n = 1:length(mu)
    [e_NLMS(:,n),w_NLMS,W_NLMS] = NLMS(s1,v,mu(n),M,a);
    AvgSE_NLMS(n) = ASE(c,e_NLMS(:,n)); %Calculate Average Squared Error for each e_NLMS
    waitbar(n/length(mu))
end
close(h)
% Find the best µ
mumin = min(AvgSE_NLMS); %µ with least error
muminInd = find(AvgSE_NLMS==min(AvgSE_NLMS)); %Index of least error
playNLMS = audioplayer(e_NLMS(:,muminInd),fs); %Save recovered signal with least error for audioplaying purposes

toc
%% % % % % % % % % %
% % % % RLS % % % %
% % % % % % % % % %
tic

[c,fs] = audioread('clear_speech.wav'); %clear speech and take the sampling frequency
v_s = audioread('noise_source.wav'); %noise
s1_s = audioread('speech_and_noise_through_room_2.wav'); %room 2

h = waitbar(0,'RLS');
% Initialization
lambda = 0.999:0.0001:1; %forgetting factor
delta = 100; %init param for P
n = 0; %reset iterator (just in case)
AvgSE_RLS = zeros(length(lambda),1);
e_RLS = [];%zeros(length(lambda),length(d));

for n = 1:length(lambda)
    [e_RLS(:,n),w_RLS] = RLS(s1_s,v_s,M,lambda(n),delta);
    AvgSE_RLS(n) = ASE(c,e_RLS(:,n));
    waitbar(n/length(lambda))
end
close(h);
% Find best forgetting factor
lambdamin = min(AvgSE_RLS)
lambdaminInd = find(AvgSE_RLS==min(AvgSE_RLS));
playRLS = audioplayer(e_RLS(:,lambdaminInd),fs);

refSig = ones(length(mu),1); %initialize reference signal
refSig = refSig*ASE(c,s1_s); %copy the average square error for the length of mu

toc
%%

% % % % % % % % % % % % % % % % % %
% % % % Segmentation and LS % % % %
% % % % % % % % % % % % % % % % % %
tic
% We segment our signal into 10 separate intervals (of equal length)
% We know that our signals are of length 32000, however to make the code to
% be more general we're gonna attempt to divide any signal into 10 equal segments
% (Note: signal segmentation is not based in time, but rather in
% datapoints)

[c,fs] = audioread('clear_speech.wav'); %clear speech and take the sampling frequency
s2 = audioread('speech_and_noise_through_room_2.wav'); %room 2
v = audioread('noise_source.wav'); %noise


M = 10; %number of segments
N_max = length(c); %"32000"
N = N_max/M; %"3200"
% A = [];

vbuffer = buffer(v,M); %segment the signal into 10 equal length vectors
% vbuffer = vbuffer';
dbuffer = buffer(s1,M);

for k = 1:M
    A = zeros(N,k); %initialize A
    for n = 1:M
        A(n:end,n) = vbuffer(k,1:end-n+1); %create matrix A
    end
    d = dbuffer(k,:)';
    w = A\d;
    hatd = A*w;
    e_LS(k,:) = d-hatd;
end
e_LS = reshape(e_LS,[N_max,1]);

AvgSE_LS = ASE(c,e_LS)
playLS = audioplayer(e_LS,fs);
toc
%%
tic
% % % % % % % % % % % % % % % % % %
% % % % % Order Selection % % % % %
% % % % % % % % % % % % % % % % % %

% s1_s = audioread('speech_and_noise_through_room_2.wav'); %room 2
% v_s = audioread('noise_source.wav'); %noise

NrOfTrials = 250;
N = length(s1_s);
N1 = N/10; %define amount of N (small or large)
N2 = N-N1+1;
nK = 10; %order count

for i = 1:NrOfTrials
    for k = 1:nK
        A = zeros(N1,k);
        for ij = 1:k
            A(:,ij) = v_s((N2-ij):(N-ij));
        end
        d = v_s(N2:N);
        w = A\d;
        hatd = A*w;
        e = d - hatd;
        
        %From lecture 8, slide 28
        FPE(k) = ((N1+k)/(N1-k))*var(e);
        AIC(k) = N1*log(var(e)) + 2*k;
        MDL(k) = N1*log(var(e)) + k*log(N1);
    end
    [~,FPE_order(i)] = min(FPE);
    [~,AIC_order(i)] = min(AIC);
    [~,MDL_order(i)] = min(MDL);
end

toc

%% Plotting

figure(1);
subplot(3,1,1); plot(c,'b');title('Clear speech signal');axis([0 length(c) -0.2 0.2]);
subplot(3,1,2); plot(s1,'r');title('Input signal (with noise)');axis([0 length(c) -0.2 0.2]);
subplot(3,1,3); plot(e_NLMS,'g');title('Recovered signal using NLMS');axis([0 length(c) -0.2 0.2]);

figure(2);
subplot(3,1,1); plot(c,'b');title('Clear speech signal');axis([0 length(c) -0.2 0.2]);
subplot(3,1,2); plot(s1,'r');title('Input signal (with noise)');axis([0 length(c) -0.2 0.2]);
subplot(3,1,3); plot(e_RLS,'g');title('Recovered signal using RLS');axis([0 length(c) -0.2 0.2]);
% Notes: notice by looking at plot (and listening) NLMS takes time to
% converge and actually remove the noise, whereas RLS does a good job from
% the beginning

figure(3);plot(mu,AvgSE_NLMS,'b'); legend('NLMS residual');
xlabel('µ - step-size');ylabel('Average Square Error'); hold on; 
plot(mu,refSig,'r'); legend('NLMS residual','e=s'),hold off;

figure(4);plot(lambda,AvgSE_RLS,'b'); 
xlabel('lambda - forgetting factor'); ylabel('Average Square Error RLS');

figure(5);
subplot(311);hist(FPE_order,(1:nK));
title('FPE criterion');xlabel('order k^*');ylabel('Trials');

subplot(312);hist(AIC_order,(1:nK));
title('AIC criterion');xlabel('order k^*');ylabel('Trials');

subplot(313);hist(MDL_order,(1:nK));
title('MDL criterion');xlabel('order k^*');ylabel('Trials');
