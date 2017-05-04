function [TSNR,HRNR] = Wiener_KIST(ns,fs,IS)

% Input Parameters :
%   ns          Noisy speech
%   fs          Sampling frequency (in Hz)
%   IS          Initial Silence (or non-speech activity) Period (in number of samples)
%
% Output Parameters : enhanced speech
%   TSNR      enhanced speech with the Two-Step Noise Reduction method
%   HNRN      enhanced speech with the Harmonic Regeneration Noise Reduction method

%% ------- input noisy speech  --------
l = length(ns);
s=ns;
wl = fix(0.020*fs);    % window length is 20 ms
NFFT=2*wl;             % FFT size is twice the window length
hanwin = hanning(wl);
if (nargin<3 || isstruct(IS))
    IS=10*wl;
end
nsum = zeros(NFFT,1);
count = 0;
for m = 0:IS-wl
    nwin = ns(m+1:m+wl).*hanwin;
    nsum = nsum + abs(fft(nwin,NFFT)).^2;
    count = count + 1;
end
d= (nsum)/count;
SP = 0.25;
normFactor=1/SP;
overlap = fix((1-SP)*wl);
offset = wl - overlap;
max_m = fix((l-NFFT)/offset);
zvector = zeros(NFFT,1);
oldmag = zeros(NFFT,1);
news = zeros(l,1);
phasea=zeros(NFFT,max_m);
xmaga=zeros(NFFT,max_m);
tsnra=zeros(NFFT,max_m);
newmags=zeros(NFFT,max_m);
alpha = 0.999;

%% --------------- TSNR ---------------------
for m = 0:max_m
    begin = m*offset+1;
    iend = m*offset+wl;
    speech = ns(begin:iend);
    winy = hanwin.*speech;
    ffty = fft(winy,NFFT);
    phasey = angle(ffty);
    phasea(:,m+1)=phasey;
    magy = abs(ffty);
    xmaga(:,m+1)= magy;
    postsnr = ((magy.^2) ./ d)-1;
    postsnr=max(postsnr,0.1);
    eta = alpha * ( (oldmag.^2)./d ) + (1-alpha) * postsnr;
    newmag = (eta./(eta+1)).*  magy;
    tsnr = (newmag.^2) ./ d;
    Gtsnr = tsnr ./ (tsnr+1);         %gain of TSNR
    tsnra(:,m+1)=Gtsnr;
    Gtsnr = gaincontrol(Gtsnr,NFFT/2);
    newmag = Gtsnr .* magy;
    newmags(:,m+1) = newmag;     %for HRNR use
    ffty = newmag.*exp(i*phasey);
    oldmag = abs(newmag);
    news(begin:begin+NFFT-1) = news(begin:begin+NFFT-1) + real(ifft(ffty,NFFT))/normFactor;
end
TSNR=news;

%% --------------- HRNR -----------------------
newharm = max(TSNR,0);
news = zeros(l,1);
%
for m = 0:max_m
    begin = m*offset+1;
    iend = m*offset+wl;
    nharm = hanwin.*newharm(begin:iend);
    ffth = abs(fft(nharm,NFFT));          %perform fast fourier transform
    snrham= ( (tsnra(:,m+1)).*(abs(newmags(:,m+1)).^2) + (1-(tsnra(:,m+1))) .* (ffth.^2) ) ./d;
    newgain= (snrham./(snrham+1));
    newgain=gaincontrol(newgain,NFFT/2);
    newmag = newgain .*  xmaga(:,m+1);
    ffty = newmag.*exp(i*phasea(:,m+1));
    news(begin:begin+NFFT-1) = news(begin:begin+NFFT-1) + real(ifft(ffty,NFFT))/normFactor;
end;
%Output
HRNR=news;

function NewGain = gaincontrol(Gain,ConstraintInLength)
%
meanGain=mean(Gain.^2);
NFFT=length(Gain);
L2=ConstraintInLength;
win=hamming(L2);
ImpulseR=real(ifft(Gain));
ImpulseR2=[ImpulseR(1:L2/2).*win(1+L2/2:L2);zeros(NFFT-L2,1);ImpulseR(NFFT-L2/2+1:NFFT).*win(1:L2/2)];
NewGain=abs(fft(ImpulseR2,NFFT));
meanNewGain=mean(NewGain.^2);
NewGain=NewGain*sqrt(meanGain/meanNewGain);





