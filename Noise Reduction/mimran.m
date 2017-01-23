function [Myout]= mimran(signal,fs,IS)

% Myout=mimran(S,FS,IS)
% S is the noisy signal
% FS is the sampling frequency
% IS is the initial silence (noise only) length in seconds (default value is .25 sec)
%
% Muhammad Imran

if (nargin<3 || isstruct(IS))
    IS=0.25; %seconds
end
W = fix(0.025*fs); %Window length is 25 ms
nfft = W;
ShiftP = 0.4; % Shift percentage is 40% (10ms) %Overlap-Add method works good with this value(.4)
windowtype = hamming(W);

if (nargin >= 3 && isstruct(IS))
    W = IS.windowsize;
    ShiftP = IS.shiftsize/W;
    nfft = IS.nfft;
    windowtype = IS.window;
    if isfield(IS,'IS')
        IS = IS.IS;
    else
        IS = 0.25;
    end
end

NIS = fix((IS*fs-W)/(ShiftP*W)+1);
% Spectral Subtractiio Rule
% 1 = Magnitude spectral subtraction
% 2 = Power spectrum subtraction
Gamma= 1;
y = segment(signal,W,ShiftP,windowtype);
Y = fft(y,nfft);
YPhase = angle(Y(1:fix(end/2)+1,:));
Y = abs(Y(1:fix(end/2)+1,:)).^Gamma;
numberOfFrames = size(Y,2);
FreqResol = size(Y,1);
N = mean(Y(:,1:NIS)')';
NRM = zeros(size(N));
NoiseCounter = 0;
NoiseLength = 20;
Beta = 0.00001;
YS = Y; % Y Magnitude Averaged

for i = 2:(numberOfFrames-1)
    YS(:,i) = (Y(:,i-1)+Y(:,i)+Y(:,i+1))/3;
end
for i=1:numberOfFrames
    [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad(Y(:,i).^(1/Gamma),N.^(1/Gamma),NoiseCounter); %Magnitude Spectrum Distance VAD
    if SpeechFlag == 0
        N = (NoiseLength*N+Y(:,i))/(NoiseLength+1);
        NRM = max(NRM,YS(:,i)-N);
        X(:,i)=Beta*Y(:,i);
    else
        DD = YS(:,i)-N;
        if i>1 && i<numberOfFrames
            for j=1:length(DD)
                if DD(j)<NRM(j)
                    DD(j)=min([DD(j) YS(j,i-1)-N(j) YS(j,i+1)-N(j)]);
                end
            end
        end
        X(:,i)=max(DD,0);
    end
end
Myout=OverlapAdd2(X.^(1/Gamma),YPhase,W,ShiftP*W);

function ReconstructedSignal = OverlapAdd2(XNEW,yphase,windowLen,ShiftLen)

if nargin<2
    yphase=angle(XNEW);
end
if nargin<3
    windowLen=size(XNEW,1)*2;
end
if nargin<4
    ShiftLen=windowLen/2;
end
if fix(ShiftLen)~=ShiftLen
    ShiftLen=fix(ShiftLen);
    disp('The shift length have to be an integer as it is the number of samples.')
    disp(['shift length is fixed to ' num2str(ShiftLen)])
end

[FreqRes FrameNum]=size(XNEW);

Spec=XNEW.*exp(j*yphase);

if mod(windowLen,2) %if FreqResol is odd
    Spec=[Spec;flipud(conj(Spec(2:end,:)))];
else
    Spec=[Spec;flipud(conj(Spec(2:end-1,:)))];
end
sig=zeros((FrameNum-1)*ShiftLen+windowLen,1);
weight=sig;
for i=1:FrameNum
    start=(i-1)*ShiftLen+1;
    spec=Spec(:,i);
    sig(start:start+windowLen-1)=sig(start:start+windowLen-1)+real(ifft(spec,windowLen));
end
ReconstructedSignal=sig;

function [NoiseFlag, SpeechFlag, NoiseCounter, Dist] = vad(signal,noise,NoiseCounter,NoiseMargin,Hangover)
if nargin<4
    NoiseMargin=3;
end
if nargin<5
    Hangover=8;
end
if nargin<3
    NoiseCounter=0;
end
FreqResol = length(signal);
SpectralDist = 20*(log10(signal)-log10(noise));
SpectralDist(find(SpectralDist<0)) = 0;
Dist=mean(SpectralDist); 
if (Dist < NoiseMargin) 
    NoiseFlag=1; 
    NoiseCounter=NoiseCounter+1;
else
    NoiseFlag=0;
    NoiseCounter=0;
end
if (NoiseCounter > Hangover) 
    SpeechFlag=0;    
else 
    SpeechFlag=1; 
end 
function Seg = segment(signal,W,SP,Window)
if nargin<3
    SP=.4;
end
if nargin<2
    W=256;
end
if nargin<4
    Window=hamming(W);
end
Window=Window(:);
L=length(signal);
SP=fix(W.*SP);
N=fix((L-W)/SP +1);
Index=(repmat(1:W,N,1)+repmat((0:(N-1))'*SP,1,W))';
hw=repmat(Window,1,N);
Seg=signal(Index).*hw;
