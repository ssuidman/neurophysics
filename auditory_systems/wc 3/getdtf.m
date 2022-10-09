function DTF = getdtf(HRTF,Fs)

% close all
% clearvars

%% Remove Directional Component
rHRTF	= rms(HRTF')';
DTF		= HRTF./rHRTF;

%% Magnitude
DTF = abs(DTF);


%% Remove echo
% DTF = remecho(DTF,Fs);

%% Smooth
% Span	= 5;
% Sweep	= locsmooth(Sweep,Span,Azimuth,Elevation);
DTF	= hsmooth(DTF,8);

%%
function hrtf = readhrtf(hrtffile)
% HRTF = READHRTF(HRTFFILE)
%
% Read hrtf-file and extract Time Series data in HRTF
% HRTF is a (k x m x n)-matrix, where
%   k = number of samples
%   m = number of trials/locations
%   n = number of microphones/channels
%
% See also READCSV, READDAT

% Copyright 2008
% Marc van Wanrooij

%% Initialization
% Check Inputs
if nargin <1
    hrtffile        = '';
end
% Check Files
hrtffile            = fcheckext(hrtffile,'.hrtf');
hrtffile            = fcheckexist(hrtffile,'.hrtf');
csvfile             = fcheckext(hrtffile,'.csv');

% Set Parameters/Constants
nbytesperfloat      = 4;

%% Read CSV-file
[exp,~,mlog]     = readcsv(csvfile);
% to determine number of trials

nrepeats			= exp(3);
ntrials             = exp(4);
% and number of channels
sel                 = ismember(mlog(:,5),[6 7]); % 6 and 7 are Inp1 and Inp2
nchan               = length(unique(mlog(sel,5)));

%% Read HRTF-file
fid             = fopen(hrtffile,'r');
% First determine number of samples
fseek(fid,0,'eof');
nbytes          = ftell(fid);
nsample         = nbytes./(nchan*ntrials*nrepeats*nbytesperfloat);
% Get HRTF data
frewind(fid)
if ispc
    hrtf         = fread(fid,inf,'float');
else %MAC
    hrtf         = fread(fid,inf,'float','l');
end
% Some Reshaping
hrtf             = reshape(hrtf,nsample,nchan,ntrials*nrepeats);
% hrtf             = permute(hrtf,[1 3 2]); % Costs memory
% And close
fclose(fid);

function [expinfo,chaninfo,cLog] = readcsv(logfile)
% Extract experimental parameters from a csv-file
%
% [EXPINFO,CHANINFO,LOG] = READCSV(CSVFILE)
%
% The EXPINFO-vector contains experimental information,
% extracted directly from the EXP-file, with column:
%               1 = 0
%               2 = Maximum Number of trials (as set in exp-file)
%               3 = Number of repeats
%               4 = Number of '=>' (trial indicator in exp-file)
%               5 = ITI start (ITI: Intertrial Interval)
%               6 = ITI stop
%               7 = Random type (0 = no, 1 = per set, 2 = all trials)
%               8 = Number of channels
%
% The CHANINFO-matrix contains data acquistion channel information,
% extracted directly from the CFG-file, with each row corresponding to a
% single channel and each column corresponding to the following parameter:
%               1 = 0
%               2 = Channel number
%               3 = Channel name
%               4 = Low-Pass frequency
%               5 = Sample rate
%               6 = Number of samples
%
% For example, to obrain the sample rate of channel 2:
%  >> chaninfo(2,5) 
%
%
% The LOG-matrix contains trial information
%
%   1) =  Trial number (1 = trial no. 1)
%   2) =  Stimulus number (e.g. Acq-Led-Led-Snd -> 1 2 3 4, but trialnr stays the same)
%   3) =  random ITI calculated
%   4) =  random ITI actual (is larger when for instance disk is too busy)
%   5) =  Modality of stimulus (0 = LED, 1 = SKY, 2 = SND1, 3 = SND2, 4 = Acquisition, 5 = Trg0, 6 = Input 1, 7 = Input 2)
%   6) =  Stimulus location x BOOG (degrees, -180..180)
%   7) =  Stimulus location y speaker (number)
%   8) =  Stimulus onset
%   9) =  Stimulus offset
%   10) = Stimulus intensity (for LED: 0 = lowest ... 7 = highest, for SND: 0 = lowest ... 100 = higest)
%   11) = Stimulus attribute (For LED: 0 = red, 1 = green; for SND: XXX (fragment of wav-name -> sndXXX.wav))
%   12) = bit (1 = ??? ... 8 = highest; For SND: XXX (fragment of wav-name-> sndXXX.wav))
%   13) = line number in EXP-file
%
%                       
%   Type of trial effects which columns are not zero. The columns used
%   are presented below:
%
%                      ITI ITI
%              Trl Stm cal act Typ  x   y  On  Off Int Atr Bit Edg
%               1   2   3   4   5   6   7   8   9   10  11  12  13
%       Led:    X   X   X   X   X   X   X   X   X   X   X       X   
%       Snd:    X   X   X   X   X   X   X   X       X   X       X
%       Acq:    X   X   X   X   X           X                   X
%       Trg0:   X   X   X   X   X           X   X           X   X
% 
%   'Cells' that are empty in the logfile are denoted with 'NaN' in the
%   cLog-matrix
%
%  See also LOADDAT

% Copyright 2006
% Author: Tomg Oct 2006
% Modified by: Marcw 2007

%% Initialization
if nargin<1
    logfile             = '';
end
chaninfopar             = 6;
logpar                  = 13;

%% check INPUT %%%
logfile                 = fcheckext(logfile,'.csv');
logfile                 = fcheckexist(logfile);
if isempty(logfile)
    disp('No Log File has been chosen (readcsv)');
    expinfo             = [];
    chaninfo            = [];
    cLog                = [];
    return
end

%% Open csv-file
fid                     = fopen(logfile);
firstLine               = fgetl(fid);

%% Get Experimental Info from first line
expinfo(1,:)            = sscanf(firstLine,'%f;')';
Nchan                   = expinfo(1,8);
Ntrials                 = expinfo(1,4);


%% Get Channel Information from next Nchannels lines
chaninfo                = NaN*ones(Nchan,chaninfopar);
ident                   = 0;
k                       = 0;
while ~ident
    A                   = fscanf(fid,'%g;',chaninfopar)';
    ident               = A(1);
    if ~ident
        k               = k+1;
        chaninfo(k,:)   = A;
    end
end
chaninfo                = chaninfo(1:k,:);
if Nchan~=k
%     disp('Uh-oh - Number of relevant channels does not correspond to number of channel configs');
end

%% load cLog-file body
frewind(fid);
count                   = 0;
cLog                    = NaN*ones(Ntrials,logpar);
while ~feof(fid)
    curLine             = fgetl(fid);
    firstCell           = sscanf(curLine,'%g;',1);
    seps                = [0 findstr(';',curLine) length(curLine)+1];
    if firstCell
        count           = count+1;
        for i           = 2:length(seps)
            icol        = i-1;
            curcell     = curLine(seps(i-1)+1:seps(i)-1);
            if sum(isletter(curcell))
                if strcmpi(strtrim(curcell),'led')
                    cLog(count,icol)    = 0;
                elseif strcmpi(strtrim(curcell),'sky')
                    cLog(count,icol)    = 1;
                elseif strcmpi(strtrim(curcell),'snd') || strcmpi(strtrim(curcell),'snd1') 
                    cLog(count,icol)    = 2;
                elseif strcmpi(strtrim(curcell),'snd2')
                    cLog(count,icol)    = 3;
                elseif strcmpi(strtrim(curcell),'acq')
                    cLog(count,icol)    = 4;
                elseif strcmpi(strtrim(curcell),'trg0')
                    cLog(count,icol)    = 5;
                elseif strcmpi(strtrim(curcell),'inp1')
                    cLog(count,icol)    = 6;
                elseif strcmpi(strtrim(curcell),'inp2')
                    cLog(count,icol)    = 7;
                elseif strcmpi(strtrim(curcell),'nan')
                    cLog(count,icol)    = NaN;
                end
            elseif isempty(str2double(curcell))
                %Skip empty cell at end of file
                warning('READCSV:emptycell','from readcsv')
                disp(['   Cell (' num2str(count) ',' num2str(icol) ') is empty ... skip'])
                disp(['   logfile= ''' logfile ''')'])
            else
                cLog(count,icol)        = str2double(curcell);
            end
        end
    end
end


%% close Log-file %%
fclose(fid);

function [AZ,EL] = fart2azel(X,Y,SKY)
% FART2AZEL Transfrom FART to double polar coordinates.
% [AZ,EL] = FART2AZEL(X,Y) transforms corresponding elements of data
%     stored in FART (Fast Auditory Rotating Target) coordinates X 
%    (Theta deg),Y (speakernumber) to double polar coordinates (azimuth AZ and 
%     elevation EL).  The arrays X and Y must be the same size (or either 
%     can be scalar). AZ and EL are returned in deg.
%
%   WHY COORDINATION TRANSFORMATION:
%       FART coordinates *does not* directly reflect Azimuth-Elevation
%       double-polar coordinate system. The elevation of targets is
%       given in speaker number, which should be transformed to deg
%       (each speaker always has the same elevation according to
%       the double polar system irrespective of FART rotation angle
%       Theta).
%       Azimuth AZ can be deduced by applying the formula:
%
%           AZ = arcsin( sin(Theta) * cos(EL) )
%
%       where Theta = FART rotation angle (deg), and EL = elevation angle
%       (deg).
%
%       Note that the double-polar coordinate system by itself does not
%       allow for any front-back distinction. Therefore, elevation is
%       defined as being absolutely larger than 90 deg when coordinates are
%       positioned in the rear hemifield.
%
% See also AZEL2CART, AZEL2POL, CART2AZEL
%
% MarcW 2007

%% Initialization
% Convert to vectors
[M,N]               = size(X);
X                   = X(:);
Y                   = Y(:);

%% Check for old bug (before 19 March 2007)
% Anything below -200 should be recalculated
% sel             = X<-200;
% X(sel)          = -X(sel)-300;

%% Define Front&Rear Speakers, Front&RearHemifields
selfront            = Y<100;                        % Front Speakers Numbers: 1-29
selrear             = Y>100;                        % Rear Speakers Numbers 101-129
selled30            = Y == 30;                      % ET-LED1 Number: 30
selled31            = Y == 31;                      % ET-LED2 Number: 31
selleds             = selled30 | selled31;          % Both ET-LED Numbers: 30&31
selrearhemi1        = (X>90 | X<-90) & selfront;    % Rear hemifield for front speakers
selrearhemi2        = (X<90 & X>-90) & selrear;   % Rear hemifield for rear speakers

%% Speaker Phi 2 Elevation
EL                  = Y;
EL(selfront)        = -55+(Y(selfront)-1)*5;
EL(selrear)         = -57.5+(Y(selrear)-100-1)*5;
EL(selleds)         = -2.5;

%% Theta 2 Azimuth

% Obtain theta coordinates for the E.T.-L.E.D.S.
X(selled30)         = X(selled30)-14;
X(selled31)         = X(selled31)+14;

% And convert
AZ                  = asind(sind(X).*cosd(EL));
% Correct for rear speakers
AZ(selrear)         = -AZ(selrear);

%% Define rear quadrant
% as having absolute Elevations larger than 90 degrees,
% i.e. the front speakers (1-31) will be positioned in the rear when FART has
% moved farther than 90 deg.
% Also, the rear speakers (101-129) will then be positioned in the front
EL(selrearhemi1)    = sign(EL(selrearhemi1)).*(180-abs(EL(selrearhemi1)));
sel                 = selrearhemi1 & EL == 0;
EL(sel)             = -180;
EL(selrearhemi2)    = sign(EL(selrearhemi2)).*(180-abs(EL(selrearhemi2)));

%% Convert back to MTX
AZ                  = reshape(AZ,M,N);
EL                  = reshape(EL,M,N);

function [Spec,Freq,Tijd] = sweep2spec(Sweep,Wav,chan,NFFT,NSweep,Fs)
% [SPEC,FREQ,TIME] = SWEEP2SPEC(SWEEP,WAV,CHAN,NFFT,NSWEEP,CHAN)
%
% See also READHRTF, GENSWEEP


%% Initialization
if nargin<3
	chan = 1;
end
if nargin<4
	NFFT    = 1024;
end
if nargin<5
	NSweep = 18;
end
if nargin<6
	Fs     = 48828.125;
end
nloc	= size(Sweep,3);
nchan	= size(Sweep,2);

%% FFT & Reshaping
nBegin  = ceil(20/1000*Fs);
Tijd    = NaN(NFFT,nloc);
Spec    = Tijd;
L       = NaN(1,nloc);

% for each location i
for i = 1:nloc
	if nchan>1
		d           = squeeze(Sweep(:,chan,i)); % obtain the measured Sweep data from channel chan
	else
		d           = squeeze(Sweep(:,i)); % obtain the measured Sweep data from channel chan
	end
% 	d           = d(101:end);               % remove the last 100 samples
	% Rough Alignment Method 1
	if length(d)==length(Wav)
		[c,lags]    = xcorr(d,Wav,'coeff');
	else
		[c,lags]    = xcorr(d,Wav,'none');
	end
	[~,indx]    = max(abs(c));
	lag         = lags(indx);
	
	% Aligmnent Method 2
	% If Sweep is too weak, the time delay might be judged incorrectly
	% Correct for this:
	if lag>250 || lag<1
		lag     = 210;
	end
	L(i)        = lag;
	nOnset      = lag+nBegin;
	indx        = nOnset + NFFT + (1:(NFFT*NSweep));
	data        = d(indx);
	data        = reshape(data, NFFT, NSweep);
	meandata    = mean(data,2)';
	meandata    = meandata - mean(meandata);
	Tijd(:,i)       = meandata;
	nfft            = 2^(nextpow2(length(meandata)));
	X               = fft(meandata,nfft);
	Spec(:,i)       = X;
end
% NumUniquePts	= ceil((NFFT+1)/2);
% Spec			= Spec(1:NumUniquePts,:);
Freq			= (0:(NFFT-1))*Fs/NFFT;

function Mag = remecho(Mag,Fs,grph)
% MAG = REMCHO(MAG,FS)
%
% Remove Echos from Magnitude response
%
% See also READHRTF, SWEEP2SPEC, FFTZERO, ALIGNIR, ALIGNIR2

%% Initialization
Flo                     = 400; % Low Frequecny
Fhi                     = 19000; % High Frequency
Fn                      = Fs/2; % Nyquist Frequency
dt                      = 1/Fs; %Time step in s
ndelay                  = 50; %original = 50
if nargin<3
    grph = 0;
end

%% Smooth
% Zero Magnitude response before Flo and after Fhi
Mag                     = fftzero(Mag, Flo/Fn, Fhi/Fn);
% Impuls respons
ir                      = ifft(Mag);
if grph
	t = 1:512;
	t = 1000*t/Fs;
	m = ir(1:512);
% 	m = 20*log10(m./mean(m));
	plot(t,m,'-','LineWidth',2);
	xlabel('Time (ms)');
	ylabel('Amplitude');
	ax = axis;
	xlim([-1 ax(2)]);
	horline;
	verline;
end
% Setting the delay to 1 ms (coarse aligning)
ir                      = alignir(ir, ndelay);
% Fine alignment
ir                      = alignir2 (ir);
% Suppress reflections beyond 3 msec
Tsup                    = 0.003; % s
nsup                    = round(Tsup/dt);
ir(nsup:end-nsup+2,:)   = 0;
% Back to the Future/Frequency Domain
Mag                     = abs(fft(ir));

function dB = hsmooth(H,Q)
% function hsmooth(H [,Q])
% Smooths the magnitude of a transfer function
% with a simple auditory filter.
% Computes the amplitude in dB for the transfer function H,
% smoothed with a Gaussian filter with constant Q.
% If H is an array, its columns are treated as transfer functions.
% Q defaults to 8.

% Notes:
%  1.  Smooths the squared magnitude
%  2.  Q = center_frequency/half-power_bandwidth
%  3.  Scaled so that input power = output power
%  4.  Tacitly assumes a periodic extension of the spectrum

if nargin < 1
  fprintf('Format: dB = hsmooth(H [,Q])\n');
  return;
end
if nargin < 2
  Q = 8;
end

mindB = -100;
min_mag = 10^(mindB/20);

colvec = 1;
numrows = size(H,1);
numcols = size(H,2);
if numrows == 1     % if necessary, convert to column vector
  colvec = 0;
  H = H';
  numrows = numcols;
  numcols = 1;
end

alpha = 2*sqrt(2*log(2));         % half-power scale factor
if rem(numrows,2) == 0
  evencase = 1;
	npts = numrows/2 + 1;           % number from DC to Nyquist
	noffset = npts - 3;
else
  evencase = 0;
	npts = (numrows+1)/2;           % number from DC to Nyquist
	noffset = npts - 2;
end

A = max(abs(H), min_mag);
A = A.*A;

if Q > 0
  tau = (-(npts-1):(npts-1));	
  Asmoothed = zeros(npts,numcols);
  Asmoothed(1,:) = A(1,:);
  for k=2:npts
    sigma = (k-1)/(alpha*Q);
    window = exp(-(0.5*tau.*tau/(sigma^2)));
		if evencase
		  window = window(2:(2*npts-1));
		end
    window = window/sum(window);
    Asmoothed(k,:) = window*cycshift(A,k+noffset);
  end
else
  Asmoothed = A(1:npts,:);
end

dB = 10*log10(Asmoothed);
if ~colvec
  dB = dB';
end

function y = cycshift(x,n)
% y = cycshift(x,n)
% shifts x cyclically down n places
% If x is a matrix, shifts columns 
% Copyright (C) 2001 The Regents of the University of California

if nargin < 2
  fprintf('Format: y = cycshift(x,n)\n');
  return;
end;

if n ~= fix(n)
  fprintf('Error in cycshift; shift not an integer.\n');
  return;
end;

col_vec = 1;         % Check and convert to column vector if necessary
if size(x,1) == 1
  col_vec = 0;
  x = x';
end

N = size(x,1);
n = rem(n,N);
if n == 0
  y = x;
elseif n > 0
  y(n+1:N,:) = x(1:N-n,:);
  y(1:n,:) = x(N-n+1:N,:);
else
  n = -n;
  y(1:N-n,:) = x(n+1:N,:);
  y(N-n+1:N,:) = x(1:n,:);
end

if ~col_vec
  y = y';
end

function fname = fcheckext(fname,fext)
% FNAME = PA_FCHECKEXT(FNAME,FEXT)
%
% Check whether extension of file FNAME corresponds to FEXT
% If not, the extension will be replaced or added.
%
% See also PA_FCHECKEXIST

% (c) 2011 Marc van Wanrooij

%% Check dot
if ~strcmp(fext(1),'.')
    fext = ['.' fext];
end

[pathstr,name,ext] = fileparts(fname);

%% Check extension
if ~strcmp(ext,fext)
    ext     = fext;
end

%% Reconstruct filename
fname       = fullfile(pathstr,[name ext]);

