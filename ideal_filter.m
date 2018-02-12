function [Data_Filtered] = ideal_filter(Data, SamplePeriod, Band)
% FORMAT    [Data_Filtered] = IdealFilter(Data, SamplePeriod, Band)
% Input:
% 	Data		    -	2D data matrix (nDimTimePoints * nTimeSeries)
% 	SamplePeriod	-   Sample period, i.e., 1/sample frequency. E.g., TR
%   Band            -   The frequency for filtering, 1*2 Array. Could be:
%                   [LowCutoff_HighPass HighCutoff_LowPass]: band pass filtering
%                   [0 HighCutoff_LowPass]: low pass filtering
%                   [LowCutoff_HighPass 0]: high pass filtering
% Output:
%	Data_Filtered       -   The data after filtering
%-----------------------------------------------------------
% Written by YAN Chao-Gan 120504 based on REST
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com
%
% Very minorly hacked by PJH to remove dependancy on the rest package - peter.hellyer10@imperial.ac.uk


sampleFreq 	 = 1/SamplePeriod;
sampleLength = size(Data,1);
paddedLength = np2(sampleLength); %%% Not entirely sure what the difference between 2^nextpow2(sampleLength) and this is here, but have included the code below for completeness.
LowCutoff_HighPass = Band(1);
HighCutoff_LowPass = Band(2);

% Get the frequency index
if (LowCutoff_HighPass >= sampleFreq/2) % All high stop
    idxLowCutoff_HighPass = paddedLength/2 + 1;
else % high pass, such as freq > 0.01 Hz
    idxLowCutoff_HighPass = ceil(LowCutoff_HighPass * paddedLength * SamplePeriod + 1);
end

if (HighCutoff_LowPass>=sampleFreq/2)||(HighCutoff_LowPass==0) % All low pass
    idxHighCutoff_LowPass = paddedLength/2 + 1;
else % Low pass, such as freq < 0.08 Hz
    idxHighCutoff_LowPass = fix(HighCutoff_LowPass * paddedLength * SamplePeriod + 1);
end

FrequencyMask = zeros(paddedLength,1);
FrequencyMask(idxLowCutoff_HighPass:idxHighCutoff_LowPass,1) = 1;
FrequencyMask(paddedLength-idxLowCutoff_HighPass+2:-1:paddedLength-idxHighCutoff_LowPass+2,1) = 1;

FrequencySetZero_Index = find(FrequencyMask==0);

%Remove the mean before zero padding
Data = Data - repmat(mean(Data),size(Data,1),1);

Data = [Data;zeros(paddedLength -sampleLength,size(Data,2))]; %padded with zero

Data = fft(Data);

Data(FrequencySetZero_Index,:) = 0;

Data = ifft(Data);

Data_Filtered = Data(1:sampleLength,:);
end

function Result = np2(n)
%Compute the min length for FFT according to AFNI's algorithm, By Xiao-Wei Song
%------------------------------------------------------------------------------------------------------------------------------
%	Copyright(c) 2007~2010
%	State Key Laboratory of Cognitive Neuroscience and Learning in Beijing Normal University
%	Written by Xiao-Wei Song 
%	http://resting-fmri.sourceforge.net
% 	<a href="Dawnwei.Song@gmail.com">Mail to Author</a>: Xiaowei Song
%	Version=1.0;
%	Release=20070903;

    if length(n)>1
        n = cast(length(n),class(n));
    end
    if n<16
        Result =2^nextpow2(n);
        return;
    end 
    
    limit =nextpow2(n);             %n=134, limit=8
    tbl=[2^(limit-1):2^limit];      %tbl =128, 129, ... , 256
    tbl =tbl(find(tbl>=n));          %tbl =134, 135, ... , 256
    for x=1:length(tbl)
        Result =tbl(x);
        [f,p]=log2(Result);
        if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
            return;
        end
        if mod(Result,3*5)==0        
            y= Result /(3*5);
            [f,p]=log2(y);
            if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
        if mod(Result,3)==0        
            y= Result /3;
            [f,p]=log2(y);
            if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
        if mod(Result,5)==0        
            y= Result /5;
            [f,p]=log2(y);
            if ~isempty(f) & f == 0.5   %Copy from nextpow2.m
                return;
            end
        end
    end
    Result =NaN;    % Should not reach, except when n=1
end