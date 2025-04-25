function [dzfwhm irout status message]=NDACC_ResolIR(dzsampling, coef, irinp)

% BELOW ARE SUGGESTED FINAL COMMENTS FOR DISTRIBUTION (PLEASE CHECK AND REVISE AS APPROPRIATE)

% This routine computes vertical resolution using the ISSI-Team-agreed standardized defintion based on the FWHM of an Impulse Response.
% For details, please refer to ISSI Team Report at http://www.issibern.ch/teams/ndacc/ISSI_Team_Report.htm
% IMPORTANT: If several filtering events occur during data processing, this routine must be called each time a filtering event occurs.
%            The values of "irout" returned during the first call must be used for "irinp" in the second call, and so on until the last call.
%            The returned vertical resolution "dzcutf" of the last call only needs to be used (no need for vertical resolution computed during previous calls)

% File history:
% Created in 2011 by GUILLAUME? FRANCK? ??????? , email address?
% Revised in Nov. 2015 by T. Leblanc: Introduce input parameter "dzsampling" and use it to compute vertical resolution in physical unit (last line of routine)

% Inputs:
%     dzsampling = Vertical sampling resolution of signal/profile being filtered (in physical unit)
%     coef       = A vector of 2*N+1 coefficients containing the values of the filter coefficients used to smooth or differentiate the signal/profile
%                  The coefficients must be either symmetric (c(k)=c(-k) for all k), or anti-symmteric (c(k)=-c(-k) for all k)
%     irinp      = A vector (suggested length = total length of your signal/profile) containing the values of the impulse.
%                  If routine first call, it is equal to a Kronecker delta function for smoothing filters, and a Heavyside step function for differentiators
%                  If not first call, it must be equal to the impulse response "irout" returned by the routine during the previous call
% Outputs:
%     status  = Informaiton flag on computaiton status: 1 = succes (smoothing) , 2 = sucess (differentiation) , 0 = failure (cancel)
%     dzfwhm  = The value of vertical resolution using the standardized defintion based on the FWHM of an impulse response
%               (same unit as the input sampling resolution)
%     irout   = A vector of same lenght as hfinp containing the values of the calculated transfer function used to compute the cutoff frequency.
%               If this routine is called again later, use this output as the input impulse "irinp" for the next call
%     message = A string containing an error message (if status=0)

% BELOW ARE THE ORIGNAL COMMENTS BY GUILLAUME? FRANCK?

% Calculation of the normalized resolution of the filter with the FWHM
% method. To obtain the real resolution, multiplie dzfwhm by the initial
% resolution.  
%
% [dzfwhm irout status message]=NDACC_ResolDF(coef) or 
% [dzfwhm irout status message]=NDACC_ResolDF(coef, []): calculate the
% normalized resolution of the filter whose coefficients are "coef" ("coef"
% must be an odd array of value)
%
% [dzfwhm irout status message]=NDACC_ResolDF(coef, irinp): irinp is the
% impulse response of the last filter. The length of irout is the same as
% irinp. 
%
% INPUTS:
%   - coef [array of number]: coefficients of the filter
%   - irinp [array of number]: [optional] transfer function of last
%   filter(s)
%
% OUTPUT:
%   - dzfwhm [array of number]: normalized resolution of the filter
%   - irout [array of number]: transfer function of system (this filter and
%   last filters)
%   - status [Value]: 
%       - 0: problem
%       - 1: low-pass filter
%       - 2: derivative filter
%   - message [String]: error message

%% Parameters
missval = -99;
message='';
dzfwhm=missval;
irout=missval;
nf=512;
nfmin=2*length(coef);
Smax=100;

%% Checking Inputs
% nomber of input
if nargin<1
    status=0;
    message='[ERROR] NDACC_ResolIR.m: Not enought inputs for this function';
    return
end
% Checking coef
if isempty(coef)
    status=0;
    message='[ERROR] NDACC_ResolIR.m: Coef must be no empty';
    return
end
if mod(coef, 2)==0
    status=0;
    message='[ERROR] NDACC_ResolIR.m: Length of coef must be odd';
    return
end
if length(coef)==1 && coef(1)~=1
    status=0;
    message='[ERROR] NDACC_ResolIR.m: If the length of coef is 1, this value must be 1';
    return
end
if size(coef, 1)==1
    coef=coef';
end
% Checking irinp
if exist('irinp', 'var')
    if ~isempty(irinp)
        if size(irinp,1)==1
            irinp=irinp';
        end
        nf=length(irinp);
        if nf<nfmin
            disp(['[WARNING] NDACC_ResolIR.m: the length of irinp is low: it ' ...
                'is better to take a length of ' num2str(nfmin) ' in minimum'])
        end
    else
        irinp=zeros(nf,1);
        irinp(nf/2,1)=Smax;
    end
else
    irinp=zeros(nf,1);
    irinp(nf/2,1)=Smax;
end
irout=ones(nf, 1).*missval;

%% Symetry of coefficients
m=(length(coef)-1)/2;
if coef(1:m)+coef(end:-1:m+2)<1e-5
    status=2;
else
    if coef(1:m)-coef(end:-1:m+2)<1e-5
        status=1;
    else
        status=0;
        message='[ERROR] NDACC_ResolIR.m: coefficients must be symmetric or anti-symmetric';
        return
    end
end

%% Derivative filter
% If derivative filter => integration of irinp
if status==2
    irinp=cumtrapz(irinp); 
end

%% Calculation of the smoothed signal
for j=m+1:length(irinp)-m
    irout(j,1)=irinp(j-m:j+m)'*coef;
end
irout(1:m,1)=0;
irout(length(irinp)-m+1:end,1)=0;

%% Identify FWHM i.e. distance between irout=0.5 (up) and irout=0.5 (down)
i_x(1,1)=find(irout>=0.5*max(irout),1)-1;
i_x(2,1)=find(irout>=0.5*max(irout), 1, 'last');
dzfwhm=dzsampling*diff(i_x-(irout(i_x)-0.5.*max(irout))./(irout(i_x+1)-irout(i_x)));
