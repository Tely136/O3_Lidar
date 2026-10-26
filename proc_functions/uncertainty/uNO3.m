function [u,u_DET,u_BKG,u_SAT] = uNO3(P_on,P_off,R_on,R_off,dz,L,uTau,dsigma_o3,C,uB_on,uB_off,rb_on,rb_off,z,options)
% Photon-count-only subtotal: DET+BKG+SAT, not full Eq.99.
% Counts (never MHz), altitude x profile; dsigma SINGLE-PATH m^2.
% uB=[SLOPE SE; INTERCEPT SE], as returned by bg_correction.
% Optional struct: nAverage=1; channelCorrelation='independent' or 'full';
% uDET_on/off, uSAT_on/off override signal SDs propagated before averaging.
% Nonlinear acquisition means must use photon_count_uncertainty.m.
% Literal Eqs.36/43/38/45 may yield signed component profiles (README).
    if nargin<15, options=struct(); end
    if ~isfield(options,'nAverage'), options.nAverage=1; end
    if ~isfield(options,'channelCorrelation'), options.channelCorrelation='independent'; end
    nprof=size(P_on,2);
    if isvector(P_on), P_on=P_on(:); P_off=P_off(:); R_on=R_on(:); R_off=R_off(:); nprof=1; end
    if ~isequal(size(uB_on),[2 nprof]) || ~isequal(size(uB_off),[2 nprof])
        error('O3:BackgroundSize','uB must be 2 x profiles, [slope; intercept].');
    end
    b_on=uP_BKG(uB_on(2,:),uB_on(1,:),reshape(rb_on,1,[]),z);
    b_off=uP_BKG(uB_off(2,:),uB_off(1,:),reshape(rb_off,1,[]),z);
    if isfield(options,'uDET_on') && isfield(options,'uDET_off')
        det_on=options.uDET_on; det_off=options.uDET_off;
    else
        det_on=uP_DET(P_on,R_on,options.nAverage); det_off=uP_DET(P_off,R_off,options.nAverage);
    end
    if isscalar(uTau)
        tau_on=uTau;
        tau_off=uTau;
    elseif numel(uTau)==2
        tau_on=uTau(1);
        tau_off=uTau(2);
    else
        error('O3:DeadTimeSize','uTau scalar or [on off] seconds required.');
    end
    if isfield(options,'uSAT_on') && isfield(options,'uSAT_off')
        sat_on=options.uSAT_on; sat_off=options.uSAT_off;
    else
        sat_on=uP_SAT(P_on,dz,L,tau_on); sat_off=uP_SAT(P_off,dz,L,tau_off);
    end
    u_DET=uNO3_P_X_NC(P_on,P_off,det_on,det_off,dsigma_o3,dz,C,'none');
    mode=validatestring(options.channelCorrelation,{'independent','full'});
    if strcmp(mode,'independent')
        u_BKG=uNO3_P_X_NC(P_on,P_off,b_on,b_off,dsigma_o3,dz,C,'full');
        u_SAT=uNO3_P_X_NC(P_on,P_off,sat_on,sat_off,dsigma_o3,dz,C,'full');
    else
        u_BKG=uNO3_P_X_FC(P_on,P_off,b_on,b_off,dsigma_o3,dz,C);
        u_SAT=uNO3_P_X_FC(P_on,P_off,sat_on,sat_off,dsigma_o3,dz,C);
    end
    u=sqrt(u_DET.^2+u_BKG.^2+u_SAT.^2);
end
