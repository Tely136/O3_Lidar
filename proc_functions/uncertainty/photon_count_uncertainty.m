function [total,components] = photon_count_uncertainty(R,P,B,uB,rb,id_on,id_off,C,dsigma,dz,nAverage,L,fd,uTau,channelCorrelation)
% Apply the numbered paper formulas to an acquisition-averaged PC retrieval.
% R: RAW acquisition counts [altitude, acquisition, channel].
% P: averaged corrected/background-subtracted counts [altitude, group, channel].
% B: background fitted to averaged dead-time-corrected counts, same shape.
% uB=[slope; intercept] SEs; rb is their correlation [group, channel].
% uTau is scalar or [on off] SECONDS, not 1/fd or percent; fd is MHz.
% Assumes equal shot count/exposure, independent detection noise, and a
% fixed dead time shared across acquisitions. Tail acquisitions are dropped
% just as in preprocess. Keep negative acquisition noise during propagation.
% BKG/SAT propagation is literal, including signed outputs (README.md).
    if nargin<15, channelCorrelation='independent'; end
    validateattributes(nAverage,{'numeric'},{'scalar','integer','positive'});
    if isscalar(uTau), uTau=[uTau uTau]; end
    if numel(uTau)~=2, error('O3:DeadTimeSize','uTau scalar or [on off] required.'); end
    nprof=size(P,2); nbin=size(P,1);
    if size(R,2)<nAverage*nprof, error('O3:ExposureSize','Not enough raw acquisitions.'); end
    ids_channel=[id_on id_off];
    det=NaN(nbin,nprof,2); sat=det; rawMean=det;
    for j=1:2
        channel=ids_channel(j);
        for i=1:nprof
            ids=(i-1)*nAverage+(1:nAverage);
            Ri=R(:,ids,channel);
            rate=scale_binary_pc(Ri,dz,L);
            Qi=dt_correction(Ri,rate,fd);
            Pi=Qi-B(:,i,channel);
            ui=uP_DET(Pi,Ri);
            det(:,i,j)=sqrt(sum(ui.^2,2))./nAverage;
            sat(:,i,j)=mean(uP_SAT(Pi,dz,L,uTau(j)),2);
            rawMean(:,i,j)=mean(Ri,2);
            invalid=any(~isfinite(rate)|rate>=fd|rate<0,2);
            det(invalid,i,j)=NaN; sat(invalid,i,j)=NaN;
        end
    end
    opts=struct('channelCorrelation',channelCorrelation,...
        'uDET_on',det(:,:,1),'uDET_off',det(:,:,2),...
        'uSAT_on',sat(:,:,1),'uSAT_off',sat(:,:,2));
    z=(1:nbin)' .* dz;
    [total,components.DET,components.BKG,components.SAT]=uNO3(...
        P(:,:,id_on),P(:,:,id_off),rawMean(:,:,1),rawMean(:,:,2),...
        dz,L,uTau,dsigma,C,uB(:,:,id_on),uB(:,:,id_off),...
        rb(:,id_on),rb(:,id_off),z,opts);
end
