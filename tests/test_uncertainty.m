function tests = test_uncertainty
    tests = functiontests(localfunctions);
end
function setupOnce(t)
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root,'proc_functions')));
end
function testRawPoisson(t)
    verifyEqual(t,uDET([0 4 100]),[0 2 10]);
    verifyEqual(t,uDET(100,25),2);
    verifyTrue(t,isnan(uDET(-1)));
end
function testSignalEquations(t)
    verifyEqual(t,uP_DET(200,100),40);
    verifyEqual(t,uP_DET(200,100,25),8);
    verifyTrue(t,isnan(uP_DET(0,0)));
    P=[10 20;30 40]; L=[100 200]; tau=[1e-9 2e-9];
    expected=299792458./(2*3.75.*L).*P.^2.*tau;
    verifyEqual(t,uP_SAT(P,3.75,L,tau),expected,'RelTol',1e-12);
    verifyEqual(t,uP_SAT(P,3.75,100,0),zeros(2));
end
function testBackgroundEquation(t)
    z=[0 2 4];
    expected=sqrt(4+0.25*z(:).^2-0.5*z(:));
    verifyEqual(t,uP_BKG(2,0.5,-0.25,z),expected,'AbsTol',1e-12);
    verifySize(t,uP_BKG([2 3],[0.5 1],[0 0],z),[3 2]);
    verifyEqual(t,uP_BKG(2,1,-1,2),0);
end
function testDetectionDerivativeMultipleProfiles(t)
    [P,C]=fixture();
    u=uNO3_P_X_NC(P,P,2,3,2,4,C);
    % Independent ON/OFF: relative variance 13/P^2, filter energy 1/2.
    expected=sqrt(13/2)./(P(2:6,:).*16);
    verifyEqual(t,u(2:6,:),expected,'RelTol',1e-12);
    verifyTrue(t,all(isnan(u([1 7],:)),'all'));
    d=(1:7)';
    v=uNO3_P_X_NC(P,P,2,3,d,4,C);
    verifyEqual(t,v(4,2),sqrt(13/2)/(P(4,2)*32),'RelTol',1e-12);
    P(4,1)=0; v=uNO3_P_X_NC(P,P,2,3,2,4,C);
    verifyTrue(t,all(isnan(v(3:5,1))));
    verifyTrue(t,isfinite(v(4,2)));
end
function testLiteralCorrelatedEquations(t)
    [P,C]=fixture(); z=(1:7)'; ua=repmat(8-z,1,2); ub=2*ua;
    a=uNO3_P_X_NC(P,P,ua,ub,2,4,C,'full');
    expected=-sqrt(5)./(16*P(4,:));
    verifyEqual(t,a(4,:),expected,'RelTol',1e-12);
    verifyLessThan(t,a(4,1),0); % Signed printed Eq.36/43, intentionally.
    b=uNO3_P_X_FC(P,P,ua,ub,2,4,C);
    verifyEqual(t,b(4,:),-1./(16*P(4,:)),'RelTol',1e-12);
    same=uNO3_P_X_FC(P,P,ua,ua,2,4,C);
    verifyEqual(t,same(4,:),[0 0]);
end
function testCrossSections(t)
    N=[100 -200]; ds=[2 4];
    verifyEqual(t,uNO3_dS_O3_R(N,ds,3,4),[250 -250]);
    verifyEqual(t,uNO3_dS_M_R([1000 2000],ds,3,4),[2500 2500]);
    % Finite difference of N=signal/(2*ds), one cross section at a time.
    s=400; on=5; off=3; h=1e-5;
    sensitivity=(s/(2*(on+h-off))-s/(2*(on-h-off)))/(2*h);
    verifyEqual(t,uNO3_dS_O3_R(100,2,1,0),abs(sensitivity),'RelTol',1e-8);
end
function testAirDensity(t)
    ds=[2 4]; dm=[1 2]; dO=[2 4]; q=[0.25 0.5]; un=[4 8];
    factor=[0.75 1];
    verifyEqual(t,uNO3_Na1(ds,dm,dO,q,un),factor.*un);
    verifyEqual(t,uNO3_Na2_NC(ds,dm,dO,q,[10 20],100,200,3,8),factor.*[10 20]*0.05,'AbsTol',1e-12);
    verifyEqual(t,uNO3_Na2_FC(ds,dm,dO,q,[10 20],100,200,3,8),factor.*[10 20]*0.01,'AbsTol',1e-12);
end
function testRangeMerge(t)
    verifyEqual(t,uNO3_MRG_NC(4,3,1),3);
    verifyEqual(t,uNO3_MRG_NC(4,3,0),4);
    verifyEqual(t,uNO3_MRG_NC(4,3,0.5),2.5);
    verifyEqual(t,uNO3_MRG_FC(4,3,0.5),3.5);
    lo=struct('DET',[3;3;NaN],'BKG',[2;2;NaN],'SAT',[1;1;NaN]);
    hi=struct('DET',[NaN;4;4],'BKG',[NaN;4;4],'SAT',[NaN;3;3]);
    [total,m]=merge_pc_uncertainty(lo,hi,[1;0.5;0]);
    verifyEqual(t,m.DET,[3;2.5;4]);
    verifyEqual(t,m.BKG,[2;3;4]); verifyEqual(t,m.SAT,[1;2;3]);
    verifyEqual(t,total,sqrt([14;19.25;41]),'RelTol',1e-12);
end
function testSubtotalAndCoefficientOrder(t)
    [P,C]=fixture(); z=(1:7)'; nprof=2;
    opts=struct('uDET_on',2*ones(7,2),'uDET_off',3*ones(7,2),...
        'uSAT_on',z+zeros(7,2),'uSAT_off',2*z+zeros(7,2));
    uB=[0.5 1;2 3]; rb=[0 0];
    [u,det,bkg,sat]=uNO3(P,P,P,P,4,100,0,2,C,uB,uB,rb,rb,z,opts);
    ub=sqrt(4+0.25*z.^2);
    expected=(sqrt(2)*ub(5)-sqrt(2)*ub(3))/(2*P(4,1)*16);
    verifyEqual(t,bkg(4,1),expected,'RelTol',1e-12);
    verifyEqual(t,u(4,:),sqrt(det(4,:).^2+bkg(4,:).^2+sat(4,:).^2),'RelTol',1e-12);
    opts.channelCorrelation='full';
    [~,~,bkg]=uNO3(P,P,P,P,4,100,0,2,C,uB,uB,rb,rb,z,opts);
    verifyEqual(t,bkg(4,:),[0 0]);
end
function testAcquisitionAveraging(t)
    nb=7; n=2; nc=4; dz=3.75; L=100; fd=270; tau=1e-9;
    R=zeros(nb,n,nc); R(:,1,:)=10; R(:,2,:)=100;
    rates=R* (1e-6*299792458/(2*dz*L));
    Q=R./(1-rates/fd); B=ones(nb,1,nc);
    P=mean(Q,2)-B; uB=zeros(2,1,nc); rb=zeros(1,nc);
    C=repmat({[ones(3,1)/3 [-0.5;0;0.5]]},1,nb);
    [~,parts]=photon_count_uncertainty(R,P,B,uB,rb,1,2,C,2,dz,n,L,fd,tau);
    Pi=squeeze(Q(4,:,1))-1; Ri=[10 100];
    acquisition_u=(Pi./Ri).^2.*sqrt(Ri);
    mean_u=sqrt(sum(acquisition_u.^2))/2;
    % Same channels, independent DET, three-point derivative energy=1/2.
    expected=mean_u/P(4,1,1)/(4*dz);
    verifyEqual(t,parts.DET(4),expected,'RelTol',1e-12);
    naive=uP_DET(P(4),mean(Ri),2)/(P(4)*4*dz);
    verifyGreaterThan(t,abs(naive-expected),1e-6);
    verifyEqual(t,parts.SAT(4),0,'AbsTol',1e-20);
end
function testSaturationAveraging(t)
    nb=7; n=2; nc=4; dz=3.75; L=100; fd=270; tau=1e-9;
    z=(1:nb)'; R=zeros(nb,n,nc);
    for j=1:nc
        R(:,1,j)=10*z; R(:,2,j)=20*z;
    end
    rate=R.*(1e-6*299792458/(2*dz*L));
    Q=R./(1-rate/fd); B=ones(nb,1,nc); P=mean(Q,2)-B;
    C=repmat({[ones(3,1)/3 [-0.5;0;0.5]]},1,nb);
    [~,parts]=photon_count_uncertainty(R,P,B,zeros(2,1,nc),zeros(1,nc),...
        1,2,C,2,dz,n,L,fd,tau);
    % Eq.34 averaged linearly because one uncertain tau is shared in time.
    perAcq=299792458/(2*dz*L).*(Q(:,:,1)-1).^2*tau;
    avgSat=mean(perAcq,2);
    rel=sqrt(2)*avgSat./P(:,:,1);
    expected=(rel(5)-rel(3))/(8*dz);
    verifyEqual(t,parts.SAT(4),expected,'RelTol',1e-12);
    % It must NOT acquire the independent-noise 1/sqrt(n) reduction.
    verifyGreaterThan(t,abs(parts.SAT(4)-expected/sqrt(2)),1e-9);
end
function testSyntheticPreprocessPipeline(t)
    % Real background fit and averaging; 7000+ bins required by preprocess.
    nb=7010; nr=4; nc=4; dz=3.75; L=1200; fd=270;
    z=(1:nb)'; R=zeros(nb,nr,nc);
    for j=1:nc
        for i=1:nr
            R(:,i,j)=20+j+100*exp(-z/500)*(1+0.1*i)+0.1*sin(z+i);
        end
    end
    times=repmat(datetime(2026,1,1)+(0:nr-1)'*minutes(1),1,3);
    [~,Q,~,~,B,uB,rb]=preprocess(R,R,times,dz,2,L,fd);
    P=Q-B; P(P<=0)=NaN;
    C=repmat({[ones(3,1)/3 [-0.5;0;0.5]]},1,nb);
    [total,parts]=photon_count_uncertainty(R,P,B,uB,rb,1,2,C,1e-22,dz,2,L,fd,0);
    verifySize(t,total,[nb 2]);
    verifyTrue(t,all(isfinite(total(50:100,:)),'all'));
    verifyEqual(t,total(50:100,:),hypot(parts.DET(50:100,:),parts.BKG(50:100,:)),'RelTol',1e-12);
end
function testInvalidInputs(t)
    [P,C]=fixture();
    verifyError(t,@()uNO3_P_X_NC(P,P,-1,1,2,4,C),'O3:NegativeUncertainty');
    verifyError(t,@()uNO3_P_X_NC(P,P,1,1,2,4,C(1:2)),'O3:FilterSize');
end
function testCorrelatedRayleighCrossSection(t)
    Na = [1000 2000; 3000 4000];
    dS_O3 = [2 4; 5 8];

    verifyEqual(t,uNO3_dS_M_S(Na,dS_O3,0.3,0.1),Na .* 0.2 ./ dS_O3,'RelTol',1e-12);
    verifyEqual(t,uNO3_dS_M_S(Na,dS_O3,0.3,0.3),zeros(2));

    % A common 2% scale error gives 2% of the Rayleigh correction.
    correction = Na .* (5 - 3) ./ dS_O3;
    verifyEqual(t,uNO3_dS_M_S(Na,dS_O3,0.02*5,0.02*3),0.02 .* correction,'RelTol',1e-12);
end
function testOzoneCrossSectionRelativeFloor(t)
    N = [100 200; 300 400];
    sigma_on = [5 6; 7 8];
    sigma_off = ones(2);
    differential = sigma_on - sigma_off;
    expected = N .* 0.02 .* sqrt(sigma_on.^2 + sigma_off.^2) ./ differential;

    verifyEqual(t,uNO3_dS_O3_R(N,differential,0.02.*sigma_on,0.02.*sigma_off),expected,'RelTol',1e-12);
    verifyEqual(t,uNO3_dS_O3_R(N,differential,zeros(2),zeros(2)),zeros(2));
end
function [P,C]=fixture()
    P=repmat([100 200],7,1);
    C=repmat({[ones(3,1)/3 [-0.5;0;0.5]]},1,7);
end
