function tests = test_uncertainty
    tests = functiontests(localfunctions);
end
function setupOnce(t)
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(genpath(fullfile(root,'proc_functions')));
end
function testRawPoisson(t)
    verifyEqual(t,uDET([0 4 100]),[0 2 10]);
    verifyEqual(t,uDET(2500),50);
    verifyTrue(t,isnan(uDET(-1)));
end
function testSignalEquations(t)
    verifyEqual(t,uP_DET(200,100),40);
    verifyEqual(t,uP_DET(25*200,25*100),25*8);
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
    lo=struct('DET',[3;3;NaN],'BKG',[2;2;NaN],'SAT',[1;1;NaN],'O3',[5;6;NaN],'RAYLEIGH',[7;7;7]);
    hi=struct('DET',[NaN;4;4],'BKG',[NaN;4;4],'SAT',[NaN;3;3],'O3',[NaN;10;10],'RAYLEIGH',[7;7;7]);
    [total,m]=merge_pc_uncertainty(lo,hi,[1;0.5;0]);
    verifyEqual(t,m.DET,[3;2.5;4]);
    verifyEqual(t,m.BKG,[2;3;4]); verifyEqual(t,m.SAT,[1;2;3]);
    verifyEqual(t,m.O3,[5;8;10]);
    verifyEqual(t,m.RAYLEIGH,[7;7;7]);
    verifyEqual(t,total,sqrt([88;132.25;190]),'RelTol',1e-12);
end
function testSubtotalAndCoefficientOrder(t)
    [P,C] = fixture();
    R = cat(3,P,P);
    P_sum = R;
    uB = repmat([0.5 1; 2 3],1,1,2);
    rb = zeros(2,2);
    options = struct('dz',4,'L',100,'uTau',0,'dsigma',2,'C',{C});
    options.sigma_o3_on = 1;
    options.sigma_o3_off = 0;
    options.sigma_m_on = 1;
    options.sigma_m_off = 0;
    options.u_sigma_o3_rel = 0;
    options.u_sigma_m_rel = 0;

    [u,parts] = photon_count_uncertainty(R,P_sum,ones(size(P)),ones(size(P)),uB,rb,1,2,options);
    z = (1:7)' * options.dz;
    ub = sqrt(4 + 0.25*z.^2);
    expected = (sqrt(2)*ub(5) - sqrt(2)*ub(3)) / (2*P(4,1)*16);
    verifyEqual(t,parts.BKG(4,1),expected,'RelTol',1e-12);
    verifyEqual(t,u(4,:),sqrt(parts.DET(4,:).^2 + parts.BKG(4,:).^2 + parts.SAT(4,:).^2),'RelTol',1e-12);

    options.channelCorrelation = 'full';
    [~,parts] = photon_count_uncertainty(R,P_sum,ones(size(P)),ones(size(P)),uB,rb,1,2,options);
    verifyEqual(t,parts.BKG(4,:),[0 0]);
end
function testAccumulatedDetectionScaling(t)
    R = 100;
    P = 200;
    n = 30;

    single = uP_DET(P,R);
    accumulated = uP_DET(n*P,n*R);
    verifyEqual(t,accumulated,sqrt(n)*single,'RelTol',1e-12);
    verifyEqual(t,accumulated/(n*P),single/P/sqrt(n),'RelTol',1e-12);
end
function testAccumulatedSaturationScaling(t)
    P = [10 20; 30 40];
    L = 100;
    n = 30;
    uTau = 1e-9;

    single = uP_SAT(P,3.75,L,uTau);
    accumulated = uP_SAT(n*P,3.75,n*L,uTau);
    verifyEqual(t,accumulated,n*single,'RelTol',1e-12);
    verifyEqual(t,accumulated./(n*P),single./P,'RelTol',1e-12);
end
function testSyntheticAccumulationPipeline(t)
    nb = 7010; nr = 5; nc = 4; dz = 3.75; L = 1200; fd = 270;
    z = (1:nb)'; R = zeros(nb,nr,nc);
    for j = 1:nc
        for i = 1:nr
            R(:,i,j) = 20+j+100*exp(-z/500)*(1+0.1*i)+0.1*sin(z+i);
        end
    end
    times = repmat(datetime(2026,1,1)+(0:nr-1)'*minutes(1),1,3);
    [an_sum,Q,~,ids_avg,~,B,uB,rb,R_sum,shots_total] = preprocess(R,R,times,dz,2,L,fd);
    P = Q-B;
    P(P<=0) = NaN;
    C = repmat({[ones(3,1)/3 [-0.5;0;0.5]]},1,nb);
    options = struct('dz',dz,'L',shots_total,'uTau',0,'dsigma',1e-22,'C',{C});
    options.sigma_o3_on = 1;
    options.sigma_o3_off = 0;
    options.sigma_m_on = 1;
    options.sigma_m_off = 0;
    options.u_sigma_o3_rel = 0;
    options.u_sigma_m_rel = 0;

    verifyEqual(t,R_sum(:,1,:),sum(R(:,1:2,:),2));
    verifyEqual(t,an_sum(:,1,:),sum(R(:,1:2,:),2));
    verifyEqual(t,shots_total,2*L*ones(1,2,nc));
    verifyEqual(t,ids_avg,[1;1;2;2;NaN]);

    rate = R_sum * (1e-6*299792458/(2*dz)) ./ shots_total;
    expected = R_sum ./ (1-rate/fd);
    verifyEqual(t,Q,expected,'RelTol',1e-12);

    [total,parts] = photon_count_uncertainty(R_sum,P,ones(nb,2),ones(nb,2),uB,rb,1,2,options);
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
function testVariableShotAccumulation(t)
    nb = 7010; nr = 5; nc = 4; dz = 3.75; fd = 270;
    shots = repmat([100;200;300;400;500],1,nc);
    rate = 2 + 10*exp(-(1:nb)'/1000);
    factor = 1e-6*299792458/(2*dz);
    R = NaN(nb,nr,nc);
    for j = 1:nc
        for i = 1:nr
            R(:,i,j) = rate .* shots(i,j) ./ factor;
        end
    end
    times = repmat(datetime(2026,1,1)+(0:nr-1)'*minutes(1),1,3);
    [~,Q,~,ids_avg,~,~,~,~,R_sum,shots_total] = preprocess(R,R,times,dz,2,shots,fd);

    verifyEqual(t,shots_total,repmat([300 700],1,1,nc));
    expected_rate = repmat(rate,1,2,nc);
    verifyEqual(t,scale_binary_pc(R_sum,dz,shots_total),expected_rate,'RelTol',1e-12);
    verifyEqual(t,scale_binary_pc(Q,dz,shots_total),expected_rate./(1-expected_rate/fd),'RelTol',1e-12);
    verifyEqual(t,ids_avg,[1;1;2;2;NaN]);
end
function testExplicitPropagationInputs(t)
    [P,C] = fixture();
    options = struct('dz',4,'dsigma',2,'C',{C});
    zero = zeros(size(P));

    [u,det,bkg,sat] = uNO3(P,P,2*ones(size(P)),3*ones(size(P)),zero,zero,zero,zero,options);
    verifyEqual(t,u(4,:),sqrt(13/2)./(P(4,:)*16),'RelTol',1e-12);
    verifyEqual(t,u,det);
    verifyEqual(t,bkg(4,:),[0 0]);
    verifyEqual(t,sat(4,:),[0 0]);
end
function testAccumulationBeforeDeadTime(t)
    nb = 7010; nc = 4; dz = 3.75; L = 100; fd = 270;
    R = zeros(nb,2,nc);
    R(:,1,:) = 10;
    R(:,2,:) = 100;
    times = repmat(datetime(2026,1,1)+(0:1)'*minutes(1),1,3);
    [~,Q,~,~,~,~,~,~,R_sum,shots_total] = preprocess(R,R,times,dz,2,L,fd);

    raw_rate = 110 * (1e-6*299792458/(2*dz)) / (2*L);
    expected = 110/(1-raw_rate/fd);
    verifyEqual(t,Q(4,1,1),expected,'RelTol',1e-12);
    verifyEqual(t,R_sum(4,1,1),110);
    verifyEqual(t,shots_total(1,1,1),2*L);

    separately_corrected = 10/(1-10*(1e-6*299792458/(2*dz*L))/fd) + 100/(1-100*(1e-6*299792458/(2*dz*L))/fd);
    verifyGreaterThan(t,abs(Q(4,1,1)-separately_corrected),1);
end
function testUnifiedCrossSectionBudget(t)
    [P,C] = fixture();
    R = cat(3,P,P);
    N = 10*ones(size(P));
    Na = 1000*ones(size(P));
    uB = ones(2,2,2);
    rb = zeros(2,2);
    options = struct('dz',4,'L',100,'uTau',1e-9,'dsigma',2,'C',{C});
    options.sigma_o3_on = 5;
    options.sigma_o3_off = 3;
    options.sigma_m_on = 0.05;
    options.sigma_m_off = 0.03;
    options.u_sigma_o3_rel = 0.02;
    options.u_sigma_m_rel = 0.02;

    [total,parts] = photon_count_uncertainty(R,R,N,Na,uB,rb,1,2,options);
    verifyEqual(t,parts.O3,N*0.02*sqrt(34)/2,'RelTol',1e-12);
    verifyEqual(t,parts.RAYLEIGH,Na*0.02*0.02/2,'RelTol',1e-12);
    expected = sqrt(parts.DET.^2 + parts.BKG.^2 + parts.SAT.^2 + parts.O3.^2 + parts.RAYLEIGH.^2);
    verifyEqual(t,total,expected);
end
function testSharedAtmosphericMerge(t)
    low = struct('DET',1,'BKG',2,'SAT',3,'O3',-6,'RAYLEIGH',7);
    high = struct('DET',1,'BKG',2,'SAT',3,'O3',2,'RAYLEIGH',7);

    [total,merged] = merge_pc_uncertainty(low,high,0.75);
    verifyEqual(t,merged.O3,-4);
    verifyEqual(t,merged.RAYLEIGH,7);
    verifyEqual(t,total,sqrt(0.75^2+0.25^2+2^2+3^2+(-4)^2+7^2),'RelTol',1e-12);
end
function [P,C]=fixture()
    P=repmat([100 200],7,1);
    C=repmat({[ones(3,1)/3 [-0.5;0;0.5]]},1,7);
end
