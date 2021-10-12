function [xd] = AK_MHIS(gfunc,tdist,boolPlot)
% calculate Pf using AK-MHIS
d=size(tdist,1);
dist=repmat({'norm',0,1},d,1); % standard normal samples

%% Algorithm Parameter
N1=d*7; % initial DoE size
deltaThres=1e-4; % threshold convergence criterion
nis=1e4; % initial sampling size
dnis=2e3; % increment sampling size
inRad=0;
outRad=3; % initial radius
dR=0.2; % radius increment
nadd=5; % limit on expansion
nhis=5; % number of Pf considered for convergence
ub=4;lb=-4; % LHS boundary
learnThres=0; % REIF2 learning threshold

%% Kriging Parameter
theta = ones(1,size(dist,1));
lob = ones(1,size(dist,1))*1e-2; upb = ones(1,size(dist,1))*1e2;
regfunc=@regpoly0; corrfunc=@corrgauss;

%% Generate Sampling Space
fprintf('Generating Sampling Space\n')
xt=ballSample(nis,d,inRad,outRad);
lay_sample={xt};
lay_fp={jointProb(dist,xt)}; % calculate joint prob density
vol=hypVolBall(d,inRad,outRad); % calculate volume of hyperball

%% Generate DoE
fprintf('Generating DoE\n')
lhsrand=lhsdesign(N1,d);
xd=lhsrand*(ub-lb)+lb;
gx=gfunc(transDist(xd,tdist));

%% Train Kriging
fprintf('Train Kriging\n')
dmodel = dacefit(xd, gx, regfunc, corrfunc, theta, lob, upb);

%% Outer Loop
Pfhis=ones(nhis,1);idxt=true(nis,1);
while true
    %% Prepare Pf Calculation
    % calculate dist of layers
    hx=distBall(lay_sample,vol);
    % rejoin all layer
    nLayer=length(lay_sample);
    hp=zeros(nis+(nLayer-1)*dnis,1);
    fp=zeros(nis+(nLayer-1)*dnis,1);
    count=0;
    for i=1:nLayer
        num=length(lay_fp{i});
        fp(count+1:count+num,:)=lay_fp{i};
        hp(count+1:count+num,:)=hx(i);
        count=count+num;
    end
    
    %% Inner Loop
    % idxt ~ samples not in the DoE
    % idxt2 ~ active samples from filtering
    nxt=size(xt,1);
    oldidxt=idxt;
    idxt=true(nxt,1); % initiate idx
    idxt(1:length(oldidxt))=oldidxt; % copy previous idxt
    idxt2=true(nxt,1); % reset all turned off samples
    inner=0; deltas=[];x=[];hisidx=1;
    while true
        inner=inner+1;
        %% Finding Best Point using Learning Function
        fprintf('Finding Best Point\n');
        tic;
        idxcomb=idxt & idxt2;
        realindex=1:nxt;realindex=realindex(idxcomb); % tracking idx
        [ymat,varmat]=predictor(xt, dmodel); 
        [xstar,learn,idxs,learns]=learnREIF2(xt(idxcomb,:),ymat(idxcomb,:),varmat(idxcomb,:),fp(idxcomb,:));
        idxr=learns<learnThres; % samples to deactivate
        
        %% Calculate Pf
        ysign=ymat<=0;
        Pf1=impSamp(ysign,fp,hp); % calculate Pf
        fprintf('Pf: %.3e\n',Pf1);
        Pfhis(hisidx)=Pf1; % update Pf history
        if hisidx==nhis
            hisidx=1;
        else
            hisidx=hisidx+1;
        end
        
        %% Stopping Condition
        % convergence Pf
        if inner>nhis+1
            a=norminv(Pf0);b=norminv(Pf1);
            delta=abs((a-b)/b);
            cond1=delta<deltaThres;
            fprintf('Delta: %.3e\n < %.3e\n',delta,deltaThres);
            if Pf1==Pf0 && Pf1==0 
                % if the pf = 0 increase radius
                break
            end
            % if the delta is a number save it
            if ~isnan(delta) && ~isinf(delta) 
                deltas=[deltas;delta];
                x=[x;inner];
            end
            if length(deltas)>2
                % fitting delta to exponential model
                f=fit(x(:),deltas(:),'exp1','StartPoint',[1,-1]);
                cond2=f(inner)<deltaThres;
                fprintf('Fit Delta: %.3e\n < %.3e\n',f(inner),deltaThres);
                if cond1 && cond2
                    break
                end
            end
        end
        Pf0=mean(Pfhis);
        % learning condition or sample filtering
        fprintf('Learning Threshold: %.3e\n of %.3e\n',learn,learnThres);
        if sum(idxcomb)==0
            break
        end
        
        %% Update DoE and Kriging Model
        fprintf('Update DoE and Kriging Model\n');
        idxt(realindex(idxs))=false; % deactivate sample already in DoE
        idxt2(realindex(idxr))=false; % deactivate sample that fulfill the stopping condition
        gstar=gfunc(transDist(xstar,tdist));
        xd=[xd;xstar];
        gx=[gx;gstar];
        dmodel = dacefit(xd, gx, regfunc, corrfunc, theta, lob, upb);
        fprintf('Added samples: %d\n',size(xd,1));
    end
    
    %% Increase Radius and Check Convergence
    nLayer=length(lay_sample);
    addcount=0;
    while true
        addcount=addcount+1;
        fprintf('Adding Layer %d | ',addcount);
        % add layer temporarily
        thisLayer=length(lay_sample)+1;
        inRad=outRad;outRad=outRad+dR;
        lay_sample{thisLayer}=ballSample(dnis,d,inRad,outRad);
        lay_fp{thisLayer}=jointProb(dist,lay_sample{thisLayer});
        vol(thisLayer,1)=hypVolBall(d,inRad,outRad);
        % update distribution
        hx=distBall(lay_sample,vol);
        count=0;
        hp=zeros(nis+(thisLayer-1)*dnis,1);
        for i=1:thisLayer
            num=length(lay_fp{i});
            hp(count+1:count+num,:)=hx(i);
            count=count+num;
        end
        count=length(ysign);
        ysign(count+1:count+dnis,:)=(predictor(lay_sample{thisLayer},dmodel)<=0);
        fp(count+1:count+dnis,:)=lay_fp{thisLayer};
        % calculate new pf
        Pf2=impSamp(ysign,fp,hp);
        fprintf('Pf: %.3e\n',Pf2);
        % check converge or not
        a=norminv(Pf1);b=norminv(Pf2);
        delta=abs((a-b)/b);
        if delta<deltaThres || addcount>nadd % convegence of Pf or limit reached
            % remove last layer because it is converged
            outRad=inRad;inRad=inRad-dR;
            lay_sample(thisLayer)=[];
            lay_fp(thisLayer)=[];
            vol(thisLayer)=[];
            ysign(count+1:count+dnis,:)=[];
            break
        end
        Pf1=Pf2; % update last Pf
    end
    % if instantly converge
    if thisLayer-nLayer == 1
        break
    end

    %% Reconfigure Candidate Samples
    nLayer=length(lay_sample);
    xt=zeros(nis+(nLayer-1)*dnis,d);count=0;
    for i=1:thisLayer-1
        num=length(lay_fp{i});
        xt(count+1:count+num,:)=lay_sample{i};
        count=count+num;
    end
end
if boolPlot
    fprintf('Ploting the result\n');
    plot_result;
end
clc;
fprintf('Final Pf              : %.3e\n',Pf1);
fprintf('Total Evaluated Sample: %d\n',size(xd,1));
end
