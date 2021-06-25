clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
 CP1(:,32:36)=[];
CP1=CP1(:,[3 5 7 9 8 11]);
 %%
% CP1(:,32:35)=[]; % includes commodity only
%%
% CP1(:,32:34)=[];includes oil only
% CP1(:,33)=[]; % includes oil only
%% Data preprocessing
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
windowSize=100;
GIRF=[];
shock=[];
t= nRow-windowSize; 
progressbar
for J = 1:t
        wend=J+windowSize-1;
    y = filt(J:wend,:); 
    % Set SVAR parameters
    p = 2;      % Order of VAR lags
    q = 10;     % Order of VMA lags      

    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1);
    nconst = size(ylag,2);

    for i = 1:p;

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end
    
    % OLS Estimates for initialization
    bar    = ylag\trimr(y,p,0);
     v      = trimr(y,p,0) - ylag*bar;
% % % % % %     omegav = v'*v/size(v,1);
% % % % % %     s      = chol(omegav)';

   
    bstart = [
                0.001855486206481
                0.068693371469487
                0.023880592273368
                0.001989809184491
                0.013648672676092
                0.000281803828058
               -0.121999669995451
               -0.002888991832387
                0.002650363820455
                0.217906436715095
                0.001195294344330
                0.008613191798914
                0.175400639477836
                1.103000246244063
                0.440791629588056
                0.495605700916872
                0.365097267442524
               -0.011572135693424
                0.650892298805934
               -1.176100342250688
               -0.342100219212525
               -1.059396487450620
               -0.545902063123713
               -0.729900334997059
               -0.402000678739595
                1.013901399027622
               -1.601300249119193
                1.018606365787122
                0.791195906556810
                0.818998254746848
                0.338395871264844
                0.465287945621782
               -1.145684521466255
            ];

  options = optimset(   'Display',              'iter', ...
                        'MaxIter',              20000,   ...
                        'MaxFunEvals',          20000,   ...
                        'TolX',                 1e-8,   ...
                        'TolFun',               1e-4 );  
            
    [ b,fval]  = fminunc( @(b) neglogl(b,v),bstart,options ); 


    disp(' Log-likelihood function ');
    disp( fval );
    disp(' structural parameter ');
    disp( b );
    disp(' VAR disturbances ');
    disp( v );
    
   % Matrix required to compute time-varying variance-covariances              
    s = getSmatrix( b,v );
    omegavv = v'*v/size(v,1);
    omegavs= s'*s/size(s,1);
    omegav=omegavs(1:6,1:6);
    omegav=topdm(omegav);

    
  % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    b = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = b(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);          
    % Construct C(L) matrices (vector moving average)  
    c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ]);

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end
       tmp = ss';
       c  = [ c  tmp(:) ];

    end

%%    
   ar0=a1;
   ar0={ar0};
   ma0=tmp;
   ma0={ma0};
   InnovCov=omegav; 
   GIRF = armairf(ar0,ma0,'Method','generalized','InnovCov',InnovCov,'NumObs',10);
   disp(J);
   
   disp ('residual')
   disp([num2str(v)]);
   shock=[shock,v];
   
   disp('GIRF')
   disp([num2str(GIRF)]);
   GIRF=[GIRF,GIRF ];
  progressbar(J/t)
end
    %% The old stuff begins here, lets see
%% TO 
impulse=GIRF;
impulse=impulse(10,:,:);
impulse1=permute(impulse, [1 2 3]);
%impulse2 =reshape(impulse1, [], size(impulse, 2),1);
impulse2=impulse(:);
impulse3= repmat(impulse1,10);
impulse3=impulse3(:);
shock1=shock';
shock2=shock1(:,1);
fixing= length(impulse3(:))-length(shock2(:));
impulse4=trimr(impulse3,0,fixing);
shock4=reshapeg(shock2,[],6);
impulse5=reshapeg(impulse4,[],6);
SHOCK=shock4;
IRF=impulse5;
GHDTO=IRF.*SHOCK;

%% FROM
impulsefrom= GIRF;
impulsefrom=impulsefrom(10,:,:);
impulsefrom1=permute(impulsefrom, [1 3 2]);
impulsefrom2=impulsefrom1(:);
impulsefrom3= repmat(impulsefrom1,10);
shockfrom=shock;
shockfrom2=shockfrom(1,:);
fixing= length(impulsefrom3(:))-length(shockfrom2(:));
impulsefrom3=impulsefrom3(:);
impulsefrom4=trimr(impulsefrom3,0,fixing);
shockfrom4=reshapeg(shockfrom2,[],6);
impulsefrom5=reshapeg(impulsefrom4,[],6);
SHOCKFROM=shockfrom4;
IRFfrom=impulsefrom5;
GHDFROM=IRFfrom.*SHOCKFROM;   
%     %% The old stuff begins here, lets see
%  
% impulse=GIRF;
% impulse1=permute(impulse, [1 3 2]);
% impulse2=reshape(impulse1, [], size(impulse, 2),1);
% impulse2=impulse(:);
% impulse3= repmat(impulse2,10);
% u=reshapeg(v,[],6);
% fixing= length(impulse3(:))-length(u(:));
% impulse3=impulse3(:);
% impulse4=trimr(impulse3,0,fixing);
% %impulse4=trimr(impulse3,0,1242);
% impulse5=reshapeg(impulse4,[],6);
% %% experiment with To and From
% ACSVDTOraw = u* impulse5'; %plot(movmean(HD(:,1),250)) works
% ACSVDFROMraw= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
% % HD=zscore(HD);
% % HD1=zscore(HD1);
% ACSVDFROMraw=ACSVDFROMraw';
save ACSVDTOraw.mat GHDTO ;
save ACSVDFROMraw.mat GHDFROM ;


clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
 CP1(:,32:36)=[];
CP1=CP1(:,[28 29 31 17 30 2]);
 %%
% CP1(:,32:35)=[]; % includes commodity only
%%
% CP1(:,32:34)=[];includes oil only
% CP1(:,33)=[]; % includes oil only
%% Data preprocessing
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
windowSize=100;
%% Inputs
y=filt;
    
    % Set SVAR parameters
    p = 2;      % Order of VAR lags
    q = 10;     % Order of VMA lags      

    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1);
    nconst = size(ylag,2);

    for i = 1:p;

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end
    
    % OLS Estimates for initialization
    bar    = ylag\trimr(y,p,0);
     v      = trimr(y,p,0) - ylag*bar;
% % % % % %     omegav = v'*v/size(v,1);
% % % % % %     s      = chol(omegav)';

   
    bstart = [
                0.001855486206481
                0.068693371469487
                0.023880592273368
                0.001989809184491
                0.013648672676092
                0.000281803828058
               -0.121999669995451
               -0.002888991832387
                0.002650363820455
                0.217906436715095
                0.001195294344330
                0.008613191798914
                0.175400639477836
                1.103000246244063
                0.440791629588056
                0.495605700916872
                0.365097267442524
               -0.011572135693424
                0.650892298805934
               -1.176100342250688
               -0.342100219212525
               -1.059396487450620
               -0.545902063123713
               -0.729900334997059
               -0.402000678739595
                1.013901399027622
               -1.601300249119193
                1.018606365787122
                0.791195906556810
                0.818998254746848
                0.338395871264844
                0.465287945621782
               -1.145684521466255
            ];

  options = optimset(   'Display',              'iter', ...
                        'MaxIter',              20000,   ...
                        'MaxFunEvals',          20000,   ...
                        'TolX',                 1e-8,   ...
                        'TolFun',               1e-4 );  
            
    [ b,fval]  = fminunc( @(b) neglogl(b,v),bstart,options ); 


    disp(' Log-likelihood function ');
    disp( fval );
    disp(' structural parameter ');
    disp( b );
    disp(' VAR disturbances ');
    disp( v );
    
   % Matrix required to compute time-varying variance-covariances              
    s = getSmatrix( b,v );
    omegav = v'*v/size(v,1);
    omegavs= s'*s/size(s,1);
    omegav=omegavs(1:6,1:6);
    omegav=topdm(omegav);
    
   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    b = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = b(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);          
    % Construct C(L) matrices (vector moving average)  
    c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ]);

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end
       tmp = ss';
       c  = [ c  tmp(:) ];

    end

%%    
   ar0=a1;
   ar0={ar0};
   ma0=tmp;
   ma0={ma0};
   InnovCov=omegav; 
   GIRF = armairf(ar0,ma0,'Method','generalized','InnovCov',InnovCov,'NumObs',10);
   
    %% The old stuff begins here, lets see
 
impulse=GIRF;
impulse1=permute(impulse, [1 3 2]);
impulse2=reshape(impulse1, [], size(impulse, 2),1);
impulse2=impulse(:);
impulse3= repmat(impulse2,10);
u=reshapeg(v,[],6);
fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
%impulse4=trimr(impulse3,0,1242);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
ECSVDTOraw = u* impulse5'; %plot(movmean(HD(:,1),250)) works
ECSVDFROMraw= u'.*impulse5';%plot(movmean(HD1(1,:),250)) works
ECSVDFROMraw= ECSVDFROMraw';
% HD=zscore(HD);
% HD1=zscore(HD1);
save ECSVDTOraw.mat ECSVDTOraw ;
save ECSVDFROMraw.mat ECSVDFROMraw ;


clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
 CP1(:,32:36)=[];
CP1=CP1(:,[20 21 22 23 24 25]);
 %%
% CP1(:,32:35)=[]; % includes commodity only
%%
% CP1(:,32:34)=[];includes oil only
% CP1(:,33)=[]; % includes oil only
%% Data preprocessing
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
windowSize=100;
%% Inputs
y=filt;
    
    % Set SVAR parameters
    p = 2;      % Order of VAR lags
    q = 10;     % Order of VMA lags      

    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1);
    nconst = size(ylag,2);

    for i = 1:p;

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end
    
     % OLS Estimates for initialization
    bar    = ylag\trimr(y,p,0);
     v      = trimr(y,p,0) - ylag*bar;
% % % % % %     omegav = v'*v/size(v,1);
% % % % % %     s      = chol(omegav)';

   
    bstart = [
                0.001855486206481
                0.068693371469487
                0.023880592273368
                0.001989809184491
                0.013648672676092
                0.000281803828058
               -0.121999669995451
               -0.002888991832387
                0.002650363820455
                0.217906436715095
                0.001195294344330
                0.008613191798914
                0.175400639477836
                1.103000246244063
                0.440791629588056
                0.495605700916872
                0.365097267442524
               -0.011572135693424
                0.650892298805934
               -1.176100342250688
               -0.342100219212525
               -1.059396487450620
               -0.545902063123713
               -0.729900334997059
               -0.402000678739595
                1.013901399027622
               -1.601300249119193
                1.018606365787122
                0.791195906556810
                0.818998254746848
                0.338395871264844
                0.465287945621782
               -1.145684521466255
            ];

  options = optimset(   'Display',              'iter', ...
                        'MaxIter',              20000,   ...
                        'MaxFunEvals',          20000,   ...
                        'TolX',                 1e-8,   ...
                        'TolFun',               1e-4 );  
            
    [ b,fval]  = fminunc( @(b) neglogl(b,v),bstart,options ); 


    disp(' Log-likelihood function ');
    disp( fval );
    disp(' structural parameter ');
    disp( b );
    disp(' VAR disturbances ');
    disp( v );
    
   % Matrix required to compute time-varying variance-covariances              
    s = getSmatrix( b,v );
    omegav = v'*v/size(v,1);
    omegavs= s'*s/size(s,1);
    omegav=omegavs(1:6,1:6); 
    omegav=topdm(omegav);
    % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    b = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = b(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);          
    % Construct C(L) matrices (vector moving average)  
    c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ]);

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end
       tmp = ss';
       c  = [ c  tmp(:) ];

    end

%%    
   ar0=a1;
   ar0={ar0};
   ma0=tmp;
   ma0={ma0};
   InnovCov=omegav; 
   GIRF = armairf(ar0,ma0,'Method','generalized','InnovCov',InnovCov,'NumObs',10);
   
    %% The old stuff begins here, lets see
 
impulse=GIRF;
impulse1=permute(impulse, [1 3 2]);
impulse2=reshape(impulse1, [], size(impulse, 2),1);
impulse2=impulse(:);
impulse3= repmat(impulse2,10);
u=reshapeg(v,[],6);
fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
%impulse4=trimr(impulse3,0,1242);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
GCSVDTOraw = u* impulse5'; %plot(movmean(HD(:,1),250)) works
GCSVDFROMraw= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
GCSVDFROMraw=GCSVDFROMraw';
% HD=zscore(HD);
% HD1=zscore(HD1);
save GCSVDTOraw.mat GCSVDTOraw ;
save GCSVDFROMraw.mat GCSVDFROMraw ;


clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
 CP1(:,32:36)=[];
CP1=CP1(:,[1 19 26 27 4 6]);
 %%
% CP1(:,32:35)=[]; % includes commodity only
%%
% CP1(:,32:34)=[];includes oil only
% CP1(:,33)=[]; % includes oil only
%% Data preprocessing
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
windowSize=100;
%% Inputs
y=filt;
    
    % Set SVAR parameters
    p = 2;      % Order of VAR lags
    q = 10;     % Order of VMA lags      

    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1);
    nconst = size(ylag,2);

    for i = 1:p;

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end
    
     % OLS Estimates for initialization
    bar    = ylag\trimr(y,p,0);
     v      = trimr(y,p,0) - ylag*bar;
% % % % % %     omegav = v'*v/size(v,1);
% % % % % %     s      = chol(omegav)';

   
    bstart = [
                0.001855486206481
                0.068693371469487
                0.023880592273368
                0.001989809184491
                0.013648672676092
                0.000281803828058
               -0.121999669995451
               -0.002888991832387
                0.002650363820455
                0.217906436715095
                0.001195294344330
                0.008613191798914
                0.175400639477836
                1.103000246244063
                0.440791629588056
                0.495605700916872
                0.365097267442524
               -0.011572135693424
                0.650892298805934
               -1.176100342250688
               -0.342100219212525
               -1.059396487450620
               -0.545902063123713
               -0.729900334997059
               -0.402000678739595
                1.013901399027622
               -1.601300249119193
                1.018606365787122
                0.791195906556810
                0.818998254746848
                0.338395871264844
                0.465287945621782
               -1.145684521466255
            ];

  options = optimset(   'Display',              'iter', ...
                        'MaxIter',              20000,   ...
                        'MaxFunEvals',          20000,   ...
                        'TolX',                 1e-8,   ...
                        'TolFun',               1e-4 );  
            
    [ b,fval]  = fminunc( @(b) neglogl(b,v),bstart,options ); 


    disp(' Log-likelihood function ');
    disp( fval );
    disp(' structural parameter ');
    disp( b );
    disp(' VAR disturbances ');
    disp( v );
    
   % Matrix required to compute time-varying variance-covariances              
    s = getSmatrix( b,v );
    omegav = v'*v/size(v,1);
    omegavs= s'*s/size(s,1);
    omegav=omegavs(1:6,1:6);
    omegav=topdm(omegav);
   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    b = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = b(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);          
    % Construct C(L) matrices (vector moving average)  
    c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ]);

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end
       tmp = ss';
       c  = [ c  tmp(:) ];

    end

%%    
   ar0=a1;
   ar0={ar0};
   ma0=tmp;
   ma0={ma0};
   InnovCov=omegav; 
   GIRF = armairf(ar0,ma0,'Method','generalized','InnovCov',InnovCov,'NumObs',10);
   
    %% The old stuff begins here, lets see
 
impulse=GIRF;
impulse1=permute(impulse, [1 3 2]);
impulse2=reshape(impulse1, [], size(impulse, 2),1);
impulse2=impulse(:);
impulse3= repmat(impulse2,10);
u=reshapeg(v,[],6);
fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
%impulse4=trimr(impulse3,0,1242);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
OEDSVDTOraw = u* impulse5'; %plot(movmean(HD(:,1),250)) works
OEDSVDFROMraw= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
OEDSVDFROMraw=OEDSVDFROMraw';
% HD=zscore(HD);
% HD1=zscore(HD1);
save OEDSVDTOraw.mat OEDSVDTOraw ;
save OEDSVDFROMraw.mat OEDSVDFROMraw ;
%%

clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
 CP1(:,32:36)=[];
CP1=CP1(:,[16 18 15 14 12 15]);
 %%
% CP1(:,32:35)=[]; % includes commodity only
%%
% CP1(:,32:34)=[];includes oil only
% CP1(:,33)=[]; % includes oil only
%% Data preprocessing
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
windowSize=100;
%% Inputs
y=filt;
    
    % Set SVAR parameters
    p = 2;      % Order of VAR lags
    q = 10;     % Order of VMA lags      

    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1);
    nconst = size(ylag,2);

    for i = 1:p;

        ylag = [ trimr(ylag,1,0) trimr(y,0,i) ];

    end
    
     % OLS Estimates for initialization
    bar    = ylag\trimr(y,p,0);
     v      = trimr(y,p,0) - ylag*bar;
% % % % % %     omegav = v'*v/size(v,1);
% % % % % %     s      = chol(omegav)';

   
    bstart = [
                0.001855486206481
                0.068693371469487
                0.023880592273368
                0.001989809184491
                0.013648672676092
                0.000281803828058
               -0.121999669995451
               -0.002888991832387
                0.002650363820455
                0.217906436715095
                0.001195294344330
                0.008613191798914
                0.175400639477836
                1.103000246244063
                0.440791629588056
                0.495605700916872
                0.365097267442524
               -0.011572135693424
                0.650892298805934
               -1.176100342250688
               -0.342100219212525
               -1.059396487450620
               -0.545902063123713
               -0.729900334997059
               -0.402000678739595
                1.013901399027622
               -1.601300249119193
                1.018606365787122
                0.791195906556810
                0.818998254746848
                0.338395871264844
                0.465287945621782
               -1.145684521466255
            ];

  options = optimset(   'Display',              'iter', ...
                        'MaxIter',              20000,   ...
                        'MaxFunEvals',          20000,   ...
                        'TolX',                 1e-8,   ...
                        'TolFun',               1e-4 );  
            
    [ b,fval]  = fminunc( @(b) neglogl(b,v),bstart,options ); 


    disp(' Log-likelihood function ');
    disp( fval );
    disp(' structural parameter ');
    disp( b );
    disp(' VAR disturbances ');
    disp( v );
    
   % Matrix required to compute time-varying variance-covariances              
    s = getSmatrix( b,v );
    omegavv = v'*v/size(v,1);
    omegavs= s'*s/size(s,1);
    omegav=omegavs(1:6,1:6);
    omegav=topdm(omegav);
   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    b = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = b(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);          
    % Construct C(L) matrices (vector moving average)  
    c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ]);

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end
       tmp = ss';
       c  = [ c  tmp(:) ];

    end
%%    
   ar0=a1;
   ar0={ar0};
   ma0=tmp;
   ma0={ma0};
   InnovCov=omegav; 
   GIRF = armairf(ar0,ma0,'Method','generalized','InnovCov',InnovCov,'NumObs',10);
   
    %% The old stuff begins here, lets see
 
impulse=GIRF;
impulse1=permute(impulse, [1 3 2]);
impulse2=reshape(impulse1, [], size(impulse, 2),1);
impulse2=impulse(:);
impulse3= repmat(impulse2,10);
u=reshapeg(v,[],6);
fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
%impulse4=trimr(impulse3,0,1242);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
OEESVDTOraw= u* impulse5'; %plot(movmean(HD(:,1),250)) works
OEESVDFROMraw= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
OEESVDFROMraw=OEESVDFROMraw';
% HD=zscore(HD);
% HD1=zscore(HD1);
save OEESVDTOraw.mat OEESVDTOraw ;
save OEESVDFROMraw.mat OEESVDFROMraw ;

load ACSVDTOraw;
load ECSVDTOraw;
load GCSVDTOraw;
load OEDSVDTOraw;
load OEESVDTOraw;
load ACSVDFROMraw;
load ECSVDFROMraw;
load GCSVDFROMraw;
load OEDSVDFROMraw;
load OEESVDFROMraw;

SVDtransmissionraw= [ACSVDTOraw(:,1:6) ECSVDTOraw(:,1:6) GCSVDTOraw(:,1:6) OEDSVDTOraw(:,1:6) OEESVDTOraw(:,1:6)];
SVDvulnerabilityraw= [ACSVDFROMraw ECSVDFROMraw GCSVDFROMraw OEDSVDFROMraw OEESVDFROMraw];
save SVDtransmissionwithoutzscore.mat SVDtransmissionraw;
save SVDvulnerabilitywithoutzscore.mat SVDvulnerabilityraw;





%
%--------------------------- Functions ----------------------------------
% 
%--------------------------------------------------------------------------
% Wrapper function for log-liklihood
%--------------------------------------------------------------------------
function lf = neglogl( b,v )

    lf = - mean( loglt( b,v ) );
end

%--------------------------------------------------------------------------
% Log-likelihood function for SVAR
%--------------------------------------------------------------------------

function lf = loglt( b,v )


    % Unpack parameter vector
    delta = b(16:21).^2;
    alpha = normcdf( b(22:27) );
    beta  = normcdf( b(28:33) );
    
    [ t,n ] = size( v );
    lf      = zeros( t,1 );
    
    
    binv  = [   1        0      0       0     0     0  
                b(1)     1      0       0     0     0
                b(2)   b(6)     1       0     0     0  
                b(3)   b(7)   b(10)     1     0     0  
                b(4)   b(8)   b(11)  b(13)    1     0  
                b(5)   b(9)   b(12)  b(14)  b(15)   1  ];

    u = ( inv(binv)*v' )';
            

    for i = 1:t;
        
        if i == 1;
            
            d = delta + alpha.*(std(u).^2)' + beta.*(std(u).^2)';
        else
            
            d = delta + alpha.*(u(i-1,:).^2)' + beta.*d;
        end
        
        omega = binv*diag(d)*binv';
        
        lf(i) = -0.5*n*log(2*pi) - 0.5*log(det(omega)) - 0.5*v(i,:)*inv(omega)*v(i,:)';
    end
    
    
end
%--------------------------------------------------------------------------
%   Return matrix to compute time-varying covariances
%--------------------------------------------------------------------------
function S = getSmatrix( b,v )


    % Unpack parameter vector
    delta = b(16:21).^2;
    alpha = normcdf( b(22:27) );
    beta  = normcdf( b(28:33) );
    
    [ t,n ] = size( v );
    lf      = zeros( t,1 );
    
    
    binv  = [   1        0      0       0     0     0  
                b(1)     1      0       0     0     0
                b(2)   b(6)     1       0     0     0  
                b(3)   b(7)   b(10)     1     0     0  
                b(4)   b(8)   b(11)  b(13)    1     0  
                b(5)   b(9)   b(12)  b(14)  b(15)   1  ];

    u = (inv(binv)*v')';
            
    for i = 1:t;
        
        if i == 1;
            
            d = delta + alpha.*(std(u).^2)' + beta.*(std(u).^2)';
        else
            
            d = delta + alpha.*(u(i-1,:).^2)' + beta.*d;
        end
       
        tmp    = binv*diag( sqrt(d) ); 
        
        S(i,:) = vec( tmp )';
              
    end
    
end

%--------------------------------------------------------------------------
% Reshape a matrix to agree with GAUSS reshape command
%--------------------------------------------------------------------------
function X = reshapeg(Y,r,c)

         tmp = reshape(Y',c,r);
         X   = tmp';

end


