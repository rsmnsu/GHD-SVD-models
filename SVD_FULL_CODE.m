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
    
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;
    omegav = v'*v/size(v,1);
    s      = chol(omegav)';


   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);           


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
                        'TolFun',               1e-4  );       
            
    [ b,fval ]  = fminunc( @(b) neglogl(b,v),bstart,options );     
  
    disp(' Log-likelihood function ');
    disp( fval );
    disp( b );
    
    % Matrix required to compute time-varying variance-covariances              
    s_store = getSmatrix( b,v );
    
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
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1);
save ACSVDTO.mat HD ;
save ACSVDFROM.mat HD1 ;


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
    
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;
    omegav = v'*v/size(v,1);
    s      = chol(omegav)';


   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);           


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
                        'TolFun',               1e-4  );       
            
    [ b,fval ]  = fminunc( @(b) neglogl(b,v),bstart,options );     
  
    disp(' Log-likelihood function ');
    disp( fval );
    disp( b );
    
    % Matrix required to compute time-varying variance-covariances              
    s_store = getSmatrix( b,v );
    
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
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1);
save ECSVDTO.mat HD ;
save ECSVDFROM.mat HD1 ;


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
    
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;
    omegav = v'*v/size(v,1);
    s      = chol(omegav)';


   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);           


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
                        'TolFun',               1e-4  );       
            
    [ b,fval ]  = fminunc( @(b) neglogl(b,v),bstart,options );     
  
    disp(' Log-likelihood function ');
    disp( fval );
    disp( b );
    
    % Matrix required to compute time-varying variance-covariances              
    s_store = getSmatrix( b,v );
    
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
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1);
save GCSVDTO.mat HD ;
save GCSVDFROM.mat HD1 ;


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
    
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;
    omegav = v'*v/size(v,1);
    s      = chol(omegav)';


   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);           


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
                        'TolFun',               1e-4  );       
            
    [ b,fval ]  = fminunc( @(b) neglogl(b,v),bstart,options );     
  
    disp(' Log-likelihood function ');
    disp( fval );
    disp( b );
    
    % Matrix required to compute time-varying variance-covariances              
    s_store = getSmatrix( b,v );
    
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
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1);
save OEDSVDTO.mat HD ;
save OEDVDFROM.mat HD1 ;
%%

clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
 CP1(:,32:36)=[];
CP1=CP1(:,[16 18 15 10 12 13]);
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
    
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;
    omegav = v'*v/size(v,1);
    s      = chol(omegav)';


   % Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(v,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
      
    % Invert A(1) matrix needed for the MLE
    a1inv = inv(a1);           


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
                        'TolFun',               1e-4  );       
            
    [ b,fval ]  = fminunc( @(b) neglogl(b,v),bstart,options );     
  
    disp(' Log-likelihood function ');
    disp( fval );
    disp( b );
    
    % Matrix required to compute time-varying variance-covariances              
    s_store = getSmatrix( b,v );
    
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
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1);
save OEESVDTO.mat HD ;
save OEESVDFROM.mat HD1 ;

load ACSVDTO.mat;
ACSVDTO= HD(:,1:6);
load ECSVDTO.mat;
ECSVDTO= HD(:,1:6);
load GCSVDTO.mat;
GCSVDTO=HD(:,1:6);
load OEDSVDTO.mat;
OEDSVDTO=HD(:,1:6);
load OEESVDTO.mat;
OEESVDTO=HD(:,1:6);

GHDtransmission= [ACSVDTO ECSVDTO GCSVDTO OEDSVDTO OEESVDTO];
save GHDTRANSMISSION.mat GHDtransmission;
load ACSVDFROM.mat;
AC=HD1';
load ECSVDFROM.mat;
EC=HD1';
load GCSVDFROM.mat;
GC=HD1';
load OEDSVDFROM.mat;
OED=HD1';
load OEESVDFROM.mat;
OEE=HD1';
GHDvulnerability= [AC EC GC OED OEE];
save GHDVULNERABILITY.mat GHDvulnerability;
%     
%        
%     % Construct variance decompositions 
%     tmp0   = reshapeg( cumsum( impulse.^2 ),q*k,k );
%     tmp1   = repmat( sum(tmp0,2),1,k );
%     decomp = reshapeg( 100*tmp0 ./ tmp1 , q , k^2 );
%  
%     format short
% %     disp('Variance decomposition (%) over time of libor in terms of factors');
% %     disp('*****************************************************************');
% % 	disp(' ' );
% %     disp('Period         y1         y2        y3        y4         y5       y6 ');   
% % 
% %     disp( [ (1:1:size(decomp,1))' decomp(:,k*k-(k-1):k*k) ] );
% %     
% %     % Graph the instantaneous effect of shocks on the libor over time  
%      rr            = size(v,1);
%      decomp_time2  = zeros( rr,k );           % Lag 2 (instantaneous)      
%      decomp_time5  = zeros( rr,k );           % Lag 5       
%      decomp_time10 = zeros( rr,k );           % Lag 20
%      for j = 1:rr
%   
%          s       = reshapeg(s_store(j,:),k,k)';
%          impulse = vec(s');
%  
%          for i = 2:q
%  
%              tmp      = ( reshapeg( c(:,i),k,k )*s )';
%              impulse  = [ impulse tmp(:) ];
%  
%          end 
%          impulse = impulse';    
%          tmp0    = reshapeg( cumsum( impulse.^2 ),q*k,k );
%          tmp1    = repmat( sum(tmp0,2),1,k );
%          decomp  = reshapeg( 100*tmp0 ./ tmp1 , q , k^2 );
%  
%          decomp_time2(j,:) = decomp(2,k*k-(k-1):k*k); 
%          decomp_time5(j,:) = decomp(5,k*k-(k-1):k*k); 
%          decomp_time10(j,:) = decomp(10,k*k-(k-1):k*k); 
%  
%      end   
%      
% %    dlmwrite('decomptime2.csv',decomp_time2);
% %    dlmwrite('decomptime10.csv',decomp_time10);
%      plot(movmean(HD(:,1),250))
%      hold on
%      plot(decomp_time10(:,1))
%      hold off
% %     
%     %**********************************************************************
%     %***
%     %***     Generate graph
%     %***
%     %**********************************************************************
% 
%     % Switch off TeX interpreter and clear figure
%     set(0,'defaulttextinterpreter','none');
%     figure(1);
%     clf;
%     
%     t = (2004:1/249:2008.61)';
% 
%     %--------------------------------------------------------%
%     % Panel (a)
%     subplot(3,3,1)
%     plot(t,decomp_time1(:,1),'-k','LineWidth',1);
%     title('Global Risk Factor');
%     ylabel('One Period Lag');
%     %ylim( [0 100] )
%     xlabel('t');
%     box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
%     
%     %--------------------------------------------------------%
%     % Panel (b)
%     subplot(3,3,2)
%     plot(t,decomp_time1(:,3),'-k','LineWidth',1);
%     title('Broad Liquidity Factor');
%     %ylabel('Oil price');
%     %ylim( [0 100] )
%     xlabel('t');
%     box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
% 
% 
%     %--------------------------------------------------------%
%     % Panel (c)
%     subplot(3,3,3)
%     plot(t,decomp_time1(:,6),'-k','LineWidth',1);
%     title('Idiosyncratic Factor');
%     %ylabel('Oil price');
%     %ylim( [0 100] )
%     xlabel('t');
%     box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
% 
% 
%     %--------------------------------------------------------%
%     % Panel (d)
%     subplot(3,3,4)
%     plot(t,decomp_time5(:,1),'-k','LineWidth',1);
%     %title('Money Shock');
%     ylabel('Five Period Lag');
%     xlabel('t');
%     box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
% 
%     %--------------------------------------------------------%
%     % Panel (e)
%     subplot(3,3,5)
%     plot(t,decomp_time5(:,3),'-k','LineWidth',1);
%     %ylabel('Output');
%     xlabel('t');
%     box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
% 
% 
%     %--------------------------------------------------------%
%     % Panel (f)
%     subplot(3,3,6)
%     plot(t,decomp_time5(:,6),'-k','LineWidth',1);
%     %ylabel('Output');
%     xlabel('t');
%     box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
% 
%     %--------------------------------------------------------%
%     % Panel (g)
%     subplot(3,3,7)
%     plot(t,decomp_time20(:,1),'-k','LineWidth',1);
%     ylabel('Twenty Period Lag');
%     xlabel('t');
%      box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
% 
%     %--------------------------------------------------------%
%     % Panel (h)
%     subplot(3,3,8)
%     plot(t,decomp_time20(:,3),'-k','LineWidth',1);
%     %ylabel('Output');
%     xlabel('t');
%     box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
% 
%     %--------------------------------------------------------%
%     % Panel (i)
%     subplot(3,3,9)
%     plot(t,decomp_time20(:,6),'-k','LineWidth',1);
%     ylabel('Price');
%     xlabel('t');
%     box off
%     axis tight
%     set(gca,'YLimMode','manual')
%     set(gca,'Ylim',[0 100]);
%     set(gca,'YTick',[0 20 40 60 80 100] );
%     set(gca,'XTick',[2004 2006 2008] );
%     
%     %--------------------------------------------------------%
%     % Print the tex file to the relevant directory
%     laprint(1,'decomp','options','factory');



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


