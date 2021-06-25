clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
CP1(:,32:36)=[];

% ind=GHDTO(:,3);
% mys=GHDTO(:,5);
% sgp=GHDTO(:,7);
% kor=GHDTO(:,9);
% phl=GHDTO(:,8);
% tha=GHDTO(:,11);4

CP1=CP1(:,[3 5 7 9 8 11]);
%%
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
y = filt; 
 p = 2;      % Order of VAR lags
    q = 2;     % Order of VMA lags  
    %%
    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1) ;
    nconst = size(ylag,2);
    for K = 1:p
        ylag = [ trimr(ylag,1,0) trimr(y,0,K) ];
    end
    %%
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;   
    omegav = v'*v/size(v,1);
%%
ncol = size(v,2);

u(:,1) = v(:,1);
a = zeros(ncol, ncol);
for K = 2 : size(v,2)
  a(1:K-1,K) = v(:,1:K-1)\v(:,K);
  u(:,K) = v(:,K) - v(:,1:K-1)*a(1:K-1,K);
end

b = diag(ones(1,ncol)) - tril(a.');

residual=u;
%%
% Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(u,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
%     s=chol(omegav);
  %%
    % Construct C(L) matrices (vector moving average) 
       c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ])

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end

       tmp = ss';
       c = [ c  tmp(:) ];

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

fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
u=reshapeg(u,[],6);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1');
save ACGHDTO.mat HD ;
save ACGHDFROM.mat HD1 ;

% GHDTO=zscore(HD);
% GHDFROM=zscore(HD1);

%%
clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
CP1(:,32:36)=[];

% ger=GHDTO(:,28);
% chl=GHDTO(:,29);
% fra=GHDTO(:,31);
% chn=GHDTO(:,17);
% uk=GHDTO(:,30);
% aus=GHDTO(:,2);

CP1=CP1(:,[28 29 31 17 30 2]);
%%
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
y = filt; 
 p = 2;      % Order of VAR lags
    q = 2;     % Order of VMA lags  
    %%
    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1) ;
    nconst = size(ylag,2);
    for K = 1:p
        ylag = [ trimr(ylag,1,0) trimr(y,0,K) ];
    end
    %%
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;   
    omegav = v'*v/size(v,1);
%%
ncol = size(v,2);

u(:,1) = v(:,1);
a = zeros(ncol, ncol);
for K = 2 : size(v,2)
  a(1:K-1,K) = v(:,1:K-1)\v(:,K);
  u(:,K) = v(:,K) - v(:,1:K-1)*a(1:K-1,K);
end

b = diag(ones(1,ncol)) - tril(a.');

residual=u;
%%
% Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(u,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
%     s=chol(omegav);
  %%
    % Construct C(L) matrices (vector moving average) 
       c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ])

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end

       tmp = ss';
       c = [ c  tmp(:) ];

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

fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
u=reshapeg(u,[],6);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1');
save ECGHDTO.mat HD ;
save ECGHDFROM.mat HD1 ;

% GHDTO=zscore(HD);
% GHDFROM=zscore(HD1);

%%
%%
clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
CP1(:,32:36)=[];

% grc=GHDTO(:,20);
% prt=GHDTO(:,21);
% irl=GHDTO(:,22);
% bel=GHDTO(:,23);
% crt=GHDTO(:,24);
% aut=GHDTO(:,25);

CP1=CP1(:,[20 21 22 23 24 25]);
%%
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
y = filt; 
 p = 2;      % Order of VAR lags
    q = 2;     % Order of VMA lags  
    %%
    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1) ;
    nconst = size(ylag,2);
    for K = 1:p
        ylag = [ trimr(ylag,1,0) trimr(y,0,K) ];
    end
    %%
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;   
    omegav = v'*v/size(v,1);
%%
ncol = size(v,2);

u(:,1) = v(:,1);
a = zeros(ncol, ncol);
for K = 2 : size(v,2)
  a(1:K-1,K) = v(:,1:K-1)\v(:,K);
  u(:,K) = v(:,K) - v(:,1:K-1)*a(1:K-1,K);
end

b = diag(ones(1,ncol)) - tril(a.');

residual=u;
%%
% Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(u,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
%     s=chol(omegav);
  %%
    % Construct C(L) matrices (vector moving average) 
       c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ])

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end

       tmp = ss';
       c = [ c  tmp(:) ];

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

fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
u=reshapeg(u,[],6);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1');
save GCGHDTO.mat HD ;
save GCGHDFROM.mat HD1 ;

% GHDTO=zscore(HD);
% GHDFROM=zscore(HD1);
%%
clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
CP1(:,32:36)=[];

% usa=GHDTO(:,1);
% cad=GHDTO(:,19);
% rus=GHDTO(:,26);
% nor=GHDTO(:,27);
% jap=GHDTO(:,4);
% nzl=GHDTO(:,6);

CP1=CP1(:,[1 19 26 27 4 6]);
%%
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
y = filt; 
 p = 2;      % Order of VAR lags
    q = 2;     % Order of VMA lags  
    %%
    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1) ;
    nconst = size(ylag,2);
    for K = 1:p
        ylag = [ trimr(ylag,1,0) trimr(y,0,K) ];
    end
    %%
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;   
    omegav = v'*v/size(v,1);
%%
ncol = size(v,2);

u(:,1) = v(:,1);
a = zeros(ncol, ncol);
for K = 2 : size(v,2)
  a(1:K-1,K) = v(:,1:K-1)\v(:,K);
  u(:,K) = v(:,K) - v(:,1:K-1)*a(1:K-1,K);
end

b = diag(ones(1,ncol)) - tril(a.');

residual=u;
%%
% Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(u,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
%     s=chol(omegav);
  %%
    % Construct C(L) matrices (vector moving average) 
       c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ])

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end

       tmp = ss';
       c = [ c  tmp(:) ];

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

fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
u=reshapeg(u,[],6);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1');
save OEDGHDTO.mat HD ;
save OEDGHDFROM.mat HD1 ;

% GHDTO=zscore(HD);
% GHDFROM=zscore(HD1);

%%
%%
clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\MyTesting.mat');
CP1=MyTesting;
CP1(:,14)=[];
CP1(:,32:36)=[];

% sau=GHDTO(:,16);
% isr=GHDTO(:,18);
% irq=GHDTO(:,15);
% slk=GHDTO(:,10);
% nig=GHDTO(:,12);
% ven=GHDTO(:,13);

CP1=CP1(:,[16 18 15 14 12 13]);
%%
CPNew = log(CP1); 
CPNew= 100*(trimr(CPNew,1,0)-trimr(CPNew,0,1));
window_Size = 2;
delta = (1/window_Size)*ones(1,window_Size);
gama=1;
filt=filter(delta,gama,CPNew);
[nRow, nCol] = size(filt);
y = filt; 
 p = 2;      % Order of VAR lags
    q = 2;     % Order of VMA lags  
    %%
    % Estimate the VAR with p lags and a constant     
    ylag   = ones(length(y),1) ;
    nconst = size(ylag,2);
    for K = 1:p
        ylag = [ trimr(ylag,1,0) trimr(y,0,K) ];
    end
    %%
    % OLS Estimates
    bar    = ylag\trimr(y,p,0);
    %mue    = bar(1:nconst,:);
    v      = trimr(y,p,0) - ylag*bar;   
    omegav = v'*v/size(v,1);
%%
ncol = size(v,2);

u(:,1) = v(:,1);
a = zeros(ncol, ncol);
for K = 2 : size(v,2)
  a(1:K-1,K) = v(:,1:K-1)\v(:,K);
  u(:,K) = v(:,K) - v(:,1:K-1)*a(1:K-1,K);
end

b = diag(ones(1,ncol)) - tril(a.');

residual=u;
%%
% Constuct A(L) (VAR) matrices and hence the C(1) long-run matrix
    bar = trimr(bar,nconst,0);
    k   = size(u,2);
   	a   = zeros(k^2,p);
    a1  = eye(k);

    for i = 1:p
        
        tmp    = bar(1+k*(i-1):i*k,:);
    	a(:,i) = tmp(:);
	    a1     = a1 - reshapeg(a(:,i),k,k);
    end
%     s=chol(omegav);
  %%
    % Construct C(L) matrices (vector moving average) 
       c = eye(k);
    c = c(:);

    for i = 1:q

       ss = zeros(k,k);
       j   = 1.0;
        
       while j <= min( [ p i ])

          ss = ss + reshapeg(a(:,j),k,k)*reshapeg(c(:,i-j+1),k,k);
          j   = j + 1;
       end

       tmp = ss';
       c = [ c  tmp(:) ];

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

fixing= length(impulse3(:))-length(u(:));
impulse3=impulse3(:);
impulse4=trimr(impulse3,0,fixing);
u=reshapeg(u,[],6);
impulse5=reshapeg(impulse4,[],6);
%% experiment with To and From
HD = u* impulse5'; %plot(movmean(HD(:,1),250)) works
HD1= u'.*impulse5'; %plot(movmean(HD1(1,:),250)) works
HD=zscore(HD);
HD1=zscore(HD1');
save OEEGHDTO.mat HD ;
save OEEGHDFROM.mat HD1 ;

load ACGHDTO.mat;
ACGHDTO= HD(:,1:6);
load ECGHDTO.mat;
ECGHDTO= HD(:,1:6);
load GCGHDTO.mat;
GCGHDTO=HD(:,1:6);
load OEDGHDTO.mat;
OEDGHDTO=HD(:,1:6);
load OEEGHDTO.mat;
OEEGHDTO=HD(:,1:6);

GHDtransmission= [ACGHDTO ECGHDTO GCGHDTO OEDGHDTO OEEGHDTO];
save GHDTRANSMISSION.mat GHDtransmission;
load ACGHDFROM.mat;
AC=HD1';
load ECGHDFROM.mat;
EC=HD1';
load GCGHDFROM.mat;
GC=HD1';
load OEDGHDFROM.mat;
OED=HD1';
load OEEGHDFROM.mat;
OEE=HD1';
GHDvulnerability= [AC EC GC OED OEE];
save GHDVULNERABILITY.mat GHDvulnerability;

% GHDTO=zscore(HD);
% GHDFROM=zscore(HD1);

