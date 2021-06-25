clc
clear
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results\GHDTRANSMISSION.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results\SVDTRANSMISSION.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\Only GHD-SVD vulnerability results\GHDVULNERABILITY.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\Only GHD-SVD vulnerability results\SVD_vulnerability\SVDVULNERABILITY.');
ghdto=movmean(GHDtransmission,250);
ghdfrom=movmean(GHDvulnerability,250);
svito=movmean(SVDtransmission,250);
svifrom=movmean(SVDvulnerability,250);
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results_withoutzscore\GHDtransmissionwithoutzscore.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results_withoutzscore\SVDtransmissionwithoutzscore.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results_withoutzscore\GHDvulnerabilitywithoutzscore.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results_withoutzscore\SVDvulnerabilitywithoutzscore.mat');
ghdtoraw=movmean(GHDtransmissionraw,250);
ghdfromraw=movmean(GHDvulnerabilityraw,250);
svitoraw=movmean(SVDtransmissionraw,250);
svifromraw=movmean(SVDvulnerabilityraw,250);
%% Transmission
%% AC
%India
subplot (2,3,1)
plot(ghdtoraw(:,1))
hold on
%yyaxis right
plot(svitoraw(:,1))
hold off
%Malaysia
subplot (2,3,2)
plot(ghdtoraw(:,2))
hold on
%yyaxis right
plot(svitoraw(:,2))
hold off
%Singapore
subplot (2,3,3)
plot(ghdtoraw(:,3))
hold on
%yyaxis right
plot(svitoraw(:,3))
hold off
%S.korea
subplot(2,3,4)
plot(ghdtoraw(:,4))
hold on
%yyaxis right
plot(svitoraw(:,4))
hold off
%Phillipines
subplot(2,3,5)
plot(ghdtoraw(:,5))
hold on
%yyaxis right
plot(svitoraw(:,5))
hold off
%Thailand
subplot(2,3,6)
plot(ghdtoraw(:,6))
hold on
%yyaxis right
plot(svitoraw(:,6))
hold off
%% EC
% Germany
subplot (2,3,1)
plot(ghdtoraw(:,7))
hold on
plot(svitoraw(:,7))
hold off
% Chile
subplot (2,3,2)
plot(ghdtoraw(:,8))
hold on
plot(svitoraw(:,8))
hold off
%France
subplot (2,3,3)
plot(ghdtoraw(:,9))
hold on
plot(svitoraw(:,9))
hold off
% China
subplot(2,3,4)
plot(ghdtoraw(:,10))
hold on
plot(svitoraw(:,10))
hold off
% UK
subplot(2,3,5)
plot(ghdtoraw(:,11))
hold on
plot(svitoraw(:,11))
hold off
% Australia
subplot(2,3,6)
plot(ghdtoraw(:,12))
hold on
plot(svitoraw(:,12))
hold off
%% GC
% Greece
subplot (2,3,1)
plot(ghdtoraw(:,13))
hold on
plot(svitoraw(:,13))
hold off
% Portugal
subplot (2,3,2)
plot(ghdtoraw(:,14))
hold on
plot(svitoraw(:,14))
hold off
% Ireland
subplot (2,3,3)
plot(ghdtoraw(:,15))
hold on
plot(svitoraw(:,15))
hold off
% Belgium
subplot(2,3,4)
plot(ghdtoraw(:,16))
hold on
plot(svitoraw(:,16))
hold off
% Croatia
subplot(2,3,5)
plot(ghdtoraw(:,17))
hold on
plot(svitoraw(:,17))
hold off
% Austria
subplot(2,3,6)
plot(ghdtoraw(:,18))
hold on
plot(svitoraw(:,18))
hold off
%% OED
% The USA
subplot (2,3,1)
plot(ghdtoraw(:,19))
hold on
plot(svitoraw(:,19))
hold off
% Canada
subplot (2,3,2)
plot(ghdtoraw(:,20))
hold on
plot(svitoraw(:,20))
hold off
% Russia
subplot (2,3,3)
plot(ghdtoraw(:,21))
hold on
plot(svitoraw(:,21))
hold off
% Norway
subplot(2,3,4)
plot(ghdtoraw(:,22))
hold on
plot(svitoraw(:,22))
hold off
% Japan
subplot(2,3,5)
plot(ghdtoraw(:,23))
hold on
plot(svitoraw(:,23))
hold off
% New Zealand
subplot(2,3,6)
plot(ghdtoraw(:,24))
hold on
plot(svitoraw(:,24))
hold off

%% OEE
% Saudi Arabia
subplot (2,3,1)
plot(ghdtoraw(:,25))
hold on
plot(svitoraw(:,25))
hold off
% Israel
subplot (2,3,2)
plot(ghdtoraw(:,26))
hold on
plot(svitoraw(:,26))
hold off
% Iraq
subplot (2,3,3)
plot(ghdtoraw(:,27))
hold on
plot(svitoraw(:,27))
hold off
% Sri Lanka
subplot(2,3,4)
plot(ghdtoraw(:,28))
hold on
plot(svitoraw(:,28))
hold off
% Nigeria
subplot(2,3,5)
plot(ghdtoraw(:,29))
hold on
plot(svitoraw(:,29))
hold off
% Venezuela
subplot(2,3,6)
plot(ghdtoraw(:,30))
hold on
plot(svitoraw(:,30))
hold off
%% Vulnerability
%% AC
%India
subplot (2,3,1)
plot(ghdfromraw(:,1))
hold on
%yyaxis right
plot(svifromraw(:,1))
hold off
%Malaysia
subplot (2,3,2)
plot(ghdfromraw(:,2))
hold on
%yyaxis right
plot(svifromraw(:,2))
hold off
%Singapore
subplot (2,3,3)
plot(ghdfromraw(:,3))
hold on
%yyaxis right
plot(svifromraw(:,3))
hold off
%S.korea
subplot(2,3,4)
plot(ghdfromraw(:,4))
hold on
%yyaxis right
plot(svifromraw(:,4))
hold off
%Phillipines
subplot(2,3,5)
plot(ghdfromraw(:,5))
hold on
%yyaxis right
plot(svifromraw(:,5))
hold off
%Thailand
subplot(2,3,6)
plot(ghdfromraw(:,6))
hold on
%yyaxis right
plot(svifromraw(:,6))
hold off
%% EC
% Germany
subplot (2,3,1)
plot(ghdfromraw(:,7))
hold on
plot(svifromraw(:,7))
hold off
% Chile
subplot (2,3,2)
plot(ghdfromraw(:,8))
hold on
plot(svifromraw(:,8))
hold off
%France
subplot (2,3,3)
plot(ghdfromraw(:,9))
hold on
plot(svifromraw(:,9))
hold off
% China
subplot(2,3,4)
plot(ghdfromraw(:,10))
hold on
plot(svifromraw(:,10))
hold off
% UK
subplot(2,3,5)
plot(ghdfromraw(:,11))
hold on
plot(svifromraw(:,11))
hold off
% Australia
subplot(2,3,6)
plot(ghdfromraw(:,12))
hold on
plot(svifromraw(:,12))
hold off
%% GC
% Greece
subplot (2,3,1)
plot(ghdfromraw(:,13))
hold on
plot(svifromraw(:,13))
hold off
% Portugal
subplot (2,3,2)
plot(ghdfromraw(:,14))
hold on
plot(svifromraw(:,14))
hold off
% Ireland
subplot (2,3,3)
plot(ghdfromraw(:,15))
hold on
plot(svifromraw(:,15))
hold off
% Belgium
subplot(2,3,4)
plot(ghdfromraw(:,16))
hold on
plot(svifromraw(:,16))
hold off
% Croatia
subplot(2,3,5)
plot(ghdfromraw(:,17))
hold on
plot(svifromraw(:,17))
hold off
% Austria
subplot(2,3,6)
plot(ghdfromraw(:,18))
hold on
plot(svifromraw(:,18))
hold off
%% OED
% The USA
subplot (2,3,1)
plot(ghdfromraw(:,19))
hold on
plot(svifromraw(:,19))
hold off
% Canada
subplot (2,3,2)
plot(ghdfromraw(:,20))
hold on
plot(svifromraw(:,20))
hold off
% Russia
subplot (2,3,3)
plot(ghdfromraw(:,21))
hold on
plot(svifromraw(:,21))
hold off
% Norway
subplot(2,3,4)
plot(ghdfromraw(:,22))
hold on
plot(svifromraw(:,22))
hold off
% Japan
subplot(2,3,5)
plot(ghdfromraw(:,23))
hold on
plot(svifromraw(:,23))
hold off
% New Zealand
subplot(2,3,6)
plot(ghdfromraw(:,24))
hold on
plot(svifromraw(:,24))
hold off

%% OEE
% Saudi Arabia
subplot (2,3,1)
plot(ghdfromraw(:,25))
hold on
plot(svifromraw(:,25))
hold off
% Israel
subplot (2,3,2)
plot(ghdfromraw(:,26))
hold on
plot(svifromraw(:,26))
hold off
% Iraq
subplot (2,3,3)
plot(ghdfromraw(:,27))
hold on
plot(svifromraw(:,27))
hold off
% Sri Lanka
subplot(2,3,4)
plot(ghdfromraw(:,28))
hold on
plot(svifromraw(:,28))
hold off
% Nigeria
subplot(2,3,5)
plot(ghdfromraw(:,29))
hold on
plot(svifromraw(:,29))
hold off
% Venezuela
subplot(2,3,6)
plot(ghdfromraw(:,30))
hold on
plot(svifromraw(:,30))
hold off