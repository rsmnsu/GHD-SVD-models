clc
clear
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results_withoutzscore\GHDtransmissionwithoutzscore.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results_withoutzscore\GHDvulnerabilitywithoutzscore.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\NEW_SVD_experimenting\SVD_revised_corrected\SVDtransmissionwithoutzscore.mat');
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\NEW_SVD_experimenting\SVD_revised_corrected\SVDvulnerabilitywithoutzscore.mat');
ghdtoraw=movmean(GHDtransmissionraw,250);
ghdfromraw=movmean(GHDvulnerabilityraw,250);
svdtoraw=movmean(SVDtransmissionraw,250);
svdfromraw=movmean(SVDvulnerabilityraw,250);


%% Transmission
%% Plot with functions
%% AC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdtoraw(:,1),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdtoraw(:,1),'Parent',subplot1,'DisplayName','Volatility');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('India');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdtoraw(:,2),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdtoraw(:,2),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Malaysia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdtoraw(:,3),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdtoraw(:,3),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Singapore');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdtoraw(:,4),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdtoraw(:,4),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('South Korea');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdtoraw(:,5),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdtoraw(:,5),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('The Phillipines');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdtoraw(:,6),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdtoraw(:,6),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Thailand');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');


savefig('ac_transmission_GHD_SVD_v2.fig');

%% EC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdtoraw(:,7),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdtoraw(:,7),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Germany');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdtoraw(:,8),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdtoraw(:,8),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Chile');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdtoraw(:,9),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdtoraw(:,9),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('France');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdtoraw(:,10),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdtoraw(:,10),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('China');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdtoraw(:,11),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdtoraw(:,11),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('UK');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdtoraw(:,12),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdtoraw(:,12),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Australia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
savefig('ec_transmission_GHD_SVD_v2.fig');
%% GC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdtoraw(:,13),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdtoraw(:,13),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Greece');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdtoraw(:,14),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdtoraw(:,14),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Portugal');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdtoraw(:,15),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdtoraw(:,15),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Ireland');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdtoraw(:,16),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdtoraw(:,16),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Belgium');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdtoraw(:,17),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdtoraw(:,17),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Croatia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdtoraw(:,18),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdtoraw(:,18),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Austria');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');


savefig('gc_transmission_GHD_SVD_v2.fig');

%% OED
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdtoraw(:,19),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdtoraw(:,19),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('The USA');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdtoraw(:,20),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdtoraw(:,20),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Canada');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdtoraw(:,21),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdtoraw(:,21),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Russia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdtoraw(:,22),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdtoraw(:,22),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Norway');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdtoraw(:,23),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdtoraw(:,23),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Japan');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdtoraw(:,24),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdtoraw(:,24),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('New Zealand');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
savefig('oed_transmission_GHD_SVD_v2.fig');

%% OEE
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdtoraw(:,25),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdtoraw(:,25),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Saudi Arabia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdtoraw(:,26),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdtoraw(:,26),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Israel');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdtoraw(:,27),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdtoraw(:,27),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Iraq');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdtoraw(:,28),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdtoraw(:,28),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Kuwait');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdtoraw(:,29),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdtoraw(:,29),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Nigeria');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdtoraw(:,30),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdtoraw(:,30),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Venezuela');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
savefig('oee_transmission_GHD_SVD_v2.fig');


%% Vulnerability
%% Plot with functions
%% AC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdfromraw(:,1),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdfromraw(:,1),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('India');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdfromraw(:,2),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdfromraw(:,2),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Malaysia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdfromraw(:,3),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdfromraw(:,3),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Singapore');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdfromraw(:,4),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdfromraw(:,4),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('South Korea');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdfromraw(:,5),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdfromraw(:,5),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('The Phillipines');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdfromraw(:,6),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdfromraw(:,6),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Thailand');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');


savefig('ac_vulnerability_GHD_SVD_v2.fig');

%% EC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdfromraw(:,7),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdfromraw(:,7),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Germany');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdfromraw(:,8),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdfromraw(:,8),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Chile');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdfromraw(:,9),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdfromraw(:,9),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('France');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdfromraw(:,10),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdfromraw(:,10),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('China');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdfromraw(:,11),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdfromraw(:,11),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('UK');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdfromraw(:,12),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdfromraw(:,12),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Australia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
savefig('ec_vulnerability_GHD_SVD_v2.fig');
%% GC
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdfromraw(:,13),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdfromraw(:,13),'Parent',subplot1,'DisplayName','Vulnerability');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Greece');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdfromraw(:,14),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdfromraw(:,14),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Portugal');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdfromraw(:,15),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdfromraw(:,15),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Ireland');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdfromraw(:,16),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdfromraw(:,16),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Belgium');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdfromraw(:,17),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdfromraw(:,17),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Croatia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdfromraw(:,18),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdfromraw(:,18),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Austria');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');


savefig('gc_vulnerability_GHD_SVD_v2.fig');

%% OED
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdfromraw(:,19),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdfromraw(:,19),'Parent',subplot1,'DisplayName','Volatility');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('The USA');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdfromraw(:,20),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdfromraw(:,20),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Canada');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdfromraw(:,21),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdfromraw(:,21),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Russia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdfromraw(:,22),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdfromraw(:,22),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Norway');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdfromraw(:,23),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdfromraw(:,23),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Japan');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdfromraw(:,24),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdfromraw(:,24),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('New Zealand');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
savefig('oed_vulnerability_GHD_SVD_v2.fig');

%% OEE
figure1 = figure;

% Create subplot
subplot1 = subplot(2,3,1,'Parent',figure1);
hold(subplot1,'on');

% Activate the left side of the axes
yyaxis(subplot1,'left');
% Create plot
plot(ghdfromraw(:,25),'Parent',subplot1,'DisplayName','Return');

% Set the remaining axes properties
set(subplot1,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot1,'right');
% Create plot
plot(svdfromraw(:,25),'Parent',subplot1,'DisplayName','Volatility');

% Set the remaining axes properties
set(subplot1,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Saudi Arabia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5140]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.00993017153793759 0.616464259596181 0.0927601793772495 0.256528881844142]);

% Create subplot
subplot2 = subplot(2,3,2,'Parent',figure1);
hold(subplot2,'on');

% Activate the left side of the axes
yyaxis(subplot2,'left');
% Create plot
plot(ghdfromraw(:,26),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot2,'right');
% Create plot
plot(svdfromraw(:,26),'Parent',subplot2);

% Set the remaining axes properties
set(subplot2,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Israel');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5140]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot3 = subplot(2,3,3,'Parent',figure1);
hold(subplot3,'on');

% Activate the left side of the axes
yyaxis(subplot3,'left');
% Create plot
plot(ghdfromraw(:,27),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot3,'right');
% Create plot
plot(svdfromraw(:,27),'Parent',subplot3);

% Set the remaining axes properties
set(subplot3,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Iraq');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5140]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot4 = subplot(2,3,4,'Parent',figure1);
hold(subplot4,'on');

% Activate the left side of the axes
yyaxis(subplot4,'left');
% Create plot
plot(ghdfromraw(:,28),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot4,'right');
% Create plot
plot(svdfromraw(:,28),'Parent',subplot4);

% Set the remaining axes properties
set(subplot4,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Kuwait');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5140]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot5 = subplot(2,3,5,'Parent',figure1);
hold(subplot5,'on');

% Activate the left side of the axes
yyaxis(subplot5,'left');
% Create plot
plot(ghdfromraw(:,29),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0 0 0],'YMinorTick','off');
% Activate the right side of the axes
yyaxis(subplot5,'right');
% Create plot
plot(svdfromraw(:,29),'Parent',subplot5);

% Set the remaining axes properties
set(subplot5,'YColor',[0.85 0.325 0.098]);
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Nigeria');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5140]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
% Create subplot
subplot6 = subplot(2,3,6,'Parent',figure1);
hold(subplot6,'on');

% Activate the left side of the axes
yyaxis(subplot6,'left');
% Create plot
plot(ghdfromraw(:,30),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0 0 0]);
% Activate the right side of the axes
yyaxis(subplot6,'right');
% Create plot
plot(svdfromraw(:,30),'Parent',subplot6);

% Set the remaining axes properties
set(subplot6,'YColor',[0.85 0.325 0.098],'YMinorTick','off');
% Create xlabel
xlabel('Year','FontSize',9);

% Create title
title('Venezuela');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5140]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'XGrid','on','XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'},'YGrid','on');
savefig('oee_vulnerability_GHD_SVD_v2.fig');