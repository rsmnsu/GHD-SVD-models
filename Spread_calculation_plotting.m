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
ghdto=zscore(ghdtoraw);
ghdfrom=zscore(ghdfromraw);
svdto=zscore(svdtoraw);
svdfrom=zscore(svdfromraw);
transmissionspread= svdto-ghdto;
vulnerabilityspread= svdfrom-ghdfrom;
%% plot
%CREATEFIGURE(Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, Y10, Y11, Y12, Y13, Y14, Y15, Y16, Y17, Y18, Y19, Y20, Y21, Y22, Y23, Y24, Y25, Y26, Y27, Y28, Y29, Y30)
% Y1= [transmissionspread(:,1) vulnerabilityspread(:,1)];


vars = {'Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Y6', 'Y7','Y8', 'Y9', 'Y10', 'Y11', 'Y12', 'Y13', 'Y14', 'Y15', 'Y16', 'Y17', 'Y18', 'Y19', 'Y20', 'Y21', 'Y22', 'Y23', 'Y24', 'Y25', 'Y26', 'Y27', 'Y28', 'Y29', 'Y30'};
for i=1:30
eval([vars{i} '=  [transmissionspread(:,i) vulnerabilityspread(:,i)]';])
end

%  Auto-generated by MATLAB on 28-Dec-2018 05:04:41

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1]);

% Create subplot
subplot1 = subplot(5,6,1);
hold(subplot1,'on');

% Create plot
plot(Y1);

% Create title
title('India');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot1,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot1,[0.9 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot1,[-1 1]);
box(subplot1,'on');
% Set the remaining axes properties
set(subplot1,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot2 = subplot(5,6,2);
hold(subplot2,'on');

% Create plot
plot(Y2);

% Create title
title('Malaysia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot2,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot2,[0.6 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot2,[-1 1]);
box(subplot2,'on');
% Set the remaining axes properties
set(subplot2,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot3 = subplot(5,6,3);
hold(subplot3,'on');

% Create plot
plot(Y3);

% Create title
title('Singapore');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot3,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot3,[0.4 0.8]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot3,[-1 1]);
box(subplot3,'on');
% Set the remaining axes properties
set(subplot3,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot4 = subplot(5,6,4);
hold(subplot4,'on');

% Create plot
plot(Y4);

% Create title
title('South Korea');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot4,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot4,[0 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot4,[-1 1]);
box(subplot4,'on');
% Set the remaining axes properties
set(subplot4,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5140],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot5 = subplot(5,6,5);
hold(subplot5,'on');

% Create plot
plot(Y5);

% Create title
title('The Phillipines');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot5,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot5,[0.6 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot5,[-1 1]);
box(subplot5,'on');
% Set the remaining axes properties
set(subplot5,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot6 = subplot(5,6,6);
hold(subplot6,'on');

% Create plot
plot(Y6);

% Create title
title('Thailand');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot6,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot6,[0.5 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot6,[-1 1]);
box(subplot6,'on');
% Set the remaining axes properties
set(subplot6,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot7 = subplot(5,6,7);
hold(subplot7,'on');

% Create plot
plot(Y7);

% Create title
title('Germany');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot7,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot7,[0.95 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot7,[-1 1]);
box(subplot7,'on');
% Set the remaining axes properties
set(subplot7,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot8 = subplot(5,6,8);
hold(subplot8,'on');

% Create plot
plot(Y8);

% Create title
title('Chile');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot8,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot8,[0.5 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot8,[-1 1]);
box(subplot8,'on');
% Set the remaining axes properties
set(subplot8,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot9 = subplot(5,6,9);
hold(subplot9,'on');

% Create plot
plot(Y9);

% Create title
title('France');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot9,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot9,[0.9 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot9,[-1 1]);
box(subplot9,'on');
% Set the remaining axes properties
set(subplot9,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot10 = subplot(5,6,10);
hold(subplot10,'on');

% Create plot
plot(Y10);

% Create title
title('China');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot10,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot10,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot10,[-1 1]);
box(subplot10,'on');
% Set the remaining axes properties
set(subplot10,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot11 = subplot(5,6,11);
hold(subplot11,'on');

% Create plot
plot(Y11);

% Create title
title('UK');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot11,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot11,[0.5 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot11,[-1 1]);
box(subplot11,'on');
% Set the remaining axes properties
set(subplot11,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot12 = subplot(5,6,12);
hold(subplot12,'on');

% Create plot
plot(Y12);

% Create title
title('Australia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot12,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot12,[0.9 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot12,[-1 1]);
box(subplot12,'on');
% Set the remaining axes properties
set(subplot12,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot13 = subplot(5,6,13);
hold(subplot13,'on');

% Create plot
plot(Y13);

% Create title
title('Greece');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot13,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot13,[0 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot13,[-1 1]);
box(subplot13,'on');
% Set the remaining axes properties
set(subplot13,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot14 = subplot(5,6,14);
hold(subplot14,'on');

% Create plot
plot(Y14);

% Create title
title('Portugal');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot14,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot14,[0.5 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot14,[-1 1]);
box(subplot14,'on');
% Set the remaining axes properties
set(subplot14,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot15 = subplot(5,6,15);
hold(subplot15,'on');

% Create plot
plot(Y15);

% Create title
title('Ireland');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot15,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot15,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot15,[-1 1]);
box(subplot15,'on');
% Set the remaining axes properties
set(subplot15,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot16 = subplot(5,6,16);
hold(subplot16,'on');

% Create plot
plot(Y16);

% Create title
title('Belgium');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot16,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot16,[0.5 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot16,[-1 1]);
box(subplot16,'on');
% Set the remaining axes properties
set(subplot16,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot17 = subplot(5,6,17);
hold(subplot17,'on');

% Create plot
plot(Y17);

% Create title
title('Croatia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot17,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot17,[0.6 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot17,[-1 1]);
box(subplot17,'on');
% Set the remaining axes properties
set(subplot17,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot18 = subplot(5,6,18);
hold(subplot18,'on');

% Create plot
plot(Y18);

% Create title
title('Austria');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot18,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot18,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot18,[-1 1]);
box(subplot18,'on');
% Set the remaining axes properties
set(subplot18,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot19 = subplot(5,6,19);
hold(subplot19,'on');

% Create plot
plot(Y19);

% Create title
title('USA');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot19,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot19,[0.6 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot19,[-1 1]);
box(subplot19,'on');
% Set the remaining axes properties
set(subplot19,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot20 = subplot(5,6,20);
hold(subplot20,'on');

% Create plot
plot(Y20);

% Create title
title('Canada');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot20,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot20,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot20,[-1 1]);
box(subplot20,'on');
% Set the remaining axes properties
set(subplot20,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot21 = subplot(5,6,21);
hold(subplot21,'on');

% Create plot
plot(Y21);

% Create title
title('Russia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot21,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot21,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot21,[-1 1]);
box(subplot21,'on');
% Set the remaining axes properties
set(subplot21,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot22 = subplot(5,6,22);
hold(subplot22,'on');

% Create plot
plot(Y22);

% Create title
title('Norway');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot22,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot22,[0.6 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot22,[-1 1]);
box(subplot22,'on');
% Set the remaining axes properties
set(subplot22,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot23 = subplot(5,6,23);
hold(subplot23,'on');

% Create plot
plot(Y23);

% Create title
title('Japan');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot23,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot23,[0.6 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot23,[-1 1]);
box(subplot23,'on');
% Set the remaining axes properties
set(subplot23,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot24 = subplot(5,6,24);
hold(subplot24,'on');

% Create plot
plot(Y24);

% Create title
title('New Zealand');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot24,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot24,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot24,[-1 1]);
box(subplot24,'on');
% Set the remaining axes properties
set(subplot24,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot25 = subplot(5,6,25);
hold(subplot25,'on');

% Create plot
plot(Y25);

% Create title
title('Saudi Arabia');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot25,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot25,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot25,[-1 1]);
box(subplot25,'on');
% Set the remaining axes properties
set(subplot25,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot26 = subplot(5,6,26);
hold(subplot26,'on');

% Create plot
plot(Y26);

% Create title
title('Israel');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot26,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot26,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot26,[-1 1]);
box(subplot26,'on');
% Set the remaining axes properties
set(subplot26,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot27 = subplot(5,6,27);
hold(subplot27,'on');

% Create plot
plot(Y27);

% Create title
title('Iraq');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot27,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot27,[0 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot27,[-1 1]);
box(subplot27,'on');
% Set the remaining axes properties
set(subplot27,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot28 = subplot(5,6,28);
hold(subplot28,'on');

% Create plot
plot(Y28);

% Create title
title('Sri Lanka');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot28,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot28,[0.5 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot28,[-1 1]);
box(subplot28,'on');
% Set the remaining axes properties
set(subplot28,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot29 = subplot(5,6,29);
hold(subplot29,'on');

% Create plot
plot(Y29);

% Create title
title('Nigeria');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot29,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot29,[0.6 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot29,[-1 1]);
box(subplot29,'on');
% Set the remaining axes properties
set(subplot29,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});
% Create subplot
subplot30 = subplot(5,6,30);
hold(subplot30,'on');

% Create plot
plot(Y30);

% Create title
title('Venezuela');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(subplot30,[0 5130]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(subplot30,[0.8 1]);
% Uncomment the following line to preserve the Z-limits of the axes
% zlim(subplot30,[-1 1]);
box(subplot30,'on');
% Set the remaining axes properties
set(subplot30,'FontSize',7,'XTick',...
    [0 514 1028 1542 2056 2570 3084 3598 4112 4626 5130],'XTickLabel',...
    {'98','00','02','04','06','08','10','12','14','16','18'});

