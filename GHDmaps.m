%load GHDTRANSMISSION.mat;
load('F:\OneDrive - University of Tasmania\Structuring_thesis\Chapter1\ONLY GHD-SVD transmission results\GHDTRANSMISSION.mat');
%% begins here is ANN

X=GHDtransmission;
N=40;
R=size(X,1);
pos=floor(linspace(1,R,N+1));
pos(end)=R+1;
sizes=diff(pos);
subsets=mat2cell(X,sizes,size(X,2));
A1=subsets{1,1};
A2=subsets{2,1};
A3=subsets{3,1};
A4=subsets{4,1};
A5=subsets{5,1};
A6=subsets{6,1};
A7=subsets{7,1};
A8=subsets{8,1};
A9=subsets{9,1};
A10=subsets{10,1};
A11=subsets{11,1};
A12=subsets{12,1};
A13=subsets{13,1};
A14=subsets{14,1};
A15=subsets{15,1};
A16=subsets{16,1};
A17=subsets{17,1};
A18=subsets{18,1};
A19=subsets{19,1};
A20=subsets{20,1};
A21=subsets{21,1};
A22=subsets{22,1};
A23=subsets{23,1};
A24=subsets{24,1};
A25=subsets{25,1};
A26=subsets{26,1};
A27=subsets{27,1};
A28=subsets{28,1};
A29=subsets{29,1};
A30=subsets{30,1};
A31=subsets{31,1};
A32=subsets{32,1};
A33=subsets{33,1};
A34=subsets{34,1};
A35=subsets{35,1};
A36=subsets{36,1};
A37=subsets{37,1};
A38=subsets{38,1};
A39=subsets{39,1};
A40=subsets{40,1};
%% FullSampleMap
inputs = GHDtransmission;

% Create a Self-Organizing Map
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);

% Train the Network
[net,tr] = train(net,inputs);

% Test the Network
outputs = net(inputs);

% View the Network
view(net)

% Plots
% Uncomment these lines to enable various plots.
% figure, plotsomtop(net)
% figure, plotsomnc(net)
fig= plotsomnd(net);
savefig('GHDtransmission.fig');
saveas (fig,'GHDtransmission.tif');
% figure, plotsomplanes(net)
% figure, plotsomhits(net,inputs)
% figure, plotsompos(net,inputs)
%% A1
inputs = A1;

% Create a Self-Organizing Map
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);

% Train the Network
[net,tr] = train(net,inputs);

% Test the Network
outputs = net(inputs);

% View the Network
view(net)

% Plots
% Uncomment these lines to enable various plots.
% figure, plotsomtop(net)
% figure, plotsomnc(net)
fig= plotsomnd(net);
savefig('A1.fig');
saveas (fig,'A1.tif');
% figure, plotsomplanes(net)
% figure, plotsomhits(net,inputs)
% figure, plotsompos(net,inputs)

%% A2
inputs = A2;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A2.fig');
saveas (fig,'A2.tif');
%% A3
inputs = A3;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A3.fig');
saveas (fig,'A3.tif');
%% A4
inputs = A4;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A4.fig');
saveas (fig,'A4.tif');
%% A5
inputs = A5;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A5.fig');
saveas (fig,'A5.tif');
%% A6
inputs = A6;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A6.fig');
saveas (fig,'A6.tif');
%% A7
inputs = A7;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A7.fig');
saveas (fig,'A7.tif');
%% A8
inputs = A8;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A8.fig');
saveas (fig,'A8.tif');
%% A9
inputs = A9;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A9.fig');
saveas (fig,'A9.tif');
%% A10
inputs = A10;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A10.fig');
saveas (fig,'A10.tif');
%% A11
inputs = A11;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A11.fig');
saveas (fig,'A11.tif');
%% A12
inputs = A12;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A12.fig');
saveas (fig,'A12.tif');
%% A13
inputs = A13;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A13.fig');
saveas (fig,'A13.tif');
%% A14
inputs = A14;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A14.fig');
saveas (fig,'A14.tif');
%% A15
inputs = A15;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A15.fig');
saveas (fig,'A15.tif');
%% A16
inputs = A16;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A16.fig');
saveas (fig,'A16.tif');
%% A17
inputs = A17;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A17.fig');
saveas (fig,'A17.tif');
%% A18
inputs = A18;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A18.fig');
saveas (fig,'A18.tif');
%% A19
inputs = A19;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A19.fig');
saveas (fig,'A19.tif');
%% A20
inputs = A20;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A20.fig');
saveas (fig,'A20.tif');
%% A21
inputs = A21;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A21.fig');
saveas (fig,'A21.tif');
%% A22
inputs = A22;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A22.fig');
saveas (fig,'A22.tif');
%% A23
inputs = A23;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A23.fig');
saveas (fig,'A23.tif');
%% A24
inputs = A24;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A24.fig');
saveas (fig,'A24.tif');
%% A25
inputs = A25;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A25.fig');
saveas (fig,'A25.tif');
%% A26
inputs = A26;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A26.fig');
saveas (fig,'A26.tif');
%% A27
inputs = A27;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A27.fig');
saveas (fig,'A27.tif');
%% A28
inputs = A28;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A28.fig');
saveas (fig,'A28.tif');
%% A29
inputs = A29;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A29.fig');
saveas(fig,'A29.tif');
%% A30
inputs = A30;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A30.fig');
saveas(fig,'A30.tif');
%% A31
inputs = A31;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A31.fig');
saveas(fig,'A31.tif');
%% A32
inputs = A32;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A32.fig');
saveas(fig,'A32.tif');
%% A33
inputs = A33;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A33.fig');
saveas(fig,'A33.tif');
%% A34
inputs = A34;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A34.fig');
saveas(fig,'A34.tif');
%% A35
inputs = A35;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A35.fig');
saveas(fig,'A35.tif');
%% A36
inputs = A36;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A36.fig');
saveas(fig,'A36.tif');
%% A37
inputs = A37;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A37.fig');
saveas(fig,'A37.tif');
%% A38
inputs = A38;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A38.fig');
saveas(fig,'A38.tif');
%% A39
inputs = A39;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A39.fig');
saveas(fig,'A39.tif');
%% A40
inputs = A40;
dimension1 = 31;
dimension2 = 31;
net = selforgmap([dimension1 dimension2]);
[net,tr] = train(net,inputs);
outputs = net(inputs);
fig= plotsomnd(net);
savefig('A40.fig');
saveas(fig,'A40.tif');
