clc
clear
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\Important_codes\IMPORTANT DY CODES\Working_on_DY_graphs\CTO.mat');
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\Important_codes\IMPORTANT DY CODES\Working_on_DY_graphs\CFO.mat');
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\Important_codes\IMPORTANT DY CODES\Working_on_DY_graphs\CTOOil.mat');
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\Important_codes\IMPORTANT DY CODES\Working_on_DY_graphs\CFOOil.mat');
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\Important_codes\IMPORTANT DY CODES\Working_on_DY_graphs\CTOCom.mat');
load('F:\OneDrive - University of Tasmania\Mardi Meetings\Meeting 22\Testing\Important_codes\IMPORTANT DY CODES\Working_on_DY_graphs\CFOCom.mat');

%% Recognise all indices

USF=     movmean(CFO(:,1),250);
AUSF=    movmean(CFO(:,2),250);
INDF=    movmean(CFO(:,3),250);
JAPF=    movmean(CFO(:,4),250);
MYSF=    movmean(CFO(:,5),250);
NZLF=    movmean(CFO(:,6),250);
SGPF=    movmean(CFO(:,7),250);
PHLF=    movmean(CFO(:,8),250);
KORF=    movmean(CFO(:,9),250);
SLKF=    movmean(CFO(:,10),250);
THAF=    movmean(CFO(:,11),250);
NGAF=    movmean(CFO(:,12),250);
VENF=    movmean(CFO(:,13),250);
KWTF=    movmean(CFO(:,14),250);
IRQF=    movmean(CFO(:,15),250);
SAUF=    movmean(CFO(:,16),250);
CHNF=    movmean(CFO(:,17),250);
ISRF=    movmean(CFO(:,18),250);
CADF=    movmean(CFO(:,19),250);
GRCF=    movmean(CFO(:,20),250);
PRTF=    movmean(CFO(:,21),250);
IRLF=    movmean(CFO(:,22),250);
BELF=    movmean(CFO(:,23),250);
CRTF=    movmean(CFO(:,24),250);
AUTF=    movmean(CFO(:,25),250);
RUSF=    movmean(CFO(:,26),250);
NORF=    movmean(CFO(:,27),250);
GERF=    movmean(CFO(:,28),250);
CHLF=    movmean(CFO(:,29),250);
UKF=     movmean(CFO(:,30),250);
FRAF=    movmean(CFO(:,31),250);


%%

UST=     movmean(CTO(:,1),250);
AUST=    movmean(CTO(:,2),250);
INDT=    movmean(CTO(:,3),250);
JAPT=    movmean(CTO(:,4),250);
MYST=    movmean(CTO(:,5),250);
NZLT=    movmean(CTO(:,6),250);
SGPT=    movmean(CTO(:,7),250);
PHLT=    movmean(CTO(:,8),250);
KORT=    movmean(CTO(:,9),250);
SLKT=    movmean(CTO(:,10),250);
THAT=    movmean(CTO(:,11),250);
NGAT=    movmean(CTO(:,12),250);
VENT=    movmean(CTO(:,13),250);
KWTT=    movmean(CTO(:,14),250);
IRQT=    movmean(CTO(:,15),250);
SAUT=    movmean(CTO(:,16),250);
CHNT=    movmean(CTO(:,17),250);
ISRT=    movmean(CTO(:,18),250);
CADT=    movmean(CTO(:,19),250);
GRCT=    movmean(CTO(:,20),250);
PRTT=    movmean(CTO(:,21),250);
IRLT=    movmean(CTO(:,22),250);
BELT=    movmean(CTO(:,23),250);
CRTT=    movmean(CTO(:,24),250);
AUTT=    movmean(CTO(:,25),250);
RUST=    movmean(CTO(:,26),250);
NORT=    movmean(CTO(:,27),250);
GERT=    movmean(CTO(:,28),250);
CHLT=    movmean(CTO(:,29),250);
UKT=     movmean(CTO(:,30),250);
FRAT=    movmean(CTO(:,31),250);

%% Oil

USFO=     movmean(CFOOil(:,1),250);
AUSFO=    movmean(CFOOil(:,2),250);
INDFO=    movmean(CFOOil(:,3),250);
JAPFO=    movmean(CFOOil(:,4),250);
MYSFO=    movmean(CFOOil(:,5),250);
NZLFO=    movmean(CFOOil(:,6),250);
SGPFO=    movmean(CFOOil(:,7),250);
PHLFO=    movmean(CFOOil(:,8),250);
KORFO=    movmean(CFOOil(:,9),250);
SLKFO=    movmean(CFOOil(:,10),250);
THAFO=    movmean(CFOOil(:,11),250);
NGAFO=    movmean(CFOOil(:,12),250);
VENFO=    movmean(CFOOil(:,13),250);
KWTFO=    movmean(CFOOil(:,14),250);
IRQFO=    movmean(CFOOil(:,15),250);
SAUFO=    movmean(CFOOil(:,16),250);
CHNFO=    movmean(CFOOil(:,17),250);
ISRFO=    movmean(CFOOil(:,18),250);
CADFO=    movmean(CFOOil(:,19),250);
GRCFO=    movmean(CFOOil(:,20),250);
PRTFO=    movmean(CFOOil(:,21),250);
IRLFO=    movmean(CFOOil(:,22),250);
BELFO=    movmean(CFOOil(:,23),250);
CRTFO=    movmean(CFOOil(:,24),250);
AUTFO=    movmean(CFOOil(:,25),250);
RUSFO=    movmean(CFOOil(:,26),250);
NORFO=    movmean(CFOOil(:,27),250);
GERFO=    movmean(CFOOil(:,28),250);
CHLFO=    movmean(CFOOil(:,29),250);
UKFO=     movmean(CFOOil(:,30),250);
FRAFO=    movmean(CFOOil(:,31),250);
%%
USTO=     movmean(CTOOil(:,1),250);
AUSTO=    movmean(CTOOil(:,2),250);
INDTO=    movmean(CTOOil(:,3),250);
JAPTO=    movmean(CTOOil(:,4),250);
MYSTO=    movmean(CTOOil(:,5),250);
NZLTO=    movmean(CTOOil(:,6),250);
SGPTO=    movmean(CTOOil(:,7),250);
PHLTO=    movmean(CTOOil(:,8),250);
KORTO=    movmean(CTOOil(:,9),250);
SLKTO=    movmean(CTOOil(:,10),250);
THATO=    movmean(CTOOil(:,11),250);
NGATO=    movmean(CTOOil(:,12),250);
VENTO=    movmean(CTOOil(:,13),250);
KWTTO=    movmean(CTOOil(:,14),250);
IRQTO=    movmean(CTOOil(:,15),250);
SAUTO=    movmean(CTOOil(:,16),250);
CHNTO=    movmean(CTOOil(:,17),250);
ISRTO=    movmean(CTOOil(:,18),250);
CADTO=    movmean(CTOOil(:,19),250);
GRCTO=    movmean(CTOOil(:,20),250);
PRTTO=    movmean(CTOOil(:,21),250);
IRLTO=    movmean(CTOOil(:,22),250);
BELTO=    movmean(CTOOil(:,23),250);
CRTTO=    movmean(CTOOil(:,24),250);
AUTTO=    movmean(CTOOil(:,25),250);
RUSTO=    movmean(CTOOil(:,26),250);
NORTO=    movmean(CTOOil(:,27),250);
GERTO=    movmean(CTOOil(:,28),250);
CHLTO=    movmean(CTOOil(:,29),250);
UKTO=     movmean(CTOOil(:,30),250);
FRATO=    movmean(CTOOil(:,31),250);
%%
USFC=     movmean(CFOCom(:,1),250);
AUSFC=    movmean(CFOCom(:,2),250);
INDFC=    movmean(CFOCom(:,3),250);
JAPFC=    movmean(CFOCom(:,4),250);
MYSFC=    movmean(CFOCom(:,5),250);
NZLFC=    movmean(CFOCom(:,6),250);
SGPFC=    movmean(CFOCom(:,7),250);
PHLFC=    movmean(CFOCom(:,8),250);
KORFC=    movmean(CFOCom(:,9),250);
SLKFC=    movmean(CFOCom(:,10),250);
THAFC=    movmean(CFOCom(:,11),250);
NGAFC=    movmean(CFOCom(:,12),250);
VENFC=    movmean(CFOCom(:,13),250);
KWTFC=    movmean(CFOCom(:,14),250);
IRQFC=    movmean(CFOCom(:,15),250);
SAUFC=    movmean(CFOCom(:,16),250);
CHNFC=    movmean(CFOCom(:,17),250);
ISRFC=    movmean(CFOCom(:,18),250);
CADFC=    movmean(CFOCom(:,19),250);
GRCFC=    movmean(CFOCom(:,20),250);
PRTFC=    movmean(CFOCom(:,21),250);
IRLFC=    movmean(CFOCom(:,22),250);
BELFC=    movmean(CFOCom(:,23),250);
CRTFC=    movmean(CFOCom(:,24),250);
AUTFC=    movmean(CFOCom(:,25),250);
RUSFC=    movmean(CFOCom(:,26),250);
NORFC=    movmean(CFOCom(:,27),250);
GERFC=    movmean(CFOCom(:,28),250);
CHLFC=    movmean(CFOCom(:,29),250);
UKFC=     movmean(CFOCom(:,30),250);
FRAFC=    movmean(CFOCom(:,31),250);
%%
USTC=     movmean(CTOCom(:,1),250);
AUSTC=    movmean(CTOCom(:,2),250);
INDTC=    movmean(CTOCom(:,3),250);
JAPTC=    movmean(CTOCom(:,4),250);
MYSTC=    movmean(CTOCom(:,5),250);
NZLTC=    movmean(CTOCom(:,6),250);
SGPTC=    movmean(CTOCom(:,7),250);
PHLTC=    movmean(CTOCom(:,8),250);
KORTC=    movmean(CTOCom(:,9),250);
SLKTC=    movmean(CTOCom(:,10),250);
THATC=    movmean(CTOCom(:,11),250);
NGATC=    movmean(CTOCom(:,12),250);
VENTC=    movmean(CTOCom(:,13),250);
KWTTC=    movmean(CTOCom(:,14),250);
IRQTC=    movmean(CTOCom(:,15),250);
SAUTC=    movmean(CTOCom(:,16),250);
CHNTC=    movmean(CTOCom(:,17),250);
ISRTC=    movmean(CTOCom(:,18),250);
CADTC=    movmean(CTOCom(:,19),250);
GRCTC=    movmean(CTOCom(:,20),250);
PRTTC=    movmean(CTOCom(:,21),250);
IRLTC=    movmean(CTOCom(:,22),250);
BELTC=    movmean(CTOCom(:,23),250);
CRTTC=    movmean(CTOCom(:,24),250);
AUTTC=    movmean(CTOCom(:,25),250);
RUSTC=    movmean(CTOCom(:,26),250);
NORTC=    movmean(CTOCom(:,27),250);
GERTC=    movmean(CTOCom(:,28),250);
CHLTC=    movmean(CTOCom(:,29),250);
UKTC=     movmean(CTOCom(:,30),250);
FRATC=    movmean(CTOCom(:,31),250);

%% Begin plotting
%% AC
subplot (2,3,1)
plot(INDT)
hold on
plot(INDTO)
plot(INDTC)
hold off

subplot (2,3,2)
plot(MYST)
hold on
plot(MYSTO)
plot(MYSTC)
hold off

subplot (2,3,3)
plot(SGPT)
hold on
plot(SGPTO)
plot(SGPTC)
hold off

subplot(2,3,4)
plot(KORT)
hold on
plot(KORTO)
plot(KORTC)
hold off

subplot(2,3,5)
plot(PHLT)
hold on
plot(PHLTO)
plot(PHLTC)
hold off

subplot (2,3,6)
plot(THAT)
hold on
plot(THATO)
plot(THATC)
hold off

%% AC
subplot (2,3,1)
plot(INDF)
hold on
plot(INDFO)
plot(INDFC)
hold off

subplot (2,3,2)
plot(MYSF)
hold on
plot(MYSFO)
plot(MYSFC)
hold off

subplot (2,3,3)
plot(SGPF)
hold on
plot(SGPFO)
plot(SGPFC)
hold off

subplot(2,3,4)
plot(KORF)
hold on
plot(KORFO)
plot(KORFC)
hold off

subplot(2,3,5)
plot(PHLF)
hold on
plot(PHLFO)
plot(PHLFC)
hold off

subplot (2,3,6)
plot(THAF)
hold on
plot(THAFO)
plot(THAFC)
hold off
%% EC
subplot (2,3,1)
plot(GERT)
hold on
plot(GERTO)
plot(GERTC)
hold off

subplot (2,3,2)
plot(CHLT)
hold on
plot(CHLTO)
plot(CHLTC)
hold off

subplot (2,3,3)
plot(FRAT)
hold on
plot(FRATO)
plot(FRATC)
hold off

subplot(2,3,4)
plot(CHNT)
hold on
plot(CHNTO)
plot(CHNTC)
hold off

subplot(2,3,5)
plot(UKT)
hold on
plot(UKTO)
plot(UKTC)
hold off

subplot (2,3,6)
plot(AUST)
hold on
plot(AUSTO)
plot(AUSTC)
hold off

%% EC
subplot (2,3,1)
plot(GERF)
hold on
plot(GERFO)
plot(GERFC)
hold off

subplot (2,3,2)
plot(CHLF)
hold on
plot(CHLFO)
plot(CHLFC)
hold off

subplot (2,3,3)
plot(FRAF)
hold on
plot(FRAFO)
plot(FRAFC)
hold off

subplot(2,3,4)
plot(CHNF)
hold on
plot(CHNFO)
plot(CHNFC)
hold off

subplot(2,3,5)
plot(UKF)
hold on
plot(UKFO)
plot(UKFC)
hold off

subplot (2,3,6)
plot(AUSF)
hold on
plot(AUSFO)
plot(AUSFC)
hold off
%% GC

subplot (2,3,1)
plot(GRCT)
hold on
plot(GRCTO)
plot(GRCTC)
hold off

subplot (2,3,2)
plot(PRTT)
hold on
plot(PRTTO)
plot(PRTTC)
hold off

subplot (2,3,3)
plot(IRLT)
hold on
plot(IRLTO)
plot(IRLTC)
hold off

subplot(2,3,4)
plot(BELT)
hold on
plot(BELTO)
plot(BELTC)
hold off

subplot(2,3,5)
plot(CRTT)
hold on
plot(CRTTO)
plot(CRTTC)
hold off

subplot (2,3,6)
plot(AUTT)
hold on
plot(AUTTO)
plot(AUTTC)
hold off

%% GC
subplot (2,3,1)
plot(GRCF)
hold on
plot(GRCFO)
plot(GRCFC)
hold off

subplot (2,3,2)
plot(PRTF)
hold on
plot(PRTFO)
plot(PRTFC)
hold off

subplot (2,3,3)
plot(IRLF)
hold on
plot(IRLFO)
plot(IRLFC)
hold off

subplot(2,3,4)
plot(BELF)
hold on
plot(BELFO)
plot(BELFC)
hold off

subplot(2,3,5)
plot(CRTF)
hold on
plot(CRTFO)
plot(CRTFC)
hold off

subplot (2,3,6)
plot(AUTF)
hold on
plot(AUTFO)
plot(AUTFC)
hold off

%% OC1

subplot (2,3,1)
plot(UST)
hold on
plot(USTO)
plot(USTC)
hold off

subplot (2,3,2)
plot(CADT)
hold on
plot(CADTO)
plot(CADTC)
hold off

subplot (2,3,3)
plot(RUST)
hold on
plot(RUSTO)
plot(RUSTC)
hold off

subplot(2,3,4)
plot(NORT)
hold on
plot(NORTO)
plot(NORTC)
hold off

subplot(2,3,5)
plot(JAPT)
hold on
plot(JAPTO)
plot(JAPTC)
hold off

subplot (2,3,6)
plot(NZLT)
hold on
plot(NZLTO)
plot(NZLTC)
hold off

%% OC1
subplot (2,3,1)
plot(USF)
hold on
plot(USFO)
plot(USFC)
hold off

subplot (2,3,2)
plot(CADF)
hold on
plot(CADFO)
plot(CADFC)
hold off

subplot (2,3,3)
plot(RUSF)
hold on
plot(RUSFO)
plot(RUSFC)
hold off

subplot(2,3,4)
plot(NORF)
hold on
plot(NORFO)
plot(NORFC)
hold off

subplot(2,3,5)
plot(JAPF)
hold on
plot(JAPFO)
plot(JAPFC)
hold off

subplot (2,3,6)
plot(NZLF)
hold on
plot(NZLFO)
plot(NZLFC)
hold off

%% OC2


subplot (2,3,1)
plot(SAUT)
hold on
plot(SAUTO)
plot(SAUTC)
hold off

subplot (2,3,2)
plot(ISRT)
hold on
plot(ISRTO)
plot(ISRTC)
hold off

subplot (2,3,3)
plot(IRQT)
hold on
plot(IRQTO)
plot(IRQTC)
hold off

subplot(2,3,4)
plot(SLKT)
hold on
plot(SLKTO)
plot(SLKTC)
hold off

subplot(2,3,5)
plot(NGAT)
hold on
plot(NGATO)
plot(NGATC)
hold off

subplot (2,3,6)
plot(VENT)
hold on
plot(VENTO)
plot(VENTC)
hold off

%% OC2
subplot (2,3,1)
plot(SAUF)
hold on
plot(SAUFO)
plot(SAUFC)
hold off

subplot (2,3,2)
plot(ISRF)
hold on
plot(ISRFO)
plot(ISRFC)
hold off

subplot (2,3,3)
plot(IRQF)
hold on
plot(IRQFO)
plot(IRQFC)
hold off

subplot(2,3,4)
plot(SLKF)
hold on
plot(SLKFO)
plot(SLKFC)
hold off

subplot(2,3,5)
plot(NGAF)
hold on
plot(NGAFO)
plot(NGAFC)
hold off

subplot (2,3,6)
plot(VENF)
hold on
plot(VENFO)
plot(VENFC)
hold off