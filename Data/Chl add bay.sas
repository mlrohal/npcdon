data t; set tx.tx_allchem; drop bay; run;

data c; 
attrib Est length=$2 format=$2. informat=$2. label='Estuary';
attrib Bay length=$2 format=$2. informat=$2. label='Bay';
attrib Sta length=$2 format=$2. informat=$2. label='Station';
set t;
if est="BB" and  sta="3" then Bay="BB";
if est="BB" and  sta="6" then Bay="BB";
if est="BB" and  sta="8" then Bay="BB";
if est="GE" and  sta="A" then Bay="US";
if est="GE" and  sta="B" then Bay="US";
if est="GE" and  sta="C" then Bay="LS";
if est="GE" and  sta="D" then Bay="LS";
if est="LC" and  sta="A" then Bay="LB";
if est="LC" and  sta="B" then Bay="LB";
if est="LC" and  sta="C" then Bay="MB";
if est="LC" and  sta="D" then Bay="MB";
if est="LC" and  sta="E" then Bay="EM";
if est="LC" and  sta="8" then Bay="EM";
if est="LC" and  sta="15" then Bay="EM";
if est="LC" and  sta="F" then Bay="EM";
if est="LC" and  sta="FD" then Bay="LB";
if est="NC" and  sta="A" then Bay="NB";
if est="NC" and  sta="B" then Bay="NB";
if est="NC" and  sta="C" then Bay="CC";
if est="NC" and  sta="D" then Bay="CC";
if est="NC" and  sta="E" then Bay="CC";
format date mmddyy10.;
Chl = Chl*8.85*0.001; /*conversion to from ug/l to mg N/L based on average N:Chl ratio in shallow marine water from Yentsch and Vaccaro */
run;
proc sort data=c; BY bay date est sta DEPTH; run;
proc means data=c noprint;
 by bay date est sta;
 var Chl;
 output out=cr mean=Chl;
 run;

proc means data=cr noprint;
by bay date sta ;
 var Chl;
 output out=cs(drop=_type_ _freq_ sta) mean=Chl;
 run;

*proc means data=cs noprint;
 *by bay date ; 
 *var DIN DON;
 *output out=N(drop=_type_ _freq_) mean=DIN DON;
 *run;

data Chl;
 set cs;
 run;

%macro ez (bay);
proc export data=Chl (where=(Bay="&bay")) 
  outfile = "C:\Users\mrohal\Desktop\NPZ Model\npcdon\Data\&bay._Chl.csv"
  dbms=csv replace;
  run;
%mend ez;

%ez (BB)
%ez (US)
%ez (LS)
%ez (LB)
%ez (MB)
%ez (EM)
%ez (NB)
%ez (CC)
