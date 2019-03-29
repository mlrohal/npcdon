data t; 
attrib Est length=$2 format=$2. informat=$2. label='Estuary';
attrib Bay length=$2 format=$2. informat=$2. label='Bay';
attrib Sta length=$2 format=$2. informat=$8. label='Station';
set tx.txzoop;
if est="BB" and  sta="1" then Bay="BB";
if est="BB" and  sta="2" then Bay="BB";
if est="BB" and  sta="3" then Bay="BB";
if est="BB" and  sta="4" then Bay="BB";
if est="BB" and  sta="5" then Bay="BB";
if est="BB" and  sta="6" then Bay="BB";
if est="GE" and  sta="A" then Bay="US";
if est="GE" and  sta="B" then Bay="US";
if est="GE" and  sta="C" then Bay="LS";
if est="GE" and  sta="D" then Bay="LS";
if est="LC" and  sta="A" then Bay="LB";
if est="LC" and  sta="B" then Bay="LB";
if est="LC" and  sta="C" then Bay="MB";
if est="LC" and  sta="D" then Bay="MB";
if est="LC" and  sta="E" then Bay="EM";
if est="LC" and  sta="F" then Bay="EM";
if est="NC" and  sta="A" then Bay="NB";
if est="NC" and  sta="B" then Bay="NB";
if est="NC" and  sta="C" then Bay="CC";
if est="NC" and  sta="D" then Bay="CC";
if est="NC" and  sta="E" then Bay="CC";
format date mmddyy10.;
Zoop=ddens*.07/1000;  /* Zooplankton biomasss as "mg N/L" using an assumed dw:N ratio of 7% (Phytoplankton 4-9% of dry weight Parson's et al. 1961), and 1000 L/m³  */
label Zoop="Zooplankton(mg N/L)";
run;

proc sort data=t; BY bay date est sta rep; run;
proc means data=t noprint;
 by bay date est sta;
 var Zoop;
 output out=tr mean=Zoop;
 run;

proc means data=tr noprint;
 by bay date est ;
 var Zoop;
 output out=ts mean=Zoop;
 run;

proc means data=ts noprint;
 by bay date ; 
 var Zoop;
 output out=Zoop(drop=_type_ _freq_) mean=Zoop;
 run;

%macro ez (bay);
proc export data=Zoop (where=(Bay="&bay")) 
  outfile = "C:\Users\mrohal\Desktop\NPZ Model\npcdon\Data\&bay._Zoop.csv"
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
