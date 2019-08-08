
data M; 
set tx.txmamg;
if est="MA" then delete;
format date mmddyy10.;
GM2=GM2*.1;   /* (Originally g/m2) MACROFAUNA biomasss as "mg/M²" * 1000, Convert to Liters 1/1000, Convert to N * .10 (using a rough assumed dw:N ratio of 10% (Chapter 12 in Higgins and Thiel, 1988)  */
label gm2="Infauna(mg N/L)";
run;

proc sort data=M; BY bay date est sta rep SEC; run;
proc means data=M noprint;
 by bay date est sta rep;
 var gm2;
 output out=M1 sum=gm2;
 run;

proc means data=M1 noprint;
 by bay date est sta;
 var gm2;
 output out=Mr mean=gm2;
 run;

*proc means data=Mr noprint;
 *by bay date est ;
 *var gm2;
 *output out=Ms mean=gm2;
 *run;

proc means data=Mr noprint;
 by bay date sta ; 
 var gm2;
 output out=MG(drop=_type_ _freq_ sta) mean=gm2;
 run;


%macro ez (bay);
proc export data=MG (where=(Bay="&bay")) 
  outfile = "C:\Users\mrohal\Desktop\NPZ Model\npcdon\Data\&bay._MG.csv"
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
