data  B2;
set B;
label gm2 = "Benthos in mg N / L";
label DATE = "Date";
run;

proc sgplot data=B2;
   title 'Upper San Antonio Bay Benthos';
   series x=Date y=gm2 / markers;
run;

proc sgplot data=B2;
   title 'Upper San Antonio Bay Benthos';
   scatter x=Date y=gm2;
run;

data N2;
set N;
label DIN = "Dissolved Inorganic Nitrogen in mg N / L";
label DON = "Dissolved Organic Nitrogen in mg N / L";
label DATE = "Date";
run;
proc sgplot data=N2;
   title 'Upper San Antonio Bay Inorganic Nitrogen Concentrations';
   scatter x=Date y=DIN;
run;
proc sgplot data=N2;
   title 'Upper San Antonio Bay Organic Nitrogen Concentrations';
   scatter x=Date y=DON;
run;

data  Z2;
set Z;
label Z = "Zooplankton in mg N / L";
label DATE = "Date";
run;

proc sgplot data=Z2;
   title 'Upper San Antonio Bay Zooplankton';
   scatter x=Date y=Z;
run;

data  P2;
set P;
label Chl = "Phytplankton in mg N / L";
label DATE = "Date";
run;

proc sgplot data=P2;
   title 'Upper San Antonio Bay Phytoplankton';
   scatter x=Date y=Chl;
run;


data  I2;
set I;
label ddin = "Dissolved Inorganic Nitrogen in mg N / L";
label ddon = "Dissolved Organic Nitrogen in mg N / L";
label DATE = "Date";
run;

proc sgplot data=I2;
   title 'Nueces River Dissolved Inorganic Nitrogen Input';
   scatter x=Date y=ddin;
run;

proc sgplot data=I2;
   title 'Nueces River Dissolved Organic Nitrogen Input';
   scatter x=Date y=ddon;
run;
data  G2;
set G;
label ddin = "Dissolved Inorganic Nitrogen in mg N / L";
label ddon = "Dissolved Organic Nitrogen in mg N / L";
label DATE = "Date";
run;

proc sgplot data=G2;
   title 'Guadalupe River Dissolved Inorganic Nitrogen Input';
   scatter x=Date y=ddin;
run;

proc sgplot data=G2;
   title 'Guadalupe River Dissolved Organic Nitrogen Input';
   scatter x=Date y=ddon;
run;
