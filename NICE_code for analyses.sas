options nofmterr;

data master;
set 'R:\Medicine\GIM\McDermott\NICE\NICE Data Analysis\Dora NICE\Final analyses\Manuscript\revision\dataset and code for submission\NICE_manuscript_dataset';

ch_dist6min_m30 = dist6min_m3 - dist6min_m0;
ch_dist6min_m60 = dist6min_m6 - dist6min_m0;

ch_time60 = time6 - time0;

ch_distscr60 = distscr6 - distscr0;
ch_speedscr60 = speedscr6 - speedscr0;
ch_climbscr60 = climbscr6 - climbscr0;

ch_PF60 = PF6 - PF0;

ch_PGC1a60 = PGC1a6 - PGC1a0;
ch_SirT160 = SirT16 - SirT10;
ch_PARP_full60 = PARP_full6 - PARP_full0;
ch_PARP_89kDa60 = PARP_89kDa6 - PARP_89kDa0;
ch_eNOS60 = eNOS6 - eNOS0;
ch_p_eNOS60 = p_eNOS6 - p_eNOS0;

ch_CPM60 = m_CPM6 - m_CPM0;
ch_tot_cnt60 = m_tot_cnt6 - m_tot_cnt0;

ch_CS60 = Citrate_Synthase_Activity6 - Citrate_Synthase_Activity0;
ch_COX60 = Cytochrome_c_Oxidase_Activity6 - Cytochrome_c_Oxidase_Activity0;

ch_NAD60 = NAD6 - NAD0;

ch_Capillary_Density60 = Capillary_Density6 - Capillary_Density0;
ch_SC_Total_Fibers60 = SC_Total_Fibers6 - SC_Total_Fibers0;
ch_SC_Type1_Fibers_pct_60 = SC_Type1_Fibers_pct_6 - SC_Type1_Fibers_pct_0;
ch_SC_Type2_Fibers_pct_60 = SC_Type2_Fibers_pct_6 - SC_Type2_Fibers_pct_0;

if group = 'NR + Resveratrol' then do;
   group2 = 'NR     ';
   resv_ind = 1;
end;
else if group = 'NR alone' then do;
   group2 = 'NR     ';
   resv_ind = 0;
end;
else if group = 'Placebo' then do;
   group2 = 'Placebo';
   resv_ind = 0;
end;

label ch_dist6min_m30 = '3-month change in six-minute walk (meters)'
      ch_dist6min_m60 = '6-month change in six-minute walk (meters)'
	  ch_time60 = '6-month change in maximal treadmill walking time (minutes)'
	  ch_distscr60 = '6-month change in WIQ distance score'
      ch_speedscr60 = '6-month change in WIQ speed score'
      ch_climbscr60 = '6-month change in WIQ stair-climbing score'
	  ch_PF60 = '6-month change in SF-36 physical functioning score'
      ch_PGC1a60 = '6-month change in PGC1a (normalized value)'
      ch_SirT160 = '6-month change in SirT1 (normalized value)'
      ch_PARP_full60 = '6-month change in PARP (full length) (normalized value)'
      ch_PARP_89kDa60 = '6-month change in PARP (89kDa) (normalized value)'
      ch_eNOS60 = '6-month change in eNOS (normalized value)'
      ch_p_eNOS60 = '6-month change in p-eNOS (normalized value)'
      ch_CPM60 = '6-month change in ActiGraph activity counts per minute (CPM)'
      ch_tot_cnt60 = '6-month change in ActiGraph total activity counts per day'
      ch_CS60 = '6-month change in Citrate Synthase Activity (nmol/min/mg protein)'
      ch_COX60 = '6-month change in Cytochrome c Oxidase Activity (nmol/min/mg protein)'
      ch_NAD60 = '6-month change in NAD+'
      ch_Capillary_Density60 = '6-month change in capillary density'
      ch_SC_Total_Fibers60 = '6-month change in Satellite Cells Total SCs/100 Fibers (Pax7)'
      ch_SC_Type1_Fibers_pct_60 = '6-month change in % Type 1 fibers'
      ch_SC_Type2_Fibers_pct_60 = '6-month change in % Type 2 fibers';
run;


**********************************************************
** Primary Analyses - MMRM for six-minute walk distance
**********************************************************;
data dist6min30(keep=id group group2 resv_ind age sex aarace dist6min_m0 ch_dist6min_m30) 
     dist6min60(keep=id group group2 resv_ind age sex aarace dist6min_m0 ch_dist6min_m60) ;
set master;
run;

data master_long;
set dist6min30(rename=(ch_dist6min_m30=ch_dist6min) in=a) 
    dist6min60(rename=(ch_dist6min_m60=ch_dist6min) in=b);
if a then visit = 'FV3';
if b then visit = 'FV6';
label ch_dist6min = 'Change in six-minute walk';
run;
proc sort data = master_long; by id visit; run;


proc means data = master n mean stddev maxdec=2;
class group;
var dist6min_m0 dist6min_m6;
run;

proc means data = master n mean stddev maxdec=2;
class group;
var dist6min_m0 dist6min_m3;
run;


**Adjusted for baseline six-minute walk + age + sex + race;
**one-sided 90% CI (use one side 80% CI to infinite);
proc mixed data=master_long;
  class id group(ref='Placebo') visit sex aarace;
  model ch_dist6min = visit group visit*group dist6min_m0 age sex aarace/solution;
  repeated visit / subject=id type=un;
  lsmeans group*visit / cl diff alpha=0.20;
run;



*******************************************
** ANCOVA for all the other outcomes
*******************************************;

%macro ancova(blv=, fv=, y=);
proc means data = master(where=(&y^=.)) n mean stddev maxdec=2;
class group;
var &blv &fv &y;
run;

title 'Adjused for baseline + age + sex + race - one-sided 90% CI';
proc glm data = master alpha=0.20;
class group(ref='Placebo') sex aarace;
model &y = group &blv age sex aarace/solution ;
lsmeans group/pdiff tdiff cl stderr;
run;

quit;
%mend;


%macro ancova_comb(blv=, fv=, y=);
proc means data = master(where=(&y^=.)) n mean stddev maxdec=2;
class group2;
var &blv &fv &y;
run;

title 'Adjused for baseline + age + sex + race + resv_ind - one-sided 90% CI';
proc glm data = master alpha=0.20;
class group2(ref='Placebo') sex aarace resv_ind(ref='0');
model &y = group2 &blv age sex aarace resv_ind/solution ;
lsmeans group2/pdiff tdiff cl stderr;
run;

title 'Adjused for baseline + age + sex + race - one-sided 90% CI';
proc glm data = master alpha=0.20;
class group2(ref='Placebo') sex aarace;
model &y = group2 &blv age sex aarace /solution ;
lsmeans group2/pdiff tdiff cl stderr;
run;

quit;
%mend;


%ancova(blv=time0, fv=time6, y=ch_time60);

%ancova(blv=distscr0, fv=distscr6, y=ch_distscr60);
%ancova(blv=speedscr0, fv=speedscr6, y=ch_speedscr60);
%ancova(blv=climbscr0, fv=climbscr6, y=ch_climbscr60);

%ancova(blv=PF0, fv=PF6, y=ch_PF60);

%ancova(blv=PGC1a0, fv=PGC1a6, y=ch_PGC1a60);
%ancova(blv=SirT10, fv=SirT16, y=ch_SirT160);

%ancova(blv=PARP_full0, fv=PARP_full6, y=ch_PARP_full60);
%ancova(blv=PARP_89kDa0, fv=PARP_89kDa6, y=ch_PARP_89kDa60);

%ancova(blv=eNOS0, fv=eNOS6, y=ch_eNOS60);
%ancova(blv=p_eNOS0, fv=p_eNOS6, y=ch_p_eNOS60);

%ancova(blv=m_tot_cnt0, fv=m_tot_cnt6, y=ch_tot_cnt60);
%ancova(blv=m_CPM0, fv=m_CPM6, y=ch_CPM60);

%ancova(blv=Citrate_Synthase_Activity0, fv=Citrate_Synthase_Activity6, y=ch_CS60);
%ancova(blv=Cytochrome_c_Oxidase_Activity0, fv=Cytochrome_c_Oxidase_Activity6, y=ch_COX60);

%ancova(blv=NAD0, fv=NAD6, y=ch_NAD60);

%ancova(blv=Capillary_Density0, fv=Capillary_Density6, y=ch_Capillary_Density60);
%ancova(blv=SC_Total_Fibers0, fv=SC_Total_Fibers6, y=ch_SC_Total_Fibers60);
%ancova(blv=SC_Type1_Fibers_pct_0, fv=SC_Type1_Fibers_pct_6, y=ch_SC_Type1_Fibers_pct_60);
%ancova(blv=SC_Type2_Fibers_pct_0, fv=SC_Type2_Fibers_pct_6, y=ch_SC_Type2_Fibers_pct_60);


%ancova_comb(blv=time0, fv=time6, y=ch_time60);

%ancova_comb(blv=distscr0, fv=distscr6, y=ch_distscr60);
%ancova_comb(blv=speedscr0, fv=speedscr6, y=ch_speedscr60);
%ancova_comb(blv=climbscr0, fv=climbscr6, y=ch_climbscr60);

%ancova_comb(blv=PF0, fv=PF6, y=ch_PF60);

%ancova_comb(blv=PGC1a0, fv=PGC1a6, y=ch_PGC1a60);
%ancova_comb(blv=SirT10, fv=SirT16, y=ch_SirT160);

%ancova_comb(blv=PARP_full0, fv=PARP_full6, y=ch_PARP_full60);
%ancova_comb(blv=PARP_89kDa0, fv=PARP_89kDa6, y=ch_PARP_89kDa60);

%ancova_comb(blv=eNOS0, fv=eNOS6, y=ch_eNOS60);
%ancova_comb(blv=p_eNOS0, fv=p_eNOS6, y=ch_p_eNOS60);

%ancova_comb(blv=m_tot_cnt0, fv=m_tot_cnt6, y=ch_tot_cnt60);
%ancova_comb(blv=m_CPM0, fv=m_CPM6, y=ch_CPM60);

%ancova_comb(blv=Citrate_Synthase_Activity0, fv=Citrate_Synthase_Activity6, y=ch_CS60);
%ancova_comb(blv=Cytochrome_c_Oxidase_Activity0, fv=Cytochrome_c_Oxidase_Activity6, y=ch_COX60);

%ancova_comb(blv=NAD0, fv=NAD6, y=ch_NAD60);

%ancova_comb(blv=Capillary_Density0, fv=Capillary_Density6, y=ch_Capillary_Density60);
%ancova_comb(blv=SC_Total_Fibers0, fv=SC_Total_Fibers6, y=ch_SC_Total_Fibers60);
%ancova_comb(blv=SC_Type1_Fibers_pct_0, fv=SC_Type1_Fibers_pct_6, y=ch_SC_Type1_Fibers_pct_60);
%ancova_comb(blv=SC_Type2_Fibers_pct_0, fv=SC_Type2_Fibers_pct_6, y=ch_SC_Type2_Fibers_pct_60);
