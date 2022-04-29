/*
Primary author: Junfeng Sun, PhD
Organizational contact information: junfeng.sun@nih.gov
Date of release: 4/29/2022
Version: 1.0
License details: public domain 
Description that clearly explains what the software/script does and for what purpose: this file documents the SAS program used for my analysis used in the paper "SARS-CoV-2 infection and persistence throughout the human body and brain among COVID-19 autopsy cases"
Usage instructions: run in SAS version 9.4
Proper attribution to others: The optimal cut-offs are selected based on Youden's J index using code similar to 
https://www.listendata.com/2015/03/sas-calculating-optimal-predicted.html
*/

PROC IMPORT OUT= WORK.pcr 
            DATAFILE= "SuppData1.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
proc print;RUN;


PROC IMPORT OUT= WORK.pats 
            DATAFILE= "pats.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
proc print;RUN;

proc sort data=pcr; by patient; run;
proc sort data=pats; by patient; run;
data dan.data1; merge pcr pats; by patient; run;
data dan.data1;set dan.data1;
ddPCR_copies_per_ng= (ddPCR_copies_per_ng1+ ddPCR_copies_per_ng2)/2; run;

PROC IMPORT OUT= WORK.tissues 
            DATAFILE= "Tissues2.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 	GUESSINGROWS=max;
data tissues; set tissues; if (Tissue^="");
proc print;RUN;

data t2; set tissues; if Tissue_Group^=.; keep Tissue_Group  Group_Name;
proc sort nodupkey; by Tissue_Group; run;
proc sort data=tissues; by Tissue_Group; run;
data dan.tissues; merge tissues(drop=Group_Name Comment) t2;by Tissue_Group; run;

data dan.data2; set dan.data1; 
keep Patient Tissue sgRNA_Cq sgRNA_copies_per_ul DOI Group ddPCR_copies_per_ng ; 
run;

proc sort data=dan.tissues; by Tissue; run;
proc sort data=dan.data2; by Tissue; run;
data dan.data2; merge dan.data2 dan.tissues; by Tissue;  run;
data dan.data2; set dan.data2; 
if (ddPCR_copies_per_ng^=.)|(sgRNA_copies_per_ul^=.); run;

data dan.data2; set dan.data2; rename Resp_vs_NonResp=Resp; run;

proc means data =dan.data2 min median max mean n;
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;
*1771, 1053;
proc means data =dan.data2 min median max mean n; where (ddPCR_copies_per_ng>0);
var ddPCR_copies_per_ng; run; *1014/1771;
proc means data =dan.data2 min median max mean n; where (sgRNA_copies_per_ul>0);
var sgRNA_copies_per_ul; run; *306/1053;


goptions reset=all;
proc gplot data =dan.data2; plot ddPCR_copies_per_ng*sgRNA_copies_per_ul; run; 

proc sort data=dan.data2; by Group; run;
ods rtf file="corr_rev.rtf";
title "Correlation between ddPCR and sgRNA";
proc corr data =dan.data2 spearman Fisher(biasadj=no); where (Resp^="Fluids");
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run; *0.759;
proc corr data =dan.data2 spearman Fisher(biasadj=no); by Group;
var ddPCR_copies_per_ng sgRNA_copies_per_ul; where (Resp^="Fluids");run;
title "When ddPCR_copies_per_ng>0";
proc corr data =dan.data2 spearman Fisher(biasadj=no); 
where (ddPCR_copies_per_ng>0)&(Resp^="Fluids");
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run; *0.762;
title "When both are non-zero";
proc corr data =dan.data2 spearman Fisher(biasadj=no); 
where (ddPCR_copies_per_ng>0)&(sgRNA_copies_per_ul>0)&(Resp^="Fluids");
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;quit; *0.762;
ods rtf close;

proc print data=dan.data2; where (Resp="Fluids"); 
var ddPCR_copies_per_ng sgRNA_copies_per_ul Resp;run; *all w/ missing;

*Compare correlation between R vs NR;
proc sort data=dan.data2; by Resp; run;
ods rtf file="corr2.rtf";
title "Correlation between ddPCR and sgRNA";
proc corr data =dan.data2 spearman Fisher(biasadj=no); 
where (Resp^="Fluids"); by Resp;
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;quit;
title "When ddPCR_copies_per_ng>0";
proc corr data =dan.data2 spearman Fisher(biasadj=no); 
where (ddPCR_copies_per_ng>0)&(Resp^="Fluids");by Resp;
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run; *0.762;
ods rtf close;

*positive pairs only;
data data2; set dan.data2; if (ddPCR_copies_per_ng>0)&(sgRNA_copies_per_ul>0); run;

ods rtf file="corr_pos.rtf";
title "Correlation between ddPCR and sgRNA (positive pairs only)";
proc corr data =data2 spearman Fisher(biasadj=no); 
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;  
proc sort data=data2; by Group; run;
proc corr data =data2 spearman Fisher(biasadj=no); by Group;
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;
proc sort data=data2; by Resp; run;
proc corr data =data2 spearman Fisher(biasadj=no); 
where (Resp^="Fluids"); by Resp;
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;quit;
ods rtf close;

*ROC curve;
data dan.roc; set dan.data2;
if (sgRNA_copies_per_ul>0) then sgRNA_pos=1;
else sgRNA_pos=sgRNA_copies_per_ul;
if sgRNA_copies_per_ul^=.; run;

proc logistic data=dan.roc plots(only)=(roc); where (Resp^="Fluids");
class sgRNA_pos(ref='0' param=ref) Resp / param=glm;
model sgRNA_pos=log10_ddPCR|Resp;
run;

ods rtf file="ROC.rtf";
title "ROC curves (Combined)";
ods select ParameterEstimates  Association ROCCurve ROCAssociation; 
proc logistic data=dan.roc plots(only)=(roc); where (Resp^="Fluids");
class sgRNA_pos(ref='0' param=ref) / param=glm;
model sgRNA_pos=log10_ddPCR;
roc "95% CI" log10_ddPCR;
run;
title "ROC curves (Resp)";
ods select ParameterEstimates  Association ROCCurve; 
proc logistic data=dan.roc plots(only)=(roc); where (Resp="Resp");
class sgRNA_pos(ref='0' param=ref) / param=glm;
model sgRNA_pos=log10_ddPCR;
run;
title "ROC curves (NonResp)";
ods select ParameterEstimates  Association ROCCurve; 
proc logistic data=dan.roc plots(only)=(roc); where (Resp="NonResp");
class sgRNA_pos(ref='0' param=ref) / param=glm;
model sgRNA_pos=log10_ddPCR;
run;
ods rtf close;


ods rtf file="cutoff.rtf";
title "log10_ddPCR: cutoff";
ods select ParameterEstimates  Association ROCCurve;
proc logistic data=dan.roc plots(only)=(roc); where (Resp^="Fluids");
class sgRNA_pos(ref='0' param=ref) / param=glm;
model sgRNA_pos=log10_ddPCR  / outroc=rocstats;
run;

data check;
set rocstats;
_SPECIF_ = (1 - _1MSPEC_);
J = _SENSIT_ + _SPECIF_ - 1;
D= Sqrt((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
*proc print data=check;run;

proc sql noprint;
create table cutoff as
select _PROB_ , J
from check
having J = max(J);
run;
proc print data=cutoff; run; *_prob_=0.29638, sens=0.93046, spec=0.91563;
proc print data=check; where _prob_>0.296 & _prob_<0.297; run;
ods rtf close;

data dan.roc; set dan.roc; 
if (ddPCR_copies_per_ng>1.469) then ddPCR_high="Y";
if (ddPCR_copies_per_ng<1.469)&(ddPCR_copies_per_ng^=.) then ddPCR_high="N";
proc freq; tables sgRNA_pos*ddPCR_high / nocol nopercent; run;
*correct;




*********************************************************************;
ods rtf file="ddPCR_LMM.rtf";
title "ddPCR";
proc print data=dan.data2; where (patient="P1");
var Patient Tissue Group ddPCR_copies_per_ng Tissue_Group Resp Group_Name; 
RUN; 
title "ddPCR: LMM";
ods select solutionf tests3;
proc mixed data=dan.data2 covtest;  
class Patient Group Resp;
model ddPCR_copies_per_ng  = Group|Resp / s residual outp=res1; 
random intercept / subject=Patient; run;
title "Residuals:";
proc univariate data=res1 noprint; qqplot PearsonResid ; run;
proc print data=res1; where ((PearsonResid<-3 | PearsonResid>3) & PearsonResid^=.); 
var Patient Tissue Group ddPCR_copies_per_ng Tissue_Group Resp Group_Name PearsonResid; 
run; 
ods rtf close;


*corrections;
proc PRINT data =dan.data2; 
where (ddPCR_copies_per_ng<0.00363)&(ddPCR_copies_per_ng>0);
var Patient Tissue ddPCR_copies_per_ng; run;
proc PRINT data =dan.data1; where (Patient="P38")&(Tissue="Testis"); run;
data dan.data1; set dan.data1;
if (Patient="P38")&(Tissue="Testis") then do;
	ddPCR_copies_per_ng2=0;
	ddPCR_copies_per_ng=0.0049;
	end;
run;

proc PRINT data =dan.data2; where (Patient="P38")&(Tissue="Testis"); run;
data dan.data2; set dan.data2;
if (Patient="P38")&(Tissue="Testis") then ddPCR_copies_per_ng=0.0049;
run;


*impute and log-transform;
data dan.data2; set dan.data2;
if (ddPCR_copies_per_ng=0) then do;
	if (Patient="P3") then log10_ddPCR=log10(rand("Uniform")*0.000217);
	else log10_ddPCR=log10(rand("Uniform")*0.00182);
	end;
else log10_ddPCR=log10(ddPCR_copies_per_ng);
run;

proc means data =dan.data2 min median max mean n; where (ddPCR_copies_per_ng=0);
var ddPCR_copies_per_ng log10_ddPCR; run; 

proc sort data=dan.data2; by Resp; run;
proc gplot data =dan.data2; plot log10_ddPCR*DOI; by Resp; run;
data dan.data2; set dan.data2; log10_DOI=log10(DOI); run;
proc gplot data =dan.data2; plot log10_ddPCR*log10_DOI; by Resp; run;
*looks much better;
proc gplot data =dan.data2; plot log10_ddPCR*log10_DOI= Resp; run; quit;

ods rtf file="log10_ddPCR_LMM.rtf";
title "log10(ddPCR): LMM";
ods select solutionf tests3;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient Group Resp(ref="NonResp") tissue;
model log10_ddPCR  = Group|Resp / s residual outp=res1; 
repeated tissue / subject=Patient type=cs;  run;
ods select solutionf estimates;
proc mixed data=dan.data2 covtest;   where (Resp^="Fluids");
class Patient Group Resp(ref="NonResp") tissue;
model log10_ddPCR  = Group Resp(Group) / s ; 
estimate "Mid vs Early" Resp(Group) 1 -1 0 0 -1 1;
estimate "Late vs Mid" Resp(Group) 0 0 -1 1 1 -1;
repeated tissue / subject=Patient type=cs;  run;
proc univariate data=res1 noprint; qqplot PearsonResid ; run;
proc print data=res1; where ((PearsonResid<-3 | PearsonResid>3) & PearsonResid^=.); 
var Patient Tissue Group ddPCR_copies_per_ng log10_ddPCR Tissue_Group Resp Group_Name PearsonResid; 
run; 

title "log10(ddPCR): LMM (log10_DOI)";
proc gplot data =dan.data2;  where (Resp^="Fluids"); plot log10_ddPCR*log10_DOI= Resp; run;quit;
ods select solutionf tests3;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient  Resp(ref="NonResp") tissue;
model log10_ddPCR  = log10_DOI|Resp / s residual outp=res1; 
repeated tissue / subject=Patient type=cs;  run;
ods select solutionf;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient  Resp(ref="NonResp") tissue;
model log10_ddPCR  = Resp log10_DOI(Resp) / s noint; 
repeated tissue / subject=Patient type=cs; run;
title "Residuals:";
proc univariate data=res1 noprint; qqplot PearsonResid ; run;
proc print data=res1; where ((PearsonResid<-3 | PearsonResid>3) & PearsonResid^=.); 
var Patient Tissue Group log10_DOI ddPCR_copies_per_ng log10_ddPCR Tissue_Group Resp Group_Name PearsonResid; 
run; 
ods rtf close;

*figure;
proc print data =dan.data2(obs=50);where (Resp^="Fluids");run;
data data2; set dan.data2; if (Resp^="Fluids");
keep Patient Resp log10_ddPCR log10_DOI Group; proc print; run;
PROC EXPORT DATA= data2 
            OUTFILE= "figure.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


ods pdf file="figure.pdf";
goptions reset=all hsize=20cm vsize=15cm;
title "log10(ddPCR): LMM (log10_DOI)";
proc gplot data =dan.data2;  where (Resp^="Fluids"); 
plot log10_ddPCR*log10_DOI= Resp; run;quit;
title "log10(ddPCR) vs (log10_DOI): Resp";
proc gplot data =data2;  where (Resp="Resp"); 
plot log10_ddPCR*log10_DOI; run;quit;
title "log10(ddPCR) vs (log10_DOI): NonResp";
proc gplot data =data2;  where (Resp="NonResp"); 
plot log10_ddPCR*log10_DOI; run;quit;
ods pdf close;



*non parametric;
proc rank data=dan.data2 out=r1; where (Resp^="Fluids");
var ddPCR_copies_per_ng; ranks r; run;
/*proc sort data=r1 out=r1; by Patient Tissue; run;*/

ods rtf file="log10_ddPCR_NP.rtf";
title "NonParametric";
ods select tests3;
proc mixed data=r1 ANOVAF METHOD=MIVQUE0;  where (Resp^="Fluids");
class Patient Group Resp(ref="NonResp") tissue;
model r  = Group|Resp / s chisq; 
repeated tissue / subject=Patient type=cs; run;*grp big difference;
ods select solutionf estimates;
proc mixed data=r1 ANOVAF METHOD=MIVQUE0;  where (Resp^="Fluids");
class Patient Group Resp(ref="NonResp") tissue;
model r  = Group Resp(Group) / s chisq; 
estimate "Mid vs Early" Resp(Group) 1 -1 0 0 -1 1;
estimate "Late vs Mid" Resp(Group) 0 0 -1 1 1 -1;
repeated tissue / subject=Patient type=cs; run;*grp big difference;

title "NonParametric(log10_DOI)";
ods select solutionf tests3;
proc mixed data=r1 ANOVAF METHOD=MIVQUE0;  where (Resp^="Fluids");
class Patient Group Resp(ref="NonResp") tissue;
model r  = log10_DOI|Resp / s chisq; 
repeated tissue / subject=Patient type=cs; run;*grp big difference;
ods select solutionf;
proc mixed data=r1 ANOVAF METHOD=MIVQUE0;  where (Resp^="Fluids");
class Patient Group Resp(ref="NonResp") tissue;
model r  = Resp log10_DOI(Resp) / s noint chisq; 
repeated tissue / subject=Patient type=cs; run;*grp big difference;
ods rtf close;

*compare tissues;
proc freq data=dan.data2; tables Group_Name*resp / nocol norow nopercent; run;
proc print data=dan.data2; run;

proc sort data=dan.data2; by Resp; run;
ods rtf file="log10_ddPCR_LMM_tissues.rtf";
title "log10(ddPCR): LMM ";
ods select tests3;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient Group Group_Name tissue;
model log10_ddPCR  = Group|Group_Name / s residual outp=res1; 
repeated tissue / subject=Patient type=cs; by Resp; run;
ods select solutionf;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient Group(ref="Mid") Group_Name tissue;
model log10_ddPCR  = Group_Name Group(Group_Name) / s ; 
repeated tissue / subject=Patient type=cs; by Resp; run;

title "log10(ddPCR): LMM (log10_DOI)";
ods select tests3;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient  Group_Name tissue;
model log10_ddPCR  = log10_DOI|Group_Name / s residual outp=res1; 
repeated tissue / subject=Patient type=cs;  by Resp; run;
ods select solutionf;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient  Group_Name tissue;
model log10_ddPCR  = Group_Name log10_DOI(Group_Name) / s noint; 
repeated tissue / subject=Patient type=cs; by Resp; run;
ods rtf close;

*Extended Data Table 3a;
proc sort data=dan.data2; by Group_Name Group; run;
ods trace off;
proc means data=dan.data2 mean std n q1 median q3 ; where (Resp^="Fluids");
var ddPCR_copies_per_ng; ods output summary=ss1;
by Group_Name Group; run;
proc print data=ss1; run;

PROC EXPORT DATA= ss1 
            OUTFILE= "EDTable3a.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


*******************************************************;
*virus isolation;
*add one value;
PROC IMPORT OUT= WORK.virus 
            DATAFILE= "Virus Isolation2.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
proc print;RUN;
proc means data=virus min median max mean; 
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;
proc freq data=virus; tables CPE; run;

data dan.virus; set virus;
log10_ddPCR=log10(ddPCR_copies_per_ng);
log10_sgRNA=log10(sgRNA_copies_per_ul);
log10_DOI=log10(DOI); 
run;

ods rtf file="corr_virus2.rtf";
title "Correlation between ddPCR and sgRNA";
proc corr data =dan.virus spearman Fisher(biasadj=no); 
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run; *0.936;
ods rtf close;


ods rtf file="ROC_virus_ddPCR2.rtf";
title "ROC curves: ddPCR ";
ods select ParameterEstimates  Association ROCCurve ROCAssociation; 
proc logistic data=dan.virus plots(only)=(roc);  
class Patient / param=glm;
model CPE(EVENT='yes')=log10_ddPCR  / outroc=rocstats;
roc "95% CI" log10_ddPCR;
run;

data check;
set rocstats;
_SPECIF_ = (1 - _1MSPEC_);
J = _SENSIT_ + _SPECIF_ - 1;
D= Sqrt((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
*proc print data=check;run;

proc sql noprint;
create table cutoff as
select _PROB_ , J
from check
having J = max(J);
run;
proc print data=cutoff; run;  
proc print data=check; where _prob_>0.53 & _prob_<0.54; run;
ods rtf close;

data dan.virus; set dan.virus; 
if (ddPCR_copies_per_ng>758) then ddPCR_high="Y";
if (ddPCR_copies_per_ng<758)&(ddPCR_copies_per_ng^=.) then ddPCR_high="N";
proc freq; tables CPE*ddPCR_high / nocol nopercent; run;
*correct;

proc template; 
   define table Stat.Logistic.ParameterEstimates; 
      dynamic NRows; 
      column Variable GenericClassValue Response DF Estimate StdErr WaldChiSq 
      ProbChiSq StandardizedEst ExpEst Label; 
      define Estimate; 
         header = "Estimate"; 
         parent = Stat.Logistic.vbest8;
         format = 20.8 ;
      end; 
   end;
run;

ods rtf file="ROC_virus_sgRNA2.rtf";
title "ROC curves: sgRNA ";
ods select ParameterEstimates  Association ROCCurve ROCAssociation; 
proc logistic data=dan.virus plots(only)=(roc);  
class Patient / param=glm;
model CPE(EVENT='yes')=log10_sgRNA  / outroc=rocstats;
roc "95% CI" log10_sgRNA;
run;

data check;
set rocstats;
_SPECIF_ = (1 - _1MSPEC_);
J = _SENSIT_ + _SPECIF_ - 1;
D= Sqrt((1-_SENSIT_)**2 + (1-_SPECIF_)**2);
*proc print data=check;run;

proc sql noprint;
create table cutoff as
select _PROB_ , J
from check
having J = max(J);
run;
proc print data=cutoff; run;  
proc print data=check; where _prob_>0.787 & _prob_<0.788; run;
ods rtf close;

proc template;
   delete Stat.Logistic.ParameterEstimates;;
run;

data dan.virus; set dan.virus; 
if (sgRNA_copies_per_ul>25069) then sgRNA_high="Y";
if (sgRNA_copies_per_ul<25069)&(sgRNA_copies_per_ul^=.) then sgRNA_high="N";
proc freq; tables CPE*sgRNA_high / nocol nopercent; run;
*rounding issue with the borderline case of 25069;



*******************************************************;
*ISH;
PROC IMPORT OUT= WORK.ish 
            DATAFILE= "ISH.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
proc print;RUN;

ods rtf file="corr_ISH.rtf";
title "Correlation between ddPCR and qISH";
proc corr data =ish spearman pearson Fisher(biasadj=no); 
var ddPCR_avg_N_copies_ng; with red_dots_median red_dots_mean ; run; *0.759;
ods rtf close;

