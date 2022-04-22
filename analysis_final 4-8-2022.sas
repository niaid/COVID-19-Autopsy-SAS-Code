libname dan "H:\Chertow\Sydney\";

****************************************************;
*data prep;

PROC IMPORT OUT= WORK.pcr 
            DATAFILE= "H:\Chertow\Sydney\used\SuppData1.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
RUN;

PROC IMPORT OUT= WORK.pats 
            DATAFILE= "H:\Chertow\Sydney\used\pats.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
RUN;

proc sort data=pcr; by patient; run;
proc sort data=pats; by patient; run;
data dan.data1; merge pcr pats; by patient; run;
data dan.data1;set dan.data1;
ddPCR_copies_per_ng= (ddPCR_copies_per_ng1+ ddPCR_copies_per_ng2)/2; run;

PROC IMPORT OUT= WORK.tissues 
            DATAFILE= "H:\Chertow\Sydney\used\Tissues2.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 	GUESSINGROWS=max;
data tissues; set tissues; if (Tissue^="");
RUN;

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

*corrections;
data dan.data2; set dan.data2;
if (Patient="P38")&(Tissue="Testis") then ddPCR_copies_per_ng=0.0049;
run;


*impute and log-transform;
proc freq data=dan.data2;tables resp;where (Resp^="Fluids");run; *1770;
proc print data=dan.data2; where (ddPCR_copies_per_ng=0)&(Resp^="Fluids"); run; *756(42.7%);
data dan.data2; set dan.data2;
if (ddPCR_copies_per_ng=0) then do;
	if (Patient="P3") then log10_ddPCR=log10(rand("Uniform")*0.000217);
	else log10_ddPCR=log10(rand("Uniform")*0.00182);
	end;
else log10_ddPCR=log10(ddPCR_copies_per_ng);
run;
data dan.data2; set dan.data2; log10_DOI=log10(DOI); run;

*ROC curve;
data dan.roc; set dan.data2;
if (sgRNA_copies_per_ul>0) then sgRNA_pos=1;
else sgRNA_pos=sgRNA_copies_per_ul;
if sgRNA_copies_per_ul^=.; run;

****************************************************;
*analysis;
ods rtf file="H:\Chertow\Sydney\used\Fig2ab.rtf";
title "log10(ddPCR): LMM";
ods select solutionf;
proc mixed data=dan.data2 covtest;   where (Resp^="Fluids");
class Patient Group Resp(ref="NonResp") tissue;
model log10_ddPCR  = Group Resp(Group) / s residual outp=res1; 
repeated tissue / subject=Patient type=cs;  run;
proc univariate data=res1 noprint; qqplot PearsonResid ; run;

title "log10(ddPCR): LMM (log10_DOI)";
ods select solutionf;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient  Resp(ref="NonResp") tissue;
model log10_ddPCR  = Resp log10_DOI(Resp) / s noint residual outp=res1; 
repeated tissue / subject=Patient type=cs; run;
ods select tests3;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient  Resp(ref="NonResp") tissue;
model log10_ddPCR  = log10_DOI|Resp / s residual outp=res1; 
repeated tissue / subject=Patient type=cs;  run;
title "Residuals:";
proc univariate data=res1 noprint; qqplot PearsonResid ; run;
ods rtf close;

*Fig 2b;
ods pdf file="H:\Chertow\Sydney\used\Fig2b.pdf";
goptions reset=all hsize=20cm vsize=15cm;
title "Fig. 2b";
proc gplot data =dan.data2;  where (Resp^="Fluids"); 
plot log10_ddPCR*log10_DOI= Resp; run;quit;
ods pdf close;


*Fig 2c;
proc sort data=dan.data2; by Resp; run;
ods rtf file="H:\Chertow\Sydney\used\Fig2c.rtf";
title "log10(ddPCR): LMM (log10_DOI)";
ods select solutionf;
proc mixed data=dan.data2 covtest;  where (Resp^="Fluids");
class Patient  Group_Name tissue;
model log10_ddPCR  = Group_Name log10_DOI(Group_Name) / s noint; 
repeated tissue / subject=Patient type=cs; by Resp; run;
ods rtf close;

*********************************************************************;

ods rtf file="H:\Chertow\Sydney\used\corr_Fig2d.rtf";
title "Correlation between ddPCR and sgRNA";
ods select FisherSpearmanCorr;
proc corr data =dan.data2 spearman Fisher(biasadj=no); where (Resp^="Fluids");
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run; *0.759;
proc sort data=dan.data2; by Group; run;
ods select FisherSpearmanCorr;
proc corr data =dan.data2 spearman Fisher(biasadj=no); by Group;
var ddPCR_copies_per_ng sgRNA_copies_per_ul; where (Resp^="Fluids");run;
proc sort data=dan.data2; by Resp; run;
ods select FisherSpearmanCorr;
proc corr data =dan.data2 spearman Fisher(biasadj=no); 
where (Resp^="Fluids"); by Resp;
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;quit;
ods rtf close;


ods rtf file="H:\Chertow\Sydney\used\corr_Fig2e.rtf";
title "Correlation between ddPCR and sgRNA";
ods select FisherSpearmanCorr;
proc corr data =dan.data2 spearman Fisher(biasadj=no); 
where (ddPCR_copies_per_ng>0)&(sgRNA_copies_per_ul>0)&(Resp^="Fluids");
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run; *0.759;
proc sort data=dan.data2; by Group; run;
ods select FisherSpearmanCorr;
proc corr data =dan.data2 spearman Fisher(biasadj=no); by Group;
where (ddPCR_copies_per_ng>0)&(sgRNA_copies_per_ul>0)&(Resp^="Fluids");
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;
proc sort data=dan.data2; by Resp; run;
ods select FisherSpearmanCorr;
proc corr data =dan.data2 spearman Fisher(biasadj=no); by Resp;
where (ddPCR_copies_per_ng>0)&(sgRNA_copies_per_ul>0)&(Resp^="Fluids");
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run;quit;
ods rtf close;


****************************************************;


ods rtf file="H:\Chertow\Sydney\used\Fig2f.rtf";
title "ROC curves: predict positive sgRNA using ddPCR";
ods select ROCCurve ROCAssociation; 
proc logistic data=dan.roc ; where (Resp^="Fluids");
class sgRNA_pos(ref='0' param=ref) / param=glm;
model sgRNA_pos=log10_ddPCR;
roc "predicting positive sgRNA using ddPCR" log10_ddPCR;
run;
ods rtf close;


ods rtf file="H:\Chertow\Sydney\used\Fig2f_cutoff.rtf";
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



*Extended Data Table 3a;
proc sort data=dan.data2; by Group_Name Group; run;
proc means data=dan.data2 mean std n q1 median q3 ; where (Resp^="Fluids");
var ddPCR_copies_per_ng; ods output summary=ss1;
by Group_Name Group; run;
PROC EXPORT DATA= ss1 
            OUTFILE= "H:\Chertow\Sydney\EDTable3a.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;




*******************************************************;
*virus isolation;
*add one value;
PROC IMPORT OUT= WORK.virus 
            DATAFILE= "H:\Chertow\Sydney\used\Virus Isolation2.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
run;

data dan.virus; set virus;
log10_ddPCR=log10(ddPCR_copies_per_ng);
log10_sgRNA=log10(sgRNA_copies_per_ul);
log10_DOI=log10(DOI); 
run;

ods rtf file="H:\Chertow\Sydney\used\Fig3b.rtf";
title "Correlation between ddPCR and sgRNA";
ods select FisherSpearmanCorr;
proc corr data =dan.virus spearman Fisher(biasadj=no); 
var ddPCR_copies_per_ng sgRNA_copies_per_ul; run; *0.936;
ods rtf close;


ods rtf file="H:\Chertow\Sydney\used\Fig3c.rtf";
title "ROC curves: ddPCR ";
ods select ParameterEstimates  ROCCurve ROCAssociation; 
proc logistic data=dan.virus plots(only)=(roc);  
class Patient / param=glm;
model CPE(EVENT='yes')=log10_ddPCR  / outroc=rocstats;
roc "predicting CPE using ddPCR" log10_ddPCR;
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

ods rtf file="H:\Chertow\Sydney\used\Fig3d.rtf";
title "ROC curves: sgRNA ";
ods select ParameterEstimates  ROCCurve ROCAssociation; 
proc logistic data=dan.virus plots(only)=(roc);  
class Patient / param=glm;
model CPE(EVENT='yes')=log10_sgRNA  / outroc=rocstats;
roc "predicting CPE using sgRNA" log10_sgRNA;
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





*******************************************************;
*ISH;
PROC IMPORT OUT= WORK.ish 
            DATAFILE= "H:\Chertow\Sydney\used\ISH.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
	 guessingrows=max;
proc print;RUN;

ods rtf file="H:\Chertow\Sydney\used\corr_ISH.rtf";
title "Correlation between ddPCR and qISH";
ods select FisherSpearmanCorr;
proc corr data =ish spearman pearson Fisher(biasadj=no); 
var ddPCR_avg_N_copies_ng; with red_dots_median red_dots_mean ; run; *0.759;
ods rtf close;
