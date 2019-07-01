proc delete data=_all_;run;

%macro ODSOff(); /* call prior to BY-group processing */
ods graphics off;
ods exclude all; /* all open destinations */
ods results off; /* no updates to tree view */
options nonotes; /* optional, but sometimes useful */
%mend;

%macro ODSOn(); /* call after BY-group processing */
ods graphics on;
ods exclude none;
ods results on;
options notes;
%mend;

%include "C:\Users\André Felipe\SkyDrive\Projeto-Beta\Scripts\score-fisher.sas";

%macro montecarlo(W=);
	%do j=1 %to 3;
	%let phi = %scan(&p., &j.);
		%do i = 1 %to 5;
		%let mu = %sysevalf(%scan(&m., &i.) / 100);
			%do n = 20 %to 100 %by 20; 
			/*Simulando dados*/
					data simulation;
						call streaminit(1299);
						do b=1 to &W;
							do i=1 to &n;
								x = rand("BETA", &mu*&phi, (1-&mu)*&phi);
								if(find(x,1) = 1) then x = 0.999999999;
								output;
							end;
						end;
					run;
			/*Obtendo estimativas.*/
					proc nlmixed data=simulation df=1e6 technique=quanew update=bfgs;
						parms mu=&mu, phi=&phi;
						loglikehood=logpdf('BETA',x, mu*phi, (1-mu)*phi);
						model x ~ general(loglikehood);
						contrast "H0" mu-&mu, phi-&phi;
						by b;
						ods output Parameters=TRV1(drop=mu phi rename=(NegLogLike=lik_H0));
						ods output IterHistory=TRV2(keep=iter negloglike rename=(NegLogLike=lik_H1));
						ods output Contrasts=Wald(keep=b ProbF rename=(probf=pvalue_wald));
					run;
			/*Ajustando teste da razão de verossimilhança.*/
					data trv2;
						set trv2 end=last;
						if not last then set trv2(firstobs=2 keep=iter rename=(iter=niter));
						if niter eq 1 or last then output;
						drop niter iter;
					run;
					data trv(keep=b pvalue_trv);
						merge trv1 trv2;
						LR = 2*(-lik_H1 + lik_H0);
						pvalue_trv = 1-CDF('CHISQUARE', LR, 2); 
					run;
			/*Ajustando teste escore.*/
					proc iml;
						load module = _all_; *lendo as funções escore e informção de fisher;
						use simulation;
						read all;
						close simulation;
						th 				= {&mu, &phi};
						p1 				= do(1,  &W*&n, &n);
						p2 				= do(&n, &W*&n, &n);
						pvalue_escore	= j(&W,1);
						do k=1 to &W;
							yy 					= x[p1[k]:p2[k]];
							T3 		 			= T(U(th,yy))*inv(If(th,yy))*U(th,yy);
							pvalue_escore[k] 	= 1 - CDF("CHISQUARE", T3, 2);
						end;
						create Escore var {"pvalue_escore"}; /* name the vars */
						append; /* write the data */
						close Escore;
					quit;
					data results;
						merge trv wald escore;
						trv    = (pvalue_trv <= 0.05);
						wald   = (pvalue_wald <= 0.05);
						escore = (pvalue_escore <= 0.05);
					run;
					proc freq data=results;
						tables trv wald escore / nocum binomial(level="1" p=0.05);
						ods output binomial=dados&n(keep = table nvalue1 label1 where=(Label1='Proportion'));
					run;
					%if &n=20 and &phi=5 and &mu=0.05 %then %do; 
						data dados(drop=label1);
							set dados&n;
							n 	= &n;
							mu  = &mu;
							phi = &phi;
						run;
					%end;
					%else %do;
						data dados&n(drop=label1);
							set dados&n;
							n = &n;
							mu  = &mu;
							phi = &phi;
						run;
						data dados;
							set dados dados&n;
						run;
					%end;
			%end;
		%end;
	%end;
%mend;

%let m = 5 25 50 75 95;
%let p = 5 15 50;

%ODSOff
%montecarlo(W=1000);
%ODSOn;


proc export data = dados 
            outfile= "C:\Users\André Felipe\SkyDrive\Projeto-Beta\simulacao2.txt"
            dbms=dlm replace;
			delimiter=",";
			putnames=yes;
run;
