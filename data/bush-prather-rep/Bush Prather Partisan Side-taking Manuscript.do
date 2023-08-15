*Load Tunisia data
use "Bush Prather Partisan Side-taking Manuscript Tunisia.dta", clear

*Tunisia tables
/*Table 2*/
/*Rows*/
ttest give_aid if treat_islamist==1, by(top2)
ttest give_aid if treat_islamist==0, by(top2)
/*Columns*/
ttest give_aid if top2==1, by(treat_islamist)
ttest give_aid if top2==0, by(treat_islamist)

/*Table 3*/
/*Rows*/
ttest trade_with if treat_islamist==1, by(top2)
ttest trade_with if treat_islamist==0, by(top2)
/*Columns*/
ttest trade_with if top2==1, by(treat_islamist)
ttest trade_with if top2==0, by(treat_islamist)

/*Table 4*/
/*Rows*/
ttest personally_benefit if treat_islamist==1, by(top2)
ttest personally_benefit if treat_islamist==0, by(top2)
/*Columns*/
ttest personally_benefit if top2==1, by(treat_islamist)
ttest personally_benefit if top2==0, by(treat_islamist)

/*Table 5*/
/*Rows*/
ttest economy_improve if treat_islamist==1, by(top2)
ttest economy_improve if treat_islamist==0, by(top2)
/*Columns*/
ttest economy_improve if top2==1, by(treat_islamist)
ttest economy_improve if top2==0, by(treat_islamist)


*Tunisia figures
/*Figure 1*/
tab influ_france
tab influ_usa
tab influ_qatar
tab influ_uae
tab influ_algeria
tab influ_saudi

/*Figure 3*/
tab trade_with
tab give_aid



*Load United States data
use "Bush Prather Partisan Side-taking Manuscript USA.dta", clear


*United States Tables
/*Table 6*/
/*Rows*/
ttest ipe_support if treat_germrus==0 & treat_ipetype==1 & treat_ipesidetaking==1, by(w2_vote_hill)
ttest ipe_support if treat_germrus==1 & treat_ipetype==1 & treat_ipesidetaking==1, by(w2_vote_hill)
/*Columns*/
ttest ipe_support if w2_vote_hill==0 & treat_ipetype==1 & treat_ipesidetaking==1, by(treat_germrus)
ttest ipe_support if w2_vote_hill==1 & treat_ipetype==1 & treat_ipesidetaking==1, by(treat_germrus)

/*Table 7*/
/*Rows*/
ttest ipe_support if treat_germrus==0 & treat_ipetype==2 & treat_ipesidetaking==1, by(w2_vote_hill)
ttest ipe_support if treat_germrus==1 & treat_ipetype==2 & treat_ipesidetaking==1, by(w2_vote_hill)
/*Columns*/ 
ttest ipe_support if w2_vote_hill==0 & treat_ipetype==2 & treat_ipesidetaking==1, by(treat_germrus)
ttest ipe_support if w2_vote_hill==1 & treat_ipetype==2 & treat_ipesidetaking==1, by(treat_germrus)

/*Table 8*/
/*Rows*/
ttest ipe_benefit_personal if treat_germrus==0 & treat_ipetype<3 & treat_ipesidetaking==1, by(w2_vote_hill)
ttest ipe_benefit_personal if treat_germrus==1 & treat_ipetype<3 & treat_ipesidetaking==1, by(w2_vote_hill)
/*Columns*/
ttest ipe_benefit_personal if w2_vote_hill==0 & treat_ipetype<3 & treat_ipesidetaking==1, by(treat_germrus)
ttest ipe_benefit_personal if w2_vote_hill==1 & treat_ipetype<3 & treat_ipesidetaking==1, by(treat_germrus)

/*Table 9*/
/*Rows*/
ttest ipe_benefit_usecon if treat_germrus==0 & treat_ipetype<3 & treat_ipesidetaking==1, by(w2_vote_hill)
ttest ipe_benefit_usecon if treat_germrus==1 & treat_ipetype<3 & treat_ipesidetaking==1, by(w2_vote_hill)
/*Columns*/
ttest ipe_benefit_usecon if w2_vote_hill==0 & treat_ipetype<3 & treat_ipesidetaking==1, by(treat_germrus)
ttest ipe_benefit_usecon if w2_vote_hill==1 & treat_ipetype<3 & treat_ipesidetaking==1, by(treat_germrus)


*United States Figures
/*Figure 2*/
tab influence_russia 
tab influence_china 
tab influence_mexico 
tab influence_germany 
tab influence_uk 
tab influence_canada 

/*Figure 4*/
tab ipe_support if treat_ipetype==1 /*investment*/
tab ipe_support if treat_ipetype==2 /*trade*/

/*Figure 5a - Clinton */
reg ipe_support treat_ipesidetaking##treat_germrus if w2_vote_hill==1 & treat_ipetype==1
margins treat_germrus#treat_ipesidetaking, asbalanced
marginsplot, x(treat_ipesidetaking)  plot2opts(connect(no)) plot1opts(connect(no)) title("") ytitle("Linear Prediction of Support for Engagement") ylabel(1(1)4,grid) yscale(range(1 4))
/*Figure 5b - Trump*/
reg ipe_support treat_ipesidetaking##treat_germrus  if w2_vote_hill==0 & treat_ipetype==1
margins treat_germrus#treat_ipesidetaking, asbalanced
marginsplot, x(treat_ipesidetaking)  plot2opts(connect(no)) plot1opts(connect(no)) title("") ytitle("Linear Prediction of Support for Engagement") ylabel(1(1)4,grid) yscale(range(1 4))

/*Figure 6a - Clinton*/
reg ipe_support treat_ipesidetaking##treat_germrus if w2_vote_hill==1 & treat_ipetype==2
margins treat_germrus#treat_ipesidetaking, asbalanced
marginsplot, x(treat_ipesidetaking)  plot2opts(connect(no)) plot1opts(connect(no)) title("") ytitle("Linear Prediction of Support for Engagement") ylabel(1(1)4,grid) yscale(range(1 4))
/*Figure 6b - Trump*/
reg ipe_support treat_ipesidetaking##treat_germrus  if w2_vote_hill==0 & treat_ipetype==2
margins treat_germrus#treat_ipesidetaking, asbalanced
marginsplot, x(treat_ipesidetaking)  plot2opts(connect(no)) plot1opts(connect(no)) title("") ytitle("Linear Prediction of Support for Engagement") ylabel(1(1)4,grid) yscale(range(1 4))

/*Figure 7a - Clinton*/
logit dembenonly treat_ipesidetaking##treat_germrus i.treat_ipetype if w2_vote_hill==1 & treat_ipetype<3
margins treat_germrus#treat_ipesidetaking, asbalanced
marginsplot, x(treat_ipesidetaking)  plot2opts(connect(no)) plot1opts(connect(no)) title("") ytitle("Probability of Viewing Democrats as Benefiting") ylabel(0(.1)0.5,grid) yscale(range(0 0.5))
/*Figure 7b - Trump*/
logit dembenonly treat_ipesidetaking##treat_germrus i.treat_ipetype if w2_vote_hill==0 & treat_ipetype<3
margins treat_germrus#treat_ipesidetaking, asbalanced
marginsplot, x(treat_ipesidetaking)  plot2opts(connect(no)) plot1opts(connect(no)) title("") ytitle("Probability of Viewing Democrats as Benefiting") ylabel(0(.1)0.5,grid) yscale(range(0 0.5))



