log using "replication appendix.txt", text replace


*Load Tunisia data
use "Bush Prather Partisan Side-taking Appendix Tunisia.dta", clear

/*Table 4*/
ologit give_aid treat_islamist##top2 age sex employ rural polint polknowledge vote_parl


/*Table 5*/
ologit trade_with treat_islamist##top2 age sex employ rural polint polknowledge vote_parl


/*Table 6*/
/*Rows*/
ttest give_aid if treat_islamist==1, by(polislam_bin)
ttest give_aid if treat_islamist==0, by(polislam_bin)
/*Columns*/
ttest give_aid if polislam_bin==0, by(treat_islamist)
ttest give_aid if polislam_bin==1, by(treat_islamist)


/*Table 7*/
/*Rows*/
ttest trade_with if treat_islamist==1, by(polislam_bin)
ttest trade_with if treat_islamist==0, by(polislam_bin)
/*Columns*/
ttest trade_with if polislam_bin==0, by(treat_islamist)
ttest trade_with if polislam_bin==1, by(treat_islamist)


/*Table 8*/
/*Rows*/
ttest engage if treat_islamist==1 & polknowledge<3, by(top2)
ttest engage if treat_islamist==0 & polknowledge<3, by(top2)
/*Columns*/
ttest engage if top2==1 & polknowledge<3, by(treat_islamist)
ttest engage if top2==0 & polknowledge<3, by(treat_islamist)


/*Table 9*/
/*Rows*/
ttest engage if treat_islamist==1 & polknowledge>=3, by(top2)
ttest engage if treat_islamist==0 & polknowledge>=3, by(top2)
/*Columns*/
ttest engage if top2==1 & polknowledge>=3, by(treat_islamist)
ttest engage if top2==0 & polknowledge>=3, by(treat_islamist)


/*Table 10*/
gen give_aidbin = cond(give_aid==1 | give_aid==2, 0, 1) if give_aid!=.

/*Rows*/
prtest give_aidbin if treat_islamist==1, by(top2)
prtest give_aidbin if treat_islamist==0, by(top2)
/*Columns*/
prtest give_aidbin if top2==1 , by(treat_islamist)
prtest give_aidbin if top2==0 , by(treat_islamist)


/*Table 11*/
gen trade_withbin = cond(trade_with==1 | trade_with==2, 0, 1) if trade_with!=.

/*Rows*/
prtest trade_withbin if treat_islamist==1, by(top2)
prtest trade_withbin if treat_islamist==0 , by(top2)
/*Columns*/
prtest trade_withbin if top2==1, by(treat_islamist)
prtest trade_withbin if top2==0, by(treat_islamist)


/*Table 12*/
gen personally_benefit2 = cond(personally_benefit>2.5,1,0) if personally_benefit!=.

/*Rows*/
ttest give_aid if treat_islamist==1 & personally_benefit2==1, by(top2)
ttest give_aid if treat_islamist==0 & personally_benefit2==1, by(top2)
/*Columns*/
ttest give_aid if top2==0 & personally_benefit2==1, by(treat_islamist)
ttest give_aid if top2==1 & personally_benefit2==1, by(treat_islamist)


/*Table 13*/
/*Rows*/
ttest trade_with if treat_islamist==1 & personally_benefit2==1, by(top2)
ttest trade_with if treat_islamist==0 & personally_benefit2==1, by(top2)
/*Columns*/
ttest trade_with if top2==0 & personally_benefit2==1, by(treat_islamist)
ttest trade_with if top2==1 & personally_benefit2==1, by(treat_islamist)


/*Table 28*/
/*Tunisia model*/
reg engage personally_benefit economy_improve top2 educ age sex employ rural polint polknowledge vote_parl i.treat_islamist


/*Table 29*/
/*Rows*/
ttest engage if treat_islamist==1 & educ<4, by(top2)
ttest engage if treat_islamist==0 & educ<4, by(top2)
/*Columns*/
ttest engage if top2==1 & educ<4, by(treat_islamist)
ttest engage if top2==0 & educ<4, by(treat_islamist)



*Load United States 2016 data
use "Bush Prather Partisan Side-taking Appendix USA 2016.dta", clear


/*Table 14*/
ologit ipe_support treat_germrus##w2_vote_hill age woman income employ polint polknowledge vote if treat_ipetype==1 & treat_ipesidetaking==1 


/*Table 15*/
ologit ipe_support treat_germrus##w2_vote_hill age woman income employ polint polknowledge vote if treat_ipetype==2 & treat_ipesidetaking==1 


/*Table 16 - FDI*/
/*Rows*/
ttest ipe_support if treat_germrus==0 & treat_ipetype==1 & treat_ipesidetaking==1 & party<=2, by(party)
ttest ipe_support if treat_germrus==1 & treat_ipetype==1 & treat_ipesidetaking==1 & party<=2, by(party)
/*Columns*/ 
ttest ipe_support if party==1 & treat_ipetype==1 & treat_ipesidetaking==1 & party<=2, by(treat_germrus)
ttest ipe_support if party==2 & treat_ipetype==1 & treat_ipesidetaking==1 & party<=2, by(treat_germrus)


/*Table 17 - Trade*/
/*Rows*/
ttest ipe_support if treat_germrus==0 & treat_ipetype==2 & treat_ipesidetaking==1 & party<=2, by(party)
ttest ipe_support if treat_germrus==1 & treat_ipetype==2 & treat_ipesidetaking==1 & party<=2, by(party)
/*Columns*/
ttest ipe_support if party==1 & treat_ipetype==2 & treat_ipesidetaking==1 & party<=2, by(treat_germrus)
ttest ipe_support if party==2 & treat_ipetype==2 & treat_ipesidetaking==1 & party<=2, by(treat_germrus)


/*Table 18*/
/*Rows*/
ttest ipe_support if treat_germrus==0 & treat_ipesidetaking==1 & polknowledge<3, by(w2_vote_hill)
ttest ipe_support if treat_germrus==1 & treat_ipesidetaking==1 & polknowledge<3, by(w2_vote_hill)
/*Columns*/
ttest ipe_support if w2_vote_hill==0 & treat_ipesidetaking==1 & polknowledge<3, by(treat_germrus)
ttest ipe_support if w2_vote_hill==1 & treat_ipesidetaking==1 & polknowledge<3, by(treat_germrus)


/*Table 19*/
/*Rows*/
ttest ipe_support if treat_germrus==0 & treat_ipesidetaking==1 & polknowledge==3, by(w2_vote_hill)
ttest ipe_support if treat_germrus==1 & treat_ipesidetaking==1 & polknowledge==3, by(w2_vote_hill)
/*Columns*/ 
ttest ipe_support if w2_vote_hill==0 & treat_ipesidetaking==1 & polknowledge==3, by(treat_germrus)
ttest ipe_support if w2_vote_hill==1 & treat_ipesidetaking==1 & polknowledge==3, by(treat_germrus)


/*Table 20*/
gen ipe_supportbin = cond(ipe_support==1 | ipe_support==2, 0, 1) if ipe_support!=.

/*Rows*/
prtest ipe_supportbin if treat_germrus==0 & treat_ipetype==1 & treat_ipesidetaking==1, by(w2_vote_hill)
prtest ipe_supportbin if treat_germrus==1 & treat_ipetype==1 & treat_ipesidetaking==1, by(w2_vote_hill)
/*Columns*/
prtest ipe_supportbin if w2_vote_hill==0 & treat_ipetype==1 & treat_ipesidetaking==1, by(treat_germrus)
prtest ipe_supportbin if w2_vote_hill==1 & treat_ipetype==1 & treat_ipesidetaking==1, by(treat_germrus)


/*Table 21*/
/*Rows*/
prtest ipe_supportbin if treat_germrus==0 & treat_ipetype==2 & treat_ipesidetaking==1, by(w2_vote_hill)
prtest ipe_supportbin if treat_germrus==1 & treat_ipetype==2 & treat_ipesidetaking==1, by(w2_vote_hill)
/*Columns*/ 
prtest ipe_supportbin if w2_vote_hill==0 & treat_ipetype==2 & treat_ipesidetaking==1, by(treat_germrus)
prtest ipe_supportbin if w2_vote_hill==1 & treat_ipetype==2 & treat_ipesidetaking==1, by(treat_germrus)


/*Table 22*/
/*Rows*/
ttest ipe_support if treat_germrus==0 & treat_ipesidetaking==1 & ipe_benefit_personal>=3, by(w2_vote_hill)
ttest ipe_support if treat_germrus==1 & treat_ipesidetaking==1 & ipe_benefit_personal>=3, by(w2_vote_hill)
/*Columns*/ 
ttest ipe_support if w2_vote_hill==0 & treat_ipesidetaking==1 & ipe_benefit_personal>=3, by(treat_germrus)
ttest ipe_support if w2_vote_hill==1 & treat_ipesidetaking==1 & ipe_benefit_personal>=3, by(treat_germrus)


/*Table 23*/
/*Rows*/
ttest ipe_support if treat_ipesidetaking==0 & treat_germrus==1 & trump_germany==0, by(w2_vote_hill)
ttest ipe_support if treat_ipesidetaking==0 & treat_germrus==0 & trump_russia==2, by(w2_vote_hill)
/*Columns*/ 
ttest ipe_support if treat_ipesidetaking==0 & w2_vote_hill==1 & ((treat_germrus==1 & trump_germany==0) | (treat_germrus==0 & trump_russia==2)), by(treat_germrus)
ttest ipe_support if treat_ipesidetaking==0 & w2_vote_hill==0 & ((treat_germrus==1 & trump_germany==0) | (treat_germrus==0 & trump_russia==2)), by(treat_germrus)


/*Figure 1*/
/*7a - Clinton*/
logit repbenonly treat_ipesidetaking##treat_germrus i.treat_ipetype if w2_vote_hill==1
margins treat_germrus#treat_ipesidetaking, asbalanced
marginsplot, x(treat_ipesidetaking)  plot2opts(connect(no)) plot1opts(connect(no)) title("") ytitle("Probability of Viewing Republicans as Benefiting") ylabel(0(.1)0.7,grid) yscale(range(0 0.7))
/*7b - Trump*/
logit repbenonly treat_ipesidetaking##treat_germrus i.treat_ipetype if w2_vote_hill==0
margins treat_germrus#treat_ipesidetaking, asbalanced
marginsplot, x(treat_ipesidetaking)  plot2opts(connect(no)) plot1opts(connect(no)) title("") ytitle("Probability of Viewing Republicans as Benefiting") ylabel(0(.1)0.7,grid) yscale(range(0 0.7))


/*Table 27*/
/*Rows*/
ttest ipe_support if treat_ipesidetaking==1 & treat_germrus==1 & ipe_contact==1, by(w2_vote_hill)
ttest ipe_support if treat_ipesidetaking==1 & treat_germrus==0 & ipe_contact==1, by(w2_vote_hill)
/*Columns*/ 
ttest ipe_support if treat_ipesidetaking==1 & ipe_contact==1 & w2_vote_hill==1, by(treat_germrus)
ttest ipe_support if treat_ipesidetaking==1 & ipe_contact==1 & w2_vote_hill==0, by(treat_germrus)

    
/*Table 28*/
/*U.S., 2016 model*/
reg ipe_support ipe_benefit_personal ipe_benefit_usecon w2_vote_hill i.trump_russia i.trump_germany educ age woman employ income  polint polknowledge vote i.treat_germrus if treat_ipesidetaking==0


/*Table 30*/
/*Rows*/
ttest ipe_support if treat_germrus==0 & treat_ipesidetaking==1 & educ>4, by(w2_vote_hill)
ttest ipe_support if treat_germrus==1 & treat_ipesidetaking==1 & educ>4, by(w2_vote_hill)
/*Columns*/
ttest ipe_support if w2_vote_hill==0 & treat_ipesidetaking==1 & educ>4, by(treat_germrus)
ttest ipe_support if w2_vote_hill==1 & treat_ipesidetaking==1 & educ>4, by(treat_germrus)



*Load United States 2018 data
use "Bush Prather Partisan Side-taking Appendix USA 2018.dta", clear

/*Table 24*/
/*Rows*/
ttest engage if treat_support==1 & manip_pass==1, by(vote_trump)
ttest engage if treat_support==0 & manip_pass==1, by(vote_trump)
/*Columns*/
ttest engage if vote_trump==1 & manip_pass==1, by(treat_support)
ttest engage if vote_trump==0 & manip_pass==1, by(treat_support)

/*Figure 2*/
reg dembenonly i.treat_support i.treat_trade 
margins treat_support, asbalanced
#delimit ;
marginsplot,
plot1opts(mcolor(black) lcolor(black) connect(none))
ci1opts(lcolor(black)) 
graphregion(fcolor(white) lcolor(white))
plotregion(fcolor(white) lstyle(none) lcolor(white) ilstyle(none))
title("", color(black))
subtitle("", color(black))
ytitle("Probability of Viewing Democrats as Benefitting")
yscale(range(0 .5))
xscale(range(-0.9 1.9))
ylabel(0(.1).5, labsize(small) glcolor(gs14) angle(horizontal))
xlabel(,valuelabel labsize(medsmall) glcolor(white) angle(horizontal)) 
xtitle("Foreign Country's Side-taking", margin(medsmall));
#delimit cr;

/*Table 25*/
/*Rows*/
ttest ipe_benefit_personal if treat_support==1 & manip_pass==1, by(vote_trump)
ttest ipe_benefit_personal if treat_support==0 & manip_pass==1, by(vote_trump)
/*Columns*/
ttest ipe_benefit_personal if vote_trump==1 & manip_pass==1, by(treat_support)
ttest ipe_benefit_personal if vote_trump==0 & manip_pass==1, by(treat_support)


/*Table 26*/
/*Rows*/
ttest ipe_benefit_usecon if treat_support==1 & manip_pass==1, by(vote_trump)
ttest ipe_benefit_usecon if treat_support==0 & manip_pass==1, by(vote_trump)
/*Columns*/
ttest ipe_benefit_usecon if vote_trump==1 & manip_pass==1, by(treat_support)
ttest ipe_benefit_usecon if vote_trump==0 & manip_pass==1, by(treat_support)

