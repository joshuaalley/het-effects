* Estimate Main Effects
* January 20, 2021
version 16.1
set more off

**************
* OPEN A LOG *
**************

capture log close
log using ../Output/2017-04-YouGov-mainfx.txt, text replace nomsg

*************
* LOAD DATA *
*************

use ../Output/2017-04-YouGov-prepped, clear

*********************
* RESCALE VARIABLES *
*********************

replace hawk = (hawk-1)/4 
replace intl = (intl-1)/4	
replace pid7 = (pid7-1)/6
label drop pid7
replace nat = (nat-1)/4
replace age = age/10 // now age in centuries
replace ed4 = (ed4-1)/3
label drop ed4

***********************
* SUMMARIZE VARIABLES *
***********************

local cs hawk intl pid7 nat male white age ed4 // controls
su strike2 strike ally ints dem cost ib4.region `cs'

********************
* GENERATE FIGURES *
********************

* Version Control for Programs
quietly adopath +AdoFiles
capture program drop main_fx
capture program drop main_fx_allt
capture program drop main_fx_byregion
capture program drop lvls
capture program drop fx_by_tmts
* set matsize 1000

* Effect of Alliances on Support for War
main_fx strike2 `cs'
quietly graph export ../Output/Figure01.pdf, replace
quietly graph export ../Output/Figure01.eps, replace
quietly graph export ../Output/A_Figure01.pdf, replace

* Effect of Alliances, Weighted
main_fx strike2 `cs' [pw=weight]
quietly graph export ../Output/A_Figure02.pdf, replace
	
* Effect of Alliances on Five-Point Dependent Variable
recode strike 0=0 25=1 50=2 75=3 100=4, gen(strike5)	
main_fx strike5 `cs', opts(xtitle("Support for War", margin(t+1)) xlabel(0(1)4, format(%1.0f))) format(%3.1f)
quietly graph export ../Output/A_Figure03.pdf, replace
	
* Effect of Alliances on Five-Point Dependent Variable, Weighted
main_fx strike5 `cs' [pw=weight], opts(xtitle("Support for War", margin(t+1)) xlabel(0(1)4, format(%1.0f))) format(%3.1f)
quietly graph export ../Output/A_Figure04.pdf, replace

* Footnote: No demographic controls or weights, no balancing other treatments
reg strike2 if ally==1, robust
reg strike2 if ally==0, robust
reg strike2 ally, robust	

* Same, but using 5-point scale
reg strike5 if ally==1, robust
reg strike5 if ally==0, robust
reg strike5 ally, robust	
	
* Effect of Alliances, by Region
main_fx_byregion strike2 `cs'
quietly graph export ../Output/A_Figure05.pdf, replace	

* Effect of Alliances on People with High Political Interest
main_fx strike2 `cs' if polint == 4 
quietly graph export ../Output/A_Figure06.pdf, replace

* Effect of Alliances on People with Low Political Interest
main_fx strike2 `cs' if polint < 4 
quietly graph export ../Output/A_Figure07.pdf, replace

* Effect of Alliances on Democrats
main_fx strike2 `cs' if pid3==1
quietly graph export ../Output/A_Figure08.pdf, replace

* Effect of Alliances on Independents
main_fx strike2 `cs' if pid3==2
quietly graph export ../Output/A_Figure09.pdf, replace

* Effect of Alliances on Republicans
main_fx strike2 `cs' if pid3==3
quietly graph export ../Output/A_Figure10.pdf, replace
	
* Effects of All Treatments on Support for War
main_fx_allt strike2 `cs'
quietly graph export ../Output/Figure02.pdf, replace
quietly graph export ../Output/Figure02.eps, replace
	 
* Support for War With and Without Alliances, by Context	
lvls strike2 `cs'
quietly graph export ../Output/Figure03.pdf, replace
quietly graph export ../Output/Figure03.eps, replace
	
* Effects of Alliances, By Context
fx_by_tmts strike2 `cs'
matrix fx = r(fx) // required for fx_by_baseline
quietly graph export ../Output/Figure04.pdf, replace
quietly graph export ../Output/Figure04.eps, replace

log close
