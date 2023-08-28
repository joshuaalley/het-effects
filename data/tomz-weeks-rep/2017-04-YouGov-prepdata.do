* Prepare Data
* January 20, 2021
version 16.1
set more off

*****************
* LOAD RAW DATA *
*****************

use ../Data/2017-04-YouGov-extracted, clear

**************
* TREATMENTS *
**************

recode alliance 1=0 2=1, gen(ally)
label define ally 0 "No" 1 "Yes"
label value ally ally
label var ally "T: Alliance?"

recode regime 1=0 2=1, gen(dem)
label define dem 0 "Autoc" 1 "Democ"
label values dem dem
label var dem "T: Democracy?"

recode stakes 1=1 2=0, gen(ints)
label define ints 0 "low" 1 "high"
label values ints ints
label var ints "T: Interests/Stakes?"

recode costs 1=0 2=1, gen(cost)
label define cost 0 "high" 1 "low"  // NOTE 1 = LOW
label values cost cost
label var cost "T: Costs Low?"

label var region "T: Region"

********************************
* DEPENDENT VARS AND MEDIATORS *
********************************

gen strike = 100*(pref-1)/4
label var strike "Y: Favor intervention (5 level)"
recode pref 4/5=100 1/3=0, gen(strike2)
label var strike2 "Y: Favor Intervention (binary)"
drop pref

* recode mediators to be on 0-100 scale
local ms r_hm r_hn m_oblig // mediators
recode `ms' (5=100) (4=75) (3=50) (2=25) (1=0)
foreach m in `ms' {
   label values `m'  // remove value labels
}

********************
* ATTENTION CHECKS *
********************

gen byte ck_alliance = cond(alliance == mck_alliance,1,0)
label var ck_alliance "Attn: Alliance"
gen byte ck_regime = cond((regime ~= mck_regime) & mck_regime~=.,1,0)
label var ck_regime "Attn: Democracy"
gen byte ck_ints = cond(stakes == mck_stakes,1,0)
label var ck_ints "Attn: Interests/Stakes"
gen byte ck_costs = cond(costs == mck_costs,1,0)
label var ck_costs "Attn: Costs"
gen byte ck_region = cond(region == mck_region,1,0)
label var ck_region "Attn: Region"
drop alliance regime stakes costs
drop mck_alliance mck_regime mck_stakes mck_costs mck_region 
	
**************
* COVARIATES *
**************

label drop hawk
recode hawk 1=5 2=4 3=3 4=2 5=1
label var hawk "H: Hawkishness (1-5)"
label values hawk // remove value labels

label var intl "I: Internationalism (1-5)"
label values intl // remove value labels

recode pid7 8=4 // recode "not sure" to indep
label var pid7 "P: Party 7 cat"

recode pid7 1/2=1 3/5=2 6/7=3, gen(pid3)
label var pid3 "P: Party 3 cat"

label var nat1 "E: US superior (1-5)"
label var nat2 "E: Rather be US Cit (1-5)"
label values nat1 // remove value labels
label values nat2
alpha nat1 nat2, gen(nat)
label var nat "E: Nat Idx (1-5)"

gen polint = newsint
recode polint 4=1 3=2 2=3 1=4 *=.
label var polint "A: Pol Int (1-4)"
drop newsint

****************
* DEMOGRAPHICS *
****************

recode gender 1=1 2=0, gen(male)
label var male "D: Male?"
drop gender

recode race 1=1 *=0, gen(white)
label var white "D: White?"
drop race

gen age = (2016-birth)/10 // age in decades
label var age "D: Age in decades"
drop birth

recode educ 1/2=1 3/4=2 5=3 6=4, gen(ed4)
label define ed4 1 "HS or Less" 2 "Some Coll" 3 "BA" 4 "Post-grad"
label values ed4 ed4
label var ed4 "D: Educ 4 Cat"
drop educ

************
* FINALIZE *
************

local vlist starttime ally dem ints cost region ///
   ck_alliance ck_regime ck_ints ck_costs ck_region ///
   strike strike2 r_hm r_hn m_oblig ///
   hawk intl pid7 pid3 nat1 nat2 nat polint male white age ed4 weight
order `vlist'	
compress
save ../Output/2017-04-YouGov-prepped, replace

