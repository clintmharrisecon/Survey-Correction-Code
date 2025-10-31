/*
This file uses our method on response group means from "SELECTION IN SURVEYS: USING RANDOMIZED INCENTIVES TO DETECT AND ACCOUNT FOR NONRESPONSE BIAS", downloaded on 2/24/2025 from NBER.
*/
clear
cls

set maxiter 100
set seed 123
version 18.0

*Generate data that matches moments in the paper (Table 4)

*DHLMTV have a final sample of about 3720 people who didn't receive any incentive (10,000 survey individuals in NCT survey, 93% invited to complete online (used by MTW). 40% randomly received no incentive).

local N = 10000*.93*.4
set obs `N'

//generate participation aversion percentiles (for illustrations)
gen U = _n/_N

*Set i for panel
gen id = _n

*Define ordering of participation aversion
gen survey_response0 = 0 //no one responses when requests have R=0
gen survey_response1 = U<=.38 //38% always-takers Table 4
gen survey_response2 = U<=.38+.07 //7% reminder compliers table 4
gen complier1 = survey_response1-survey_response0
gen complier2 = survey_response2-survey_response1
gen survey_ordered = survey_response1+survey_response2

*Enter data from Table 4 for outcomes for Always Takers (people who respond with no incentive and no reminder)
*We want to get comparable SEs as DHLMTV et al. for each outcome.
/*
SD = SE*sqrt(N-1)

Where SE is known from DHLMTV et al, and N is set.
*/


*Request 1 compliers:
*Rounding and odd number of people in each complier group can complicate matching DHLMTV numbers for continuous variables.

*Get equally sized groups of people to give opposite signed SD
bysort complier1 : egen medc1 = median(U) if complier1==1
gen add_sd = U < medc1 & complier1==1
gen sub_sd = U > medc1 & complier1==1
count if complier1==1 & add_sd==0 & sub_sd==0 //count people who should just get the mean
local samp_fix = r(N)
count if complier1 == 1
local N_1 = r(N) //sample of request 1 compliers ("always-takers")

sort id

*Continuous
local SD = 116*sqrt(`N_1'-1)*(`N_1'+`samp_fix')/`N_1' //adjust for SD for some people being given the mean
gen earnings_before = 3746 if complier1==1
replace earnings_before = 3746+`SD' if complier1==1 & add_sd==1 //add SD for half
replace earnings_before = 3746-`SD' if complier1==1 & sub_sd==1 //subtract SD for half
local SD = 107*sqrt(`N_1'-1)*(`N_1'+`samp_fix')/`N_1'
gen earnings_after = 3783 if complier1==1
replace earnings_after = 3783+`SD' if complier1==1 & add_sd==1 //add SD for half
replace earnings_after = 3783-`SD' if complier1==1 & sub_sd==1 //subtract SD for 
*Binary (make the appropriate percentages 1 and 0 for the relevant group)
gen earnings_large_loss = 1 if _n<=round(.38*_N)*.13
replace earnings_large_loss = 0 if _n<=round(.38*_N) & earnings_large_loss==.
gen employment_before = 1 if _n<=round(.38*_N)*.65
replace employment_before = 0 if _n<=round(.38*_N) & employment_before==.
gen employment_after = 1 if _n<=round(.38*_N)*.64
replace employment_after = 0 if _n<=round(.38*_N) & employment_after==.
gen employment_loss = 1 if _n<=round(.38*_N)*.03
replace employment_loss = 0 if _n<=round(.38*_N) & employment_loss==.

*Enter Data from Table 4 for outcomes for Reminder Compliers (people who respond with no incentive after reminders)
*Get equally sized groups of people to give opposite signed SD
bysort complier2 : egen medc2 = median(U) if complier2==1
replace add_sd = U < medc2 & complier2==1
replace sub_sd = U > medc2 & complier2==1
count if complier2==1 & add_sd==0 & sub_sd==0 //count people who should just get the mean
local samp_fix = r(N)
count if complier2 == 1
local N_2 = r(N) //sample of request 1 compliers ("always-takers")

sort id
*Continuous
local SD = 256*sqrt(`N_2'-1)*(`N_2'+`samp_fix')/`N_2'
replace earnings_before = 3244 if complier2==1
replace earnings_before = 3244+`SD' if complier2==1 & add_sd==1
replace earnings_before = 3244-`SD' if complier2==1 & sub_sd==1
local SD = 251*sqrt(`N_2'-1)*(`N_2'+`samp_fix')/`N_2'
replace earnings_after = 3257 if complier2==1
replace earnings_after = 3257+`SD' if complier2==1 & add_sd==1
replace earnings_after = 3257-`SD' if complier2==1 & sub_sd==1
*Binary (make the appropriate percentages 1 and 0 for the relevant group)
replace earnings_large_loss = 1 if _n<=round(.38*_N)+.07*_N*.12 & earnings_large_loss==.
replace earnings_large_loss = 0 if _n<=round(.38*_N)+.07*_N & earnings_large_loss==.
replace employment_before = 1 if _n<=round(.38*_N)+.07*_N*.55 & employment_before==.
replace employment_before = 0 if _n<=round(.38*_N)+.07*_N & employment_before==.
replace employment_after = 1 if _n<=round(.38*_N)+.07*_N*.55 & employment_after==.
replace employment_after = 0 if _n<=round(.38*_N)+.07*_N & employment_after==.
replace employment_loss = 1 if _n<=round(.38*_N)+.07*_N*.03 & employment_loss==.
replace employment_loss = 0 if _n<=round(.38*_N)+.07*_N & employment_loss==.



*Reshape to panel
reshape long survey_response , i(id) j(requests)
drop if survey_response == . //if individuals did not receive a rth reminder, drop that obs.


sort id requests
gen complier = survey_response[_n]-survey_response[_n-1] if id[_n] == id[_n-1]

foreach var of varlist earnings_before earnings_after earnings_large_loss employment_before employment_after employment_loss {
	replace `var' = `var'*survey_response
	replace `var' = 0 if `var' == . //in accordance with theory, be careful to explicitly code the Heckman correction with (survey_response = Z) rather than letting it default to when Y is missing.

}


*First stage
reg survey_response i.requests if requests<=1
reg survey_response i.requests if requests>=1


*Confirm local average responses match Table 4 of MTW
foreach var of varlist earnings_before earnings_after earnings_large_loss employment_before employment_after employment_loss {
	ivregress 2sls `var' (survey_response = i.requests) if requests<=1 , cluster(id) //always takers
	local mean_`var'1 = _b[survey_response]
	ivregress 2sls `var' (survey_response = i.requests) if requests>=1, cluster(id) //reminder compliers
	local mean_`var'2 = _b[survey_response]
	ivregress 2sls `var' (survey_response = i.requests) if requests!=1 , cluster(id) // pool them
}


*Heckman command gets mad over the R=0 rows
drop if requests==0
matrix outcomes = J(2,6,.)
local j = 1
*Use heckman correction
foreach outcome in "earnings_before" "earnings_after" {
	heckman `outcome' `controls', select(survey_response = i.requests) vce(cluster id)
	mat outcomes[1,`j'] = _b[_cons]
	mat outcomes[2,`j'] = _se[_cons]
	local ++j
}



foreach outcome in "earnings_large_loss" "employment_before" "employment_after" "employment_loss" {
	heckprobit `outcome' `controls', select(survey_response = i.requests) vce(cluster id)
	margins, expression(normal(xb(`outcome')))
	mat outcomes[1,`j'] = r(b)
	mat outcomes[2,`j'] = sqrt(r(V)[1,1])
	local ++j
}


matlist outcomes


forv c = 1/6 {
	if outcomes[1,`c'] < 1 {
	local B : di %5.3fc outcomes[1,`c']
	}
	else {
		local B : di %9.0fc outcomes[1,`c']
	}
	
}

forv c = 1/6 {
	if outcomes[1,`c'] < 1 {
	local se : di %5.3fc outcomes[2,`c']
	}
	else {
		local se : di %3.0fc outcomes[2,`c']
	}
	
}

esttab matrix(outcomes)

*2 step w/ bootstrap SEs (heckman command won't cluster SEs for 2 step)
cap prog drop heg_survey_bootstrap
prog heg_survey_bootstrap, eclass
syntax varlist
qui probit survey_response i.requests
qui predict imr_inner, xb
qui replace imr_inner = 0 if imr_inner==. & survey_ordered>0
gen imr = (normalden(imr_inner))/(normal(imr_inner))
reg `varlist' imr if survey_response==1
cap drop imr_inner imr 
end

local j=1
foreach outcome in "earnings_before" "earnings_after" "earnings_large_loss" "employment_before" "employment_after" "employment_loss" {
	bootstrap _b , reps(50) cluster(id) nowarn nodots nodrop seed(123): heg_survey_bootstrap `outcome'

	*Just overwrite same matrix for table
	mat outcomes[1,`j'] = _b[_cons]
	mat outcomes[2,`j'] = _se[_cons]
	local ++j
}


forv c = 1/6 {
	if outcomes[1,`c'] < 1 {
	local B : di %5.3fc outcomes[1,`c']
	}
	else {
		local B : di %9.0fc outcomes[1,`c']
	}

}

forv c = 1/6 {
	if outcomes[1,`c'] < 1 {
	local se : di %5.3fc outcomes[2,`c']
	}
	else {
		local se : di %3.0fc outcomes[2,`c']
	}
	

}


esttab matrix(outcomes)

