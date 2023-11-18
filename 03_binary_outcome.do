/*
This file estimates the survey selection correction for a binary outcome. 
*/

local gender_pop = 10000

clear
cls
set obs `gender_pop'
gen Female=1

set seed 123

*Draw unobserved entrepreneurship and willingness to fill out related survey from joint distribution for women
mat Sig = (1, .25, .1 \ .25, 1, .1 \ .1, .1, 1)
drawnorm u_ent u_survey STEM, cov(Sig)
*Note: STEM is some predictor of both survey response and entrepreneurship, e.g. STEM major

*Do the same for men
set obs `=`gender_pop'*2'
replace Female=0 if Female==.
mat Sig = (1, .75, .2 \ .75, 1, .1 \ .2, .1, 1)
drawnorm m_u_ent m_u_survey m_STEM, cov(Sig)

replace u_ent = m_u_ent if Female==0
replace u_survey = m_u_survey if Female==0 
replace STEM = m_STEM if Female==0

drop m_*

*Change base rates to get stuff that sorta looks like reality
replace u_ent = u_ent-0.5-0.5*Female //Women less entrepreneurial than men
replace u_survey = u_survey-1 //most people don't respond to survey
replace STEM = STEM-0.5-0.5*Female //Women less STEMy


replace STEM = STEM>0

*Generate linear ent measure
rename u_ent ent_linear
*Generate binary ent measure
gen ent_binary = ent_linear>0

*Generate 2 reminders, sufficient to perform method and checks on method
gen survey_1_notification = u_survey>0
gen survey_2_notification = u_survey+0.3>0
gen survey_3_notification = u_survey+0.5>0
gen survey_ordered = survey_1_notification+survey_2_notification+survey_3_notification


bysort Female: sum survey* ent_binary ent_linear

*Generate responses to everyone getting survey reminders version
gen ent_linear_observed = ent_linear if survey_ordered>0
gen ent_binary_observed = ent_binary if survey_ordered>0


*Randomize intensity of follow-up, following DiNardo et al 2021. We approximate this ex post, by throwing away responses received after reminders for random half of subjects.
*1. Generate binary "intensity" exclusion restriction.
gen rem_exclusion = rnormal()>0
*2. Throw away late responses from those not assigned to "intense follow-up"
gen rem_exclusion_response = survey_3_notification if rem_exclusion==1
replace rem_exclusion_response = survey_1_notification if rem_exclusion==0

gen ent_linear_observed_rem = ent_linear_observed if rem_exclusion_response==1 //keep responses received by group-specific cutoff time
gen ent_binary_observed_rem = ent_binary_observed if rem_exclusion_response==1 //keep responses received by group-specific cutoff time


*I want vars demeaned so constants can be compared across models
foreach var of varlist STEM {
qui sum `var'
replace `var' = `var'-r(mean)
}

keep if Female==0 //I coded more extreme selection for men so that gaps male-female gaps will be overstated without the correction.

*local controls = "STEM" //if you want
local controls = "" //simpler for graphs

*-------------------------------------------------------------------------------
* 1. Raw population average 
*-------------------------------------------------------------------------------
*Compare results to:
*1. Infeasible case with no selection
reg ent_binary
sum ent_binary

*Code up MLE for binary outcomes
capture program drop ml_stacked_heckman_bivariate
program ml_stacked_heckman_bivariate
    version 16
    args todo b lnf
	
	tempvar xb_main
	mleval `xb_main' = `b', eq(1)
	forv k=1/$notification_count { //try to pass in the number of notifications from outside or back it out from parameters inside
	tempvar xbs`k'
	mleval `xbs`k'' = `b', eq(`=`k'+1')
	}
	tempvar tau
	mleval `tau' = `b', eq(`=$notification_count+2') //replace first 3 with number of notifications

	tempname rho
	scalar `rho' = (expm1(2*`tau'))  /  (exp(2*`tau')+1)

	*Loop over K bivariate Heckman selection models (only 1 rho)
	quietly replace `lnf' = 0 //prime likelihood
	*First notification equation (always includes all controls)
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs1')/1)) if ${ML_y`=1+1'}==0   // Don't respond after 1st notice
	quietly replace `lnf' =  `lnf'+log(binormal(`xb_main',`xbs1',`rho')) if ${ML_y`=1+1'}==1 & $ML_y1==1 //Respond and have response=1
	quietly replace `lnf' =  `lnf'+log(binormal(-`xb_main',`xbs1',-`rho')) if ${ML_y`=1+1'}==1 & $ML_y1==0 //Respond and have response=0
	
	*All other reminders, may or may not include controls (covariate-specific responses to reminders)
	forv k=2/$notification_count {
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs1'+`xbs`k'')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+log(binormal(`xb_main',`xbs1'+`xbs`k'',`rho')) if ${ML_y`=`k'+1'}==1 & $ML_y1==1 //Respond and have response=1
	quietly replace `lnf' =  `lnf'+log(binormal(-`xb_main',`xbs1'+`xbs`k'',-`rho')) if ${ML_y`=`k'+1'}==1 & $ML_y1==0 //Respond and have response=0	
	}
		
end

sum survey_ordered
local notification_count = r(max)
global notification_count = r(max)
di $notification_count
*Always have (add controls you want here)
local Equations "`Equations' (Selection1: survey_1_notification = `controls')"
forv k = 2/`notification_count' {
	local Equations "`Equations' (Selection`k': survey_`k'_notification = `controls')"
}
di "`Equations'"


ml model lf0 ml_stacked_heckman_bivariate (Binary_Outcome: ent_binary = `controls') `Equations' /athrho,  missing vce(robust)
*ml init initial, copy
ml maximize

estimates store surv_corr

margins, expression(normal(xb(Binary_Outcome)))


sum ent_binary


*Make figure showing means by notification group w/ predicted values.

*Generate necessary objects for graph
local N_graphs = 1000 //set number of data points for illustrative graphs, doesn't need to equal N
capture set obs `N_graphs'

capture drop u
gen u = (_n)/`N_graphs' if _n<=`N_graphs'

capture drop rem_*_*
local ko = 0
local share_atleast_4 = 0 //1 more than highest value never happens
forv k=3(-1)1 {
	local ++ko
	di `ko'
count if survey_ordered>=`k'
local N_`k' = r(N)
count
local N = r(N)
local share_atleast_`k' = `N_`k''/`N'


reg ent_binary_observed if survey_ordered==`k'
gen rem_`k'_mean = _b[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''
gen rem_`k'_ciup = _b[_cons]+1.96*_se[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''
gen rem_`k'_cidown = _b[_cons]-1.96*_se[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''

local sample_str = "`sample_str' (line rem_`k'_mean u, lw(*1) lcolor(black)  xaxis(1 2)) (line rem_`k'_ciup u, lpattern(dash) lcolor(black)) (line rem_`k'_cidown u, lpattern(dash) lcolor(black))"

local notice_str = "`notice_str' `share_atleast_`k''"
local notice_lab_str = "`notice_lab_str' `share_atleast_`k'' `""`ko'"""
}
reg ent_binary_observed
capture drop rem_mean rem_ci*
gen rem_mean = _b[_cons]
gen rem_ciup = _b[_cons]+1.96*_se[_cons] 
gen rem_cidown = _b[_cons]-1.96*_se[_cons]
local uncorrection_str = "`uncorrection_str' (line rem_mean u, lw(*1) lcolor(gs8)) (line rem_ciup u, lpattern(dash) lcolor(gs8)) (line rem_cidown u, lpattern(dash) lcolor(gs8))"


*Add corrected prediction
estimates restore surv_corr
capture drop selhaz
gen selhaz = invnormal(1-u)
capture drop y_hat
capture drop y_corr_ciup y_corr_cidown
predictnl y_hat_corr = normal(selhaz*(expm1(2*[/]athrho))/(exp(2*[/]athrho)+1) + _b[_cons]), ci(y_corr_ciup y_corr_cidown)
local correction_str = "`correction_str' (line y_hat_corr u, lw(*1) lcolor(blue) xaxis(1 2)) (line y_corr_ciup u, lpattern(dash) lcolor(blue)) (line y_corr_cidown u, lpattern(dash) lcolor(blue))"


twoway `sample_str' `uncorrection_str' `correction_str', ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10, nogrid) yscale(range(0, 1)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Within-Notification Mean" 10 "Uncorrected Extrapolation" 13 "Corrected Extrapolation") ring(0) position(2) region(lstyle(1))) xla(0 "0" `share_atleast_3' "1" `share_atleast_2' "2" `share_atleast_1' "3" , axis(2)) xtitle("Notifications Received", axis(2)) xline(`notice_str', lp(solid))






*-------------------------------------------------------------------------------
* 2. Test using multiple reminders (FIML and/or twostep)
*-------------------------------------------------------------------------------
*For MLE, this involves a main rho identified off of the first and second notice (1st reminder), with rho+rho_k given for all parts of the likelihood for people who respond after 3 or more notices (k>=2). 

*Ordered selection Heckman Correction stacked estimator (one standard heckman correction for each reminder, with 1 set of parameters)
*Code up MLE for binary outcomes
capture program drop ml_stacked_heckman_biv_test
program ml_stacked_heckman_biv_test
    version 16
    args todo b lnf
	
	tempvar xb_main
	mleval `xb_main' = `b', eq(1)
	forv k=1/$notification_count { //try to pass in the number of notifications from outside or back it out from parameters inside
	tempvar xbs`k'
	mleval `xbs`k'' = `b', eq(`=`k'+1')
	}
	tempvar tau
	mleval `tau' = `b', eq(`=$notification_count+2') //replace first 3 with number of notifications
	forv k=3/$notification_count { //this should start at 3, comparison of group 1 (first responders) and k identifies these correlations
	tempvar tau`k'
	mleval `tau`k'' = `b', eq(`=$notification_count+`k'') //replace first 3 with number of notifications
	}

	tempname rho
	scalar `rho' = (expm1(2*`tau'))  /  (exp(2*`tau')+1)
	forv k=3/$notification_count {
	tempname rho_`k'
	scalar `rho_`k'' = (expm1(2*(`tau'+`tau`k'')))  /  (exp(2*(`tau'+`tau`k''))+1)
	}

	*Loop over K bivariate Heckman selection models (many rhos)
	quietly replace `lnf' = 0 //prime likelihood
	*First notification equation (always includes all controls)
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs1')/1)) if ${ML_y`=1+1'}==0   // Don't respond after 1st notice
	quietly replace `lnf' =  `lnf'+log(binormal(`xb_main',`xbs1',`rho')) if ${ML_y`=1+1'}==1 & $ML_y1==1 //Respond and have response=1
	quietly replace `lnf' =  `lnf'+log(binormal(-`xb_main',`xbs1',-`rho')) if ${ML_y`=1+1'}==1 & $ML_y1==0 //Respond and have response=0
	
	*Notification special because identifies main rho
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs1'+`xbs2')/1)) if ${ML_y`=2+1'}==0   // Don't respond after 2nd notice
	quietly replace `lnf' =  `lnf'+log(binormal(`xb_main',`xbs1'+`xbs2',`rho')) if ${ML_y`=2+1'}==1 & $ML_y1==1 //Respond and have response=1
	quietly replace `lnf' =  `lnf'+log(binormal(-`xb_main',`xbs1'+`xbs2',-`rho')) if ${ML_y`=2+1'}==1 & $ML_y1==0 //Respond and have response=0	
	
	*All other reminders, may or may not include controls (covariate-specific responses to reminders)
	forv k=3/$notification_count {
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs1'+`xbs`k'')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+log(binormal(`xb_main',`xbs1'+`xbs`k'',`rho_`k'')) if ${ML_y`=`k'+1'}==1 & $ML_y1==1 //Respond and have response=1
	quietly replace `lnf' =  `lnf'+log(binormal(-`xb_main',`xbs1'+`xbs`k'',-(`rho_`k''))) if ${ML_y`=`k'+1'}==1 & $ML_y1==0 //Respond and have response=0	
	}
		
end

sum survey_ordered
local notification_count = r(max)
global notification_count = r(max)
di $notification_count
*Always have (add controls you want here)
local Equations " (Selection1: survey_1_notification = `controls')"
forv k = 2/`notification_count' {
	local Equations "`Equations' (Selection`k': survey_`k'_notification = `controls')"
}
di "`Equations'"

local athrhos ""
forv k = 3/$notification_count { //this should always be 3! The difference between groups 1 and 0 can't identify anything (0 unobserved), the difference between 1 and 2 identifies the "main rho" and the difference between 1 and 3 (the 3rd difference) identifies a new rho parameter. Subsequent notifications allow comparisons between group 1 and group k, which identify rho_k.
	local athrhos "`athrhos' /athrho`k'"
}
di "`athrhos'"


ml model lf0 ml_stacked_heckman_biv_test (Binary_Outcome: ent_binary = `controls') `Equations' /athrho `athrhos',  missing vce(robust)
*ml init initial, copy
ml maximize

estimates store rho_test

margins, expression(normal(xb(Binary_Outcome)))

sum ent_binary

test [/]athrho3


*Add corrected prediction for each rho
estimates restore rho_test
capture drop selhaz
gen selhaz = invnormal(1-u)

predictnl y_hat_corr12 = normal(selhaz*(expm1(2*[/]athrho))/(exp(2*[/]athrho)+1) + _b[_cons]), ci(y_corr_ciup12 y_corr_cidown12)
local correction_str = "`correction_str' (line y_hat_corr12 u, lw(*1) lcolor(red) xaxis(1 2)) (line y_corr_ciup12 u, lpattern(dash) lcolor(red)) (line y_corr_cidown12 u, lpattern(dash) lcolor(red))"

predictnl y_hat_corr13 = normal(selhaz*(expm1(2*([/]athrho+[/]athrho3)))/(exp(2*([/]athrho+[/]athrho3))+1) + _b[_cons]), ci(y_corr_ciup13 y_corr_cidown13)
local correction_str = "`correction_str' (line y_hat_corr13 u, lw(*1) lcolor(purple) xaxis(1 2)) (line y_corr_ciup13 u, lpattern(dash) lcolor(purple)) (line y_corr_cidown13 u, lpattern(dash) lcolor(purple))"


twoway `sample_str' `uncorrection_str' `correction_str', ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10, nogrid) yscale(range(0, 1)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Within-Notification Mean" 10 "Uncorrected Extrapolation" 13 "Corrected Extrapolation" 16 "Corrected Extrapolation 1-2" 19 "Corrected Extrapolation 1-3") ring(0) position(2) region(lstyle(1))) xla(0 "0" `share_atleast_3' "1" `share_atleast_2' "2" `share_atleast_1' "3" , axis(2)) xtitle("Reminders Received", axis(2)) xline(`notice_str', lp(solid))


cd "C:\Users\charris\Dropbox\Research\2021_Survey_Reminders\Code\Simulations"
*Make alternatives for illustrative figures
*1. Dinardo et al
* notice 1 mean graph (Figure)
* notice 1-2 mean graph (Figure)
* notice 1, and 2 mean graph (Monotonicity/index function)
* notice 1, and 2 mean with extrapolation (Figure)

*2. Harris et al
* Notice 1 mean graph (Figure)
* Notice 2 only mean graph (Figure)
* Both above with extrapolation (Figure)

*Notice 1 mean
local sample_str = ""
local notice_str = ""
local notice_lab_str = ""

capture drop rem_*_*
local ko = 0
local share_atleast_4 = 0 //1 more than highest value never happens
forv k=3(-1)3 {
	local ++ko
	di `ko'
count if survey_ordered>=`k'
local N_`k' = r(N)
count
local N = r(N)
local share_atleast_`k' = `N_`k''/`N'


reg ent_binary_observed if survey_ordered>=`k'
gen rem_`k'_mean = _b[_cons] if  u<=`=`share_atleast_`k'''
*gen rem_`k'_ciup = _b[_cons]+1.96*_se[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''
*gen rem_`k'_cidown = _b[_cons]-1.96*_se[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''

local sample_str = "`sample_str' (line rem_`k'_mean u, lw(*1) lcolor(black)  xaxis(1 2))"

local notice_str = "`notice_str' `share_atleast_`k''"
local notice_lab_str = "`notice_lab_str' `share_atleast_`k'' `""`ko'"""
}

twoway `sample_str' , ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10, nogrid) yscale(range(0, 1)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Within-Notification Mean") ring(0) position(2) region(lstyle(1))) xla(0 "0" `share_atleast_3' "1", axis(2)) xtitle("Notifications Received", axis(2)) xline(`notice_str', lp(solid)) legend(on)
graph export F1_group1.png, replace

*Notice 1-2 mean
local sample_str = ""
local notice_str = ""
local notice_lab_str = ""

capture drop rem_*_*
local ko = 0
local share_atleast_4 = 0 //1 more than highest value never happens
forv k=2(-1)2 {
	local ++ko
	di `ko'
count if survey_ordered>=`k'
local N_`k' = r(N)
count
local N = r(N)
local share_atleast_`k' = `N_`k''/`N'


reg ent_binary_observed if survey_ordered>=`k'
gen rem_`k'_mean = _b[_cons] if  u<=`=`share_atleast_`k'''
*gen rem_`k'_ciup = _b[_cons]+1.96*_se[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''
*gen rem_`k'_cidown = _b[_cons]-1.96*_se[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''

local sample_str = "`sample_str' (line rem_`k'_mean u, lw(*1) lcolor(black)  xaxis(1 2))"

local notice_str = "`notice_str' `share_atleast_`k''"
local notice_lab_str = "`notice_lab_str' `share_atleast_`k'' `""`ko'"""
}

twoway `sample_str' , ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10, nogrid) yscale(range(0, 1)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Within-Notification Mean") ring(0) position(2) region(lstyle(1))) xla(0 "0" `share_atleast_2' "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(`notice_str', lp(solid)) legend(on)
graph export F2_group2.png, replace

*Notice 1 and 2 mean separately
local sample_str = ""
local notice_str = ""
local notice_lab_str = ""

capture drop rem_*_*
local ko = 0
local share_atleast_4 = 0 //1 more than highest value never happens
forv k=3(-1)2 {
	local ++ko
	di `ko'
count if survey_ordered>=`k'
local N_`k' = r(N)
count
local N = r(N)
local share_atleast_`k' = `N_`k''/`N'


reg ent_binary_observed if survey_ordered==`k'
gen rem_`k'_mean = _b[_cons] if  u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''
*gen rem_`k'_ciup = _b[_cons]+1.96*_se[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''
*gen rem_`k'_cidown = _b[_cons]-1.96*_se[_cons] if u>`=`share_atleast_`=`k'+1''' & u<=`=`share_atleast_`k'''

local sample_str = "`sample_str' (line rem_`k'_mean u, lw(*1) lcolor(black)  xaxis(1 2))"

local notice_str = "`notice_str' `share_atleast_`k''"
local notice_lab_str = "`notice_lab_str' `share_atleast_`k'' `""`ko'"""
}

twoway `sample_str', ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10, nogrid) yscale(range(0, 1)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Within-Notification Mean") ring(0) position(2) region(lstyle(1))) xla(0 "0" `share_atleast_3' "1" `share_atleast_2' "2" , axis(2)) xtitle("Notifications Received", axis(2)) xline(`notice_str', lp(solid)) legend(on)
graph export F3_both_groups.png, replace



*Add extrapolation
twoway `sample_str' (line y_hat_corr u, lw(*1) lcolor(blue) xaxis(1 2)), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10, nogrid) yscale(range(0, 1)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Within-Notification Mean" 3 "Corrected Extrapolation") ring(0) position(2) region(lstyle(1))) xla(0 "0" `share_atleast_3' "1" `share_atleast_2' "2" , axis(2)) xtitle("Notifications Received", axis(2)) xline(`notice_str', lp(solid)) legend(on)
graph export F4_extrapolation.png, replace





