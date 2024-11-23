/*
This file does Harris, Eckhardt, and Goldfarb survey selection correction method.
It models the outcome as continuous. 
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

*Graphs look cleaner with no controls, and this lets you check rho against the male-only value above
keep if Female==0 //I coded more extreme selection for men so that gaps male-female gaps will be overstated without the correction.

*local controls = "Female STEM" //if you want
local controls = "" //simpler for graphs
*-------------------------------------------------------------------------------
* 1. Raw population average 
*-------------------------------------------------------------------------------
*Compare results to:
*1. Infeasible case with no selection
reg ent_linear
sum ent_linear
local trueval = r(mean)
*2. Ignoring selection
reg ent_linear_observed
*3. Heckman correction assigning random half of surveyed pop to no-followup group
program drop _all

heckman ent_linear_observed_rem STEM, select(rem_exclusion_response = STEM  rem_exclusion) twostep

//Confirm equivalent to ordered selection correction (binary is special case of ordered)
*Stata command in progress
*hegsurvey ent_linear_observed_rem STEM, select(rem_exclusion_response = STEM rem_exclusion) twostep


*4. Now perform preferred technique
*Two step version
*Simple ordered probit (common shift for everyone)

* if you just use 1 reminder
gen survey_ordered_simple = survey_ordered
replace survey_ordered_simple=0 if survey_ordered<2 //only use first reminder

oprobit survey_ordered_simple
gen imr_simple = (normalden([/]cut1)-normalden([/]cut2))/(normal([/]cut2)-normal([/]cut1)) if survey_ordered_simple ==2
replace imr_simple = (normalden([/]cut2))/(1-normal([/]cut2)) if survey_ordered_simple ==3

*This uses the first reminder
reg ent_linear_observed imr_simple

oprobit survey_ordered
*Generate ordered probit IMR (cut1 has opposite sign of constant in probit for selection>0)
predict imr_inner, xb
replace imr_inner = 0 if imr_inner==. & survey_ordered>0

*Generalized inverse mills ratio from ordered probit
gen imr = (normalden(imr_inner-[/]cut2)-normalden(imr_inner-[/]cut1))/(normal(imr_inner-[/]cut2)-normal(imr_inner-[/]cut1)) if survey_ordered ==1
replace imr = (normalden(imr_inner-[/]cut3)-normalden(imr_inner-[/]cut2))/(normal(imr_inner-[/]cut3)-normal(imr_inner-[/]cut2)) if survey_ordered ==2
replace imr = (normalden(imr_inner-[/]cut3))/(normal(imr_inner-[/]cut3)) if survey_ordered ==3

*This uses all the reminders, more efficient
/*
*Stata command in progress
program drop _all
hegsurvey ent_linear_observed, select(survey_ordered =) twostep
hegsurvey ent_linear_observed_rem, select(rem_exclusion_response = rem_exclusion) twostep mills(gimr)
*/

reg ent_linear_observed imr 

local twostep = _b[_cons]

gen u_hat = imr*_b[imr]


*reg ent_linear_observed imr c.imr#i.survey_ordered

reg ent_linear

di `twostep'

*Infeasible Heckman correction controlling for people's actual preference for survey
sort u_survey
gen p_true = _n/_N

gen imr_infeasible = (normalden(invnormal(p_true-1/_N))-normalden(invnormal(p_true)))/(1/_N)

reg ent_linear_observed imr_infeasible

*Code up the above in MLE for correct SEs and generalizeability (arbitrary number of reminders, allow for control variables)

*Code up linear 2nd stage - MLE
*Ordered selection Heckman Correction stacked estimator (one standard heckman correction for each reminder, with 1 set of parameters)
capture program drop ml_stacked_ordered_heckman
program ml_stacked_ordered_heckman
    version 16
    args todo b lnf
	

	
	tempvar xb_main
	mleval `xb_main' = `b', eq(1)
	forv k=1/$notification_count { //try to pass in the number of notifications from outside or back it out from parameters inside
	tempvar xbs`k'
	mleval `xbs`k'' = `b', eq(`=`k'+1')
	}
	tempvar tau lns
	mleval `tau' = `b', eq(`=$notification_count+2') //replace first 3 with number of notifications
	mleval `lns' = `b', eq(`=$notification_count+3') //replace first 3 with number of notifications
	
	tempname sig lambda
	tempvar rho
	gen double `rho' = (expm1(2*`tau'))  /  (exp(2*`tau')+1)
	scalar `sig' = exp(`lns')
	scalar `lambda' = `rho'*`sig'

	*Loop over K Heckman selection models (only 1 sigma and rho)
	quietly replace `lnf' = 0 //prime likelihood
	*First notification equation (always includes all controls)
		quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs1')/1)) if ${ML_y`=1+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+lnnormalden($ML_y1, `xb_main', `sig') + lnnormal((`xbs1' + `rho'/`sig'*($ML_y1 - `xb_main'))/sqrt(1 - `rho'^2)) if ${ML_y`=1+1'}==1
	*All other reminders, may or may not include controls (covariate-specific responses to reminders)
	forv k=2/$notification_count {
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs1'+`xbs`k'')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+lnnormalden($ML_y1, `xb_main', `sig') + lnnormal((`xbs1'+`xbs`k'' + `rho'/`sig'*($ML_y1 - `xb_main'))/sqrt(1 - `rho'^2)) if ${ML_y`=`k'+1'}==1 
	
	}	
end

/*I want an arbitrary number of reminders where individuals are allowed to be more/less responsive to reminders based on observables - XB_k for all K reminders.
1. Looping from 1 to K where K is the total number of notices (k-1 is number of reminders)
2. Adding equation info to a string for each k
3. Sticking that string into the ml model command

This should allow for an arbitrary number of selection equations, with notification-specific coefficients
 */

sum survey_ordered
global notification_count = r(max)
di $notification_count
*Always have 
local Equations "`Equations' (Selection1: survey_1_notification = `controls')"
forv k = 2/$notification_count {
	local Equations "`Equations' (Selection`k': survey_`k'_notification = `controls')"
}
di "`Equations'"

local athrho_eq "(athrho_XB:  = )" //set up an entire linear equation for rho - current is constant only, but you can input a treatment variable for select/outcome correlation to vary e.g. by gender

ml model lf0 ml_stacked_ordered_heckman (Linear_s2: ent_linear_observed = `controls') `Equations' `athrho_eq' /lnsigma,  missing vce(robust)
*ml init initial, copy
ml maximize

test [athrho_XB]_cons

predict Y_hat_stacked_mle 


local rho = (expm1(2*[athrho_XB]_cons))  /  (exp(2*[athrho_XB]_cons)+1)
local sigma = exp(/lns)
di `rho'
di `sigma'

local cons_mle = _b[_cons]


test [Linear_s2]_cons = `trueval'

reg ent_linear_observed imr 
test _cons = `trueval'

di `twostep'
di `cons_mle'
di `trueval'

*Make figures

local N_graphs = 1000 //set number of data points for illustrative graphs, doesn't need to equal N
capture set obs `N_graphs'

*Make graph for actual observed data, showing avg outcome on Y axis, x is percentile of preference for answering survey
sum ent_linear_observed if survey_ordered==3
local rem_0_mean = r(mean)
sum ent_linear_observed if survey_ordered==2
local rem_1_mean = r(mean)
sum ent_linear_observed if survey_ordered==1
local rem_2_mean = r(mean)

sum survey_1_notification
local share_rem_0 = r(mean)
sum survey_2_notification
local share_rem_1 = r(mean)
sum survey_3_notification
local share_rem_2 = r(mean)

/*right to left
di `=`share_rem_0'+`share_rem_1'+`share_rem_2''
capture drop Y_bar_data
gen Y_bar_data = .
replace Y_bar_data = 0 if _n<=100 & _n<`=100-`share_rem_2'*100'
replace Y_bar_data = `rem_0_mean' if _n>=`=100-`share_rem_0'*100'
replace Y_bar_data = `rem_1_mean' if _n>=`=100-`share_rem_1'*100' & _n<`=100-`share_rem_0'*100'
replace Y_bar_data = `rem_2_mean' if _n>=`=100-`share_rem_2'*100' & _n<`=100-`share_rem_1'*100'
*/

*left to right
di `=`share_rem_0'+`share_rem_1'+`share_rem_2''
capture drop Y_bar_data
gen Y_bar_data = .
*replace Y_bar_data = 0 if _n<=100 & _n>`=`share_rem_2'*100'
replace Y_bar_data = `rem_0_mean' if _n<=`=`share_rem_0'*`N_graphs''
replace Y_bar_data = `rem_1_mean' if _n<=`=`share_rem_1'*`N_graphs'' & _n>`=`share_rem_0'*`N_graphs''
replace Y_bar_data = `rem_2_mean' if _n<=`=`share_rem_2'*`N_graphs'' & _n>`=`share_rem_1'*`N_graphs''

*Same thing for percentile rank derived from imr
sum imr if survey_ordered==3
local rem_0_imr = r(mean)
sum imr if survey_ordered==2
local rem_1_imr = r(mean)
sum imr if survey_ordered==1
local rem_2_imr = r(mean)

capture drop p_imr_data
gen p_imr_data = .
*replace Y_bar_data = 0 if _n<=100 & _n>`=`share_rem_2'*100'
replace p_imr_data = normal(`=-1*`rem_0_imr'')*100 if _n<=`=`share_rem_0'*`N_graphs''
replace p_imr_data = normal(`=-1*`rem_1_imr'')*100 if _n<=`=`share_rem_1'*`N_graphs'' & _n>`=`share_rem_0'*`N_graphs''
replace p_imr_data = normal(`=-1*`rem_2_imr'')*100 if _n<=`=`share_rem_2'*`N_graphs'' & _n>`=`share_rem_1'*`N_graphs''


capture drop p_data
gen p_data = .
replace p_data = `share_rem_0'/2*100 if _n<=`=`share_rem_0'*`N_graphs''
replace p_data = (`share_rem_0'+(`share_rem_1'-`share_rem_0')/2)*100 if _n<=`=`share_rem_1'*`N_graphs'' & _n>`=`share_rem_0'*`N_graphs''
replace p_data = (`share_rem_1'+(`share_rem_2'-`share_rem_1')/2)*100 if _n<=`=`share_rem_2'*`N_graphs'' & _n>`=`share_rem_1'*`N_graphs''





capture drop p
gen p = _n*100/`N_graphs' if _n<=`N_graphs'



*Generate guess for unobserved preference for surveys (this code only works correctly with no controls, use for illustrations!)
capture drop u 
gen u = .
replace u = runiform(0,`share_rem_0')*100 if survey_ordered==3
replace u = runiform(`share_rem_0',`share_rem_1')*100 if survey_ordered==2
replace u = runiform(`share_rem_1',`share_rem_2')*100 if survey_ordered==1
replace u = runiform(`share_rem_2',1)*100 if survey_ordered==0
reg ent_linear_observed u
local b_up = _b[u]
local Y_cons = _b[_cons]


di `share_rem_0'
capture drop p_mean
gen p_mean = .
replace p_mean = `share_rem_0'/2*100 if survey_ordered==3
sum p_mean if survey_ordered==3
replace p_mean = (`share_rem_0'+(`share_rem_1'-`share_rem_0')/2)*100 if survey_ordered==2
replace p_mean = (`share_rem_1'+(`share_rem_2'-`share_rem_1')/2)*100 if survey_ordered==1
replace p_mean = (`share_rem_2'+(1-`share_rem_2')/2)*100 if survey_ordered==0




*twoway (line Y_bar_data p if p<`=`share_rem_0'*100') ///
*(line Y_bar_data p_imr if p<=`=`share_rem_2'*100' & p>`=`share_rem_1'*100') ///
*(line Y_bar_data p if p<=`=`share_rem_1'*100' & p>`=`share_rem_0'*100') ///
*(line Y_bar_data p if p<=`=`share_rem_0'*100') ///
*(line vert1 p_max_rem_0)

*gen imrp = invnormal(p*0.01)

*Within Sample mean
reg ent_linear_observed //for SEs
capture drop Y_bar_mean
gen Y_bar_mean = _b[_cons] if p<.

*Normal (Heckman) Correction extrapolation 
capture drop imr_p
gen imr_p = invnormal(1-p/100) //the inverse mills ratio is just a guess at a person's t-stat - it is more complicated usually because we DON'T know people's unobserved preferences - for the graph we do.

capture drop Y_norm_corr
gen Y_norm_corr = imr_p*`rho'*`sigma'+`cons_mle'

*Generate linear extrapolation
reg ent_linear_observed p_mean
local b_p_mean = _b[p_mean]
local Y_cons = _b[_cons]
capture drop Y_lin_corr
gen Y_lin_corr = `b_p_mean'*p+`Y_cons'

*Get y_range from max/min of objects in all graphs
qui sum Y_norm_corr
local y_min_norm = r(min)
local y_max_norm = r(max)
qui sum Y_lin_corr
local y_min_lin = r(min)
local y_max_lin = r(max)

local y_min = min(`y_min_norm', `y_min_lin')
local y_max = max(`y_max_norm', `y_max_lin')



*Show data for average survey preference percentile with average outcome for each bin.
twoway (scatter Y_bar_data p_data if Y_bar_data<.), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,100)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile"))

*Add Within-sample mean
twoway (scatter Y_bar_data p_data if Y_bar_data<.) (line Y_bar_mean p), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,100)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome"))

*Add linear extrapolation
twoway (scatter Y_bar_data p_data if Y_bar_data<.) (line Y_bar_mean p) (line Y_lin_corr p), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,100)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome" 3 "Linear Extrapolation"))

*Add percentile rank of average imr for individuals in each bin
twoway (scatter Y_bar_data p_data if Y_bar_data<.) (line Y_bar_mean p) (line Y_lin_corr p) (scatter Y_bar_data p_imr_data) , ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,100)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome" 3 "Linear Extrapolation" 4 "Percentile, Mean(IMR)"))

*Add correction under normality assumption
twoway (scatter Y_bar_data p_data if Y_bar_data<.) (line Y_bar_mean p) (line Y_lin_corr p) (scatter Y_bar_data p_imr_data) (line Y_norm_corr p), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,100)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome" 3 "Linear Extrapolation" 4 "Percentile, Mean(IMR)" 5 "Normal Extrapolation"))



*Have only sample mean and correction
twoway (scatter Y_bar_data p_data if Y_bar_data<.) (line Y_bar_mean p) (line Y_norm_corr p), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,100)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome" 3 "Normal Extrapolation"))

*Have only sample mean and both corrections
twoway (scatter Y_bar_data p_data if Y_bar_data<.) (line Y_bar_mean p) (line Y_norm_corr p) (line Y_lin_corr p), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,100)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome" 3 "Normal Extrapolation" 4 "Linear Extrapolation"))


*-------------------------------------------------------------------------------
* 2. Test using multiple reminders (FIML)
*-------------------------------------------------------------------------------
*We jointly estimate the model as above with another model that includes fixed effects for all notification groups other than the first. Jointly testing that all FE's are 0 is a test of whether the model fits the data well, very similar to a pretrends test.

*Ordered selection Heckman Correction stacked estimator (one standard heckman correction for each reminder, with 1 set of parameters)
capture program drop ml_stacked_ordered_heckman_test
program ml_stacked_ordered_heckman_test
    version 16
    args todo b lnf
	
	
	
	tempvar xb_main
	mleval `xb_main' = `b', eq(1)
	forv k=1/$notification_count { //try to pass in the number of notifications from outside or back it out from parameters inside
	tempvar xbs`k'
	mleval `xbs`k'' = `b', eq(`=`k'+1')
	}
	
	tempvar tau lns
	mleval `tau' = `b', eq(`=$notification_count+2') //replace first 3 with number of notifications
	mleval `lns' = `b', eq(`=$notification_count+3')
	

	tempname sig lambda
	tempvar rho
	gen double `rho' = (expm1(2*`tau'))  /  (exp(2*`tau')+1)
	scalar `sig' = exp(`lns')
	scalar `lambda' = `rho'*`sig'

	*Overidentified model
	forv k=2/$notification_count { //notifcation fixed effects to rationalize average outcome taking rho as given, 0 if model is perfect
	tempvar FE`k'
	mleval `FE`k'' = `b', eq(`=$notification_count+3+`k'-1')
	}
	
	tempvar xb_main_oid
	mleval `xb_main_oid' = `b', eq(`=$notification_count+3+$notification_count-1+1')
	forv k=1/$notification_count { //try to pass in the number of notifications from outside or back it out from parameters inside
	tempvar xbs`k'_oid
	mleval `xbs`k'_oid' = `b', eq(`=$notification_count+3+$notification_count-1+1+`k'')
	}

	tempvar lns_oid

	mleval `lns_oid' = `b', eq(`=$notification_count+3+$notification_count-1+1+$notification_count+1')
	
	tempname sig_oid lambda_oid
	scalar `sig_oid' = exp(`lns_oid')
	scalar `lambda_oid' = `rho'*`sig_oid'

	
	
	
*-------Estimate base model
	*Loop over K Heckman selection models (only 1 sigma and rho)
	quietly replace `lnf' = 0 //prime likelihood

	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs1')/1)) if ${ML_y`=1+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+lnnormalden($ML_y1, `xb_main', `sig') + lnnormal((`xbs1' + `rho'/`sig'*($ML_y1 - `xb_main'))/sqrt(1 - `rho'^2)) if ${ML_y`=1+1'}==1 

	forv k=1/$notification_count { //these identify deviations from main rho
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs`k'')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+lnnormalden($ML_y1, `xb_main', `sig') + lnnormal((`xbs`k'' + `rho'/`sig'*($ML_y1 - `xb_main'))/sqrt(1 - `rho'^2)) if ${ML_y`=`k'+1'}==1 
	}
	
	
*------Jointly estimate model with FEs
	*Loop over K Heckman selection models (only 1 sigma and rho)
	forv k=1/1 { //these identify main rho
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs`k'_oid')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+lnnormalden($ML_y1, `xb_main_oid', `sig_oid') + lnnormal((`xbs`k'_oid' + `rho'/`sig_oid'*($ML_y1 - `xb_main_oid'))/sqrt(1 - `rho'^2)) if ${ML_y`=`k'+1'}==1 
	}	
	forv k=2/$notification_count { //these identify deviations from main rho
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs`k'_oid')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+lnnormalden($ML_y1, `xb_main_oid'+`FE`k'', `sig_oid') + lnnormal((`xbs`k'_oid' + `rho'/`sig_oid'*($ML_y1 - (`xb_main_oid'+`FE`k'')))/sqrt(1 - `rho'^2)) if ${ML_y`=`k'+1'}==1 
	}
	
end

/*I want an arbitrary number of reminders where individuals are allowed to be more/less responsive to reminders based on observables - XB_k for all K reminders. 
1. Looping from 1 to K where K is the total number of notices (k-1 is number of reminders)
2. Adding equation info to a string for each k
3. Sticking that string into the ml model command

This should allow for an arbitrary number of selection equations, with notification-specific coefficients
 */

sum survey_ordered
local notification_count = r(max)
global notification_count = r(max)

*Always have (add controls you want here)
local Equations " (Selection1: survey_1_notification = `controls')"
local Equations_oid " (Selection1_oid: survey_1_notification = `controls')"
forv k = 2/`notification_count' {
	local Equations "`Equations' (Selection`k': survey_`k'_notification = `controls')"
	local Equations_oid "`Equations_oid' (Selection`k'_oid: survey_`k'_notification = `controls')"
}
di "`Equations'"

local FE_notification ""
local FE_test
forv k = 2/$notification_count { 
	local FE_notification "`FE_notification' /FE_notification`k'"
	local FE_test "`FE_test' [/]FE_notification`k'"
}
di "`FE_notification'"



local athrho_eq "(athrho_XB:  = )" //set up an entire linear equation for rho - current is constant only, but you can input a treatment variable for select/outcome correlation to vary e.g. by gender


ml model lf0 ml_stacked_ordered_heckman_test (Linear_s2: ent_linear_observed = `controls') `Equations' `athrho_eq' /lnsigma `FE_notification' (Binary_Outcome_oid: ent_binary = `controls') `Equations_oid' /lnsigma_oid, vce(robust) missing
*ml check
*ml init initial, copy
ml maximize

estimate store overid

test `FE_test'






