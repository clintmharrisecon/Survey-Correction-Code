/*
*Binary outcome
This file demonstrates the Harris, Eckhardt, and Goldfarb survey selection correction method. It creates simulated data on gender entrepreneurship gaps, with a STEM binary indicator as a "control".
This is the same approach we would use for treatment effects by replacing the Female indicator with a (randomly assigned) treatment indicator.

1. DGP generates
	a. A continuous and binary outcome in a DGP
	b. Selection in nonresponse
	c. Observed response timing
	
2. Estimates nonresponse-corrected MLE
	a. Allows for subgroup-specific nonresponse bias parameters (rho)
	b. Test response_preference-outcome correlation within groups, that each is zero, that they are jointly zero, and that difference between groups is zero. 
	
3. Estimates overidentification test
	a. Overidentification test using "predictable trends" test.
	b. Produces "predictable trends" overidentification figure.
	
*/

clear
cls
version 18.0
*-------------------------------------------------------------------------------
* 0. Choose controls/etc
*-------------------------------------------------------------------------------
set seed 12345
*Outcome
*local Y = "ent_linear_observed" //
local Y = "ent_binary_observed"
*Covariates in main equation
local X = "Female i.Female#c.STEM"
*Covariates in survey respone equation (fine if the same as X)
local Z = "Female i.Female#c.STEM"
*Covariates to model heterogeneity in rho (subgroup selection bias)
local X_rho = "Female i.Female#c.STEM"
*Covariates to model heterogeneity in effects of later notifications (overidentification test variables)
local X_oid = "Female i.Female#c.STEM"

/*
local X = "Female"
*Covariates in survey respone equation (fine if the same as X)
local Z = "Female"
*Covariates to model heterogeneity in rho (subgroup selection bias)
local X_rho = "Female"
*Covariates to model heterogeneity in effects of later notifications (overidentification test variables)
local X_oid = "Female"
*/

*Choose whether you want an illustrative graph (one group, no controls)
local graphs = 0
if `graphs' == 1 {
	local X = ""
	local Z = ""
	local X_rho = ""
}


*Group indicator - used to test whether rho is equal between groups and whether overidentification test passes within each group
global S_ML_groupvar = "Female" //define for group-specific overidentification test code

*-------------------------------------------------------------------------------
* 1. DGP
*-------------------------------------------------------------------------------

*Setup
local gender_pop = 2000
clear
cls
set obs `gender_pop'
gen Female=1


*Draw unobserved entrepreneurship and willingness to fill out related survey from joint distribution for women
mat Sig = (1, -.3 \ -.3, 1)
drawnorm u_ent u_survey, cov(Sig)

*Do the same for men
set obs `=`gender_pop'*2'
replace Female=0 if Female==.
mat Sig = (1, .3 \ .3, 1)
drawnorm m_u_ent m_u_survey, cov(Sig)

gen STEM = rnormal()-.5*Female>0
*Note: STEM is some predictor of both survey response and entrepreneurship, e.g. STEM major

replace u_ent = m_u_ent if Female==0
replace u_survey = m_u_survey if Female==0 

drop m_*

*Create latent variables
replace u_ent = u_ent-0.5-0.5*Female+STEM //Women less entrepreneurial than men
replace u_survey = u_survey-1+STEM //most people don't respond to survey


*Generate linear ent measure
rename u_ent ent_linear
*Generate binary ent measure
gen ent_binary = ent_linear>0
replace ent_linear = ent_linear+5

*Generate 2 reminders, sufficient to perform method and overid test.
gen survey_1_response = u_survey>0
gen survey_2_response = u_survey+0.5>0
gen survey_3_response = u_survey+.75>0
gen survey_4_response = u_survey+1>0
gen survey_ordered = survey_1_response+survey_2_response+survey_3_response+survey_4_response




*Generate responses to everyone getting survey reminders version
gen ent_linear_observed = ent_linear if survey_ordered>0
gen ent_binary_observed = ent_binary if survey_ordered>0


*Randomize intensity of follow-up, following DiNardo et al 2021. We approximate this ex post, by throwing away responses received after reminders for random half of subjects.
*1. Generate binary "intensity" exclusion restriction.
gen rem_exclusion = rnormal()>0
*2. Throw away late responses from those not assigned to "intense follow-up"
gen rem_exclusion_response = survey_4_response if rem_exclusion==1
replace rem_exclusion_response = survey_1_response if rem_exclusion==0

gen ent_linear_observed_rem = ent_linear_observed if rem_exclusion_response==1 //keep responses received by group-specific cutoff time
gen ent_binary_observed_rem = ent_binary_observed if rem_exclusion_response==1 //keep responses received by group-specific cutoff time


if `graphs' == 1 {
	keep if Female==0
}

*-------------------------------------------------------------------------------
* 2. Estimate corrected model
*-------------------------------------------------------------------------------
*Compare results to:
*1. Infeasible case with no selection
probit ent_binary `X', r
*2. Ignoring selection
probit `Y' `X', r
*3. Feasible, but inefficient: Heckman correction assigning random half of surveyed pop to no-followup group - bysort is to give it gender-specific selection corrections
bysort $S_ML_groupvar: heckprobit ent_binary_observed_rem `X', select(rem_exclusion_response = `Z' rem_exclusion) r

*------------------------------------
* 2.a Estimate model
*------------------------------------
*Code up MLE for binary outcomes
capture program drop ml_hegsurvey_binary
program ml_hegsurvey_binary
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

	tempvar rho
	gen double `rho' = (expm1(2*`tau'))  /  (exp(2*`tau')+1)

	*Loop over K bivariate Heckman selection models (only 1 rho)
	quietly replace `lnf' = 0 //prime likelihood
		forv k=1/$notification_count {
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs`k'')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
	quietly replace `lnf' =  `lnf'+log(binormal(`xb_main',`xbs`k'',`rho')) if ${ML_y`=`k'+1'}==1 & $ML_y1==1 //Respond and have response=1
	quietly replace `lnf' =  `lnf'+log(binormal(-(`xb_main'),`xbs`k'',-(`rho'))) if ${ML_y`=`k'+1'}==1 & $ML_y1==0 //Respond and have response=0	
	}
		
end







sum survey_ordered
global notification_count = r(max)
*Always have (add controls you want here)
local Equations "`Equations' (Selection1: survey_1_response = `Z')"
forv k = 2/$notification_count {
	local Equations "`Equations' (Selection`k': survey_`k'_response = `Z')"
}
di "`Equations'"

ml model lf0 ml_hegsurvey_binary (Binary_Outcome: `Y' = `X') `Equations' (athrho_XB: = `X_rho'),  missing vce(robust)
*ml init initial, copy
ml maximize

estimates store corrected_estimates
mat initial = e(b)

*------------------------------------
* 2.b Test rho 
*------------------------------------
*Test rho = 0 and rho 
*This errors out with no subgroups
if `graphs' == 0 {
estimates restore corrected_estimates
*Test men and women
margins , expression((expm1(2*(xb(athrho_XB))))  /  (exp(2*(xb(athrho_XB)))+1)) over($S_ML_groupvar) pwcompare post //post is required for joint test
*Men and women's rho, SEs and p-values
matlist r(table)
*Gender difference in rho
matlist r(table_vs)
*Women and men rhos are jointly 0
test (_b[0bn.$S_ML_groupvar]=0) (_b[1.$S_ML_groupvar]=0)
di r(p)
}
estimates restore corrected_estimates



*Produce figure (only nice with no controls)
if `graphs' == 1 {
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
estimates restore corrected_estimates
capture drop selhaz
gen selhaz = invnormal(1-u)
capture drop y_hat
capture drop y_corr_ciup y_corr_cidown
predictnl y_hat_corr = normal(selhaz*(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1) + _b[_cons]), ci(y_corr_ciup y_corr_cidown)
local correction_str = "`correction_str' (line y_hat_corr u, lw(*1) lcolor(blue) xaxis(1 2)) (line y_corr_ciup u, lpattern(dash) lcolor(blue)) (line y_corr_cidown u, lpattern(dash) lcolor(blue))"


twoway `sample_str' `uncorrection_str' `correction_str', ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10, nogrid) yscale(range(0, 1)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Within-Notification Mean" 10 "Uncorrected Extrapolation" 13 "Corrected Extrapolation") ring(0) position(2) region(lstyle(1))) xla(0 "0" `share_atleast_3' "1" `share_atleast_2' "2" `share_atleast_1' "3" , axis(2)) xtitle("Notifications Received", axis(2)) xline(`notice_str', lp(solid))
}




*-------------------------------------------------------------------------------
* 3. Overidentification test using multiple reminders
*-------------------------------------------------------------------------------
capture program drop ml_hegsurvey_bin_test
program ml_hegsurvey_bin_test
    version 16
    args todo b lnf

	tempvar xb_main
	mleval `xb_main' = `b', eq(1)
	forv k=1/$notification_count { //pass in the number of notifications from outside
	tempvar xbs`k'
	mleval `xbs`k'' = `b', eq(`=`k'+1')
	}
	
	tempvar tau
	mleval `tau' = `b', eq(`=$notification_count+2') //replace first 3 with number of notifications

	tempvar rho
	gen double `rho' = (expm1(2*`tau'))  /  (exp(2*`tau')+1)

	*Overidentified model
	forv k = 3/$notification_count {
	tempvar oid_xb`k'
	mleval `oid_xb`k'' = `b', eq(`=$notification_count+`k'')	
	}
	
	
*-------Estimate base model
	*Loop over K Heckman selection models (only 1 rho)
	quietly replace `lnf' = 0 //prime likelihood
			   
	forv k=1/2 {
               quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs`k'')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
               quietly replace `lnf' =  `lnf'+log(binormal(`xb_main',`xbs`k'',`rho')) if ${ML_y`=`k'+1'}==1 & $ML_y1==1 //Respond and have response=1
               quietly replace `lnf' =  `lnf'+log(binormal(-(`xb_main'),`xbs`k'',-(`rho'))) if ${ML_y`=`k'+1'}==1 & $ML_y1==0 //Respond and have response=0            
               }
			   
	forv k=3/$notification_count {
               quietly replace `lnf' =  `lnf'+log(1 - normal((`xbs`k'')/1)) if ${ML_y`=`k'+1'}==0   // Don't respond after kth notice
               quietly replace `lnf' =  `lnf'+log(binormal(`xb_main'+`oid_xb`k'',`xbs`k'',`rho')) if ${ML_y`=`k'+1'}==1 & $ML_y1==1 //Respond and have response=1
               quietly replace `lnf' =  `lnf'+log(binormal(-(`xb_main'+`oid_xb`k''),`xbs`k'',-`rho')) if ${ML_y`=`k'+1'}==1 & $ML_y1==0 //Respond and have response=0            
               }

end









sum survey_ordered
global notification_count = r(max)

*Always have (add controls you want here)
local Equations " (Selection1: survey_1_response = `Z')"
forv k = 2/$notification_count {
	local Equations "`Equations' (Selection`k': survey_`k'_response = `Z')"
}
di "`Equations'"

*Add equation for each overidentified notification
local overid ""
forv k = 3/$notification_count {
	local overid "`overid' (overid`k': `Y' = `X_oid')"
}


gen r = rnormal()
qui reg r `X_oid' //get the right number of zeros
mat initial_oid = initial
forv k = 3/$notification_count {
mat initial_oid = initial_oid, e(b)
}


ml model lf0 ml_hegsurvey_bin_test (Binary_Outcome: `Y' = `X') `Equations' (athrho_XB: = `X_rho') `overid', vce(robust) missing
*ml check
ml init initial_oid, copy skip
ml maximize

estimates store overid


/*
*Confirm equivalent estimates to dropping all late responders (true if you fully-flexibly specify heterogeneous effects of notifications across X)
preserve
*Drop last reminder
replace `Y' = . if survey_ordered<$notification_count-1
forv k = 3/$notification_count {
	replace survey_`k'_response = survey_2_response
}
ml model lf0 ml_hegsurvey_binary (Binary_Outcome: `Y' = `X') `Equations' (athrho_XB: = `X_rho'),  missing vce(robust)
ml init initial, copy
ml maximize
restore
*/


*------------------------------------
* 3.a Test overidentification (null of zero effects in terms of probit coefficients - eta in the paper) 
*------------------------------------
local over_id_predict_str  "predict(xb equation(overid3))"
local test_0 "(_b[1bn._predict#0bn.Female]=0)"
local test_1 "(_b[1bn._predict#1.Female]=0)"
forv k = 4/$notification_count {
	local over_id_predict_str  "`over_id_predict_str' predict(xb equation(overid`k'))"
	local test_0 "`test_0' (_b[`=`k'-2'._predict#0bn.Female]=0)"
	local test_1 "`test_1' (_b[`=`k'-2'._predict#1.Female]=0)"
}
estimates restore overid
margins , `over_id_predict_str' post over($S_ML_groupvar)


*Test men 
test `test_0'
*Test women
test `test_1'
*Test both jointly
test `test_0' `test_1'


*------------------------------------
* 3.b Predictable trends figure (in effects by notification group for interpretability) 
*------------------------------------
*Check Y=Y_hat for each response group, for both men and women
*Create indicators for survey response that are binary and mutually exclusive
local timing = 0
forv k = $notification_count(-1)1 {
	local ++timing
	capture gen response_group_`timing' = survey_ordered==`k'
}
capture gen response_group = $notification_count-survey_ordered+1

local Y_hat_str = "(response_group_1==1)*binormal(xb(Binary_Outcome),xb(Selection1),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))/normal(xb(Selection1))"
forv k = 2/`=$notification_count' {
	local Y_hat_str = "`Y_hat_str' + (response_group_`k'==1)* 	(binormal(xb(Binary_Outcome),xb(Selection`k'),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))-binormal(xb(Binary_Outcome),xb(Selection`=`k'-1'),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))) / 	(normal(xb(Selection`k'))-normal(xb(Selection`=`k'-1')))"
}

local Y_hat_oid_str = "(response_group_1==1)*binormal(xb(Binary_Outcome),xb(Selection1),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))/normal(xb(Selection1))"
forv k = 2/2 { //2 has no overidentification shifts
	local Y_hat_oid_str = "`Y_hat_oid_str' + (response_group_`k'==1)* 	(binormal(xb(Binary_Outcome),xb(Selection`k'),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))-binormal(xb(Binary_Outcome),xb(Selection`=`k'-1'),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))) / 	(normal(xb(Selection`k'))-normal(xb(Selection`=`k'-1')))"
}
forv k = 3/3 { //3 has 1 overidentification shift
	local Y_hat_oid_str = "`Y_hat_oid_str' + (response_group_`k'==1)* 	(binormal(xb(Binary_Outcome)+xb(overid`k'),xb(Selection`k'),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1)) -binormal(xb(Binary_Outcome),xb(Selection`=`k'-1'),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))) / 	(normal(xb(Selection`k'))-normal(xb(Selection`=`k'-1')))"
}
forv k = 4/`=$notification_count' { //4+ have the current and prior overid shift
	local Y_hat_oid_str = "`Y_hat_oid_str' + (response_group_`k'==1)* 	(binormal(xb(Binary_Outcome)+xb(overid`k'),xb(Selection`k'),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1)) -binormal(xb(Binary_Outcome)+xb(overid`=`k'-1'),xb(Selection`=`k'-1'),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))) / 	(normal(xb(Selection`k'))-normal(xb(Selection`=`k'-1')))"
}

*Y_hat with and without overid parameters
/*
estimates restore overid
margins , expression(`Y_hat_str') over(response_group $S_ML_groupvar)
margins , expression(`Y_hat_oid_str') over(response_group $S_ML_groupvar)
bysort response_group Female : sum ent_binary_observed
*/

*Overidentification test with effects of notifications on outcome in terms of observed outcome variable
estimates restore overid
margins if response_group>=3 & response_group<=$notification_count, expression(`Y_hat_oid_str'-(`Y_hat_str')) over(response_group $S_ML_groupvar) post

mat r_table = r(table)

estimates store event_study

*Create test and terms for event study figure
local over_id_predict_str  "predict(xb equation(overid3))"
local test_0 "(_b[3bn.response_group#0bn.Female]=0)"
local test_1 "(_b[3bn.response_group#1.Female]=0)"
forv k = 4/$notification_count {
	local over_id_predict_str  "`over_id_predict_str' predict(xb equation(overid`k'))"
	local test_0 "`test_0' (_b[`k'.response_group#0bn.Female]=0)"
	local test_1 "`test_1' (_b[`k'.response_group#1.Female]=0)"
}

*Test men 
test `test_0'
*Test women
test `test_1'
*Test both jointly
test `test_0' `test_1'



estimates restore event_study
*Produce figure
matrix coef_B0 = J($notification_count,3,0)
matrix coef_B1 = J($notification_count,3,0)

local es_label_string = ""
forv k = 3/$notification_count {
	*get B, lb, and ub for 95% confidence
	*Group 0
	matrix coef_B0[`k',1] = r_table[1,(`k'-2)*2-1], r_table[5,(`k'-2)*2-1], r_table[6,(`k'-2)*2-1]
	*Group 1
	matrix coef_B1[`k',1] = r_table[1,(`k'-2)*2], r_table[5,(`k'-2)*2], r_table[6,(`k'-2)*2]
	
	local es_label_string "`es_label_string' `=`k'' "`k'""
}

coefplot (matrix(coef_B0[,1]), ci((coef_B0[,2] coef_B0[,3]))) ///
(matrix(coef_B1[,1]), ci((coef_B1[,2] coef_B1[,3]))) ///
, vertical ytitle("Effects of Notifications on Outcome") xtitle("Notification") yline(0) ciopts(recast(rcaps)) xlab(1 "1" 2 "2" `es_label_string') legend(order(1 "Group 0" 3 "Group 1"))



