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
* 0.a. Set options for which things to see
*-------------------------------------------------------------------------------
*Choose whether you want an illustrative graph (one group, plots MSR function). The MSR function is X-specific, so this uses no controls to make things clean.
local graphs = 0

*Choose whether you want to see differences between model Y_hat and E[Y|S(r)-S(r-1)=1] for the overidentified model. This increases computation time minimally if graphs==1 (no controls), but substantially if graphs==0 (takes about 1 minute with N=200, NT=1000)
local Y_oid = 0

*-------------------------------------------------------------------------------
* 0.b. Choose controls/etc
*-------------------------------------------------------------------------------
set seed 12345
*Outcome
*local Y = "ent_continuous_observed" //
local Y = "ent_binary_observed"
*Covariates in main equation
local X = "c.Female i.Female#c.STEM"
*Extra covariates in survey response equation (reminder indicators)
local Z = "i.requests i.requests#c.Female i.requests#i.Female#c.STEM"
*Covariates to model heterogeneity in rho (subgroup selection bias)
local X_rho = "c.Female i.Female#c.STEM"
*Covariates to model heterogeneity in effects of later requests (overidentification interactions with request indicators) - leave empty if you only want requests 3+ indicators, to see if Y at those values is predicted by model estimated off of first 2 requests on average
local X_oid = "Female STEM"
*local X_oid = ""

/*
local X = "Female"
*Covariates in survey respone equation (fine if the same as X)
local Z = "Female"
*Covariates to model heterogeneity in rho (subgroup selection bias)
local X_rho = "Female"
*Covariates to model heterogeneity in effects of later Requests (overidentification test variables)
local X_oid = "Female"
*/


if `graphs' == 1 {
	local X = ""
	local Z = "i.requests"
	local X_rho = ""
	local X_oid = ""
}


*Group indicator - used to test whether rho is equal between groups and whether overidentification test passes within each group
global S_ML_groupvar = "Female" //define for group-specific overidentification test code

*-------------------------------------------------------------------------------
* 1. DGP
*-------------------------------------------------------------------------------

*Setup
local gender_pop = 500
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


*Generate continuous ent measure
rename u_ent ent_continuous
*Generate binary ent measure
gen ent_binary = ent_continuous>0
replace ent_continuous = ent_continuous

*Generate 4 requests, sufficient to perform method and overid test.
gen survey_response1 = u_survey>0
gen survey_response2 = u_survey+0.5>0
gen survey_response3 = u_survey+.75>0
gen survey_response4 = u_survey+1>0
egen survey_ordered = rowtotal(survey_response*)




*Generate responses to everyone getting survey requests version
gen ent_continuous_observed = ent_continuous if survey_ordered>0
gen ent_binary_observed = ent_binary if survey_ordered>0


*Randomize intensity of follow-up, following DiNardo et al 2021. We approximate this ex post, by throwing away responses received after requests for random half of subjects.
*1. Generate binary "intensity" exclusion restriction.
gen rem_exclusion = rnormal()>0
*2. Throw away late responses from those not assigned to "intense follow-up"
gen rem_exclusion_response = survey_response4 if rem_exclusion==1
replace rem_exclusion_response = survey_response1 if rem_exclusion==0

gen ent_continuous_observed_rem = ent_continuous_observed if rem_exclusion_response==1 //keep responses received by group-specific cutoff time
gen ent_binary_observed_rem = ent_binary_observed if rem_exclusion_response==1 //keep responses received by group-specific cutoff time


if `graphs' == 1 {
	keep if Female==0
}

gen id = _n

reshape long survey_response , i(id) j(requests)
drop if survey_response == . //if individuals did not receive a kth reminder, drop them.

replace ent_binary_observed = ent_binary_observed*survey_response
replace ent_binary_observed = 0 if ent_binary_observed == . //in accordance with LAR proof, be careful to explicitly code the Heckman correction with (survey_response = Z) rather than letting it default to when Y is missing.
*-------------------------------------------------------------------------------
* 2. Estimate corrected model
*-------------------------------------------------------------------------------
*Compare results to those that discard panel variation
preserve 
bysort id: egen T = max(requests)
keep if (requests == T & rem_exclusion==1) | (requests == 1 & rem_exclusion==0)
*1. Infeasible case with no selection
probit ent_binary `X', r
*2. Ignoring selection
probit `Y' `X', r
*3. Feasible, but inefficient: Heckman correction assigning random half of surveyed pop to no-followup group - bysort is to give it gender-specific selection corrections
bysort $S_ML_groupvar: heckprobit `Y' `X', select(survey_response = `X' `Z') cluster(id)
restore

*4. Leverage panel variation (cluster SEs on i!) - we recommend interacting requests with all X in selection equation following Blandhol et al. (2022).
bysort $S_ML_groupvar: heckprobit `Y' `X', select(survey_response = `X' `Z') cluster(id)


heckprobit `Y' `X', select(survey_response = `X' `Z') cluster(id)
mat initial = e(b) //initialize preferred version below with single-rho estimates 

*We also recommend allowing for the possibility of different corr(Y,S) parameters between groups - can be done via bysort for group means, but differences between groups (e.g. treatment effects or gaps) are easier to calculate SEs for if you estimate everything jointly.
*------------------------------------
* 2.a Preferred estimator
*------------------------------------
*Code up MLE for binary outcomes - standard heckman correction with equations rather than scalars - any item (XB, Za, or rho_i = Xi*rho) can have heterogeneity by X. Note that convergence on rho will be touchy if far from zero.
capture program drop ml_hegsurvey_binary
program ml_hegsurvey_binary
    version 16
    args todo b lnf
    
    tempvar xb_main xb_s tau rho
    
    mleval `xb_main' = `b', eq(1)
    mleval `xb_s' = `b', eq(2)
    mleval `tau' = `b', eq(3)
    
    gen `rho' = (expm1(2*`tau')) / (exp(2*`tau')+1)
    
    *Code Likelihood
    quietly replace `lnf' = 0 //prime likelihood
    quietly replace `lnf' = `lnf'+log(1 - normal((`xb_s'))) if ${ML_y2}==0 //no response
    quietly replace `lnf' = `lnf'+log(binormal(`xb_main',`xb_s',`rho')) if ${ML_y2}==1 & $ML_y1==1 //Respond and have response=1
    quietly replace `lnf' = `lnf'+log(binormal(-(`xb_main'),`xb_s',-(`rho'))) if ${ML_y2}==1 & $ML_y1==0 //Respond and have response=0
end

ml model lf0 ml_hegsurvey_binary (Binary_Outcome: `Y' = `X') (Selection: survey_response = `X' `Z') (athrho_XB: = `X_rho'),  missing cluster(id)
ml maximize, iterate(0) //get size of beta vector

*get as many zeros as you need to initialize with parameters you have with copy option
if colsof(e(b))-colsof(initial)!=0 {
mat zeros = J(1,colsof(e(b))-colsof(initial),0)
	mat initial_now = initial, zeros
}
else {
	mat initial_now = initial
}
ml model lf0 ml_hegsurvey_binary (Binary_Outcome: `Y' = `X') (Selection: survey_response = `X' `Z') (athrho_XB: = `X_rho'),  missing cluster(id)
ml init initial_now, copy
ml maximize


estimates store corrected_estimates
mat initial = e(b)
mat XB_initial = e(b)[1, "Binary_Outcome:"] //main equation coefficients


/*Compare to uncorrected version if you want
constraint define 1 _b[athrho_XB:_cons]=0
ml model lf0 ml_hegsurvey_binary (Binary_Outcome: `Y' = `X') (Selection: survey_response = `X' `Z') (athrho_XB: = ),  missing cluster(id) constraints(1)
ml init initial, skip
ml maximize

probit `Y'  `X' if survey_response==1
*/

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
*Make figure showing means by Request group w/ predicted values.

*Prime strings
local sample_str = ""
local notice_str = ""

*Produce R=0 row
expand 2 if requests==1 , gen(R0)
replace requests = 0 if R0==1
replace survey_response = 0 if R0==1
replace `Y' = 0 if R0==1
drop R0

*Generate necessary objects for graph
local N_graphs = 1000 //set number of data points for illustrative graphs, doesn't need to equal N
capture set obs `N_graphs'

capture drop u
gen u = (_n)/`N_graphs' if _n<=`N_graphs'


*Estimate uncorrected model
reg `Y' if requests==4 & survey_response==1
predictnl `Y'_uncorr = _b[_cons], ci(y_uncorr_ciup y_uncorr_cidown)

*uncorrected graph string
local uncorrection_str = "(line `Y'_uncorr u, lw(*1) lcolor(gs8) xaxis(1 2)) (line y_uncorr_ciup u, lpattern(dash) lcolor(gs8)) (line y_uncorr_cidown u, lpattern(dash) lcolor(gs8))"


*Get response probabilities by R
reg survey_response i.requests, cluster(id)

local p_0 = 0
forv r = 1/4 {
	local p_`r' = _b[`r'.requests]
}

*Code Yhat = YS
replace `Y' = 0 if survey_response==0

*LARs
forv r = 1/4 {
	ivregress 2sls `Y' (survey_response = i.requests) if requests>=`=`r'-1' & requests<=`r'

	local LAR`r' = _b[survey_response]
	di `p`=`r'-1''
	di `p`=`r'''
	gen LAR_`r' = _b[survey_response] if u >`p_`=`r'-1'' & u <=`p_`=`r'''
	
	local sample_str = "`sample_str' (line LAR_`r' u, lw(*2) lcolor(black) xaxis(1 2))"
	
	
	local notice_str = "`notice_str' `p_`r''"
}


*Corrected model
estimates restore corrected_estimates
capture drop selhaz
gen selhaz = invnormal(1-u)
capture drop y_hat
capture drop y_corr_ciup y_corr_cidown
predictnl y_hat_corr = (selhaz*((expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1)) + _b[_cons])/(sqrt(1-((expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))^2)), ci(y_corr_ciup y_corr_cidown)

replace y_hat_corr = normal(y_hat_corr)
replace y_corr_ciup = normal(y_corr_ciup)
replace y_corr_cidown = normal(y_corr_cidown)



local correction_str = "`correction_str' (line y_hat_corr u, lw(*1) lcolor(blue) xaxis(1 2)) (line y_corr_ciup u, lpattern(dash) lcolor(blue)) (line y_corr_cidown u, lpattern(dash) lcolor(blue))"


twoway `sample_str' `uncorrection_str' `correction_str', ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10, nogrid) yscale(range(0, 1)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Within-Request Mean" 5 "Uncorrected Extrapolation" 8 "Corrected Extrapolation") ring(0) position(2) region(lstyle(1))) xla(0 "0" `p_1' "1" `p_2' "2" `p_3' "3"  `p_4' "4" , axis(2)) xtitle("Requests Received", axis(2)) xline(`notice_str', lp(solid))


drop if requests==0 //unnecessary for anything other than convenient LAR for r=1
}




*-------------------------------------------------------------------------------
* 3. Overidentification test using multiple requests
*-------------------------------------------------------------------------------
*Make some extra vars - careful to avoid collinearity, can cause convergence issues. Going through and creating c. type interactions rather than i., which produce a few redundant/collinear interations in different equations.
sum survey_ordered
local K = r(max)
tab requests, gen(requests_)
drop requests_1 requests_2 //omitted from main equation for overid test

local oid_var_count = 0
local X_oid_all = ""
forv k = 3/`K' { //for each request past 2
	local X_oid_all = "`X_oid_all' requests_`k'"
	local ++oid_var_count
	if "`X_oid'" == "" {
	}
	else {
	foreach var of varlist `X_oid' {
		if "`var'" == "$S_ML_groupvar" {
			local X_oid_all = "`X_oid_all' c.requests_`k'#c.$S_ML_groupvar"
		}
		else {
			local X_oid_all = "`X_oid_all' c.requests_`k'#i.$S_ML_groupvar#c.`var'" //the c. instead i. is used because all basic X variables are already defined in the X local
		}
		local ++oid_var_count
	}
	}
}


*------------------------------------
* 3.a Simple case using built in commands
*------------------------------------
heckprobit `Y' `X' io(1 2).requests if $S_ML_groupvar==0, select(survey_response = `X' `Z') cluster(id)
*produce p-value on test that all excess requests have no effect on outcome
test (_b[ent_binary_observed:3.requests]=0) (_b[ent_binary_observed:4.requests]=0)

if `graphs' == 0 {
heckprobit `Y' `X' io(1 2).requests if $S_ML_groupvar==1, select(survey_response = `X' `Z') cluster(id)
*produce p-value on test that all excess requests have no effect on outcome
test (_b[ent_binary_observed:3.requests]=0) (_b[ent_binary_observed:4.requests]=0)
}



*------------------------------------
* 3.b Test overidentification using jointly estimated model with heterogeneous rho and overid parameters (null of zero effects in terms of probit coefficients - beta_R in the paper). Estimate model.
*------------------------------------

*Jointly estimating men and women with gender/etc-specific parameters (different results from simple version above because rho is allowed to vary with all of X
*Initialize

*Tweaking program is convenient for testing below because we can put overid vars in their own "equation" in Stata MLE syntax, though the actual covariates and their coefficients will be in the main equation.
capture program drop ml_hegsurvey_binary_oid
program ml_hegsurvey_binary_oid
    version 16
    args todo b lnf
	
	tempvar xb_main xb_s tau xb_oid rho
	
	mleval `xb_main' = `b', eq(1)

	mleval `xb_s' = `b', eq(2)
	

	mleval `tau' = `b', eq(3) //replace first 3 with number of Requests
	
	*Simply add these to xb_main in main equation, we'll test that this whole predicted value is zero
	mleval `xb_oid' = `b', eq(4) //overid XB

	gen `rho' = (expm1(2*`tau'))  /  (exp(2*`tau')+1)
	
	
	*Code Likelihood
	quietly replace `lnf' = 0 //prime likelihood
	quietly replace `lnf' =  `lnf'+log(1 - normal((`xb_s'))) if ${ML_y2}==0 //no response
	quietly replace `lnf' =  `lnf'+log(binormal(`xb_main'+`xb_oid',`xb_s',`rho')) if ${ML_y2}==1 & $ML_y1==1 //Respond and have response=1
	quietly replace `lnf' =  `lnf'+log(binormal(-(`xb_main'+`xb_oid'),`xb_s',-(`rho'))) if ${ML_y2}==1 & $ML_y1==0 //Respond and have response=0	
		
end

ml model lf0 ml_hegsurvey_binary_oid (Binary_Outcome: `Y' = `X') (Selection: survey_response = `X' `Z') (athrho_XB: = `X_rho') (xb_oid: = `X_oid_all' , nocons),   missing cluster(id)
ml init initial, skip
*ml check
ml maximize

estimates store overid

*Confirm equivalent estimates to adding overid vars to main equation or dropping all late responders (true if you fully-flexibly specify heterogeneous effects of Requests across X)
/*
ml model lf0 ml_hegsurvey_binary (Binary_Outcome: `Y' = `X' `X_oid_all' ) (Selection: survey_response = `X' `Z') (athrho_XB: = `X_rho'),  missing cluster(id)
ml init initial, skip
*ml check
ml maximize
mat initial_oid = e(b)
mat XB_initial_plus_oid = e(b)[1, "Binary_Outcome:"] //main equation coefficients
mat XB_initial_sel = e(b)[1, "Selection:"] //selection equation coefficients
mat XB_initial_rho = e(b)[1, "athrho_XB:"] //rho equation coefficients

mat initial_oid = XB_initial_plus_oid[1,1..`=colsof(XB_initial)-1'] ///main X nocons
, XB_initial_plus_oid[1,colsof(XB_initial_plus_oid)] /// constant
, XB_initial_sel /// selection 
, XB_initial_rho /// rho
, XB_initial_plus_oid[1,colsof(XB_initial)..`=colsof(XB_initial_plus_oid)-1'] // overid coefficients in main eq


ml model lf0 ml_hegsurvey_binary (Binary_Outcome: `Y' = `X') (Selection: survey_response = `X' `Z') `Equations' (athrho_XB: = `X_rho') if requests <=2,  missing cluster(id)
*ml init initial, copy
ml maximize
*/

*------------------------------------
* 3.c Test coefficients = 0. There are a few ways to do this depending on whether you want ALL effects of requests to be zero or you want them to be zero on average. Standard heckman test is that rho=0 on average is a weak test of the model because LARS could vary wildly around m(u).
*------------------------------------
*This tests that all coefficients are jointly zero. Putting a bunch of low-info interactions seems likely to kill power for this, a bit arbitrary and at discretion of researcher, not our favorite.
estimates restore overid
test [xb_oid]

*Test overall average overidentification effect is zero - good if you lack N to estimate a model with a lot of heterogeneity and you don't think rho varies between groups. This assumption is violated in this DGP, so not recommended here but shown for reference.
estimates restore overid
margins if requests>2, predict(xb equation(xb_oid))

*Test that overall *average* overidentification effect is zero *for each group* - reasonable if you are concerned about a group-specific rho, as we have done here to match the DGP. If deviations of model predictions from true responses are large but average to zero, this will fail to reject despite bad misspecification.
estimates restore overid
margins if requests>2, predict(xb equation(xb_oid)) over($S_ML_groupvar) post

*The above tests for each group that their individual overid predicted value is zero for all requests on average - joint test is better
if `graphs' == 0 {
test (_b[0bn.$S_ML_groupvar]=0) (_b[1.$S_ML_groupvar]=0)
}

*This test is whether all deviations from predictions are zero for predetermined groups that are important, akin to parallel trends test. Our preferred test.
estimates restore overid
margins , predict(xb equation(xb_oid)) post over($S_ML_groupvar requests) coefleg
mat r_table = r(table)

local test_0 "(_b[0bn.$S_ML_groupvar#3.requests]=0)"
local test_1 "(_b[1.$S_ML_groupvar#3.requests]=0)"

qui sum requests
local max_R = `r(max)'
forv k = 4/`max_R' {
	local test_0 "`test_0' (_b[0bn.$S_ML_groupvar#`k'.requests]=0)"
	local test_1 "`test_1' (_b[1.$S_ML_groupvar#`k'.requests]=0)"
}

*Test men 
test `test_0'
*Test women
if `graphs' == 0 {
test `test_1'
*Test both jointly - preferred
test `test_0' `test_1'
}

*Graph event study of coefficients that are deviations from model predictions for subjects who respond to or prior to request r for all r.
matrix coef_B0 = J(`max_R',3,0)
if `graphs' == 0 {
matrix coef_B1 = J(`max_R',3,0)
}
local es_label_string = ""
forv k = 3/`max_R' {
	*get B, lb, and ub for 95% confidence
	*Group 0 men
	matrix coef_B0[`k',1] = r_table[1,`k'], r_table[5,`k'], r_table[6,`k']
	*Group 1 women
	if `graphs' == 0 {
	matrix coef_B1[`k',1] = r_table[1,`max_R'+`k'], r_table[5,`max_R'+`k'], r_table[6,`max_R'+`k']
	}
	
	local es_label_string "`es_label_string' `=`k'' "`k'""
}

if `graphs' == 0 {
	coefplot (matrix(coef_B0[,1]), ci((coef_B0[,2] coef_B0[,3]))) ///
(matrix(coef_B1[,1]), ci((coef_B1[,2] coef_B1[,3]))) ///
, vertical ytitle("Effects of Requests on Outcome") xtitle("Request") yline(0) ciopts(recast(rcaps)) xlab(1 "1" 2 "2" `es_label_string') legend(order(1 "Group 0" 3 "Group 1"))
}
else {
	coefplot (matrix(coef_B0[,1]), ci((coef_B0[,2] coef_B0[,3]))) ///
, vertical ytitle("Effects of Requests on Outcome") xtitle("Request") yline(0) ciopts(recast(rcaps)) xlab(1 "1" 2 "2" `es_label_string') legend(order(1 "Group 0"))
}


if `Y_oid' == 1 {
/*------------------------------------
* 3.d Test overidentification for joint model and produce figure in terms of predicted LAR (from model estimated on first 2 requests) vs actual LAR (from model with requests 3+ in main equation)

This makes a slightly more interpretable event study figure, but is much more complicated conceptually and computationally than just testing whether beta=0 on requests in the main equation. It is computationally prohibitive for the real application given the lack of benefit relative to testing beta. We show it here, but it is slow (It takes 1 minute on my Dell Precision 7560 Laptop with N=200 and NT=1000).


------------------------------------
*/
/*
*we're going to expand the data and generate an indicator A=1 for the duplicate observations. 
Then our command
>margins, expression (binormal(XB,Za,rho)*A+normal(Za)*)1-A)) over(A $S_ML_groupvar Requests) post
is going to give us Pr(YS=1) AND Pr(S=1) depending on the value of A. We can then play games with these averages to construct IV-estimates for LARs that look like
PR(YS|R=r)-PR(YS|R=r-1) /
PR(S|R=r) - P(S|R=r-1)
*/


/*
First, confirm we can match true conditional probabilities with just-id estimates (request dummies)
estimates restore overid
margins, expression(binormal(xb(Binary_Outcome)+xb(xb_oid),xb(Selection),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))) over(requests $S_ML_groupvar) //this is just the joint probabilities, check to make sure the sketchy trick you do below matches this

*While we're at it, I also want to match the marginal probability of S
margins, expression(normal(xb(Selection))) over(requests $S_ML_groupvar)

expand 2, gen(A) //indicator for cheesing stata into doing both joint and marginal probabilities in a single expression

*It really is this easy.
margins, expression(binormal(xb(Binary_Outcome)+xb(xb_oid),xb(Selection),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))*(1-A) ///
+ normal(xb(Selection))*A ///
) over(A requests $S_ML_groupvar) noesample post coefleg

estimates restore overid
*Now we can get predicted local average responses from the model
margins, expression(binormal(xb(Binary_Outcome)+xb(xb_oid),xb(Selection),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))*(1-A) ///
+ normal(xb(Selection))*A ///
) over(A $S_ML_groupvar requests) noesample post coefleg

*Request 1 compliers are different unless you explicitly code requests=0 observations (minus 0 added to show this is actually the same formula, S=0 when R=0)
nlcom (_b[0bn.A#0bn.$S_ML_groupvar#1bn.requests]-0) ///
/(_b[1.A#0bn.$S_ML_groupvar#1bn.requests]-0)

*Request 2 different because Stata has bn etc prefixes
nlcom (_b[0bn.A#0bn.$S_ML_groupvar#2.requests]-_b[0bn.A#0bn.$S_ML_groupvar#1bn.requests]) ///
/(_b[1.A#0bn.$S_ML_groupvar#2.requests]-_b[1.A#0bn.$S_ML_groupvar#1bn.requests])

qui sum requests
local R_max = r(max)
*All later requests follow the same pattern
forv r=3/`R_max' {
	nlcom (_b[0bn.A#0bn.$S_ML_groupvar#`r'.requests]-_b[0bn.A#0bn.$S_ML_groupvar#`=`r'-1'.requests]) ///
/(_b[1.A#0bn.$S_ML_groupvar#`r'.requests]-_b[1.A#0bn.$S_ML_groupvar#`=`r'-1'.requests])
}


*Check that it matches the average values for each complier group - yes it does.
bysort $S_ML_groupvar survey_ordered: sum ent_binary if survey_response==1
*/

*Now, code up the thing you actually want.

*Cheese Stata into doing P(YS=1) without requests, P(YS=1) with request effects, and P(S), for each r.

expand 2, gen(A) //indicator for P(r)

expand 2, gen(B) //indicator for over-id models

estimates restore overid
margins, expression( ///
binormal(xb(Binary_Outcome),xb(Selection),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))*(1-A)*(1-B) ///
+binormal(xb(Binary_Outcome)+xb(xb_oid),xb(Selection),(expm1(2*xb(athrho_XB)))/(exp(2*xb(athrho_XB))+1))*(1-A)*B ///
+ normal(xb(Selection))*A*B ///
) over(B A requests $S_ML_groupvar) noesample post coefleg

estimates store oid_margins


*Code difference in predicted LARs vs actual LARs.

*Numerator is the same, so we have
/*
For model 1 identified off r=1 and r=2,
and model 2 having request indicators in main equation
E[YS_1|r]-E[YS_1|r']
/(P(r)-P(r')
-
E[YS_2|r]-E[YS_2|r']
/(P(r)-P(r'))
=
(E[YS_1|r]-E[YS_2|r']-(E[YS_2|r]-E[YS_2|r']))
/(P(r)-P(r'))

*/

estimates restore oid_margins
local nlcom_str = ""

*For group 0 (Men)

*Request 1 (different because no request 0) - should be 0 difference because model fits first 2 requests
nlcom (((_b[1.B#0bn.A#1bn.requests#0bn.$S_ML_groupvar]-0) ///
-(_b[0bn.B#0bn.A#1bn.requests#0bn.$S_ML_groupvar]-0)) ///
/(_b[1.B#1.A#1bn.requests#0bn.$S_ML_groupvar]-0))
local nlcom_str = "`nlcom_str' (((_b[1.B#0bn.A#1bn.requests#0bn.$S_ML_groupvar]-0) -(_b[0bn.B#0bn.A#1bn.requests#0bn.$S_ML_groupvar]-0)) /(_b[1.B#1.A#1bn.requests#0bn.$S_ML_groupvar]-0))"

*Request 2 (different because of stata syntax) - should be 0 difference because model fits first 2 requests
nlcom (((_b[1.B#0bn.A#2.requests#0bn.$S_ML_groupvar]-_b[1.B#0bn.A#1bn.requests#0bn.$S_ML_groupvar]) ///
-(_b[0bn.B#0bn.A#2.requests#0bn.$S_ML_groupvar]-_b[0bn.B#0bn.A#1bn.requests#0bn.$S_ML_groupvar])) ///
/(_b[1.B#1.A#2.requests#0bn.$S_ML_groupvar]-_b[1.B#1.A#1bn.requests#0bn.$S_ML_groupvar]))

local nlcom_str = "`nlcom_str' (((_b[1.B#0bn.A#2.requests#0bn.$S_ML_groupvar]-_b[1.B#0bn.A#1bn.requests#0bn.$S_ML_groupvar]) -(_b[0bn.B#0bn.A#2.requests#0bn.$S_ML_groupvar]-_b[0bn.B#0bn.A#1bn.requests#0bn.$S_ML_groupvar])) /(_b[1.B#1.A#2.requests#0bn.$S_ML_groupvar]-_b[1.B#1.A#1bn.requests#0bn.$S_ML_groupvar])) "

qui sum requests
local R_max = r(max)
*Request 3+
forv r = 3/`R_max' {
	nlcom (((_b[1.B#0bn.A#`r'.requests#0bn.$S_ML_groupvar]-_b[1.B#0bn.A#`=`r'-1'.requests#0bn.$S_ML_groupvar]) ///
-(_b[0bn.B#0bn.A#`r'.requests#0bn.$S_ML_groupvar]-_b[0bn.B#0bn.A#`=`r'-1'.requests#0bn.$S_ML_groupvar])) ///
/(_b[1.B#1.A#`r'.requests#0bn.$S_ML_groupvar]-_b[1.B#1.A#`=`r'-1'.requests#0bn.$S_ML_groupvar]))

local nlcom_str = "`nlcom_str' (((_b[1.B#0bn.A#`r'.requests#0bn.$S_ML_groupvar]-_b[1.B#0bn.A#`=`r'-1'.requests#0bn.$S_ML_groupvar]) -(_b[0bn.B#0bn.A#`r'.requests#0bn.$S_ML_groupvar]-_b[0bn.B#0bn.A#`=`r'-1'.requests#0bn.$S_ML_groupvar])) /(_b[1.B#1.A#`r'.requests#0bn.$S_ML_groupvar]-_b[1.B#1.A#`=`r'-1'.requests#0bn.$S_ML_groupvar])) "
}

nlcom `nlcom_str' 

if `graphs' == 0 {
*Now for group 1 (Female)
*Request 1 (different because no request 0) - should be 0 difference because model fits first 2 requests
nlcom (((_b[1.B#0bn.A#1bn.requests#1.$S_ML_groupvar]-0) ///
-(_b[0bn.B#0bn.A#1bn.requests#1.$S_ML_groupvar]-0)) ///
/(_b[1.B#1.A#1bn.requests#1.$S_ML_groupvar]-0))
local nlcom_str = "`nlcom_str' (((_b[1.B#0bn.A#1bn.requests#1.$S_ML_groupvar]-0) -(_b[0bn.B#0bn.A#1bn.requests#1.$S_ML_groupvar]-0)) /(_b[1.B#1.A#1bn.requests#1.$S_ML_groupvar]-0))"

*Request 2 (different because of stata syntax) - should be 0 difference because model fits first 2 requests
nlcom (((_b[1.B#0bn.A#2.requests#1.$S_ML_groupvar]-_b[1.B#0bn.A#1bn.requests#1.$S_ML_groupvar]) ///
-(_b[0bn.B#0bn.A#2.requests#1.$S_ML_groupvar]-_b[0bn.B#0bn.A#1bn.requests#1.$S_ML_groupvar])) ///
/(_b[1.B#1.A#2.requests#1.$S_ML_groupvar]-_b[1.B#1.A#1bn.requests#1.$S_ML_groupvar]))

local nlcom_str = "`nlcom_str' (((_b[1.B#0bn.A#2.requests#1.$S_ML_groupvar]-_b[1.B#0bn.A#1bn.requests#1.$S_ML_groupvar]) -(_b[0bn.B#0bn.A#2.requests#1.$S_ML_groupvar]-_b[0bn.B#0bn.A#1bn.requests#1.$S_ML_groupvar])) /(_b[1.B#1.A#2.requests#1.$S_ML_groupvar]-_b[1.B#1.A#1bn.requests#1.$S_ML_groupvar])) "

qui sum requests
local R_max = r(max)
*Request 3+
forv r = 3/`R_max' {
	nlcom (((_b[1.B#0bn.A#`r'.requests#1.$S_ML_groupvar]-_b[1.B#0bn.A#`=`r'-1'.requests#1.$S_ML_groupvar]) ///
-(_b[0bn.B#0bn.A#`r'.requests#1.$S_ML_groupvar]-_b[0bn.B#0bn.A#`=`r'-1'.requests#1.$S_ML_groupvar])) ///
/(_b[1.B#1.A#`r'.requests#1.$S_ML_groupvar]-_b[1.B#1.A#`=`r'-1'.requests#1.$S_ML_groupvar]))

local nlcom_str = "`nlcom_str' (((_b[1.B#0bn.A#`r'.requests#1.$S_ML_groupvar]-_b[1.B#0bn.A#`=`r'-1'.requests#1.$S_ML_groupvar]) -(_b[0bn.B#0bn.A#`r'.requests#1.$S_ML_groupvar]-_b[0bn.B#0bn.A#`=`r'-1'.requests#1.$S_ML_groupvar])) /(_b[1.B#1.A#`r'.requests#1.$S_ML_groupvar]-_b[1.B#1.A#`=`r'-1'.requests#1.$S_ML_groupvar])) "
}

}

nlcom `nlcom_str', post
mat r_table = r(table)
matlist r_table

local test_nl_str_0 = ""
local test_nl_str_1 = ""
forv r = 3/`R_max' {
	local test_nl_str_0 = "`test_nl_str_0' (_b[_nl_`r'] = 0)"
	local test_nl_str_1 = "`test_nl_str_1' (_b[_nl_`=`R_max'+`r''] = 0)"
}

test (_b[_nl_3] = 0)

test `test_nl_str_0' 
if `graphs' == 0 {
test `test_nl_str_1' 
test `test_nl_str_0' `test_nl_str_1' 
}



/* This checks out
*Confirm that result is consistent with predicted LAR...
nlcom (_b[0bn.B#0bn.A#3.requests#0bn.Female]-_b[0bn.B#0bn.A#2.requests#0bn.Female]) ///
/(_b[1.B#1.A#3.requests#0bn.Female]-_b[1.B#1.A#2.requests#0bn.Female])

nlcom (_b[0bn.B#0bn.A#4.requests#0bn.Female]-_b[0bn.B#0bn.A#3.requests#0bn.Female]) ///
/(_b[1.B#1.A#4.requests#0bn.Female]-_b[1.B#1.A#3.requests#0bn.Female])

*And actual LAR
bysort $S_ML_groupvar survey_ordered: sum ent_binary if survey_response==1
*/



*Produce event study figure

local es_label_string = ""
forv k = 1/`max_R' {
	if r_table[5,`k']==. {
		mat r_table[5,`k'] = 0
		mat r_table[6,`k'] = 0
		if `graphs' == 0 {
		mat r_table[5,`max_R'+`k'] = 0
		mat r_table[6,`max_R'+`k'] = 0
		}
	}
	*get B, lb, and ub for 95% confidence
	*Group 0 men
	matrix coef_B0[`k',1] = r_table[1,`k'], r_table[5,`k'], r_table[6,`k']
	*Group 1 women
	if `graphs' == 0 {
	matrix coef_B1[`k',1] = r_table[1,`max_R'+`k'], r_table[5,`max_R'+`k'], r_table[6,`max_R'+`k']
	}
	
	local es_label_string "`es_label_string' `=`k'' "`k'""
}


if `graphs' == 0 {
	coefplot (matrix(coef_B0[,1]), ci((coef_B0[,2] coef_B0[,3]))) ///
(matrix(coef_B1[,1]), ci((coef_B1[,2] coef_B1[,3]))) ///
, vertical ytitle("Effects of Requests on Outcome") xtitle("Request") yline(0) ciopts(recast(rcaps)) xlab(1 "1" 2 "2" `es_label_string') legend(order(1 "Group 0" 3 "Group 1"))
}
else {
	coefplot (matrix(coef_B0[,1]), ci((coef_B0[,2] coef_B0[,3]))) ///
, vertical ytitle("Effects of Requests on Outcome") xtitle("Request") yline(0) ciopts(recast(rcaps)) xlab(1 "1" 2 "2" `es_label_string') legend(order(1 "Group 0" ))
}
}







*------------------------------------
