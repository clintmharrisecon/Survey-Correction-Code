/*
This file uses our method on response group means from "SELECTION IN SURVEYS: USING RANDOMIZED INCENTIVES TO DETECT AND ACCOUNT FOR NONRESPONSE BIAS", downloaded on 2/24/2025 from NBER.
*/
clear
cls

set maxiter 100
set seed 123
version 18.0

*Generate data that matches moments in the paper (Table 4)

*10,000 survey individuals in NCT survey.
set obs 10000

//generate participation aversion percentiles (for illustrations)
gen U = _n/_N


*Define ordering of participation aversion
gen survey_1_notification = U<=.38
gen survey_2_notification = U<=.38+.07
gen survey_ordered = survey_1_notification+survey_2_notification

*Enter data from Table 4 for outcomes for Always Takers (people who respond with no incentive and no reminder)
*We want to get comparable SEs as DHLMTV et al. for each outcome.
/*
SD = SE*sqrt(N-1)

Where SE is known from DHLMTV et al, and N=10000.
*/
*Continuous
local SD = 116*sqrt(.38*_N-1)
gen earnings_before = 3746+`SD' if _n<=.38*_N*0.5 //add SD for half
replace earnings_before = 3746-`SD' if _n>.38*_N*0.5 & _n<=.38*_N //subtract SD for half
local SD = 107*sqrt(.38*_N-1)
gen earnings_after = 3783+`SD' if _n<=.38*_N*0.5 //add SD for half
replace earnings_after = 3783-`SD' if _n>.38*_N*0.5 & _n<=.38*_N //subtract SD for 
*Binary (make the appropriate percentages 1 and 0 for the relevant group)
gen earnings_large_loss = 1 if _n<=.38*_N*.13
replace earnings_large_loss = 0 if _n<=.38*_N & earnings_large_loss==.
gen employment_before = 1 if _n<=.38*_N*.65
replace employment_before = 0 if _n<=.38*_N & employment_before==.
gen employment_after = 1 if _n<=.38*_N*.64
replace employment_after = 0 if _n<=.38*_N & employment_after==.
gen employment_loss = 1 if _n<=.38*_N*.03
replace employment_loss = 0 if _n<=.38*_N & employment_loss==.

*Enter Data from Table 4 for outcomes for Reminder Compliers (people who respond with no incentive after reminders)
*Continuous
local SD = 256*sqrt(.07*_N-1)
replace earnings_before = 3244+`SD' if _n<=(.38+.07*0.5)*_N & earnings_before==.
replace earnings_before = 3244-`SD' if _n<=(.38+.07)*_N & earnings_before==.
local SD = 251*sqrt(.07*_N-1)
replace earnings_after = 3257+`SD' if _n<=(.38+.07*0.5)*_N & earnings_after==.
replace earnings_after = 3257-`SD' if _n<=(.38+.07)*_N & earnings_after==.
*Binary (make the appropriate percentages 1 and 0 for the relevant group)
replace earnings_large_loss = 1 if _n<=.38*_N+.07*_N*.12 & earnings_large_loss==.
replace earnings_large_loss = 0 if _n<=.38*_N+.07*_N & earnings_large_loss==.
replace employment_before = 1 if _n<=.38*_N+.07*_N*.55 & employment_before==.
replace employment_before = 0 if _n<=.38*_N+.07*_N & employment_before==.
replace employment_after = 1 if _n<=.38*_N+.07*_N*.55 & employment_after==.
replace employment_after = 0 if _n<=.38*_N+.07*_N & employment_after==.
replace employment_loss = 1 if _n<=.38*_N+.07*_N*.03 & employment_loss==.
replace employment_loss = 0 if _n<=.38*_N+.07*_N & employment_loss==.

*Code linear correction (one standard heckman correction for each reminder, with 1 set of outcome/correlation parameters)
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
	
	tempname rho sig lambda
	scalar `rho' = (expm1(2*`tau'))  /  (exp(2*`tau')+1)
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

sum survey_ordered
global notification_count = r(max)
di $notification_count
*Always have 
local Equations "`Equations' (Selection1: survey_1_notification = `controls')"
forv k = 2/$notification_count {
	local Equations "`Equations' (Selection`k': survey_`k'_notification = `controls')"
}
di "`Equations'"

matrix outcomes = J(2,6,.)
local j = 1

foreach outcome in "earnings_before" "earnings_after" {
ml model lf0 ml_stacked_ordered_heckman (Linear_s2: `outcome' = `controls') `Equations' /athrho /lnsigma,  missing vce(robust)
*ml init initial, copy
ml check
ml search
ml maximize

mat outcomes[1,`j'] = _b[_cons]
mat outcomes[2,`j'] = _se[_cons]
local ++j
}



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
local Equations ""
*Always have (add controls you want here)
local Equations "`Equations' (Selection1: survey_1_notification = `controls')"
forv k = 2/`notification_count' {
               local Equations "`Equations' (Selection`k': survey_`k'_notification = `controls')"
}
di "`Equations'"


foreach outcome in "earnings_large_loss" "employment_before" "employment_after" "employment_loss" {
ml model lf0 ml_stacked_heckman_bivariate (Binary_Outcome: `outcome' = ) `Equations' /athrho,  missing vce(robust)
*ml init initial, copy
ml maximize

margins, expression(normal(xb(Binary_Outcome)))
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

*2 step w/ bootstrap SEs
cap prog drop heg_survey_bootstrap
prog heg_survey_bootstrap, eclass
syntax varlist
qui oprobit survey_ordered 
qui predict imr_inner, xb
qui replace imr_inner = 0 if imr_inner==. & survey_ordered>0
qui gen double gimr = (normalden(imr_inner-[/]cut2)-normalden(imr_inner-[/]cut1))/(normal(imr_inner-[/]cut2)-normal(imr_inner-[/]cut1)) if survey_ordered ==1
replace gimr = (normalden(imr_inner-[/]cut2))/(normal(imr_inner-[/]cut2)) if survey_ordered ==2
reg `varlist' gimr 
cap drop imr_inner gimr 
end

local j=1
foreach outcome in "earnings_before" "earnings_after" "earnings_large_loss" "employment_before" "employment_after" "employment_loss" {
bootstrap _b , reps(50) nowarn nodots nodrop seed(123): heg_survey_bootstrap `outcome'

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


*Make visual illustrations - use first outcome (earnings_before)
oprobit survey_ordered 
*Generate ordered probit IMR
predict imr_inner, xb
replace imr_inner = 0 if imr_inner==. & survey_ordered>0

gen imr = (normalden(imr_inner-[/]cut2)-normalden(imr_inner-[/]cut1))/(normal(imr_inner-[/]cut2)-normal(imr_inner-[/]cut1)) if survey_ordered ==1
replace imr = (normalden(imr_inner-[/]cut2))/(normal(imr_inner-[/]cut2)) if survey_ordered ==2

foreach outcome in "earnings_before" "earnings_after" "earnings_large_loss" "employment_before" "employment_after" "employment_loss" {
reg `outcome' imr 
}
reg earnings_before imr 

local N_graphs = 10000 //set number of data points for illustrative graphs, doesn't need to equal N
capture set obs `N_graphs'

forv k = 1/`=$notification_count' {
               di `k'
               
               sum earnings_before if survey_ordered == `k'
               local notification_`k'_mean = r(mean)
               di `notification_`k'_mean'
               
               count if survey_ordered == `k'
               local share_group_`k' = r(N)
               local share_group_`k' = `share_group_`k''/`N_graphs'
               di `share_group_`k''
               
               
               *sum survey_`k'_notification
               *local share_notification_`=$notification_count-`k'' = r(mean)
}



*left to right
capture drop Y_bar_data
gen Y_bar_data = .
*replace Y_bar_data = 0 if _n<=100 & _n>`=`share_rem_2'*100'
replace Y_bar_data = `notification_2_mean' if _n<=`=`share_group_2'*`N_graphs''
replace Y_bar_data = `notification_1_mean' if _n<=`=(`share_group_1'+`share_group_2')*`N_graphs'' & _n>`=`share_group_2'*`N_graphs''


*Same thing for percentile rank derived from imr
sum imr if survey_ordered==2
local group_2_imr = r(mean)
sum imr if survey_ordered==1
local group_1_imr = r(mean)


capture drop p_imr_data
gen p_imr_data = .
*replace Y_bar_data = 0 if _n<=100 & _n>`=`share_rem_2'*100'
replace p_imr_data = normal(`=-1*`group_2_imr'')*100 if _n<=`=`share_group_2'*`N_graphs''
replace p_imr_data = normal(`=-1*`group_1_imr'')*100 if _n<=`=(`share_group_1'+`share_group_2')*`N_graphs'' & _n>`=`share_group_2'*`N_graphs''


capture drop U_imr_data
gen U_imr_data = p_imr_data/100


capture drop p_data
gen p_data = .
replace p_data = `share_group_2'/2*100 if _n<=`=`share_group_2'*`N_graphs''
replace p_data = (`share_group_2'+(`share_group_1')/2)*100 if _n<=`=(`share_group_1'+`share_group_2')*`N_graphs'' & _n>`=`share_group_2'*`N_graphs''

capture drop U_data
gen U_data = p_data/100


*Within Sample mean
reg earnings_before //for SEs
capture drop Y_bar_mean
gen Y_bar_mean = _b[_cons] if U<.

*Normal (Heckman) Correction extrapolation 
capture drop imr_U
gen imr_U = invnormal(1-U) //the inverse mills ratio is just a guess at a person's t-stat - it is more complicated usually because we DON'T know people's unobserved preferences - for the graph we do.

*For illustrative graph
reg earnings_before imr
capture drop Y_norm_corr
gen Y_norm_corr = imr_U*_b[imr]+_b[_cons]

*Generate linear extrapolation
reg earnings_before U_data
local b_U_mean = _b[U_data]
local Y_cons = _b[_cons]
capture drop Y_lin_corr
gen Y_lin_corr = _b[U_data]*U+_b[_cons]

*Get y_range from max/min of objects in all graphs
qui sum Y_norm_corr
local y_min_norm = r(min)
local y_max_norm = r(max)
qui sum Y_lin_corr
local y_min_lin = r(min)
local y_max_lin = r(max)

local y_min = min(`y_min_norm', `y_min_lin')
local y_max = max(`y_max_norm', `y_max_lin')


capture gen pop_mean = 3095 //Ground Truth from Table 5

*Graphs (some other plausibly interesting graphs are commented out below)
*Add percentile rank of average imr for individuals in each bin
twoway (line pop_mean U, xaxis(1 2)) (scatter Y_bar_data U_imr_data, xaxis(1 2)) (line Y_bar_mean U, lp(dash)) , ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Earnings Before") legend(order(1 "Population Mean" 2 "Mean by Response Time" 3 "Sample Mean") position(6) region(lcolor(black))) xla(.38 "1" .45 "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(.38 .45, lp(solid))

*Normal Extrapolation
twoway (line pop_mean U, xaxis(1 2)) (scatter Y_bar_data U_imr_data, xaxis(1 2)) (line Y_bar_mean U, lp(dash)) (line Y_norm_corr U, lp(dash)) , ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10) yscale(range(`y_min', `y_max')) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Earnings Before") legend(order(1 "Population Mean" 2 "Mean by Response Time" 3 "Sample Mean" 4 "Normal Extrapolation") position(6) region(lcolor(black)))  xla(.38 "1" .45 "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(.38 .45, lp(solid))

/*
*Show data for average survey preference percentile with average outcome for each bin.
twoway (scatter Y_bar_data U_data if Y_bar_data<., xaxis(1 2)), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10) yscale(range(2500, 4500)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile")) xla(.38 "1" .45 "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(.38 .45)

*Add Within-sample mean
twoway (scatter Y_bar_data U_data if Y_bar_data<., xaxis(1 2)) (line Y_bar_mean U), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10) yscale(range(2500, 4500)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome")) xla(.38 "1" .45 "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(.38 .45)

*Add linear extrapolation
twoway (scatter Y_bar_data U_data if Y_bar_data<., xaxis(1 2)) (line Y_bar_mean U) (line Y_lin_corr U), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10) yscale(range(2500, 4500)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome" 3 "Linear Extrapolation")) xla(.38 "1" .45 "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(.38 .45)

*Add percentile rank of average imr for individuals in each bin
twoway (scatter Y_bar_data U_data if Y_bar_data<., xaxis(1 2)) (line Y_bar_mean U) (line Y_lin_corr U) (scatter Y_bar_data U_imr_data) , ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10) yscale(range(2500, 4500)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome" 3 "Linear Extrapolation" 4 "Percentile, Mean(IMR)")) xla(.38 "1" .45 "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(.38 .45)

*Add correction under normality assumption
twoway (scatter Y_bar_data U_data if Y_bar_data<., xaxis(1 2)) (line Y_bar_mean U) (line Y_lin_corr U) (scatter Y_bar_data U_imr_data) (line Y_norm_corr U), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) xlabel(#10) yscale(range(2500, 4500)) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Mean Percentile" 2 "Sample Mean Outcome" 3 "Linear Extrapolation" 4 "Percentile, Mean(IMR)" 5 "Normal Extrapolation")) xla(.38 "1" .45 "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(.38 .45)


*Add population mean
capture gen pop_mean = 3095.8 //from Table A.3
twoway (line pop_mean U, xaxis(1 2)) (line Y_bar_mean U) (line Y_lin_corr U) (line Y_norm_corr U), ylabel(, nogrid) graphregion(color(white)) plotregion(color(white)) xscale(range(0,1)) yscale(range(2500, 4500))  xlabel(#10) ylabel(#10) xtitle("Survey Aversion Percentile") ytitle("Outcome") legend(order(1 "Population Mean" 2 "Sample Mean" 3 "Linear Extrapolation" 4 "Normal Extrapolation")) xla(.38 "1" .45 "2", axis(2)) xtitle("Notifications Received", axis(2)) xline(.38 .45)
*/


