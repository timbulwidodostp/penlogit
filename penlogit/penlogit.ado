/*
pllf2 is an adapted code of the -pllf- command by Royston
*/

*! 1.0.1 Discacciati A, Orsini N 21feb2014
*! 1.0.2 Discacciati A, Orsini N 01sep2014
*! 1.1.0 Discacciati A, Orsini N 19jul2017

capture program drop penlogit
program penlogit, eclass properties(mi) byable(onecall)
version 12 

if _by() {
	local BY `"by `_byvars'`_byrc0':"'
}
	
`version' `BY' _vce_parserun penlogit, mark(OFFset CLuster) : `0'

if "`s(exit)'" != "" {
	version 12: ereturn local cmdline `"penlogit `0'"'
	exit
}

if replay() {
	if ("`e(cmd)'"!="penlogit")  error 301  
	Replay `0'
}
else `version' `BY' Estimate `0'
end

capture program drop Estimate
program Estimate, eclass byable(recall)
version 12
syntax varlist [if] [in] [fweight],  ///
[ ///	
	NPrior(string) ///
	LFPrior(string) ///
	ppl(string)  ///
	BINomial(string) ///
	Level(integer $S_level) ///
	or ///
	NOList ///
	nppl(integer 100) ///
	NOCONstant ///
]


local cmdline : copy local 0
marksample touse 

gettoken depv indepv : varlist	
noi _rmcoll `indepv' if `touse' , forcedrop `=cond("`noconstant'"=="noconstant","noconstant"," ")'
local indepv `r(varlist)' 

local nr_profile : word count `ppl'


*Normal priors
local nprior_elements : word count `nprior'
if mod(`nprior_elements', 3) != 0 {
	di as err "incorrect number of elements in nprior()"
	exit 198
}
local npriors = `nprior_elements' / 3

local npriornames ""
local priorm ""
local priorv ""
tokenize "`nprior'"

forv j = 1/3 {
	forv k =`j'(3)`nprior_elements' {
		if (`j'==1) {
			confirm numeric variable ``k''
			local npriornames "`npriornames' ``k''"
		}
		if (`j'==2) {
			local priorm "`priorm' `=``k'''"
		}
		if (`j'==3) {
			local priorv "`priorv' `=``k'''"
		}
	}
}

unab npriornames : `npriornames' , min(`npriors') max(`npriors')

foreach v of local npriornames {
	if (`: list v in indepv' != 1) {
		di as err "prior `v' is not among [indepvars]"
		exit 198
	}
}
*

*Generalized Log-f priors
local lfprior_elements: word count `lfprior'
if mod(`lfprior_elements', 5) != 0 {
	di as err "incorrect number of elements in lfprior()"
	exit 198
}
local lfpriors = `lfprior_elements' / 5

local lfpriornames ""
local mm ""
local d0 ""
local d1 ""
local ss ""
tokenize "`lfprior'"

forv j = 1/5 {
	forv k =`j'(5)`lfprior_elements' {
		if (`j'==1) {
			confirm numeric variable ``k''
			local lfpriornames "`lfpriornames' ``k''" 
		}
		if (`j'==2) {
			local mm "`mm' `=``k'''"
		}
		if (`j'==3) {
			local d1 "`d1' `=``k'''"
		}
		if (`j'==4) {
			local d0 "`d0' `=``k'''"
		}
		if (`j'==5) {
			local ss "`ss' `=``k'''"
		}
	}
}

unab lfpriornames : `lfpriornames', min(`lfpriors') max(`lfpriors')

foreach v of local lfpriornames {
	if (`: list v in indepv' != 1) {
		di as err "prior `v' is not among [indepvars]"
		exit 198
	}
}
*

*Check for repeated priors
local names "`lfpriornames' `npriornames'"
local dup: list dups names
if "`dup'" != "" {
	local rep: list uniq dup
	di as err "prior(s) for `dup' specified more than once"
	exit 198
}
*

*Check penalized profile likelihood variables
if "`ppl'" != "" {
	local bign = c(N)

	capture unab ppl : `ppl' , min(`nr_profile') max(`nr_profile')

	foreach v of local ppl {
		if (`: list v in indepv' != 1) {
			di as err "var `v' (option ppl) is not among [indepvars]"
			exit 198
		}
	}
}
*

tempvar constant H n tag
tempname A prior_v H prior_logrr lb ub s lfc lfn lfexpoff nc nexpoff retnorm retlogf subjn

gen `constant' = 1
gen `H' = 0

if "`binomial'" == "" {
	gen `n' = 1
}
else {
	confirm variable `binomial'
	gen `n' = `binomial'
}


if "`binomial'" != "" {
	capture  assert `n' >= `depv' if `touse' & !missing(`depv'), fast
	if _rc != 0 {
		di as err "`depv' > `binomial' in some cases"
		exit 499
	}
}
else {
	_rmcoll `depv' if `touse', logit
}


tempvar wgt
if "`weight'" == ""  {
     gen `wgt' = 1
}
else  {
	local weights [`weight'`exp']
 	gen `wgt' `exp'
}

gen `tag' = 0

* Normal priors
if `npriors' != 0 {
	mat `retnorm' = J(2,`npriors',.)
	mat rownames `retnorm' = m v
	mat colnames `retnorm' = `npriornames'
}
local c = 1
qui foreach v of local npriornames {
	qui set obs `=_N+1'
	local prior_logrr : word  `c' of `priorm'
	local prior_v : word  `c' of `priorv'
	
	mat `retnorm'[1,`c'] = `prior_logrr'
	mat `retnorm'[2,`c'] = `prior_v'
	
	scalar `A' = 900
	scalar `s' = (`prior_v'/(2/scalar(`A')))^0.5
	scalar `H' = - (`prior_logrr'/scalar(`s'))
	
	replace `depv' = `A' in l
	replace `v' = `=1/`s'' in l
	replace `H' = scalar(`H') in l
	replace `constant' = 0 in l 
	replace  `n' = 2*scalar(`A') in l
	replace `tag' = 1 in l
	
	if "`nolist'" == "" {
		scalar `lb' = `=exp(`prior_logrr')'*(invF(`=2*scalar(`A')', /*
			*/ `=2*scalar(`A')', `=(1-(`level'/100))/2'))^`s'
		scalar `ub' = `=exp(`prior_logrr')'*(invF(`=2*scalar(`A')', /*
			*/ `=2*scalar(`A')', `=((1-(`level'/100))/2)+(`level'/100)'))^`s'
		local normal`c' `: di " exact prior median OR (`level'% PL): " %4.2fc  `=exp(`prior_logrr')*(invF(2*scalar(`A'), /*
			*/ 2*scalar(`A'), 0.5))^`s'' /*
			*/ " (" %4.2f scalar(`lb') ", " %4.2f scalar(`ub') ")"'
			
		get_df_invF `=invF(`=2*scalar(`A')', `=2*scalar(`A')', 0.975)^`s''
		scalar `subjn' = r(df)/2

		local normalz`c' `: di "cases=" %-5.2f scalar(`subjn') " noncases=" %-5.2f scalar(`subjn') " exp(offset)=" /*
		*/ %-5.0g `=exp(scalar(`H'))' '
	
	}
		
	local `c++'
}
* 

* Generalized Log-f priors
if `lfpriors' != 0 {
	mat `retlogf' = J(4,`lfpriors',.)
	mat rownames `retlogf' = m df1 df2 s
	mat colnames `retlogf' = `lfpriornames'
}
local c = 1
qui foreach v of local lfpriornames {
	qui set obs `=_N+1'
	local mmx : word  `c' of `mm'
	local d1x : word  `c' of `d1'
	local d0x : word  `c' of `d0'
	local ssx : word  `c' of `ss'
	
	mat `retlogf'[1,`c'] = `mmx'
	mat `retlogf'[2,`c'] = `d1x'
	mat `retlogf'[3,`c'] = `d0x'
	mat `retlogf'[4,`c'] = `ssx'


	scalar `H' = ln(`d1x' / `d0x') - (`mmx' / `ssx')
	
	replace `depv' = `=(`d1x'/2)' in l
	replace `v' = `=1/`ssx'' in l
	replace `H' = scalar(`H') in l
	replace `constant' = 0 in l 
	replace  `n' = `=(`d1x'+`d0x')/2' in l
	replace `tag' = 1 in l
	
	if "`nolist'" == "" {
		scalar `lb' = `=exp(`mmx')'*(invF(`d1x', /*
			*/ `d0x', `=(1-(`level'/100))/2'))^`ssx'
		scalar `ub' = `=exp(`mmx')'*(invF(`d1x', /*
			*/ `d0x', `=((1-(`level'/100))/2)+(`level'/100)'))^`ssx'
		local logf`c' `: di " exact prior median OR (`level'% PL): " %4.2fc  `=exp(`mmx')*(invF(`d1x', /*
			*/ `d0x', 0.5))^`ssx'' /*
			*/ " (" %4.2f scalar(`lb') ", " %4.2f scalar(`ub') ")"'
	
		scalar `lfc' = `d1x'/(2 * `ssx'^2)
		scalar `lfn' = `d0x'/(2 * `ssx'^2)
		scalar `lfexpoff' = exp(scalar(`H'))
		local logfz`c' `: di "cases=" %-5.2f scalar(`lfc') " noncases=" %-5.2f scalar(`lfn') " exp(offset)=" /*
		*/ %-5.0g `=exp(scalar(`H'))' '
	
	}
	
	local `c++'
}
*

qui foreach v of local indepv {
	replace `v' = 0 if `tag' == 1 & `v' == .
}
qui replace `wgt' = 1 if `tag'


local c = 1
if "`ppl'" != "" {
	local bign = c(N)
				
	foreach x of local ppl {
		qui pllf2  glm `depv' `indepv' `=cond("`noconstant'"=="noconstant"," ","`constant'")' /*
			*/ if `touse' [fw = `wgt'], /*
			*/ fam(binomial `n') offzet(`H') nocons /*
			*/ profile(`x')  nodots nograph level(`level') n(`nppl')
		tempname lb_`x' ub_`x'
		scalar `lb_`x'' =  r(l_llci)
		scalar `ub_`x'' =  r(l_ulci)
		local `c++'
	}
	qui drop if _n > `bign'
}	

		
tempname coefs VCE ilog


qui glm `depv' `indepv' `=cond("`noconstant'"=="noconstant"," ","`constant'")' if `touse' [fw = `wgt'], /*
	*/  offset(`H') nocons  nolog  fam(binomial `n')
	
	
local ll = e(ll)
local ic = e(ic)
local k  = e(k)
local converged = e(converged)
local nobs = e(N)-`npriors'-`lfpriors'
local nobsda = e(N)
local conams "`indepv' `=cond("`noconstant'"=="noconstant", "", "_cons")'"

mat `coefs' = e(b)
mat `VCE' = e(V)
mat `ilog' = e(ilog)

mat colnames `coefs' = `conams'
mat rownames `VCE' = `conams'
mat colnames `VCE' = `conams'

mat coleq `coefs' = ""
mat roweq `VCE' = ""
mat coleq `VCE' =  "" 
 
ereturn post `coefs' `VCE',  obs(`nobs') depn(`depv')

ereturn matrix ilog = `ilog'
if `npriors' != 0 {
	ereturn matrix nprior = `retnorm'
}
if `lfpriors' != 0 {
	ereturn matrix lfprior = `retlogf'
}
ereturn scalar N = `nobs'
ereturn scalar N_da = `nobsda'
ereturn scalar pll = `ll'
ereturn scalar ic = `ic'
ereturn scalar k = `k'
ereturn scalar converged = `converged'

ereturn local cons = "`noconstant'"
ereturn local m = cond("`binomial'"!="", "`binomial'", "1")
ereturn local varfunc = "glim_v2"
ereturn local link = "glim_l02"
ereturn local opt = "moptimize"
ereturn local predict "glim_p"	
if "`weight'" != "" {
	ereturn local wexp = "`exp'"
	ereturn local wtype = "fweight" 
}
ereturn local indepvars "`indepv'"
ereturn local cmd "penlogit"
ereturn local cmdline `"penlogit `cmdline'"'

if  `=`nprior_elements' + `lfprior_elements'' == 0 {
	local title "Logistic regression"
	local lltitle "Log likelihood"
}
else {
	local title "Penalized logistic regression"
	local lltitle "Penalized log likelihood"
}

di _n in gr "`title'"  _col(55) "No. of obs  = " in ye %10.0g e(N) _n

if "`nolist'" == "" {
	local c = 1
	foreach v of local npriornames {
		di in gr "Normal prior for `=abbrev("`v'", 12)': " in y "`normal`c''"
		di in gr "Data approx. equivalent to prior: " in y "`normalz`c''" _n
		local `c++'
	}

	local c = 1
	foreach v of local lfpriornames {
		di in gr "Log-F prior for `=abbrev("`v'", 21)': " in y "`logf`c''"
		di in gr "Data approx. equivalent to prior: " in y "`logfz`c''" _n
		local `c++'
	}
}


di in gr "`lltitle' = " in ye `ll' _n

Replay, level(`level')  `or'

local c = 1
if "`ppl'" != "" {
	tempname plm
	mat define `plm' = J(`nr_profile', 2, .)
	mat colnames `plm' = lb ub
	mat rownames `plm' = `ppl'
	local ytitle "`depv'"
	if "`or'" == "" local tt "Coef."
	else local tt "Odds Ratio"
	di as txt   "{hline 13}{c TT}{hline 24}"
	di as txt %12s abbrev("`ytitle'",12) _col(14)"{c |}" /*
		*/ %1s  " [`level'% PL Conf. Interval]"
	di as txt "{hline 13}{c +}{hline 24}"
	foreach x of local ppl {
		mat `plm'[`c',1] = `lb_`x''
		mat `plm'[`c',2] = `ub_`x''
		if "`or'" == "" di in gr %12s abbrev("`x'",12) _col(14) "{c |}" /*
			*/ as res _col(18) %9.0g `lb_`x'' _col(30) %9.0g `ub_`x''
		else di in gr %12s abbrev("`x'",12) _col(14) "{c |}" /* 
			*/ as res _col(18) %9.0g exp(`lb_`x'') _col(30) %9.0g exp(`ub_`x'')
		local `c++'
	}
	di as txt "{hline 13}{c BT}{hline 24}"
	ereturn matrix ppl = `plm'
}
	

qui drop if `tag' == 1	

ereturn repost, esample(`touse')	
end

capture program drop Replay
program Replay
	syntax [, Level(cilevel) or ]
	if "`or'"== "" ereturn display, level(`level')
	else ereturn display, level(`level') eform("Odds Ratio")
end





*------- pllf2 -------*

capture program drop pllf2
program define pllf2, rclass sortpreserve
version 9.0
/*
	Currently supported commands include at least the following:

	clogit cnreg cox ereg fit glm gnbreg heckman logistic logit	///
		   mlogit nbreg ologit oprobit poisson probit regress reg3	///
		   streg stcox stpm weibull
*/
gettoken cmd 0 : 0
if "`cmd'"=="stpm" local eqxb [xb]
else if ("`cmd'"=="fit") | ("`cmd'"=="reg") | (substr("`cmd'",1,4)=="regr") local cmd regress

syntax anything [if] [in] [aweight fweight pweight iweight], ///
 [DEViance FORMula(string) gen(string) DIFFerence LEVel(cilevel) ///
 PLaceholder(string) PROfile(string) range(string) ///
 MAXCost(int -1) n(integer 100) noci noDOTs nograph noCONStant gropt(string asis)		///
 LEVLINe(string asis) CILINes(string asis) offzet(string) *]

if `maxcost'<0 local maxcost = int(`n'/2)

if "`placeholder'"!="" {
	if wordcount("`placeholder'")!=1 {
		di as err "invalid placeholder()"
		exit 198
	}
}
else local placeholder X

if "`formula'"!="" {
	if "`profile'"!="" {
		di as txt "[profile() ignored]"
		local profile
	}
	if "`range'"=="" {
		di as err "range() required"
		exit 198
	}
	// Check for `placeholder' in `anything'
	local result: subinstr local anything "`placeholder'" "", count(local nat)
	if `nat'==0 {
		di as err `"`placeholder' not found in regression_cmd_stuff ( `anything' )"'
		exit 198
	}
}
else if "`profile'"!="" {
	local varlist "`anything'"
	// Disentangle profile; extract eq from it, if present
	tokenize `profile', parse("[]")
	if "`5'"!="" {
		di as err "syntax error in profile(`profile'), invalid equation name"
		exit 198
	}
	if "`2'"!="" {
		// Rudimentary check that user has entered the eq correctly
		if "`1'"!="[" | "`3'"!="]" {
			di as err "syntax error in profile(`profile'), invalid equation name"
			exit 198
		}
		local profile `4'
		local eq [`2']
		constraint free
		local cuse `r(free)'
	}
	if "`profile'"!="_cons"{
		unab profile: `profile'
	}
	else {
		if "`eq'"=="" {
			di as err "profile log likelihood for _cons not supported"
			exit 198
		}
	}
}
else {
	di as err "must supply either formula() or profile()"
	exit 198
}
/*
if "`gen'"=="" {
	local gen1 _beta
	local gen2 _pll
}
else gettoken gen1 gen2: gen
if "`gen2'"=="" local gen2 _pll
*/
if "`range'"!="" {
	gettoken from to: range
	confirm num `from'
	confirm num `to'
	if `from' > `to' {
		local temp `from'
		local from `to'
		local to `temp'
		local temp
	}
}
if "`weight'" != "" local wt [`weight'`exp']
if "`profile'" != "" { // ------------ begin linear profiling --------

	// Fit model and get level% ci. Program terminates if invalid cmd attempted.
	if "`eq'"==""  & "`constant'"!="noconstant" {
		qui _rmcoll `varlist' `profile'	// strips `profile' if already mentioned in `varlist'
		local tempvl `r(varlist)'
	}
	else local tempvl `varlist'
	capture noisily `cmd' `tempvl' `if' `in' `wt', `options' `constant' offset(`offzet')
	local ytitle `e(depvar)'
	quietly {
		// Check that alleged parameter exists. One or both of `eqxb' or `eq' will always be null.
		capture local b0 = `eq'`eqxb'_b[`profile']
		if "`b0'"=="" local b0 .
		if "`cmd'"=="regress" local z = invttail(e(df_r), (100-`level')/200)
		else local z = -invnorm((100-`level')/200)
		local nobs = e(N)
		local ll0  = e(ll)
		capture local se = `eq'`eqxb'_se[`profile']
		if _rc==0 {
			local llci = `b0'-`z'*`se'
			local ulci = `b0'+`z'*`se'
		}
		else {
			local se .
			if "`range'"=="" {
				noi di as err "could not select range for `profile' (could not estimate MLE). try supplying range()"
				exit 198
			}
		}
		if "`range'"=="" {
			// previously used default range as Wald-based `level' CI; now using +/-(z*1.2)SE for this.
			local from = `b0'-`z'*1.2*`se'
			local to =   `b0'+`z'*1.2*`se'
		}

		// If no equation specified, unabbreviate varlist and remove `profile' from it
		if "`eq'"=="" {
			unab varlist: `varlist'
			local varlist: list varlist - profile
		}

		if "`cmd'"=="regress" {
			constraint free
			local cuse `r(free)'

			// For regress, need to identify yvar; otherwise, not required
			gettoken yvar varlist: varlist
		}
		if `n'>_N {
			di as txt "[Increasing the dataset size to `n'] " _cont
			set obs `n'
		}
		tempvar X Y order // values of X = regression coefficient, Y = profile likelihood for `profile'
		gen `X' = .
		gen `Y' = .
		gen long `order' = _n
		local stepsize = (`to'-`from')/(`n'-1)
		forvalues i=1/`n' {
			local b = `from'+(`i'-1)*`stepsize'
			if "`eq'"!="" {
				// Use constrained regression
				constraint define `cuse' `eq'`profile'=`b'
				`cmd' `varlist' `if' `in' `wt', `options' constraint(`cuse') `constant'
			}
			else {
				if `i'==1 {
					tempvar offset
					gen `offset' = .
				}
				if "`cmd'"=="regress" {
					replace `offset' = `yvar'-`b'*`profile'
					regress `offset' `varlist' `if' `in' `wt', `options'
				}
				else {
					replace `offset' = `b'*`profile'+`offzet'
					`cmd' `varlist' `if' `in' `wt', `options' offset(`offset') `constant'
				}
			}
			sort `order'
			replace `X' = `b' in `i'
			replace `Y' = e(ll) in `i'
			if "`dots'"!="nodots" noi di "." _c
		}
		local cost 0	// "cost": number of extra evaluations of pll needed to find likelihood based CI
		local left_limit .
		local right_limit .
		if "`ci'"!="noci" {
/*
			Search for likelihood based CI bounds: VERY crude!
			Left limit first.
*/
			local target = `ll0'-`z'^2/2	// pll value for computing likelihood based CI
			if `Y'[1]<`target' {
/*
			First evaluated pll is below target pll on left of mle.
			Bracket target from already known values of pll and interpolate.
*/
				forvalues i=2/`n' {
					if `Y'[`i']>=`target' {
						local left_limit = `X'[`i'-1]+`stepsize'*(`target'-`Y'[`i'-1])/(`Y'[`i']-`Y'[`i'-1])
						continue, break	// exit forvalues loop
					}
				}
			}
			else {
/*
			Search for left ll-based confidence limit to the left of first value of b.
			Requires new evaluations of ll - no more than `maxcost' allowed.
*/
				local Yold = `Y'[1]
				local bold `from'
				local i 1
				while missing(`left_limit') & `i'<=`maxcost' {
					local b = `from'-`i'*`stepsize'
					// evaluate pll
					if "`eq'"!="" {
						// Use constrained regression
						constraint define `cuse' `eq'`profile'=`b'
						`cmd' `varlist' `if' `in' `wt', `options' constraint(`cuse') `constant'
					}
					else {
						if "`cmd'"=="regress" {
							replace `offset' = `yvar'-`b'*`profile'
							regress `offset' `varlist' `if' `in' `wt', `options'
						}
						else {
							replace `offset' = `b'*`profile'+`offzet'
							`cmd' `varlist' `if' `in' `wt', `options' offset(`offset') `constant'
						}
					}
					local Ynew = e(ll)
					if `Ynew'<`target' {	// now bracketing target
						local left_limit = `bold'-`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `b'
						local ++i
					}
					if "`dots'"!="nodots" noi di "." _c
				}
				local cost `i'
				if missing(`left_limit') noi di as txt _n "[note: failed to find left-hand confidence limit]"
			}
/*
		Now right limit
*/
			if `Y'[`n']>`target' {	// search for right_limit to the right of last value of b
				local Yold = `Y'[`n']
				local bold `to'
				local i 1
				while missing(`right_limit') & `i'<=`maxcost' {
					local b = `to'+`i'*`stepsize'
					// evaluate pll
					if "`eq'"!="" {
						// Use constrained regression
						constraint define `cuse' `eq'`profile'=`b'
						`cmd' `varlist' `if' `in' `wt', `options' constraint(`cuse') `constant'
					}
					else {
						if "`cmd'"=="regress" {
							replace `offset' = `yvar'-`b'*`profile'
							regress `offset' `varlist' `if' `in' `wt', `options'
						}
						else {
							replace `offset' = `b'*`profile'+`offzet'
							`cmd' `varlist' `if' `in' `wt', `options' offset(`offset') `constant'
						}
					}
					local Ynew = e(ll)
					if `Ynew'<`target' {	// now bracketing target
						local right_limit = `bold'+`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `b'
						local ++i
					}
					if "`dots'"!="nodots" noi di "." _c
				}
				local cost = `cost'+`i'
			}
			else {	// pll is below target on right of mle, bracket target and interpolate
				local n1 = `n'-1
				forvalues i=`n1'(-1)1 {
					if `Y'[`i']>=`target' {
						local right_limit = `X'[`i']+`stepsize'*(`target'-`Y'[`i'])/(`Y'[`i'+1]-`Y'[`i'])
						continue, break	// exit forvalues loop
					}
				}
			}
			if missing(`right_limit') noi di as txt _n "[note: failed to find right-hand confidence limit]"
		}
		lab var `X' "`eq'_b[`profile']"
	}
}
else {	// --------------- begin non-linear profiling ---------------
	tempvar xx
	qui gen `xx' = .
	// Trial fit of model with central value of param. Program terminates if invalid cmd attempted.
	local A = (`to'-`from')/2
	parsat `"`anything'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder')
	cap `cmd' `r(result)' `if' `in' `wt', `options' `constant'
	local rc = _rc
	if `rc' error `rc'

	local ytitle `e(depvar)'
	cap local ll = e(ll)
	if _rc | ("`ll'"==".") {
		di as err "valid log likelihood not returned in e(ll)"
		exit 198
	}

	quietly {
		if "`cmd'"=="regress" local z = invttail(e(df_r), (100-`level')/200)
		else local z = -invnorm((100-`level')/200)

		local nobs = e(N)
		replace `touse' = e(sample)

		if `n'>_N {
			di as txt "[Increasing the dataset size to `n'] " _cont
			set obs `n'
		}
		tempvar X Y order // values of X = regression coefficient, Y = profile likelihood for `profile'
		gen `X' = .
		gen `Y' = .
		gen long `order' = _n
		local stepsize = (`to'-`from')/(`n'-1)
		local y1 .
		local y2 .
		local y3 .
		local b0 .
		local ll0 .
		local done 0
		forvalues i=1/`n' {
			local A = `from'+(`i'-1)*`stepsize'
			// Substitute A in `formula' and fit model
			parsat `"`anything'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder')
			cap `cmd' `r(result)' `if' `in' `wt', `options' `constant'
			if _rc {
				noi di as err "model fit failed at parameter = " `A'
				exit 198
			}
			local y3 = e(ll)
			if !missing(`y1') & !missing(`y2') {
				if `y1'<`y2' & `y2'>`y3' {
					// Solve quadratic y = a + b*x + c*x^2 through 3 points to get MLE b0 and ll at MLE ll0
					local c = (2*(`y3'-`y2')-(`y3'-`y1'))/(2*`stepsize'^2)
					local b = (`y3'-`y2')/`stepsize'-`c'*(2*`A'-`stepsize')
					local a = `y3'-`A'*(`b'+`c'*`A')
					local b0 = -`b'/(2*`c')
					local ll0 = `a'+`b0'*(`b'+`c'*`b0')
					local done 1
				}
			}
			local y1 `y2'
			local y2 `y3'
			sort `order'
			replace `X' = `A' in `i'
			replace `Y' = e(ll) in `i'
			if "`dots'"!="nodots" noi di "." _c
		}
		if !`done' {
			// MLE not straddled. Attempt quadratic solution using terminals of range and a midpoint. Maintain equal spacing.
			if mod(`n',2)==0 {	// even number of evaluations
				local mid = `n'/2
				local last `n'-1
			}
			else {	// odd number of evaluations
				local mid = (`n'+1)/2
				local last `n'
			}
			local A = `X'[`last']
			local s = `X'[`mid']-`X'[1]	// stepsize for this exercise
			local y1 = `Y'[1]
			local y2 = `Y'[`mid']
			local y3 = `Y'[`last']
			local c = (2*(`y3'-`y2')-(`y3'-`y1'))/(2*`s'^2)
			local b = (`y3'-`y2')/`s'-`c'*(2*`A'-`s')
			local a = `y3'-`A'*(`b'+`c'*`A')
			local b0 = -`b'/(2*`c')
			local ll0 = `a'+`b0'*(`b'+`c'*`b0')
			noi di as txt _n "Note: range does not include MLE of parameter - estimate may be inaccurate."
		}
		local cost 0	// "cost": number of extra evaluations of pll needed to find likelihood based CI
		local left_limit .
		local right_limit .
		if "`ci'"!="noci" {
/*
		Search for likelihood based CI bounds: VERY crude!
		Left limit first.
*/
			local target = `ll0'-`z'^2/2	// pll value for computing likelihood based CI
			if `Y'[1]<`target' {
/*
			First evaluated pll is below target pll on left of mle.
			Bracket target from already known values of pll and interpolate.
*/
				forvalues i=2/`n' {
					if `Y'[`i']>=`target' {
						local left_limit = `X'[`i'-1]+`stepsize'*(`target'-`Y'[`i'-1])/(`Y'[`i']-`Y'[`i'-1])
						continue, break	// exit forvalues loop
					}
				}
			}
			else {
/*
			Search for left ll-based confidence limit to the left of first value of b.
			Requires new evaluations of ll - no more than `maxcost' allowed.
*/
				local Yold = `Y'[1]
				local bold `from'
				local i 1
				while missing(`left_limit') & `i'<=`maxcost' {
					local A = `from'-`i'*`stepsize'
					// evaluate pll
					parsat `"`anything'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder')
					cap `cmd' `r(result)' `if' `in' `wt', `options' `constant'
					if _rc {
						noi di as err "model fit failed at parameter = " `A'
						exit 198
					}
					local Ynew = e(ll)
					if `Ynew'<`target' {	// now bracketing target
						local left_limit = `bold'-`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `A'
						local ++i
					}
					if "`dots'"!="nodots" noi di "." _c
				}
				local cost `i'
				if missing(`left_limit') noi di as txt _n "[note: failed to find left-hand confidence limit]"
			}
/*
		Now right limit
*/
			if `Y'[`n']>`target' {	// search for right_limit to the right of last value of b
				local Yold = `Y'[`n']
				local bold `to'
				local i 1
				while missing(`right_limit') & `i'<=`maxcost' {
					local A = `to'+`i'*`stepsize'
					// evaluate pll
					parsat `"`anything'"' `if' `in', formula(`formula') var(`xx') value(`A') placeholder(`placeholder')
					cap `cmd' `r(result)' `if' `in' `wt', `options' `constant'
					if _rc {
						noi di as err "model fit failed at parameter = " `A'
						exit 198
					}
					local Ynew = e(ll)
					if `Ynew'<`target' {	// now bracketing target
						local right_limit = `bold'+`stepsize'*(`target'-`Yold')/(`Ynew'-`Yold')
					}
					else {
						local Yold `Ynew'
						local bold `A'
						local ++i
					}
					if "`dots'"!="nodots" noi di "." _c
				}
				local cost = `cost'+`i'
			}
			else {	// pll is below target on right of mle, bracket target and interpolate
				local n1 = `n'-1
				forvalues i=`n1'(-1)1 {
					if `Y'[`i']>=`target' {
						local right_limit = `X'[`i']+`stepsize'*(`target'-`Y'[`i'])/(`Y'[`i'+1]-`Y'[`i'])
						continue, break	// exit forvalues loop
					}
				}
			}
			if missing(`right_limit') noi di as txt _n "[note: failed to find right-hand confidence limit]"
		}
		local se .
		local llci .
		local ulci .
		lab var `X' "`placeholder' in `formula'"
	}
}
// Pseudo-SE
*local pse = (`right_limit'-`left_limit')/(2*`z')
*capture drop `gen1'
*capture drop `gen2'
*rename `X' `gen1'
*rename `Y' `gen2'
*lab var `Y' "profile log likelihood function"
local ll_limit = `ll0'-`z'^2/2
local limit `ll_limit'
if "`difference'"!="" {
	// compute difference, subtract ll0
	qui replace `Y' = `Y'-`ll0'
	*lab var `gen2' "profile log likelihood difference function"
	local limit = -`z'^2/2
}
if "`deviance'"!="" {
	qui replace `Y' = -2*`Y'
	if "`difference'"!="" lab var `Y' "profile deviance difference function"
	else lab var `Y' "profile deviance function"
	local limit = -2*`limit'
}


/*
local asym = 100*((`right_limit'-`b0')-(`b0'-`left_limit'))/(`right_limit'-`left_limit')
if "`graph'"!="nograph" {
	local Asym: display %4.1f `asym'

	// Extract title if present - default involves asymmetry
	_get_gropts , graphopts(`gropt') getallowed(title)
	local gropt `s(graphopts)'
	if `"`s(title)'"'=="" local title title("Asymmetry = `Asym'%")
	else if `"`s(title)'"'==`""""' local title
	else local title title(`s(title)')

	if !missing(`left_limit')  local lll `left_limit'
	if !missing(`right_limit') local rrr `right_limit'
	if ("`lll'`rrr'"!="")					///
		local xl xline(`lll' `rrr', lstyle(ci) `cilines')
	graph twoway line `gen2' `gen1', `gropt' `title' ///
	    `xl' yline(`limit', lstyle(refline) `levline')
}
if "`dots'"!="nodots" di
local tt "Coef."
di as txt _n "{hline 13}{c TT}{hline 47}"
di as txt %12s abbrev("`ytitle'",12) _col(14)"{c |}" ///
  %10s "Coef." "   Std. Err.     [`level'% PLL Conf. Int.]"
di as txt "{hline 13}{c +}{hline 47}"
if "`eq'"!="" di as res %-12s abbrev("`eq'",12) _col(14) as txt "{c |}"
if "`formula'"!="" di as txt %12s "`placeholder'" _cont
else di as txt %12s abbrev("`profile'",12) _cont
di _col(14) "{c |}" as res ///
  _col(16) %9.0g `b0' ///
  _col(28) %9.0g `pse'  ///
  _col(41) %9.0g `left_limit' ///
  _col(53) %9.0g `right_limit'
di as txt "{hline 13}{c BT}{hline 47}"
di as txt "Note: Std. Err. is pseudo standard error, derived from PLL CI"
*/
// MLE of beta
return scalar b = `b0'
return scalar se = `se'
*return scalar pse = `pse'

// Upper confidence limits
return scalar n_ulci = `ulci'
return scalar l_ulci = `right_limit'

// Lower confidence limits
return scalar n_llci = `llci'
return scalar l_llci = `left_limit'

// maximised likelihood
return scalar ll = `ll0'

// target ll for computing likelihood-based confidence limits
return scalar ll_limit = `ll_limit'
return scalar nobs = `nobs'

// Cost of searching for pll-based confidence interval
return scalar cost = `cost'

// Percentage asymmetry
*return scalar asym = `asym'


end

capture program drop parsat
program define parsat, rclass
version 9.0
syntax anything [if], Formula(string) var(varname) VALue(string) Placeholder(string)

// Replace placeholder with var and update var
local result: subinstr local anything "`placeholder'" " `var' ", all

// Substitute `placeholder' with `value' in `formula'
local f = subinstr(`"`formula'"', `"`placeholder'"' , "`value'", .)
quietly replace `var' = `f' `if'
return local result `result'

end
 
*Calculate approximate number of equivalent subjects 
capture program drop get_df_invF
program get_df_invF, rclass
args ub
local p1 = 0.01
local p2 = 1000
while  abs(`p1'-`p2') >  1e-3  {
	local t = (`p1'+`p2')/2
	local diff = `ub' - invF(`t' , `t', .975)    
	if `diff' < 0 local p1 = `t'
	else local p2 = `t'
}
return scalar df = `t'
end


exit

clear
input y x n  
173 1 307
602 0 1265
end

penlogit y x , binomial(n) nprior(x .25 4) pl(x)
