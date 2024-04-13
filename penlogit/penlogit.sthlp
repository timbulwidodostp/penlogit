{smcl}
{* *! version 1.1.0 19jul2017}{...}
{cmd:help penlogit}{right: ({browse "http://www.stata-journal.com/article.html?article=st0400":SJ15-3: st0400})}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{hi:penlogit} {hline 2}}Penalized logistic regression{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 13 2}
{cmd:penlogit} {depvar} [{indepvars}] {ifin} {weight}
[{cmd:,} {it:options}]

{synoptset 20}{...}
{synopthdr}
{synoptline}
{synopt :{opt np:rior()}}impose normal priors{p_end}
{synopt :{opt lfp:rior()}}impose generalized log-F priors{p_end}
{synopt :{opth ppl(varlist)}}request penalized profile-likelihood limits{p_end}
{synopt :{opt nppl(#)}}evaluate penalized profile-likelihood limits at # points; default is {cmd:nppl(100)}{p_end}
{synopt :{opth bin:omial(varname)}}specify the variable containing the binomial denominator{p_end}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt or}}report odds ratios{p_end}
{synopt :{opt nol:ist}}suppress the summary of prior distributions{p_end}
{synopt :{opt nocon:stant}}suppress the constant term{p_end}
{synoptline}
{p2colreset}{...}
{pstd}The full specification for {cmd:nprior()} is {cmd:nprior(}{it:varname m v} [{it:varname m v} ...]{cmd:)} and for {opt lfprior()} is
{cmd:lfprior(}{it:varname m df1 df2 s} [{it:varname m df1 df2 s}
...]{cmd:)}.{p_end}
{pstd}{cmd:by}, {cmd:statsby}, and {cmd:xi} are allowed; see {help prefix}.{p_end}
{pstd}{cmd:fweight}s are allowed; see {help weight}.{p_end}
{pstd}{cmd:penlogit} allows postestimation commands such as
{cmd:testparm}, {cmd:test}, {cmd:lincom}, {cmd:predictnl}, and
{cmd:predict}; see {helpb test}, {helpb lincom}, {helpb predictnl}, and
{helpb glm_postestimation##predict:glm postestimation}.{p_end}


{title:Description}

{pstd}
{cmd:penlogit} estimates penalized logistic regression models for a binary
response via data augmentation.  It allows the user to impose Normal and
generalized log-F prior distributions on one or more model parameters (log
odds-ratios).

{pstd}
By default, no priors are imposed on the model coefficients.  If no priors are
imposed by the user (that is, if neither option {opt nprior()} nor option
{cmd:lfprior()} is used), {cmd:penlogit} reproduces the results obtained using
{helpb logit}.


{title:Options}

{phang}
{cmd:nprior(}{it:varname m v} [{it:varname m v} ...]{cmd:)} imposes a normal
prior with mean=mode=median={it:m} and variance {it:v} on the
desired model parameter (log odds-ratio).

{phang}
{cmd:lfprior(}{it:varname m df1 df2 s} [{it:varname m df1 df2 s} ...]{cmd:)}
imposes a generalized log-F prior distribution with mode {it:m}, degrees of
freedom {it:df1} and {it:df2}, and scale factor {it:s} on the desired model
parameter (log odds-ratio).

{phang}
{opth ppl(varlist)} specifies the variables for which penalized
profile-likelihood limits are required.  It calls an adapted version of the
{cmd:pllf} command (Royston 2007).

{phang}
{opt nppl(#)} evaluates penalized profile-likelihood at {it:#} equally spaced
points.  The default is {cmd:nppl(100)}.

{phang}
{opth binomial(varname)} specifies the variable containing the binomial
denominator when the data are grouped (that is, when {it:depvar} contains the
total number of successes or failures).

{phang}
{opt level(#)}; see {helpb estimation options##level():[R] estimation options}.

{phang}
{opt or} displays the exponentiated coefficients (odds ratios) and
corresponding standard errors and confidence intervals.

{phang}
{opt nolist} suppresses the summary of prior distributions in terms of exact
prior percentiles (50th, 2.5th, and 97.5th) and data approximately equivalent
to priors.

{phang}
{opt noconstant} suppresses the constant term.
 

{title:Examples}

    {title:Penalized logistic regression on grouped data (binomial data)}

{pstd}
Input data on no fetal monitoring and neonatal death (Neutra et al. 1978){p_end}
{phang2}{cmd:. input nomonit deaths n}{p_end}
{phang2}{cmd:0 3 694}{p_end}
{phang2}{cmd:1 14 2298}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}
Estimate penalized maximum likelihood odds-ratio for "no monitoring" status
(Normal prior with {it:m}={cmd:log(2)} and {it:v}={cmd:0.5}){p_end}
{phang2}{cmd:. penlogit deaths nomonit, or binomial(n) nprior(nomonit log(2) 0.5)}{p_end}

{pstd}
Estimate penalized maximum-likelihood odds-ratio for "no monitoring" status
(log-F prior with {it:m}={cmd:log(2)}, {it:df1}={cmd:2000}, {it:df2}={cmd:2},
and {it:s}={cmd:0.5}) and 95% penalized profile-likelihood limits{p_end}
{phang2}{cmd:. penlogit deaths nomonit, or binomial(n) lfprior(nomonit log(2) 2000 2 0.5) ppl(nomonit)}

    {title:Penalized logistic regression on grouped data (frequency-weighted data)}

{pstd}
Input data on no fetal monitoring and neonatal death (Neutra et al. 1978){p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. input nomonit death wgt}{p_end}
{phang2}{cmd:0 0 691}{p_end}
{phang2}{cmd:0 1 3}{p_end}
{phang2}{cmd:1 0 2284}{p_end}
{phang2}{cmd:1 1 14}{p_end}
{phang2}{cmd:end}{p_end}

{pstd}
Estimate penalized maximum-likelihood odds-ratio for "no monitoring" status
(log-F prior with {it:m}={cmd:log(2)}, {it:df1}={cmd:2000}, {it:df2}={cmd:2},
and {it:s}={cmd:0.5}) and 95% penalized profile-likelihood limits{p_end}
{phang2}{cmd:. penlogit death nomonit [fw=wgt], or lfprior(nomonit log(2) 2000 2 0.5) ppl(nomonit)}

    {title:Penalized logistic regression on individual data}

{pstd}
Load the full dataset on neonatal mortality (Neutra et al. 1978){p_end}
{phang2}{cmd:. use http://www.imm.ki.se/biostatistics/data/neutra1978.dta, clear}{p_end}

{pstd}
Estimate penalized maximum-likelihood odds-ratio for "no monitoring" status
and age at delivery{p_end}
{phang2}{cmd:. xi: penlogit death nomonit i.teenages, or lfprior(nomonit log(2) 2000 2 0.5) nprior(_Iteenages_1 log(2) 0.5 _Iteenages_2 log(4) 0.5)}
{p_end}


{title:Stored results}

{pstd}
{cmd:penlogit} stores the following in {cmd:e()}:

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_da)}}number of observations including the augmented data{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}
{synopt:{cmd:e(pll)}}penalized log-likelihood{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:penlogit}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of independent variables{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(link)}}program to calculate link function{p_end}
{synopt:{cmd:e(varfunc)}}program to calculate variance function{p_end}
{synopt:{cmd:e(m)}}number of binomial trials{p_end}
{synopt:{cmd:e(cons)}}set if {cmd:noconstant} specified{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(ppl)}}penalized profile-likelihood limits{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(nprior)}}{it:m} and {it:v} of the normal priors{p_end}
{synopt:{cmd:e(lfprior)}}{it:m}, {it:df1}, {it:df2}, and {it:s} of the log-F priors{p_end}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{title:References}

{phang}Neutra, R. R., S. E. Fienberg, S. Greenland, and E. A. Friedman. 1978.
Effect of fetal monitoring on neonatal death rates.
{it:New England Journal of Medicine} 299: 324-326.

{phang}Royston, P. 2007.
{browse "http://www.stata-journal.com/article.html?article=st0132":Profile likelihood for estimation and confidence intervals}.
{it:Stata Journal} 7: 376-387.


{title:Authors}

{pstd}Andrea Discacciati{p_end}
{pstd}{browse "http://www.imm.ki.se/biostatistics/":Unit of Biostatistics} and {browse "http://ki.se/imm/nutrition-en":Unit of Nutritional Epidemiology}{p_end}
{pstd}{browse "http://ki.se/imm":Institute of Environmental Medicine}{p_end}
{pstd}Karolinska Institutet{p_end}
{pstd}Stockholm, Sweden{p_end}
{pstd}andrea.discacciati@ki.se

{pstd}Nicola Orsini{p_end}
{pstd}{browse "http://www.imm.ki.se/biostatistics/":Unit of Biostatistics} and {browse "http://ki.se/imm/nutrition-en":Unit of Nutritional Epidemiology}{p_end}
{pstd}{browse "http://ki.se/imm":Institute of Environmental Medicine}{p_end}
{pstd}Karolinska Institutet{p_end}
{pstd}Stockholm, Sweden{p_end}
{pstd}nicola.orsini@ki.se


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 15, number 3: {browse "http://www.stata-journal.com/article.html?article=st0400":st0400}{p_end}
