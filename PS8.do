cd "/Users/zhouyeyang/Library/Mobile Documents/com~apple~CloudDocs/U of Arizona/Y1S2/ECON 522A/PS/PS8 3.28"
import excel "Engel.xls", clear
rename A foodexp
rename B income
gen log_foodexp = log(foodexp)
gen log_income = log(income)

* Log-linear model, homoskedasticity-only standard errors
reg log_foodexp log_income

* Log-linear model, heteroskedasticity-robust standard errors
reg log_foodexp log_income, robust

* Nonlinear model, homoskedasticity-only standard errors
nl (log_foodexp=exp({a=0})*log_income^{b=1})
test /b=1
* Nonlinear model, heteroskedasticity-robust standard errors
nl (log_foodexp=exp({a=0})*log_income^{b=1}), robust
test /b=1


/* 
Nonlinear Regression 
https://www.stata.com/features/overview/nonlinear-regression/
command: nl
*/