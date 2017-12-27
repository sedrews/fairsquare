from z3 import *

#currently I'm only using the real-valued attributes for ease, and the models were generated accordingly

age=Const("age", RealSort())
fnlwgt=Const("fnlwgt", RealSort())
education_num=Const("education-num", RealSort())
capital_gain=Const("capital-gain", RealSort())
capital_loss=Const("capital-loss", RealSort())
hours_per_week=Const("hours-per-week", RealSort())

"""
SVM
        -0.463  * (normalized) age
 +      -0.1386 * (normalized) fnlwgt
 +      -0.8136 * (normalized) education-num
 +     -21.1418 * (normalized) capital-gain
 +      -3.2893 * (normalized) capital-loss
 +      -0.6579 * (normalized) hours-per-week
 +       1.9875

Normalization functions:
age	N(x) = (x-17.0)/73.0
fnlwgt	N(x) = (x-12285.0)/1472420.0
education-num	N(x) = (x-1.0)/15.0
capital-gain	N(x) = (x-0.0)/99999.0
capital-loss	N(x) = (x-0.0)/4356.0
hours-per-week	N(x) = (x-1.0)/98.0
"""
#svm sat iff <=50k
svm=((-0.463   * (age-17)/73) + \
     (-0.1386  * (fnlwgt-12285)/1472420) + \
     (-0.8136  * (education_num-1)/15) + \
     (-21.1418 * (capital_gain)/99999) + \
     (-3.2893  * (capital_loss)/4356) + \
     (-0.6579  * (hours_per_week-1)/98) + \
       1.9875 ) > 0

"""
DT (depth was set to have a maximum of 4)
capital-gain < 7073.5
|   age < 28.5 : <=50K 
|   age >= 28.5
|   |   education-num < 12.5 : <=50K 
|   |   education-num >= 12.5
|   |   |   capital-loss < 1881.5 : <=50K 
|   |   |   capital-loss >= 1881.5 : >50K 
capital-gain >= 7073.5
|   age < 20 : <=50K 
|   age >= 20 : >50K 
"""
#dt sat iff <=50k
dt=Or( And( capital_gain < 7073.5, \
            Or( age < 28.5, \
                And( Not(age < 28.5), \
                     Or(education_num < 12.5, \
                        And(Not(education_num < 12.5), \
                            capital_loss < 1881.5))))), \
       And(Not(capital_gain < 7073.5), age < 20))

print("svm examples")
s=Solver()
print("<=50k example")
print("age==59, fnlwgt==148844, education_num==9, capital_gain==0, capital_loss==0, hours_per_week==40")
s.add(svm, age==59, fnlwgt==148844, education_num==9, capital_gain==0, capital_loss==0, hours_per_week==40)
print(s.check()) 
print("sat -> <=50k")
s.reset()
print(">50k example")
print("age==54, fnlwgt==166459, education_num==15, capital_gain==99999, capital_loss==0, hours_per_week==60")
s.add(svm, age==54, fnlwgt==166459, education_num==15, capital_gain==99999, capital_loss==0, hours_per_week==60)
print(s.check())
print("unsat -> >50k")

print("dt examples")
s=Solver()
print("<=50k example")
print("age==59, fnlwgt==148844, education_num==9, capital_gain==0, capital_loss==0, hours_per_week==40")
s.add(dt, age==59, fnlwgt==148844, education_num==9, capital_gain==0, capital_loss==0, hours_per_week==40)
print(s.check()) 
print("sat -> <=50k")
s.reset()
print(">50k example")
print("age==54, fnlwgt==166459, education_num==15, capital_gain==99999, capital_loss==0, hours_per_week==60")
s.add(dt, age==54, fnlwgt==166459, education_num==15, capital_gain==99999, capital_loss==0, hours_per_week==60)
print(s.check())
print("unsat -> >50k")
