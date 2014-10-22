using MixedModels, SelectMixedModels, RDatasets

ds = dataset("lme4", "Dyestuff");
fm1 = fit(lmm(Yield ~ 1|Batch, ds))
@time conAIC(fm1)
marAIC(fm1)
marBIC(fm1)

sleep = dataset("lme4", "sleepstudy");
fm2 = fit(lmm(Reaction ~ Days + (Days | Subject), sleep))
@time conAIC(fm2)
marAIC(fm2)
marBIC(fm2)

inst = dataset("lme4", "InstEval");
fm3 = fit(lmm(Y ~ Dept*Service + (1|S) , inst))
@time conAIC(fm3)
marAIC(fm3)
marBIC(fm3)

chem = dataset("mlmRev", "Chem97");
fm4 = fit(lmm(Score ~  (1|Lea), chem))
@time conAIC(fm4)
marAIC(fm4)
marBIC(fm4)

pen = dataset("lme4", "Penicillin");
fm5 = fit(lmm(Diameter ~ (1|Plate) + (1|Sample), pen))
@time conAIC(fm5)
marAIC(fm5)
marBIC(fm5)

psts = dataset("lme4", "Pastes");
fm6 = fit(lmm(Strength ~ (1|Sample) + (1|Batch), psts))
@time conAIC(fm6)
marAIC(fm6)
marBIC(fm6)

fm7 = fit(lmm(Score ~  (1|School) + (1|Lea), chem))
@time conAIC(fm7)
marAIC(fm7)
marBIC(fm7)

fm8 = fit(lmm(Y ~ Dept*Service + (1|S) + (1|D), inst))
@time conAIC(fm8)
marAIC(fm8)
marBIC(fm8)
