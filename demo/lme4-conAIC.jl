using MixedModels, SelectMixedModels, RDatasets

ds = dataset("lme4", "Dyestuff");
fm1 = fit(lmmg(Yield ~ 1|Batch, ds))
@time conAIC(fm1)

sleep = dataset("lme4", "sleepstudy");
fm2 = fit(lmmg(Reaction ~ Days + (Days | Subject), sleep))
@time conAIC(fm2)

inst = dataset("lme4", "InstEval");
fm3 = fit(lmmg(Y ~ Dept*Service + (1|S) , inst))
@time conAIC(fm3)

chem = dataset("mlmRev", "Chem97");
fm4 = fit(lmmg(Score ~  (1|Lea), chem))
@time conAIC(fm4)

pen = dataset("lme4", "Penicillin");
fm5 = fit(lmmg(Diameter ~ (1|Plate) + (1|Sample), pen))
@time conAIC(fm5)

psts = dataset("lme4", "Pastes");
fm6 = fit(lmmg(Strength ~ (1|Sample) + (1|Batch), psts))
@time conAIC(fm6)

fm7 = fit(lmmg(Score ~  (1|School) + (1|Lea), chem))
@time conAIC(fm7)

fm8 = fit(lmmg(Y ~ Dept*Service + (1|S) + (1|D), inst))
@time conAIC(fm8)
