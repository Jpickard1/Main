infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
VectorField:= [
    +x2*x3+x4*x5+x2*x5
    +x1*x3+x3*x4+x1*x5
    +x1*x2+x2*x4+x4*x5
    +x2*x3+x3*x5+x1*x5
    +x3*x4+x1*x4+x1*x2
]:
OutputSystem := [
    x1
    x2
x1
] :
OutputsVariables:= [y] :
Inputs := [] :
Parameters 	:= [] :
Variables 	:= [x1,x2,x3,x4,x5] :
libname := cat(currentdir(),"/release"),libname :
readlib(observabilityTest) :
NonObservable := observabilityTest(	VectorField	,
					Variables	,
					OutputSystem	,
					Parameters	,
					Inputs			) :
print(GroupInfGen := observabilitySymmetries(	VectorField	,
					Variables	,
					OutputSystem	,
					Parameters	,
					Inputs		,
					NonObservable		) :
print(quit :
