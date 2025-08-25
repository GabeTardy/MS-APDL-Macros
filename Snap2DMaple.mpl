
# Snap2DMaple.mw
# - by Gabriel Tardy 
# - using some code from Dr. Y. Jane Liu
# - written and compiled July 2025
restart:

# Requirements
with(plots):
with(plottools):
with(ColorTools):
with(LinearAlgebra):
# Inputs

(* --- Required nondimensionalized parameters for Maple --- *)
n := 4; # Nondimensionalized Maximum Load
tau__0 := 80; # Nondimensionalized Timestep Length
m := 0.35; # Ratio of Arch Rise to Radius of Gyration
zeta := 1/5; # Damping Ratio (*2)

(* --- Required nondimensionalized parameters for APDL --- *)
lambda := 1/50; # Ratio of a__0 to L
k := 2; # Ratio of h to b

(* --- Required dimensional parameters for APDL --- *)
a__0 := 1; # Rise of the arch
E := 10000e3; # Young's Modulus
rho := 0.0179; # Density
nu := 0; # Poisson's Ratio

(* --- Required parameters for Bimodulus Solution --- *)
isBimodulus := 1;
# 0 - Unimodulus
# 1 - Bimodulus (Liu and Peddieson Step Function)
# 2 - Bimodulus (Kuznetsov Arctangent Function)
# 3 - Reserved for future implementations

# For bimodulus solutions, you must specify any two of the following four constants:
E__1; # or E__t
E__2; # or E__c

E__r := 12000e3; # = 4*E__t*E__c/(sqrt(E__t) + sqrt(E__c))^2
eta := 2; # = E__c/E__t

# If you are using the Kuznetsov Arctangent Function, you must also specify:
beta__Kuznetsov := 0.1;

# Finally, you must specify which E will be used for nondimensionalization and calculation of the critical load
E__loading := E__r;
E__nondimensionalization := E__r;

(* --- Required parameters for Maple solution --- *)
dataRes := 0.05; # A solution is generated every (dataRes) units of tau__0. 
paletteName := "Spring";
generateAnimation := false;
plotSmoothness := 5000; # How many steps to show in the explored plot?

(* --- Required for simultaneous Ansys solution --- *)
# Location of AnsysNNN.exe (where NNN is the current version number (last two digits of year, then digit of revision number))
ansysExec := "C:\\Program Files\\ANSYS Inc\\ANSYS Student\\v252\\ansys\\bin\\winx64\\ANSYS252.exe";

(* --- Types of Analyses to Complete --- *)
# Perform static analysis?
ansysStaticAnalysis := 1;

# Perform transient analysis?
ansysTransientAnalysis := 1;

(* --- Compare Plots to Previous Solutions --- *)
# Name of current analysis?
analysisName := "NoSnapBackTest; m = 0.35, eta = 2";

# Directories in the Snap2DResults subdirectory to search for comparison plots
# Compare Ers
#compareTo := ["Er Test; Er = 5000e3, eta = 1.5", "Er Test; Er = 10000e3, eta = 1.5", "Er Test; Er = 20000e3, eta = 1.5"]; 

# Compare etas
#compareTo := ["Constant Er; eta = 0.125", "Constant Er; eta = 0.25", "Constant Er; eta = 0.5", "Constant Er; eta = 0.67", "Constant Er; eta = 1", "Constant Er; eta = 1.5", "Constant Er; eta = 2", "Constant Er; eta = 4", "Constant Er; eta = 8"];

# Compare things NOT snapping through
# compareTo := ["NoSnapThroughTest; m = 1.67, eta = 2", "NoSnapThroughTest; m = 1.33, eta = 1.5", "NoSnapThroughTest; m = 0.75, eta = 0.5", "NoSnapThroughTest; m = 0.67, eta = 0.5", "NoSnapThroughTest; m = 0.6, eta = 0.5", "NoSnapThroughTest; m = 0.83, eta = 0.75"];

# Compare things not snapping BACK
compareTo := ["NoSnapBackTest; m = 0.35, eta = 2", "NoSnapBackTest; m = 0.36, eta = 2", "NoSnapBackTest; m = 0.16, eta = 0.25", "NoSnapBackTest; m = 0.125, eta = 0.25","NoSnapBackTest; m = 0.1, eta = 0.25", "NoSnapBackTest; m = 0.15, eta = 0.25", "NoSnapBackTest; m = 0.185, eta = 0.5", "NoSnapBackTest; m = 0.225, eta = 0.5"];

compareAgainst := "eta";

compareStatic := 1;

compareTransient := 1;
(* --- Generate loading and load functions --- *)

# The loading matrix contains two columns;
# The first column is the nondimensionalized time (*n*tau__0) at which the load occurs
# The second column is the load (*n) itself. 
# For example, [1, 1] suggests that at a time of 1*n*tau__0 there is a load of 1*n. 
# Several sample loading cases are presented below. 

# 1-cycle no pull-back
loading := Matrix([
[1, 1],
[2, 0]
]);

# 2-cycle Pull-back
#loading := Matrix([
#[1, 1],
#[5/2, -1/2],
#[4, 1],
#[11/2, -1/2],
#[6,0]]);

#loading := Matrix([
#[1, 1],
#[2, 0],
#[3, 1],
#[4, 0]]);
(* --- DO NOT EDIT ANYTHING FOUND BELOW! THIS MAY BREAK THE SCRIPT. --- *)
thesisStyle := proc(plotFunc):

return display(plotFunc, thickness=3, font=["Times New Roman", roman, 16], labelfont=["Cambria Math", italic, 16], axis=[thickness=2, tickmarks=[thickness=0],gridlines=[thickness=0, colour=lightgray]], gridlines, legendstyle=[font=["Times New Roman", roman, 16]]);

end proc:

linestyles := ['solid', 'dash', 'dashdot', 'longdash', 'spacedash']: # 'dot', 
linesymbols := ['asterisk', 'box', 'circle', 'cross', 'diagonalcross', 'diamond', 'point', 'solidbox', 'solidcircle', 'soliddiamond']:
lsymb := proc(it):
return linesymbols[it mod 10 + 1];
end proc:
ls := proc(it):
return linestyles[it mod 5 + 1];
end proc:
col := proc(idx):
return cat(paletteName, " ", idx);
end proc:

# Collect current time so that total time taken to run can be analyzed
currentTime := time[real]();
# Generate the load functions (assumes ramped loading)
loadingFunctions := Vector(RowDimension(loading)):

loadingConditions := Vector(RowDimension(loading)):
# Additionally generate loading conditions if n=1, tau__0 = 1 (used for plotting)
nondimLoadingConditions := Vector(RowDimension(loading)):

for i from 1 to RowDimension(loading) do:

# Special case since the first load is always 0 at time = 0, so the first loading function has no constant intercept
if i = 1 then

# Use the point-point form of a line to generate the nondimensionalized loading functions from the load "endpoints"
loadingFunctions[i] := uneval(simplify('n'*(loading[i,2])/((loading[i,1])*'n*tau__0')*('tau')));
nondimLoadingConditions[i] := 'tau' > 0 and 'tau' <= loading[i,1];
loadingConditions[i] := 'tau' > 0 and 'tau' <= loading[i,1]*n*tau__0;

else

# Use the point-point form of a line to generate the nondimensionalized loading functions from the load "endpoints"
loadingFunctions[i] := uneval(simplify('n'*(loading[i,2] - loading[i-1,2])/((loading[i,1] - loading[i-1,1])*'n*tau__0')*('tau'-loading[i,1]*'n*tau__0') + loading[i,2]*'n'));
nondimLoadingConditions[i] := 'tau' > loading[i-1,1] and 'tau' <= loading[i,1];
loadingConditions[i] := 'tau' > n*tau__0*loading[i-1,1] and 'tau' <= n*tau__0*loading[i,1];

fi;


end do;

# Display partially-evaluated loading functions as generated by the above code.
'u' = evalr~(loadingFunctions);

# Convert loading functions to a single piecewise function for use later
uTau := subs('n'=1, 'tau__0'=1, piecewise(seq(seq~(zip(`[]`, nondimLoadingConditions, loadingFunctions))))):
uTauTrue := subs('n'=n, 'tau__0'=tau__0, piecewise(seq(seq~(zip(`[]`, loadingConditions, loadingFunctions))))):
# Plot the load functions

# Generate "nondimensionalized" versions of the loading functions wherein tau__0 and n are both 1 (i.e. plot is multiplied by n and tau__0 where necessary)
nondimensionalizedPlots := subs('tau__0' = 1, 'n' = 1, loadingFunctions):

# Generate sequence of plots using bounds listed in loading (except special case for first line, which has a left bound of 0 always).
nondimPlots := Vector(RowDimension(loading)):
firstNonDimPlot := plot(nondimensionalizedPlots[1], tau=0..loading[1,1], legend=typeset("Load Phase 1: ", u=uneval(uneval(loadingFunctions[1])), "; ", 'tau'=uneval(uneval([0, loading[1,1]*'n*tau__0']))), color = cat(paletteName, " 1")):

# Using the end-first vector notation to assign to only the last entries, then generate the plots with their respective bounds:
nondimPlots := [firstNonDimPlot, seq(plot(nondimensionalizedPlots[i], tau=loading[i-1,1]..loading[i,1], legend=typeset("Load Phase ", (RowDimension(loading)+i+1), ": ", u=uneval(loadingFunctions[i]), "; ", 'tau'=uneval([loading[i-1,1]*'n*tau__0', loading[i,1]*'n*tau__0'])), color = cat(paletteName, " ", (RowDimension(loading)+i+1))), i=(-RowDimension(loading)+1)..-1)]:
loadingFunctionPlots := display(nondimPlots, scaling=constrained, labels=['tau/(n*tau__0)', 'u/n']);
# Extra calculations in order to turn nondimensionalized problem into dimensional problem for Ansys
# Perform bimodulus work
if isBimodulus > 0 then
    bimodEq1 := E__r = 4*E__1*E__2/(sqrt(E__1) + sqrt(E__2))^2;
    bimodEq2 := eta = E__2/E__1;

    # Solve only for missing quantities (for use in analysis later)
    bimodSolu := solve({bimodEq1, bimodEq2}, select~(type, [E__r, E__2, E__1, eta], 'name'));
    assign(seq(bimodSolu));
else
E__2 := E;
E__1 := E;
E__loading := E;
fi;
# Statics
L := a__0/lambda;
h := a__0*sqrt(3*m);
b := h*k;
I__zz := b*(h^3)/12;
A := b*h;
q__0 := n*E__loading*I__zz*a__0*(Pi/L)^4;

# Dynamics
omega := (Pi/L)^2*sqrt(E__loading*I__zz/rho/A);
alpha__d := zeta*omega;
t__0 := tau__0/omega;

# Calculate "true" loading using column operations: the "time" values must all be multiplied by n*t__0, while the loading values must be multiplied by q__0.
loadStep := Vector(RowDimension(loading),i->i):
trueLoading := Copy(loading):
trueLoading[1..RowDimension(loading),1] *= n*t__0:
trueLoading[1..RowDimension(loading),2] *= q__0:

# Force load step id to be first value in matrix
trueLoading := Matrix([loadStep, trueLoading]):

# Also force column numbers (APDL requirement, hack solution to get this working because I am tired of trying to figure out why APDL is being finicky about this)
trueLoading2 := Matrix(RowDimension(trueLoading)+1, 3, [0,1,2]):
trueLoading2[2..-1,1..-1] := trueLoading:
trueLoading := trueLoading2;

# Write the loading file to a .csv for use in APDL
Export("Snap2DLoading.csv", trueLoading, base=worksheetdir):

# Use bitwise operators to merge all of the flags for the analysis options (which analyses to generate, which parameters to export, etc.)
analysisFlags := Bits:-Join([ansysStaticAnalysis, ansysTransientAnalysis]):

# Write the geometry and material data to a .csv for use in APDL.
# First column - item index in APDL (required by APDL)
# Second column - actual data value
geomMaterialData := evalf(Matrix([
[1, L],
[2, h],
[3, b],
[4, a__0],
[5, E__1],
[6, rho],
[7, nu],
[8, omega],
[9, alpha__d],
[10, isBimodulus],
[11, E__2],
[12, beta__Kuznetsov],
[13, q__0],
[14, analysisFlags]
])):

Export("Snap2DGeom.csv", geomMaterialData, base=worksheetdir):

# Get current time so that a solution file can be created
if analysisName = "timestamp" then:
ct := Date():
analysisName := cat(convert(Month(ct),'string'),"-",convert(DayOfMonth(ct),'string'),"-",convert(HourOfDay(ct),'string'),"-",convert(Minute(ct),'string')):
fi:

# Make a directory into which our solution files will go:
if not FileTools:-Exists(cat("Snap2DResults\\",analysisName)) then
mkdir(cat(interface(worksheetdir),"\\Snap2DResults\\",analysisName));
fi;

# Make sure all Maple variables are exported (so that we can read them later as metadata for comparison purposes):
# NOTE: THE ORDER SHOULD NEVER CHANGE. THIS WILL CAUSE OLD DATA TO BE INCOMPATIBLE WITH NEW SYSTEM
MapleMetadata := Transpose(Matrix([convert~('[
n,
tau__0,
m,
zeta,
lambda,
k,
a__0,
E,
rho,
nu,
isBimodulus,
E__1,
E__2,
E__r,
eta,
beta__Kuznetsov,
E__loading,
E__nondimensionalization,
dataRes,
paletteName,
generateAnimation,
plotSmoothness,
ansysExec,
ansysStaticAnalysis,
ansysTransientAnalysis,
analysisName,
compareTo
]', 'string'),[
n,
tau__0,
m,
zeta,
lambda,
k,
a__0,
E,
rho,
nu,
isBimodulus,
E__1,
E__2,
E__r,
eta,
beta__Kuznetsov,
E__loading,
E__nondimensionalization,
dataRes,
paletteName,
generateAnimation,
plotSmoothness,
ansysExec,
ansysStaticAnalysis,
ansysTransientAnalysis,
analysisName,
compareTo
]]));

Export(cat("Snap2DResults\\",analysisName,"\\Metadata.csv"), MapleMetadata, base=worksheetdir):

# Remove the output file so that we make sure to load a new one
if FileTools:-Exists("Snap2DOutput.csv") then
FileTools:-Remove("Snap2DOutput.csv");
printf("Output file removed successfully.");
else

printf("Output file does not exist to be deleted.");

fi;

printf("\n\nMaple will now stop working until APDL is complete! This is normal and APDL may take anywhere from 2-5 minutes to complete.");
# Launch APDL
#threadID := system[launch](ansysExec,"-b","-i","Snap2DMaple.mac","-o","maple.out");

APDLStopcode := ssystem(cat(ansysExec," -b -i Snap2DMaple.mac -o maple.out -j snapLastRun -smp"));
# Process WILL hang maple until APDL is complete. This is because I couldn't figure out how to implement multithreading.
# Wait until return code is [8, ""]; otherwise APDL failed (see APDL Operations Guide)
if APDLStopcode[1] <> 8 then
printf(cat("APDL run resulted in an error. (Task ended with stopcode ",APDLStopcode[1],".) Please visit the file maple.out for more information."));
fi;
# Solve for the exact solution via Maple while APDL works.
govEqs := Vector([seq((diff(r(tau),tau$2)) + zeta*(diff(r(tau),tau)) + (r(tau)^2-(1-m))*r(tau)/m = 1 - loadFunc, loadFunc in loadingFunctions)]);
ics := Vector(RowDimension(govEqs)):
solutionTaus := Vector(RowDimension(govEqs)): # Points at which we want the solution
solutionRs := Vector(RowDimension(govEqs)): # The solution itself
solutionUs := Vector(RowDimension(govEqs)): # The taus transformed to u (loading)
solutionPlots := Vector(RowDimension(govEqs)): # Plots of each solution


# First loading phase is special again since endpoint/ICs is/are load=0 at time=0.
# Specify initial conditions
ics[1] := D(r)(0) = 0, r(0) = 1:

# Generate points at which solution will occur with a frequency of dataRes
solutionTaus[1] := [seq(j, j = 0 .. tau__0*n*loading[1,1], dataRes)]:

# Determine the numerical solution at each point
# The indices at the end select the second element of the solution (the first element is an array showing where the solution is, where tau is, etc.) ----> [2]
# From that second element, the first element is selected (which is just the solution array itself) -----------------------------------------------------> [1]
# From the solution array, every row (from row 1 to the last row, which has the shorthand -1) in the second and third columns are selected. -------------------------> [1..-1,2..3]
solutionRs[1] := dsolve({ics[1], govEqs[1]}, type=numeric, output = Array(solutionTaus[1]), maxfun=0)[2][1][1..-1,2..3]:

# Note that the above array contains both r(tau) and dr/d(tau), both of which are required for the initial conditions of the next load step. 

# Generate a plot of the solution we found
solutionUs[1] := [seq(subs(tau=solutionTaus[1][j], eval(loadingFunctions[1])), j=1..ColumnDimension(solutionTaus[1]))]:
solutionPlots[1] := plot(solutionUs[1], solutionRs[1][1..-1,1], style = line, scaling = unconstrained, legend=typeset("Load Phase 1: ", u=uneval(loadingFunctions[1]), "; ", 'tau'=uneval([0, loading[1,1]*'n*tau[0]'])), color = cat(paletteName, " 1")):

# Now do the same thing for the next loading phases, letting the final conditions of the last solution be the initial conditions of the new solution
for i from 2 to RowDimension(govEqs) do:

# Specify initial conditions:
ics[i] := D(r)(solutionTaus[i-1][-1]) = solutionRs[i-1][-1,2], r(solutionTaus[i-1][-1]) = solutionRs[i-1][-1,1]:

# Generate points at which solution will occur with a frequency of dataRes
solutionTaus[i] := [seq(j, j = tau__0*n*loading[i-1,1] .. tau__0*n*loading[i,1], dataRes)]:

# Determine the numerical solution at each point
# The indices at the end select the second element of the solution (the first element is an array showing where the solution is, where tau is, etc.) ----> [2]
# From that second element, the first element is selected (which is just the solution array itself) -----------------------------------------------------> [1]
# From the solution array, every row (from row 1 to the last row, which has the shorthand -1) in the second and third columns are selected. -------------------------> [1..-1,2..3]
solutionRs[i] := dsolve({ics[i], govEqs[i]}, type=numeric, output = Array(solutionTaus[i]), maxfun=0)[2][1][1..-1,2..3]:

# Note that the above array contains both r(tau) and dr/d(tau), both of which are required for the initial conditions of the next load step. 

# Generate a plot of the solution we found
# Use the loading function with the taus used to generate the loading at each time value
solutionUs[i] := [seq(subs(tau=solutionTaus[i][j], eval(loadingFunctions[i])), j=1..ColumnDimension(solutionTaus[i]))];
solutionPlots[i] := plot(solutionUs[i], solutionRs[i][1..-1,1], style = line, scaling = unconstrained, legend=typeset("Load Phase ", i, ": ", u=uneval(loadingFunctions[i]), "; ", 'tau'=uneval([loading[i-1,1]*'n*tau[0]', loading[i,1]*'n*tau[0]'])), color = cat(paletteName, " ", i));

end do:

# Make a correlation between all solutions for the Explore command:
allTaus := [seq(seq(i), i in solutionTaus)]/(n*tau__0):
allUs := [seq(seq(i), i in solutionUs)]:
allRs := [seq(seq(i[1..-1,1]), i in solutionRs)]:

# Solve the static solution
staticGovEq := r__0^3 - (1 - m)*r__0 - m*(1 - u__0);
staticSolution:=solve(staticGovEq,u__0);
staticPlot := plot([staticSolution, r__0, r__0 = 1.1...-1.1], color=gray, style = point, thickness = 3, scaling = unconstrained, legend="Static Solution"):

# Plot all solutions together
allSolutions := display(staticPlot, seq(solutionPlots), scaling=constrained, labels=['u','r']):

# By default "allSolutions" does not actually evaluate the expressions in the legend and they are two-levels-unevaluated. In order to evaluate the terms (to prevent ugliness, repeat the last statement to evaluate everything once:
allSolutions;
tau__0__backup := tau__0:
n__backup := n:
L__backup := L:
tau__0 := 'tau__0':
n := 'n':
L := 'L':

dataIndex := 'dataIndex':
if generateAnimation then
Explore(display(<display(loadingFunctionPlots, pointplot([allTaus[dataIndex], eval(uTau, {tau=allTaus[dataIndex]})], color=black, symbol=solidcircle))|
display(allSolutions, pointplot([allUs[dataIndex], allRs[dataIndex]], color=black, symbol=solidcircle), scaling=constrained, labels=['u','r'], title=typeset('tau', " = ", allTaus[dataIndex], " ",'n*tau[0]'))|
display(plot(sin(Pi*x), x=0..1, color=black, linestyle=dash, legend="Undeformed Configuration"), plot(allRs[dataIndex]*sin(Pi*x), x=0..1, color=black, legend=typeset("Deformed Shape at ",'tau'," = ",allTaus[dataIndex]," ",'n*tau[0]')), labels=['x/L','r*sin(Pi*x/L)'], view=[0..1, -1..1])>), parameters=[dataIndex=[seq(floor(i), i=1..nops(allUs), nops(allUs)/plotSmoothness)]], animate=true, loop=true, size=[1500, 500]);
fi;
# Ansys component
# This code was previously required because the process was multithreaded. Unfortunately, the ssystem command (used to make sure that APDL actually had an output and not an error) forces Maple to hang until the APDL subprocess is complete, so this code is no longer necessary.
# do:
# until FileTools:-Exists("Snap2DOutputTransient.csv") or APDLStopcode[1] <> 8;
# Substitute constants back
tau__0 := tau__0__backup:
n := n__backup:
L := L__backup:

# Switch E to the one used for nondimensionalization
if isBimodulus > 0 then:
    E := E__nondimensionalization;
fi:

AnsysPlots := []:
currentIndex := 1;

# Process transient solution if it exists;
if ansysTransientAnalysis = 1 then:
    AnsysDataTrans := Import("Snap2DOutputTransient.csv",output=Array):

    # Multiply all values by omega to get taus
    AnsysDataTrans[1..-1,1] *= omega:

    # Use loading function to get u from tau
    for i from 1 to RowDimension(AnsysDataTrans) by 1 do:
        AnsysDataTrans[i,1] := eval(uTauTrue, tau=AnsysDataTrans[i,1]);
        AnsysDataTrans[i,2] := 1 - AnsysDataTrans[i,2]/a__0;
    end do:

    # Add transient plot to list
    AnsysPlots := [seq(AnsysPlots), listplot(AnsysDataTrans, labels=['u', 'r'], linestyle=dashdot, legend="Ansys Transient Results")]:
fi:

# Process static solution if it exists;
if ansysStaticAnalysis = 1 then:
    AnsysDataStatic := Import("Snap2DOutputStatic.csv",output=Array):
    AnsysDataStatic[1..-1,1] *= q__0/(E*I__zz*(Pi/L)^4*a__0):

    for i from 1 to RowDimension(AnsysDataStatic) by 1 do:
        AnsysDataStatic[i,2] := 1 - AnsysDataStatic[i,2]/a__0:
    end do:

    AnsysPlots := [seq(AnsysPlots), listplot(AnsysDataStatic, labels=['u', 'r'], linestyle=solid, legend="Ansys Static Results")]:
fi:

# Export static solution file with current timestamp/analysis name
Export(cat("Snap2DResults\\",analysisName,"\\AnsysStatic.csv"), AnsysDataStatic, base=worksheetdir):
Export(cat("Snap2DResults\\",analysisName,"\\AnsysTrans.csv"), AnsysDataTrans, base=worksheetdir):

# For comparison plots, we can't have any evaluated names. Therefore, unset all input variables and rely on the values given in the Metadata csv.
# As I need comparisons, I will put the unset values here:
E__r := 'E__r':
eta := 'eta':

# Process comparison plots if they exist:
ComparePlotsStatic := []:
ComparePlotsStaticDiff := []:
ComparePlotsTransient := []:

compareStatic := Vector(ColumnDimension(compareTo));
compareStaticDiff := Vector(ColumnDimension(compareTo));
compareTrans := Vector(ColumnDimension(compareTo));
compareMetadata := Vector(ColumnDimension(compareTo));

idx := 1:
for compareEach in compareTo do:
compareStatic[idx] := Import(cat("Snap2DResults\\",compareEach,"\\AnsysStatic.csv"), base=worksheetdir, output=Array):
compareTrans[idx] := Import(cat("Snap2DResults\\",compareEach,"\\AnsysTrans.csv"), base=worksheetdir, output=Array):
compareMetadata[idx] := Import(cat("Snap2DResults\\",compareEach,"\\Metadata.csv"), base=worksheetdir, output=Array):


# Find comparison item so that we can create plot legend and title
# First attempt to use cached index to find value to avoid looping (is optimization in the studio with us today???)
if currentIndex > 0 and currentIndex < RowDimension(compareMetadata[idx]) and compareMetadata[idx][currentIndex,1] = compareAgainst then:

comparisonIndex := currentIndex;

else
# Second attempt; loop to find index
comparisonIndex := -1;
currentIndex := 1;
while currentIndex < RowDimension(compareMetadata[idx]) and comparisonIndex = -1 do:

if compareMetadata[idx][currentIndex,1] = compareAgainst then:

comparisonIndex := currentIndex;

fi;

currentIndex += 1;

end do:
fi:

# Give up :\ just use the name of the analysis
if comparisonIndex = -1 then:

compareInfo := compareEach;

else

# we found it
compareInfo := convert(compareAgainst,'symbol')," = ",compareMetadata[idx][comparisonIndex,2];

fi:

ComparePlotsStatic := [seq(ComparePlotsStatic), listplot(compareStatic[idx], labels=['u', 'r'], linestyle=solid, legend=typeset(compareInfo, ": Static Results"), color=col(idx), linestyle=ls(idx))];
ComparePlotsTransient := [seq(ComparePlotsTransient), listplot(compareTrans[idx], labels=['u', 'r'], linestyle=solid, legend=typeset(compareInfo, ": Transient Results"), color=col(idx), linestyle=ls(idx))];

# also for kicks and giggles mainly, compute the differential of the static plot using finite difference
compareStaticDiff[idx] := [-1];
for j from 1 to RowDimension(compareStatic[idx])-1 do:
compareStaticDiff[idx] := [seq(compareStaticDiff[idx]), (compareStatic[idx][j+1,1] - compareStatic[idx][j,1])/(compareStatic[idx][j+1,2] - compareStatic[idx][j,2])];
end do:

ComparePlotsStaticDiff := [seq(ComparePlotsStaticDiff), pointplot(Vector(compareStaticDiff[idx][4..-1]), Vector(compareStatic[idx][4..-1,2]), style=line, labels=[typeset(Diff(u,r)), 'r'], linestyle=solid, legend=typeset(compareInfo, ": Static Results"), color=col(idx), linestyle=ls(idx))];

idx += 1;
end do:

# Plot final version of total comparison plot
thesisStyle(display(allSolutions, seq(AnsysPlots), seq(ComparePlotsStatic), seq(ComparePlotsTransient)));

thesisStyle(display(seq(ComparePlotsStatic)));
thesisStyle(display(seq(ComparePlotsTransient)));
thesisStyle(display(seq(ComparePlotsStatic), seq(ComparePlotsTransient)));
thesisStyle(display(seq(ComparePlotsStaticDiff)));
# Curve fitting testing
idx := 5;
bestFitTest := Statistics:-PolynomialFit(3, compareStatic[idx][1..-1,2], compareStatic[idx][1..-1,1], 'r');
uFit := r -> bestFitTest;
testPlot := plot([seq(uFit(r), r=-1..1, 0.05)], [seq(r, r=-1..1, 0.05)], legend="Cubic Fit"):
display(testPlot,ComparePlotsStatic[idx]);
staticSolCompare := [seq(coeff(sort(subs(r__0=r, staticSolution), r, ascending), 'r', k), k=0..3)];
bestFitTestCompare := [seq(coeff(sort(bestFitTest, r, ascending), 'r', k), k=0..3)];

for k from 1 to 4 by 1 do:
if k <> 3 then
print(bestFitTestCompare[k] / staticSolCompare[k]);
fi;
end do;


pointplot([seq(x, x=0..1, 0.05)], [seq(x^3, x=0..1, 0.05)], symbol=box, linestyle=solid, style=pointline);
# Determine elapsed analysis time (in seconds)
elapsedTime := time[real]() - currentTime;

