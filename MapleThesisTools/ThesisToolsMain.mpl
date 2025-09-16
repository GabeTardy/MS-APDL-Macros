ThesisTools := module()
    description "Tools for use in Gabriel Tardy's thesis.";
    option package;
    export thesisStyle, lsymb, ls, col, ThesisStyle, SetPalette, TimeProfile, SetLoadingType, AnsysExe, AnsysExeString, RunScenario;

    # Local variables to the ThesisTools module.
    local linestyles, linesymbols, currentTime, paletteName;

    # Snap-buckling utility variables
    local loadingTypes, loading, defaultAnalysis, ansysExec;
    
    # Location of AnsysNNN.exe (where NNN is the current version number (last two digits of year, then digit of revision number))
    ansysExec := "C:\\Program Files\\ANSYS Inc\\ANSYS Student\\v252\\ansys\\bin\\winx64\\ANSYS252.exe";

    # Comparison utility variables
    local compareTo := ["TimeTest; Hold; eta = 0.5, tau0 = 0.1", "TimeTest; Hold; eta = 0.5, tau0 = 1", "TimeTest; Hold; eta = 0.5, tau0 = 10", "TimeTest; Hold; eta = 0.5, tau0 = 50", "TimeTest; Hold; eta = 0.5, tau0 = 100", "TimeTest; Hold; eta = 2, tau0 = 0.1", "TimeTest; Hold; eta = 2, tau0 = 1", "TimeTest; Hold; eta = 2, tau0 = 10", "TimeTest; Hold; eta = 2, tau0 = 50", "TimeTest; Hold; eta = 2, tau0 = 100"];
    local compareAgainst := ["eta", "tau__0"];


    # Set current time to -1 (implies not set)
    currentTime := -1;

    (* --- Styling options --- *)
    # Plots a displayed plot using the thesis-specific style.
    thesisStyle := proc(plotFunc):
        return plots:-display(plotFunc, thickness=3, font=["Times New Roman", roman, 16], labelfont=["Cambria Math", italic, 16], axis=[thickness=2, tickmarks=[thickness=0],gridlines=[thickness=0, colour=lightgray]], gridlines, legendstyle=[font=["Times New Roman", roman, 16]]);
    end proc:
    ThesisStyle := thesisStyle: # The camelCase is a fallback name to maintain backwards compatibility.

    # All possible line styles and symbols.
    paletteName := "Spring": # Preferred palette name.
    linestyles := ['solid', 'dash', 'dashdot', 'longdash', 'spacedash']: # 'dot' is not included
    linesymbols := ['asterisk', 'box', 'circle', 'cross', 'diagonalcross', 'diamond', 'point', 'solidbox', 'solidcircle', 'soliddiamond']:

    # Proc for setting the color palette.
    SetPalette := proc(paletteStr::string)
        paletteName := paletteStr;
    end proc:

    # Proc to get a certain line symbol by index.
    lsymb := proc(it):
        return linesymbols[it mod 10 + 1];
    end proc:

    # Proc to get a certain line style by index.
    ls := proc(it):
        return linestyles[it mod 5 + 1];
    end proc:

    # Proc to get a certain color from the palette by index.
    col := proc(idx):
        return cat(paletteName, " ", idx);
    end proc:

    TimeProfile := proc()
        description "Toggles a time profiler on/off. Call once to start timing, call again to stop and return elapsed time.";
        local dt;

        if currentTime = -1 then
            currentTime := time[real]():
            return -1;
        else
            dt := time[real]() - currentTime;
            currentTime := -1:
            return dt;
        end if;
    end proc:

(* ---- SNAP BUCKLING ---- *)

RunScenario := proc(aName::string, nameSubs::uneval, {
    n::numeric := 4, 
    tau__0::numeric := 80,
    m::numeric := 1/2,
    zeta::numeric := 1/5,
    lambda::numeric := 1/50,
    k::numeric := 2,
    a__0::numeric := 1,
    E::numeric := .12000e8,
    rho::numeric := .179e-1,
    nu::numeric := 0,
    isBimodulus::numeric := 1,
    E__r::numeric := .12000e8,
    eta::numeric := .5,
    beta__Kuznetsov::numeric := .1,
    critical__E::numeric := .12000e8,
    primary__E::numeric := .12000e8,
    dataRes::numeric := .5e-1,
    paletteName::numeric := "Spring",
    generateAnimation::numeric := false,
    plotSmoothness::numeric := 5000,
    ansysExec::numeric := "C:\\Program Files\\ANSYS Inc\\ANSYS Student\\v252\\ansys\\bin\\winx64\\ANSYS252.exe",
    ansysStaticAnalysis::numeric := 1,
    ansysTransientAnalysis::numeric := 1,
    compareTo::numeric := ["Constant Er; eta = 0.25", "Constant Er; eta = 0.5", "Constant Er; eta = 0.75", "Constant Er; eta = 1", "Constant Er; eta = 1.5", "Constant Er; eta = 2", "Constant Er; eta = 4"]
}) uses LinearAlgebra;
    if assigned(loading) = false then
        error("You must set a loading type before running this script.");
    fi;

    local analysisName := aName:
    for i from 1 to nops(nameSubs) by 1 do:
        analysisName := StringTools:-Substitute(analysisName, cat("$", i), eval(eval((nameSubs[i]))));
    end do:
    print(analysisName);

    return;

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

    printf("\n\nMaple will now stop working until APDL is complete! This is normal and APDL may take anywhere from 2-5 minutes to complete.");
    APDLStopcode := ssystem(cat(ansysExec," -b -i Snap2DMaple.mac -o maple.out -j snapLastRun -smp"));

    # Wait until return code is [8, ""]; otherwise APDL failed (see APDL Operations Guide)
    if APDLStopcode[1] <> 8 then
        error(cat("APDL run resulted in an error. (Task ended with stopcode ",APDLStopcode[1],".) Please visit the file maple.out for more information."));
        return APDLStopcode[1];
    fi;

    # Switch E to the one used for nondimensionalization
    if isBimodulus > 0 then:
        E := E__nondimensionalization;
    fi:

    AnsysPlots := []:


    # Process transient solution if it exists;
    if ansysTransientAnalysis = 1 then:
        AnsysDataTrans := Import("Snap2DOutputTransient.csv",output=Array):

        # Multiply all values by omega to get taus
        AnsysDataTrans[1..-1,1] *= omega:

        # Send this version to a new variable which just includes time values
        AnsysDataTransTime := copy(AnsysDataTrans, deep);

        # Use loading function to get u from tau
        for i from 1 to RowDimension(AnsysDataTrans) by 1 do:
            AnsysDataTrans[i,1] := eval(uTauTrue, tau=AnsysDataTrans[i,1]);
            AnsysDataTrans[i,2] := 1 - AnsysDataTrans[i,2]/a__0;
            AnsysDataTransTime[i,2] := 1 - AnsysDataTransTime[i,2]/a__0;
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
    Export(cat("Snap2DResults\\",analysisName,"\\AnsysTransTime.csv"), AnsysDataTransTime, base=worksheetdir):
    return APDLStopcode[1];
end proc:

(* --- Loading --- *)
loadingTypes := table([
    onecycle = Matrix([
        [1, 1],
        [2, 0]
    ]), # 1-cycle generic
    onecyclepull = Matrix([
        [1, 1],
        [5/2, -1/2],
        [4, 1]
    ]), # 1-cycle pullback
    onecyclehold = Matrix([
        [1, 1],
        [3, 1]
    ]), # 1-cycle double-time hold (load twice as long)
    step = Matrix([
        [0, 1],
        [3, 1]
    ]), # Step load (not sure if this one works properly yet)
    twocycle = Matrix([
        [1, 1],
        [2, 0],
        [3, 1],
        [4, 0]
    ]), # 2-cycle generic
    twocyclepull = Matrix([
        [1, 1],
        [5/2, -1/2],
        [4, 1],
        [11/2, -1/2],
        [6,0]
    ]) # 2-cycle pullback
]);

SetLoadingType := proc({ type::symbol := step })
    description "Sets the loading type for snap-buckling analyses. Available types are: onecycle, onecyclepull, onecyclehold, step, twocycle, twocyclepull.";
    local possibleTypes;

    if assigned(loadingTypes[type]) = false then
        possibleTypes := [indices(loadingTypes, 'nolist')];
        error cat("Loading type '", type, "' is not recognized. Available types are: ", convert~(possibleTypes, 'string'));
        return loading;
    fi;

    loading := loadingTypes[type];
    return loading;
end proc:

# @TODO: Implement:
# - PlotLoading
# - ConvertLoading

# Snap-Buckling Utility Procs
AnsysExe := proc({version::integer := 252, edition::string := "Student"})
    description "Sets the path to the ANSYS executable by version number.";
    local path, pathStudent;

    path := cat("C:\\Program Files\\ANSYS Inc\\v", version, "\\ansys\\bin\\winx64\\ANSYS", version, ".exe");
    pathStudent := cat("C:\\Program Files\\ANSYS Inc\\ANSYS Student\\v", version, "\\ansys\\bin\\winx64\\ANSYS", version, ".exe");

    if FileTools:-Exists(path) and edition <> "Student" then
        ansysExec := path;
        return ansysExec;
    fi;

    if FileTools:-Exists(pathStudent) and edition = "Student" then
        ansysExec := pathStudent;
        return ansysExec;
    fi;

    error cat("APDL v", version, " not detected with these settings. Are you sure you are using the right version?");
    return ansysExec;
end proc:
AnsysExeString := proc(path::string)
    description "Sets the path to the ANSYS executable by string.";
    if path <> "" and FileTools:-Exists(path) then
        ansysExec := path;
        return ansysExec;
    fi;

    error "Invalid path specified.";
    return ansysExec;
end proc:
# proc(etc) uses LinearAlgebra

end module: