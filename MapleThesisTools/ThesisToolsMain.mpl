ThesisTools := module()
    description "Tools for use in Gabriel Tardy's thesis.";
    option package;
    export thesisStyle, lsymb, ls, col, ThesisStyle, SetPalette, TimeProfile;

    # Local variables to the ThesisTools module.
    local linestyles, linesymbols, currentTime, paletteName;

    # Set current time to -1 (implies not set)
    currentTime := -1;

    (* --- Styling options --- *)
    # Plots a displayed plot using the thesis-specific style.
    thesisStyle := proc(plotFunc):
        return display(plotFunc, thickness=3, font=["Times New Roman", roman, 16], labelfont=["Cambria Math", italic, 16], axis=[thickness=2, tickmarks=[thickness=0],gridlines=[thickness=0, colour=lightgray]], gridlines, legendstyle=[font=["Times New Roman", roman, 16]]);
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

# Snap buckling code is contained in its own file.
#$include "SnapBuckling.mm"

end module: