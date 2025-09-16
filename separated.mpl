for idex from 1 to nops(compareAgainst) do:
    # First attempt to find the index: try to exploit a cached value
    if currentIndex[idex] > 0 and currentIndex[idex] < RowDimension(compareMetadata[idx]) and compareMetadata[idx][currentIndex[idex],1] = compareAgainst then:
        comparisonIndex := currentIndex[idex]:
    else
        # Second attempt; loop to find index
        comparisonIndex := -1:
        currentIndex[idex] := 1:
        while currentIndex[idex] < RowDimension(compareMetadata[idx]) and comparisonIndex = -1 do:
            if compareMetadata[idx][currentIndex[idex],1] = compareAgainst then:
                comparisonIndex := currentIndex[idex];
            fi;

            currentIndex[idex] += 1;
        end do:
    fi:

    # No third attempt: give up :\ just use the name of the analysis and break the loop.
    if comparisonIndex = -1 then:
        compareInfo := compareEach;
        WARNING("Comparison data is irregular. Are you sure you are using the right symbol?");
        break; # exit loop
    else
        # we found it
        compareInfo := cat(compareInfo, convert(compareAgainst,'symbol')," = ",compareMetadata[idx][comparisonIndex,2]);
    fi:
end do: