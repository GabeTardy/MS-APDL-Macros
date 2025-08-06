# PowerShell script to set APDL environment variables for USERMAT development; also patches ANSUPF.bat to update compiler references.
# --- Requires elevation to fix ANSUPF.bat; but does not require elevation to set user environment variables.
# --- Compiled by Gabriel Tardy, but written primarily by ChatGPT and Copilot because quite honestly I did not want to do it myself.

# *** Why does this script exist? (Various notes) ***
# In this script, I modify the user environment variables to reflect the proper setup as documented in the Programmer's Reference.
# I also do a (in my opinion) dangerous manipulation of ANSUPF.bat, which is not ideal.
# This could be avoided if ANSYS updated ANSUPF.bat to reflect the latest oneAPI compiler. 
# In 2025R2, they correctly updated the Visual Studio version, but did not update the ifort compiler or bump the version, and so I do that in this script.
# Additionally, they reference an incorrect vars.bat file, which I also fix.
# I am not sure why this was not done already and I assume it must be an oversight. 
# This would not really be a problem under normal circumstances, but I cannot figure out how to legally and safely obtain old versions of the Fortran compiler.
# I have tested my modifications using the USERMAT I developed for my thesis, and it works exactly as expected.

# Define base install directory
$basePath = "C:\Program Files\ANSYS Inc\ANSYS Student"

# Check if ANSYS Student folder exists
if (-not (Test-Path $basePath)) {
    $basePath = "C:\Program Files\ANSYS Inc"

    if (-not (Test-Path $basePath)) {
        Write-Error "Neither ANSYS Student nor ANSYS base path found."
        exit 1
    } else {
        Write-Warning "ANSYS Student folder not found. Using ANSYS base path: $basePath. Make sure customization files are installed."
    }
}

# Get list of versioned folders (v###)
$versionFolders = Get-ChildItem -Path $basePath -Directory | Where-Object {
    $_.Name -match '^v\d+$'
}

if ($versionFolders.Count -eq 0) {
    Write-Error "No ANSYS version folders found in: $basePath"
    exit 1
}

# Extract version info and find highest version
$versions = $versionFolders | ForEach-Object {
    [PSCustomObject]@{
        Name = $_.Name
        Version = [int]($_.Name -replace '^v', '')
        FullPath = "$($_.FullName)\ansys\bin\winx64"
    }
}

$latest = $versions | Sort-Object -Property Version -Descending | Select-Object -First 1
$latestPath = $latest.FullPath
$latestVersion = $latest.Version

# Ensure latest path exists
if (-not (Test-Path $latestPath)) {
    Write-Error "Latest ANSYS path does not exist: $latestPath"
    exit 1
}

# Get current PATH for the user
$currentPath = [Environment]::GetEnvironmentVariable("Path", "User")
$currentPathItems = $currentPath -split ';'

# Remove old ANSYS paths (vXXX < latest)
$filteredPathItems = $currentPathItems | Where-Object {
    ($_ -notmatch '\\ansys\\bin\\winx64') -or
    ($_ -match 'v(\d+)' -and [int]($Matches[1]) -ge $latestVersion)
}

# Add the latest ANSYS path if not already present
if (-not ($filteredPathItems -contains $latestPath)) {
    $filteredPathItems += $latestPath
    Write-Host "Added ANSYS path to PATH: $latestPath"
} else {
    Write-Host "Latest ANSYS path already in PATH: $latestPath"
}

# Join back into a string and update PATH
$newPath = ($filteredPathItems -join ';').TrimEnd(';')
[Environment]::SetEnvironmentVariable("Path", $newPath, "User")

# Prompt user to select folder for ANS_USER_PATH
Add-Type -AssemblyName System.Windows.Forms
$dialog = New-Object System.Windows.Forms.FolderBrowserDialog
$dialog.Description = "Select folder for ANS_USER_PATH"

if ($dialog.ShowDialog() -eq [System.Windows.Forms.DialogResult]::OK) {
    [Environment]::SetEnvironmentVariable("ANS_USER_PATH", $dialog.SelectedPath, "User")
    Write-Host "Set ANS_USER_PATH to: $($dialog.SelectedPath)"
} else {
    Write-Warning "No folder selected. ANS_USER_PATH not set."
}

# Set ANS_USE_UPF to true
[Environment]::SetEnvironmentVariable("ANS_USE_UPF", "TRUE", "User")
Write-Host "Set ANS_USE_UPF to true."

# Get the root ANSYS version directory from $latestPath
$versionRoot = $latestPath | Split-Path | Split-Path  # navigates to ...\v252\ansys

# Define path to ANSUPF.bat using relative path from $latestPath
$batchFilePath = Join-Path -Path $versionRoot -ChildPath "custom\user\winx64\ANSUPF.bat"

# Validate that the file exists
if (-not (Test-Path $batchFilePath)) {
    Write-Error "ANSUPF.bat not found at: $batchFilePath"
    exit 1
}

# Find all system environment variables that match IFORT_COMPILERXX
$compilerVars = Get-ChildItem Env: | Where-Object {
    $_.Name -match '^IFORT_COMPILER(\d+)$'
} | ForEach-Object {
    [PSCustomObject]@{
        Name = $_.Name
        Version = [int]($_.Name -replace '^IFORT_COMPILER', '')
    }
}

if ($compilerVars.Count -eq 0) {
    Write-Error "No IFORT_COMPILERXX environment variables found."
    exit 1
}

# Get the highest IFORT_COMPILER version
$latestCompiler = $compilerVars | Sort-Object Version -Descending | Select-Object -First 1
$compilerEnvVar = "%IFORT_COMPILER$($latestCompiler.Version)%"

# Read contents of ANSUPF.bat
$originalContent = Get-Content -Path $batchFilePath

# Replace all occurrences of %IFORT_COMPILER##% and "ifort" with the latest values
$modifiedContent = $originalContent | ForEach-Object {
    $_ -replace '%IFORT_COMPILER\d+%', $compilerEnvVar -replace '\bifort\b', 'ifx' -replace '"%IFORT_HOME%\\..\\env\\vars.bat"', '"%IFORT_HOME%\..\..\setvars.bat"'
}

# Create a backup before overwriting
$backupPath = "$batchFilePath.bak"
Copy-Item -Path $batchFilePath -Destination $backupPath -Force
Write-Host "Backup of ANSUPF.bat created at: $backupPath"

# Write updated content back to the original file
Set-Content -Path $batchFilePath -Value $modifiedContent -Encoding ASCII
Write-Host "ANSUPF.bat updated successfully:"
Write-Host "- Replaced IFORT_COMPILER## with $compilerEnvVar"
Write-Host "- Replaced 'ifort' with 'ifx'"