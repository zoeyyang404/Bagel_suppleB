$Rscript = (Get-Command Rscript -ErrorAction Stop).Source

$scripts = @(
    "simulation/example1_varknown_calibration.R",
    "simulation/example1_varunknown_calibration.R",
    "simulation/example2_calibration_dis.R",
    "simulation/example2_calibration_conti.R",
    "simulation/example_seasonality_calibration.R"
)

foreach ($script in $scripts) {
    Write-Host "========================================"
    Write-Host "Running: $script"
    Write-Host "========================================"

    & $Rscript -e ".libPaths(c(normalizePath('Rlib'), .libPaths())); CALIB_NO_RUN <- TRUE; source('$script'); SKIP_EXISTING <- FALSE; run_calibration()"

    if ($LASTEXITCODE -ne 0) {
        Write-Error "$script failed with exit code $LASTEXITCODE"
        break
    }

    Write-Host "Finished: $script"
}