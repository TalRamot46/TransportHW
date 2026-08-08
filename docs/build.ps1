<#
.SYNOPSIS
    Builds a LaTeX document in docs/<name>/ with pdflatex.

.DESCRIPTION
    latexmk does not work in this repository. On this machine a file that has
    just been written cannot immediately be overwritten -- the write fails with
    a Windows EINVAL, which pdflatex reports as

        ! I can't write on file `<name>.pdf'.

    latexmk runs its pdflatex passes back to back, so its second pass always
    collides with the file its first pass just produced. Running latexmk through
    perl directly fails identically, because latexmk is not the cause.

    Deleting the output before each pass avoids the collision entirely. The same
    workaround is applied to the figures, in the _savefig helper of
    src/homework1/plots.py.

    A single pdflatex pass succeeds, so the problem only appears when a document
    is rebuilt -- which is also why it went unnoticed until the figures were
    regenerated.

.PARAMETER Name
    Document name without extension, e.g. "homework1" for docs/homework1/homework1.tex.

.PARAMETER Passes
    Number of pdflatex passes. Three is enough to resolve cross-references.

.EXAMPLE
    .\docs\build.ps1 homework1
#>
param(
    [Parameter(Mandatory = $true)][string]$Name,
    [int]$Passes = 3
)

$ErrorActionPreference = 'Stop'

$docDir = Join-Path $PSScriptRoot $Name
$texFile = Join-Path $docDir "$Name.tex"

if (-not (Test-Path $texFile)) {
    throw "No such document: $texFile"
}

Push-Location $docDir
try {
    for ($i = 1; $i -le $Passes; $i++) {
        # Delete the output first; overwriting it in place is what fails.
        Remove-Item "$Name.pdf" -Force -ErrorAction SilentlyContinue

        $output = & pdflatex -interaction=nonstopmode "$Name.tex" 2>&1
        $blocked = $output | Select-String -Pattern "can't write on file"
        if ($blocked) {
            $output | Select-String -Pattern '^!' | Select-Object -First 5
            throw "pass ${i}: pdflatex could not write its output."
        }

        $fatal = $output | Select-String -Pattern '^! '
        if ($fatal -and -not (Test-Path "$Name.pdf")) {
            $fatal | Select-Object -First 5
            throw "pass ${i}: pdflatex failed."
        }

        Write-Host "pass ${i}: ok"
    }

    $unresolved = Select-String -Path "$Name.log" -Pattern 'Reference .* undefined|Rerun to get'
    if ($unresolved) {
        Write-Warning "Unresolved references remain; try more passes:"
        $unresolved | Select-Object -First 5
    }

    $pdf = Get-Item "$Name.pdf"
    Write-Host "built $($pdf.FullName) ($($pdf.Length) bytes)"
}
finally {
    Pop-Location
}
