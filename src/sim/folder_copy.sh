sourceDir="/orozco/projects/MontecarloJurgen/data/chromatin/celegans/long"
targetDir="/orozco/projects/MontecarloJurgen/data/chromatin/2000bp_linker/long"

find "$sourceDir" -type d | sed -e "s?$sourceDir?$targetDir?" | xargs mkdir -p
