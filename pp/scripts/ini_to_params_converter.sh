#!/usr/bin/env bash

# Check if exactly two arguments are passed
if [ "$#" -ne 2 ]; then
    echo "Error: Missing arguments"
    echo "Usage: $0 ini_file params_file"
    exit 1  # Exit with a non-zero status to indicate an error
fi

ini_file="$1"
params_file="$2"
rm "$params_file" || echo "" # Overwrite previous params file

sed '/^\[/d' "$ini_file" > "$params_file"

replace_parameter_int() {
  parameter="$1"
  replacement="$2"
  params_file="$3"
  sed -i -E "s/^(${parameter}[[:space:]]*)= -1/\1= ${replacement}/" "$params_file"
}
replace_parameter_str() {
  parameter="$1"
  replacement="$2"
  params_file="$3"
  sed -i -E "s/^(${parameter}[[:space:]]*)= '-1'/\1= '${replacement}'/" "$params_file"
}
replace_parameter_real() {
  parameter="$1"
  replacement="$2"
  params_file="$3"
  sed -i -E "s/^(${parameter}[[:space:]]*)= -1.0/\1= ${replacement}/" "$params_file"
}
enclose_str_in_quotes() {
  parameter="$1"
  params_file="$2"
  # Wrap the first word
  #sed -i -E "/^${parameter}/s/=(\s*)(\S+)/= '\2'/" "${params_file}"
  # Wrap the entire line
  sed -i -E "/^${parameter}/s/=(\s*)(.*)/= '\2'/" "${params_file}"
}
replace_parameter_int "nlx" "ntile" "$params_file"
replace_parameter_int "nly" "ntile" "$params_file"
replace_parameter_int "nlz" "ntile" "$params_file"
replace_parameter_int "n1" "nmesh" "$params_file"
replace_parameter_int "n2" "nmesh" "$params_file"
replace_parameter_int "n3" "nmesh" "$params_file"
replace_parameter_int "ntasks" "nnodes * tpnode" "$params_file"

replace_parameter_str "densfilein" "run_name" "$params_file"
replace_parameter_str "densfileout" "run_name" "$params_file"
replace_parameter_str "filein" "run_name" "$params_file"
replace_parameter_str "filein" "run_name" "$params_file"
replace_parameter_str "fileout" "run_name" "$params_file"

replace_parameter_int "iZeld" "ilpt" "$params_file"
replace_parameter_int "ntilemerge" "ntile" "$params_file"
replace_parameter_int "ntasksmerge" "ntasks" "$params_file"
replace_parameter_int "nsub" "nmesh - 2 * nbuff" "$params_file"
replace_parameter_int "next" "nsub * ntile + 2 * nbuff" "$params_file"

replace_parameter_real "Omvac" "1.0 - Omx - OmB" "$params_file"
replace_parameter_real "dcore_box" "boxsize \/ ntile" "$params_file"
replace_parameter_real "cellsize" "dcore_box \/ nsub" "$params_file"
replace_parameter_real "buffersize" "cellsize * nbuff" "$params_file"
replace_parameter_real "dL_box" "dcore_box + 2 * buffersize" "$params_file"
replace_parameter_real "mlatt" "2.775e11 * (Omx + OmB) * h**2 * (cellsize)**3" "$params_file"

enclose_str_in_quotes "machine" "$params_file"
enclose_str_in_quotes "submit_command" "$params_file"
enclose_str_in_quotes "run_name" "$params_file"
enclose_str_in_quotes "short_name" "$params_file"
enclose_str_in_quotes "runtype" "$params_file"
enclose_str_in_quotes "tlimit" "$params_file"
enclose_str_in_quotes "fielddir" "$params_file"
enclose_str_in_quotes "maps" "$params_file"
enclose_str_in_quotes "tabfile_map" "$params_file"
enclose_str_in_quotes "tabfile_sfr" "$params_file"
enclose_str_in_quotes "pkfile" "$params_file"
enclose_str_in_quotes "rapi" "$params_file"
enclose_str_in_quotes "TabInterpFile" "$params_file"
enclose_str_in_quotes "filterfile" "$params_file"
