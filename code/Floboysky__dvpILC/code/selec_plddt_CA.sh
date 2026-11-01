#!/bin/bash

# Script to retrieve the plDDT scores for alpha carbon atoms (CA) from a .cif file and save them to a .json file.
# Usage: ./selec_plDDT_CA.sh /path/to/folder_with_cif_files
# NOTE: You may need to run the "dos2unix" command before you can run this script.


extract_plDDT_to_json() {
    local input_file="$1"
    local output_file="$2"

    if [[ ! -f "$input_file" ]]; then
        echo "The input file does not exist."
        return 1
    fi

    # Initializing the JSON file
    echo "{" > "$output_file"

    # Preparing arrays for each key
    declare -A keys
    keys=(
        #[ATOM]=""
        [number_atoms]=""
        #[type_atoms]=""
        [atoms]=""
        #[point]=""
        [name_res]=""
        [chain]=""
        #[nbr_chain]=""
        [nbr_res]=""
        #[intero]=""
        #[X]=""
        #[Y]=""
        #[Z]=""
        #[un]=""
        [plDDT]=""
    )

    while IFS= read -r line; do
        fields=($(echo "$line" | awk '{for (i=1; i<=NF; i++) print $i}'))

        #keys[ATOM]+="\"${fields[0]}\","      #  ATOM column
        keys[number_atoms]+="${fields[1]},"      # Column atom number
        #keys[type_atoms]+="\"${fields[2]}\"," # Column atom type
        keys[atoms]+="\"${fields[3]}\","      # Column atom name
        #keys[point]+="\"{${fields[4]}\","       # .
        keys[name_res]+="\"${fields[5]}\","   # Residue name
        keys[chain]+="\"${fields[6]}\","      # Chain
        #keys[nbr_chain]+="${fields[7]},"         # Chain number
        keys[nbr_res]+="${fields[8]},"           # Residue number
        #keys[intero]+="\"${fields[9]}\","      # ?
        #keys[X]+="${fields[10]},"                 # X coordinate
        #keys[Y]+="${fields[11]},"                 # Y coordinate
        #keys[Z]+="${fields[12]},"                # Z coordinate
        #keys[un]+="${fields[13]},"                # 1
        keys[plDDT]+="${fields[14]},"            # plDDT column
    done < <(grep "ATOM" "$input_file" | grep "CA")

    # Writing the results to the JSON file
    for key in "${!keys[@]}"; do
        # Cleaning the final commas
        cleaned_values=$(echo ${keys[$key]} | sed 's/,$//')
        echo "  \"$key\": [$cleaned_values]," >> "$output_file"
    done

    # Closing the JSON
    sed -i '$ s/,$//' "$output_file" # Removes the last comma
    echo "}" >> "$output_file"

    echo "Extraction completed. The results are in: $output_file"
}


#input_dir="$(pwd)"
input_dir="$1"

for input_file in "$input_dir"/*.cif; do
        if [[ -f "$input_file" ]]; then
            output_file="${input_file%.cif}_result.json"
            extract_plDDT_to_json "$input_file" "$output_file"
        fi
done
