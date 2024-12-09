localrules: download_intact, tfp_intact


rule download_intact:
    output:
        temp("dbs/hg38/tfp/intact/raw/intact.txt")
    shell:
        """
        wget --no-verbose https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip -O {output}.zip && \
        unzip -o {output}.zip -d $( dirname {output} ) && \
        rm {output}.zip
        """


rule tfp_intact:
    input:
        inc=rules.download_intact.output,
        lmb=rules.gen_tfs_lambert.output,
        pid=rules.gen_pid_uniprot.output,
    output: 'dbs/hg38/tfp/intact/intact.tsv'
    shell:
        """
        python workflow/scripts/dbs/tfp/intact.py \
        {input.inc} {input.lmb} {input.pid} {output}
        """


# Pubmed
# Rule to generate all unique TF pairs and save them to a file
rule generate_pairs:
    input:
        input_file = "data/lambert.csv"
    output:
        pair_file = "data/tf_pairs.txt"
    shell:
        """
        # Initialize output file
        > {output[0]}

        # Read TFs into an array
        tfs=()
        while IFS= read -r line; do
            tfs+=("$line")
        done < {input[0]}

        # Generate unique pairs and save to output file
        for ((i = 0; i < ${{#tfs[@]}}; i++)); do
            for ((j = i + 1; j < ${{#tfs[@]}}; j++)); do
                echo "${{tfs[i]}},${{tfs[j]}}" >> {output[0]}
            done
        done
        """

# Rule to perform PubMed searches for each TF in the list and save results to an output file
rule search_tf_abstracts:
    input:
        input_file = "data/lambert.csv"
    output:
        output_file = "results/tf_abstract_counts.txt" # abstracts for each tf
    shell:
        """
        # Initialize output file with a header
        echo -e "TF\tAbstracts" > {output}

        # Read TFs into an array, cleaning any extra whitespace or control characters
        tfs=()
        while IFS= read -r line; do
            clean_line=$(echo "$line" | tr -d '\r' | xargs)
            [[ -n "$clean_line" ]] && tfs+=("$clean_line")
        done < {input}

        # Loop through each TF in the array and perform the PubMed search
        for tf in "${{tfs[@]}}"; do
            echo "Searching for $tf in PubMed with Title/Abstract constraint..."

            # Perform the search with title/abstract constraint and get the count of abstracts
            count=$(esearch -db pubmed -query "$tf[Title/Abstract]" | xtract -pattern ENTREZ_DIRECT -element Count 2>/dev/null)

            # If count is empty, set it to 0
            if [[ -z "$count" ]]; then
                count=0
            fi

            # Append the TF and count to the output file
            echo -e "$tf\t$count" >> {output}
        done

        echo "Results saved to {output}"
        """

# Define the input and output files
chunk_dir = "chunks/"               # Directory to store chunk files
num_chunks = 30                    # Number of chunks (one for each core)

# Rule to split the pair file into chunks
rule split_pairs:
    input:
        input_file = "data/tf_pairs.txt"
    output:
        expand(chunk_dir + "chunk_{i}.txt", i=range(num_chunks))
    run:
        import os
        os.makedirs(chunk_dir, exist_ok=True)

        # Read all lines as full pairs, ensuring no split occurs in pairs
        with open(input_file, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        # Calculate the approximate number of pairs per chunk
        total_pairs = len(lines)
        chunk_size = total_pairs // num_chunks + (1 if total_pairs % num_chunks != 0 else 0)

        # Split the lines into chunks and write each chunk to a file
        for i in range(num_chunks):
            chunk_lines = lines[i * chunk_size : (i + 1) * chunk_size]
            with open(f"{chunk_dir}/chunk_{i}.txt", "w") as chunk_file:
                chunk_file.write("\n".join(chunk_lines) + "\n")



# Calculate the abstracts for all the possible pairs

num_chunks = 32  # Adjust this based on the actual number of chunks

# Defineix la llista completa de fitxers de sortida
#rule all:
#    input:
#        expand("results/chunk_{i}_results.txt", i=range(32))


rule pubmed_search_chunk:
    input:
        chunk_file = "chunks/chunk_{i}.txt"
    output:
        result_file = "results/chunk_{i}_results.txt"
    shell:
        """
        # Check if input file exists and is not empty
        input_file="{input.chunk_file}"
        output_file="{output.result_file}"

        if [[ ! -f "$input_file" ]]; then
            echo "Input file $input_file does not exist."
            exit 1
        fi

        if [[ ! -s "$input_file" ]]; then
            echo "Input file $input_file is empty."
            echo -e "TF_Pair\tAbstracts" > "$output_file"
            exit 0
        fi

        # Initialize an array for TF pairs
        pairs=()
        while IFS= read -r line; do
            pairs+=("$line")
        done < "$input_file"

        # Write a header to the output file
        echo -e "TF_Pair\tAbstracts" > "$output_file"
 # Check if esearch and xtract are available
        if ! command -v esearch &>/dev/null || ! command -v xtract &>/dev/null; then
            echo "Error: esearch or xtract commands are not available."
            exit 1
        fi

        # Loop through each TF pair
        for pair in "${{pairs[@]}}"; do
            IFS=',' read -r tf1 tf2 <<< "$pair"
            count=""
            attempts=0
            max_attempts=4

            while [[ -z "$count" && $attempts -lt $max_attempts ]]; do
                count=$(esearch -db pubmed -query "$tf1[Title/Abstract] AND $tf2[Title/Abstract]" | xtract -pattern ENTREZ_DIRECT -element Count 2>/dev/null)

                if [[ $? -ne 0 || -z "$count" ]]; then
                    count=""
                    ((attempts++))
                    sleep 0.125
                    echo "Retrying $tf1 and $tf2, attempt $attempts..."
                else
                    break
                fi
            done

            if [[ -z "$count" ]]; then
                count=0
            fi

            echo -e "$tf1,$tf2\t$count" >> "$output_file"
            sleep 0.125
        done
        """




# Once you have everything from pubmed, process data.
rule pubmed_final:
    input:
        pairs_path="data/pair_search_results.txt",  
        tf_path="data/tf_abstract_counts.txt"   
    output:
        out_path="data/pubmed_final.csv" 
    shell:
        """
        python scripts/pubmed_final.py -p {input.pairs_path} -t {input.tf_path} -o {output.out_path}
        """

