#to run open the terminal and type "bash ./run_file.sh 12"
id="$1"
Rscript pre_asr_formatting.R "$id"
./BayesTraitsV4 ./resolved_no_root_l1000_lineagegrp_"$id".nexus.trees ./lineagegrp_"$id"_data_file.txt < ./asr_lineagegrp_"$id"_input_file.txt > /dev/null
data_start_line=$(cat *.Log.txt | grep -n Lh | cut -d : -f 1)
tail -n +"$((data_start_line))"  *.Log.txt > dataframe.txt
python ./transpose.py dataframe.txt > asr_lineagegrp_"$id"_dataframe.txt
rm ./dataframe.txt
Rscript post_asr_formatting.R "$id"
