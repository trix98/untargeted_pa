import csv
txt_file = "mummichog/input_data/mummichog_input_ttest_rsd_1.txt"
csv_file = "mummichog/input_data/mummichog_input_ttest_rsd_1.csv"
with open(txt_file, "w") as output:
    with open(csv_file, "r") as input:
        [ output.write("\t".join(row)+'\n') for row in csv.reader(input)]
    output.close()