import csv
txt_file = "mummichog_input_ttest.txt"
csv_file = "differentially_abundant_tests/mummichog_input_ttest.csv"
with open(txt_file, "w") as output:
    with open(csv_file, "r") as input:
        [ output.write("\t".join(row)+'\n') for row in csv.reader(input)]
    output.close()