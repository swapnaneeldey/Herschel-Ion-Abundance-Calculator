import csv

from tabulate import tabulate
import abundance_calculator
import source_file as sf

source_list = []
source_rep = ""
for name in sf.source_name():
    source = name.split("_")[0]
    if source_rep != source:
        source_list.append(source)
        source_rep = source

transition_list = ["52", "57", "88", "122"]
spaxel = [('4', '0'), ('4', '1'), ('4', '2'), ('4', '3'), ('4', '4'), ('3', '0'), ('3', '1'), ('3', '2'),
          ('3', '3'), ('3', '4'), ('2', '0'), ('2', '1'), ('2', '2'), ('2', '3'), ('2', '4'), ('1', '0'),
          ('1', '1'), ('1', '2'), ('1', '3'), ('1', '4'), ('0', '0'), ('0', '1'), ('0', '2'), ('0', '3'),
          ('0', '4')]
source = "W3"
abundance_dictionary = {}
for spax in spaxel:
    abundance_dictionary[spax] = []
for transition in transition_list:
    abundance = abundance_calculator.source_abundance(f"{source}_{transition}")
    abundance_error = abundance_calculator.source_abundance_error(f"{source}_{transition}")
    for i, spax in enumerate(spaxel):
        abundance_dictionary[spax].append(abundance[i])
        abundance_dictionary[spax].append(abundance_error[i])

expression = r"spaxel, O++/H+ (52), error, N++/H+ (57), error, O++/H+ (88), error, N+/H+ (122)"

headers = ["spaxel", "O++/H+ (52)", "error", r"N++/H+ (57)", "error",
           r"O++/H+ (88)", "error", r"N+/H+ (122)", "error"]
table_data = [[key] + value for key, value in abundance_dictionary.items()]

output = expression + "\n" + "\n".join(",".join(str(item) for item in row) for row in table_data)

#Save the output to a CSV file
with open("/users/sdey/DATA/W3/table.csv", "w", newline='') as file:
    writer = csv.writer(file)
    writer.writerow(headers)
    writer.writerows(table_data)

print("Output saved to output.csv")
