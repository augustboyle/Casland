import pandas as pd
import argparse

complement = str.maketrans("ATGC","TACG")
transition = str.maketrans("ATGC","GCAT")
opposing = str.maketrans("ATGC","CGTA")

parser = argparse.ArgumentParser(description='Create a ddG_perturbation reference table for a custom sequence')
parser.add_argument('target',
	help='Perfect target sequence (nontemplate strand)')
parser.add_argument('modality', choices = ["binding","scission"],
	help='Select "binding" or "scission" to calculate ddG_pertubations')
parser.add_argument('--ddG_table', default = "S17_F7C_ddG_perturbation.tsv",
	help='Path to ddG_perturbation table')

args = parser.parse_args()
ddG_data = pd.read_table(args.ddG_table)

if args.modality == "binding":
	ddG_select = ddG_data[ddG_data['assay'] == "Productive binding"]
else:
	ddG_select = ddG_data[ddG_data['assay'] == "Scission"]

print("target_name\tddG_perturbation\tassay\ttarget_sequence")
for index, row in ddG_select.iterrows():
	name_details = row["name"].split("|")
	if name_details[0] == "perfect_target":
		continue
	if name_details[1] == "1":
		pos = int(name_details[2]) + 20
		base = args.target[pos]
		proximal_base = args.target[pos+1]
		complement_base = base.translate(complement)
		transition_base = base.translate(transition)
		opposing_base = base.translate(opposing)
		if base == proximal_base:
			insert_base = proximal_base.translate(complement)
		else:
			insert_base = proximal_base.translate(complement).translate(transition)

		if name_details[0] == "del_1":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos] + args.target[pos+1:])
		elif name_details[0] == "del_2":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos] + args.target[pos+2:])
		elif name_details[0] == "del_3":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos] + args.target[pos+3:])
		elif name_details[0] == "ins_1":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos + 1] + insert_base + args.target[pos+1:])
		elif name_details[0] == "ins_2":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos + 1] + insert_base*2 + args.target[pos+1:])
		elif name_details[0] == "ins_3":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos + 1] + insert_base*3 + args.target[pos+1:])
		elif name_details[0] == "sub_complement":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos] + complement_base + args.target[pos+1:])
		elif name_details[0] == "sub_transition":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos] + transition_base + args.target[pos+1:])
		elif name_details[0] == "sub_opposing":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos] + opposing_base + args.target[pos+1:])
	
	elif name_details[1] == "2":
		pos_1 = int(name_details[2]) + 20
		pos_2 = int(name_details[3]) + 20
		base_1 = args.target[pos_1]
		base_2 = args.target[pos_2]
		complement_base_1 = base_1.translate(complement)
		complement_base_2 = base_2.translate(complement)
		transition_base_1 = base_1.translate(transition)
		transition_base_2 = base_2.translate(transition)
		if name_details[0] == "complement_complement":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos_1] + complement_base_1 + args.target[pos_1+1:pos_2] + complement_base_2 + args.target[pos_2+1:])		
		elif name_details[0] == "transition_complement":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos_1] + complement_base_1 + args.target[pos_1+1:pos_2] + transition_base_2 + args.target[pos_2+1:])
		elif name_details[0] == "transition_transition":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos_1] + transition_base_1 + args.target[pos_1+1:pos_2] + transition_base_2 + args.target[pos_2+1:])
		elif name_details[0] == "grid_complement":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos_1] + args.target[pos_1:pos_2+1].translate(complement) + args.target[pos_2+1:])
		elif name_details[0] == "rgrid_complement":
			print(row["name"] + "\t" + str(row["ddG_perturbation"]) + "\t" + args.modality + "\t" + args.target[:pos_1+1].translate(complement) + args.target[pos_1+1:pos_2+1] + args.target[pos_2+1:20].translate(complement) + args.target[20:23])



