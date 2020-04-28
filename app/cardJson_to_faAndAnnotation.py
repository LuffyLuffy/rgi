#!/mnt/ilustre/users/huimin.zhang/newmdt/bin/anaconda/bin/python
import os
import csv
import json
working_directory = os.getcwd()

with open('card.json', 'r') as jfile:
    data = json.load(jfile)

version = data["_version"]

ncbi=True

annotations = []
with open(os.path.join(working_directory, "card_database_v{}.fasta".format(version)), 'w') as fout:
    for i in data:
        if i.isdigit():
            # use homolog models only
            if data[i]['model_type_id'] in ['40292']:

                drug_class = []
                mechanism = []
                group = []

                if "ARO_category" in data[i]:
                    for c in data[i]["ARO_category"]:
                        if "category_aro_class_name" in data[i]["ARO_category"][c]: 
                            if data[i]["ARO_category"][c]["category_aro_class_name"] in ["Drug Class"]:
                                drug_class.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))
                            if data[i]["ARO_category"][c]["category_aro_class_name"] in ["Resistance Mechanism"]:
                                mechanism.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))
                            if data[i]["ARO_category"][c]["category_aro_class_name"] in ["AMR Gene Family"]:
                                group.append(("{}".format(data[i]["ARO_category"][c]["category_aro_name"])))
                try:
                    for seq in data[i]['model_sequences']['sequence']:
                        # print(data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"])

                        if ncbi == True:
                            # header used to be able to validate CARD sequences with genbank sequences
                                header = ("{ncbi}|ARO:{ARO_accession}|ID:{model_id}|Name:{ARO_name}".format(
                                ARO_accession=data[i]['ARO_accession'],
                                model_id=data[i]['model_id'],
                                ARO_name=(data[i]['ARO_name']).replace(" ", "_"),
                                ncbi=data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"]
                                ))
                        else:
                            header = ("ARO:{}|ID:{}|Name:{}|NCBI:{}".format(
                                data[i]['ARO_accession'],
                                data[i]['model_id'],
                                (data[i]['ARO_name']).replace(" ", "_"),
                                data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"]
                                ))

                        sequence = data[i]['model_sequences']['sequence'][seq]["dna_sequence"]["sequence"]
                        fout.write(">{}\n".format(header))
                        fout.write("{}\n".format(sequence))
                        annotations.append([header, "; ".join(drug_class), "; ".join(mechanism), "; ".join(group)])

                except Exception as e:
                    print("No model sequences for model ({}, {}). Omitting this model and keep running.".format(data[i]['model_id'], data[i]['model_name']))
with open(os.path.join(working_directory, "card_annotations_v{}.txt".format(version)), "w") as af:
    writer = csv.writer(af, delimiter='\t', dialect='excel')
    writer.writerow(["header","class", "mechanism", "group"])
    for row in annotations:
        writer.writerow(row) 
