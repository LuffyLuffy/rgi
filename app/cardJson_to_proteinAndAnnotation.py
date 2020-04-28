#!/mnt/ilustre/users/huimin.zhang/newmdt/bin/anaconda/bin/python
import os
import csv
import json
working_directory = os.getcwd()
with open('card.json','r') as jsfile:
    j = json.load(jsfile)
version = j["_version"]
annotations = []
with open('proteindb2.fsa','w') as fout:
    for i in j:
        if i.isdigit():
            # model_type: protein homolog model
            if j[i]['model_type_id'] == '40292':
                drug_class = []
                mechanism = []
                group = []
                if "ARO_category" in j[i]:
                    for c in j[i]["ARO_category"]:
                        if "category_aro_class_name" in j[i]["ARO_category"][c]: 
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["Drug Class"]:
                                drug_class.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["Resistance Mechanism"]:
                                mechanism.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["AMR Gene Family"]:
                                group.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))                

                try:
                    pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
                except KeyError:
                    logger.warning("No bitscore for model (%s, %s). RGI will omit this model and keep running." \
                        % (j[i]['model_id'], j[i]['model_name']))
                    logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")
                else:
                    try:
                        for seq in j[i]['model_sequences']['sequence']:
                            fout.write('>%s_%s | model_type_id: 40292 | pass_bitscore: %s | %s\n' % (i, seq, pass_bit_score, j[i]['ARO_name']))
                            fout.write('%s\n' %(j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))
                            # header used to be able to validate CARD sequences with genbank sequences
                            header = ("{iseq}|{ncbi}|ARO:{ARO_accession}|ID:{model_id}|Name:{ARO_name}".format(
                            ARO_accession=j[i]['ARO_accession'],
                            model_id=j[i]['model_id'],
                            ARO_name=(j[i]['ARO_name']).replace(" ", "_"),
                            ncbi=j[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"],
                            iseq=str(i)+"_"+str(seq)
                            ))
                            annotations.append([header, "; ".join(drug_class), "; ".join(mechanism), "; ".join(group)])
                    except Exception as e:
                        logger.warning("No model sequences for model (%s, %s). RGI will omit this model and keep running." \
                            % (j[i]['model_id'], j[i]['model_name']))
                        logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")
                


            # model_type: protein variant model
            elif j[i]["model_type_id"] == "40293":
                drug_class = []
                mechanism = []
                group = []
                if "ARO_category" in j[i]:
                    for c in j[i]["ARO_category"]:
                        if "category_aro_class_name" in j[i]["ARO_category"][c]: 
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["Drug Class"]:
                                drug_class.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["Resistance Mechanism"]:
                                mechanism.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["AMR Gene Family"]:
                                group.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))                    
                try:
                    pass_bit_score = j[i]['model_param']['blastp_bit_score']['param_value']
                except KeyError:
                    logger.warning("No bitscore for model (%s, %s). RGI will omit this model and keep running." \
                        % (j[i]['model_id'], j[i]['model_name']))
                    logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")
                else:
                    try:
                        snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
                    except Exception as e:
                        logger.warning("No snp for model (%s, %s). RGI will omit this model and keep running." \
                            % (j[i]['model_id'], j[i]['model_name']))
                        logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")

                    try:
                        for seq in j[i]['model_sequences']['sequence']:
                            fout.write('>%s_%s | model_type_id: 40293 | pass_bit_score: %s | SNP: %s | %s\n' \
                                % (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
                            fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))
                            header = ("{iseq}|{ncbi}|ARO:{ARO_accession}|ID:{model_id}|Name:{ARO_name}".format(
                            ARO_accession=j[i]['ARO_accession'],
                            model_id=j[i]['model_id'],
                            ARO_name=(j[i]['ARO_name']).replace(" ", "_"),
                            ncbi=j[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"],
                            iseq=str(i)+"_"+str(seq)
                            ))
                            annotations.append([header, "; ".join(drug_class), "; ".join(mechanism), "; ".join(group)])
                    except Exception as e:
                        logger.warning("No model sequences for model (%s, %s). RGI will omit this model and keep running." \
                            % (j[i]['model_id'], j[i]['model_name']))
                        logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")

            # model_type: protein overexpression model
            elif j[i]["model_type_id"] == "41091":
                drug_class = []
                mechanism = []
                group = []
                if "ARO_category" in j[i]:
                    for c in j[i]["ARO_category"]:
                        if "category_aro_class_name" in j[i]["ARO_category"][c]: 
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["Drug Class"]:
                                drug_class.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["Resistance Mechanism"]:
                                mechanism.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))
                            if j[i]["ARO_category"][c]["category_aro_class_name"] in ["AMR Gene Family"]:
                                group.append(("{}".format(j[i]["ARO_category"][c]["category_aro_name"])))    
                try:
                    pass_bit_score = j[i]["model_param"]["blastp_bit_score"]["param_value"]
                except KeyError:
                    logger.warning("No bitscore for model (%s, %s). RGI will omit this model and keep running." \
                        % (j[i]['model_id'], j[i]['model_name']))
                    logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")
                else:
                    try:
                        snpList = [j[i]['model_param']['snp']['param_value'][k] for k in j[i]['model_param']['snp']['param_value']]
                    except Exception as e:
                        logger.warning("No snp for model (%s, %s). RGI will omit this model and keep running." \
                            % (j[i]['model_id'], j[i]['model_name']))
                        logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")

                    try:
                        for seq in j[i]['model_sequences']['sequence']:
                            fout.write('>%s_%s | model_type_id: 41091 | pass_bit_score: %s | SNP: %s | %s\n' \
                                % (i, seq, pass_bit_score, ','.join(snpList), j[i]['ARO_name']))
                            fout.write('%s\n' % (j[i]['model_sequences']['sequence'][seq]['protein_sequence']['sequence']))
                            header = ("{iseq}|{ncbi}|ARO:{ARO_accession}|ID:{model_id}|Name:{ARO_name}".format(
                            ARO_accession=j[i]['ARO_accession'],
                            model_id=j[i]['model_id'],
                            ARO_name=(j[i]['ARO_name']).replace(" ", "_"),
                            ncbi=j[i]['model_sequences']['sequence'][seq]["dna_sequence"]["accession"],
                            iseq=str(i)+"_"+str(seq)
                            ))
                            annotations.append([header, "; ".join(drug_class), "; ".join(mechanism), "; ".join(group)])

                    except Exception as e:
                        logger.warning("No model sequences for model (%s, %s). RGI will omit this model and keep running." \
                            % (j[i]['model_id'], j[i]['model_name']))
                        logger.info("Please let the CARD Admins know! Email: card@mcmaster.ca")

with open(os.path.join(working_directory, "card_protein_annotations_v{}.txt".format(version)), "w") as af:
    writer = csv.writer(af, delimiter='\t', dialect='excel')
    writer.writerow(["header","class", "mechanism", "group"])
    for row in annotations:
        writer.writerow(row) 
