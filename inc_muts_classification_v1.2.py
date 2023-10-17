import argparse

# Several arguments are needed to run the python script
parser = argparse.ArgumentParser()
# Arguments are added to the parser
parser.add_argument("file", help="Oncokb output csv file", type=str)
parser.add_argument("outfile", help="annotated outfile", type=str)
parser.add_argument("gene_database", help="TSG/oncogene database path", type=str)
args = parser.parse_args()


######## CRITERIOS DE ANOTACION #########
#Grupo I
#	-Variantes con tratamiento dirigido
#Grupo II
#	-Variantes oncogenic/likely oncogenic revisadas por Oncokb
#Grupo III
#	-Variantes en genes clasificados como TSG/oncogen con AF poblacional <1% (no LOW, no MODIFIER)
#	-Variantes en genes clasificados como TSG/oncogen con valor "-" en AF poblacional (no hay valor, no se pueden filtrar) (no LOW, no MODIFIER)
#Grupo IV
#	-Variantes oncogenicas en genes no clasificados como TSG/oncogen (no LOW, no MODIFIER) <1% AF poblacional
#	-Variantes oncogenic/likely oncogenic LOW, pero en zonas de splicing <1% AF poblacional
#GRUPO V
#	-Resto variantes oncogenicas AF >1%
#GRUPO VI
#	-Resto de variantes (no cumplen anteriores requisitos)


#file = open("/nfs/home/references/annotations/human/GRCh38/oncogene_TSG_merged_revisado.tsv", 'rb')
file = open(args.gene_database, 'rb')
f = file.readlines()
file.close()

genes = {}

for line in f:

    line = line.rstrip()
    values = line.split("\t")
    genes[values[0]] = values[1]

file = open(args.file, 'rb')
f = file.readlines()
file.close()


outfile = open(args.outfile, 'wb')
outfile.write(f[0].rstrip() + ",INC_DRIVER\n")


index = {}
fields = f[0].rstrip().split(",")

for i in range(0, len(fields)):
    if fields[i] == "ONCOGENIC":
        index["ONCOGENIC"] = i
    elif fields[i] == "TREATMENT":
        index["TREATMENT"] = i
    elif fields[i] == "MAX_AF":
        index["MAX_AF"] = i
    elif fields[i] == "REVIEWED_VARIANT":
        index["REVIEWED_VARIANT"] = i
    elif fields[i] == "GENE_SYMBOL":
        index["GENE_SYMBOL"] = i
    elif fields[i] == "IMPACT":
        index["IMPACT"] = i
    elif fields[i] == "CONSEQUENCE":
        index["CONSEQUENCE"] = i


for line in f[1:]:

    line = line.rstrip()
    values = line.split(",")
    if "synonymous" in values[index["CONSEQUENCE"]] or "upstream" in values[index["CONSEQUENCE"]] or\
            "downstream" in values[index["CONSEQUENCE"]] or "Likely_Neutral" in values[index["ONCOGENIC"]] or\
            "Inconclusive" in values[index["ONCOGENIC"]] or "5_prime_UTR_variant" in values[index["CONSEQUENCE"]] or\
            "3_prime_UTR_variant" in values[index["CONSEQUENCE"]] or "coding_sequence_variant" in values[index["CONSEQUENCE"]] or\
            "NMD_transcript" in values[index["CONSEQUENCE"]] or "intergenic_variant" in values[index["CONSEQUENCE"]] or\
            "intron_variant" in values[index["CONSEQUENCE"]] or "mature_miRNA_variant" in values[index["CONSEQUENCE"]] or\
            "protein_altering_variant" in values[index["CONSEQUENCE"]]  or "non_coding_transcript" in values[index["CONSEQUENCE"]] or\
            "transcript_ablation" in values[index["CONSEQUENCE"]]:
        outfile.write(line + ",VI\n")
    else:
        if values[index["TREATMENT"]] != "" and values[index["TREATMENT"]] != "-":
            if values[index["REVIEWED_VARIANT"]] == "True":
                outfile.write(line + ",Ia\n")
            else:
                outfile.write(line + ",Ib\n")
        elif values[index["REVIEWED_VARIANT"]] == "True" and values[index["ONCOGENIC"]] == "Oncogenic" or\
                values[index["ONCOGENIC"]] == "Likely_Oncogenic" or values[index["ONCOGENIC"]] == "Predicted_Oncogenic":
            outfile.write(line + ",II\n")
        elif values[index["ONCOGENIC"]] == "Oncogenic" or\
            values[index["ONCOGENIC"]] == "Likely_Oncogenic" or values[index["ONCOGENIC"]] == "Predicted_Oncogenic":
            if values[index["IMPACT"]] != "LOW" and values[index["IMPACT"]] != "MODIFIER":
                if values[index["GENE_SYMBOL"]] in genes:
                    if values[index["MAX_AF"]] != "-":
                        if float(values[index["MAX_AF"]]) < 0.05:
                            outfile.write(line + ",III\n")
                        else:
                            outfile.write(line + ",VI\n")
                    else:
                        outfile.write(line + ",III\n")

                else:
                    if values[index["MAX_AF"]] != "-":
                        if float(values[index["MAX_AF"]]) < 0.05:
                            outfile.write(line + ",IV\n")
                        else:
                            outfile.write(line + ",VI\n")
                    else:
                        outfile.write(line + ",V\n")

            else:
                if "splice" in values[index["CONSEQUENCE"]]:
                    if values[index["MAX_AF"]] != "-":
                        if float(values[index["MAX_AF"]]) < 0.05 :
                            outfile.write(line + ",IV\n")
                        else:
                            outfile.write(line + ",VI\n")
                    else:
                        outfile.write(line + ",V\n")

                else:
                    outfile.write(line + ",V\n")
        elif values[index["MAX_AF"]] != "-":
            if float(values[index["MAX_AF"]]) < 0.05:
                outfile.write(line + ",V\n")
            else:
                outfile.write(line + ",VI\n")
        else:
                outfile.write(line + ",VI\n")

outfile.close()
