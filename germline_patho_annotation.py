import gzip
import argparse

# Several arguments are needed to run the python script
parser = argparse.ArgumentParser()
# Arguments are added to the parser
parser.add_argument("infile", help="Input VCF", type=str)
parser.add_argument("database", help="Potential pathogenic germline variant database", type=str)
parser.add_argument("outfile", help="Output annotated VCF", type=str)
parser.add_argument("type", help="Cancer type (TCGA code)", type=str)
args = parser.parse_args()

#file = open(
 #           "/nfs/backupp1/TCGA_patho_germline/TCGA_mutaciones_germinales/PCA_pathVar_integrated_filtered_adjusted_hg38.vcf", 'rb')
file = open(args.database, 'rb')
f = file.readlines()
file.close()

#type = "COAD"
type = args.type

patho = []
total = []

for line in f:
    if line.startswith("#"):
        continue

    line = line.rstrip()
    values = line.split("\t")
    info = values[7].split(";")
    cancer = info[1].split("=")[1]

    if cancer == type:
        id = values[0] + "_" + values[1] + "_" + values[3] + "_" + values[4]
        patho.append(id)
    else:
        id = values[0] + "_" + values[1] + "_" + values[3] + "_" + values[4]
        total.append(id)

#file = gzip.open("/media/scratch2/FIS_exomas/TANDA3/Sheila/tmp_tejido/ET189_germline.vcf.gz", 'rb')
file = gzip.open(args.infile, 'rb')
f = file.readlines()
file.close()

#outfile = gzip.open("/media/scratch2/FIS_exomas/TANDA3/Sheila/germline_patho//ET189_germline_patho.vcf.gz", 'wb')
outfile = gzip.open(args.outfile, 'wb')

for line in f:
    if line.startswith("#"):
        outfile.write(line)
        if line.startswith("##INFO=<ID=BAM,"):
            outfile.write("##INFO=<ID=GERMLINE_PATHO,Number=1,Type=String,Description=" +
                          '"Described as a pathogenic germline variant">' + "\n")
        continue

    line = line.rstrip()
    values = line.split("\t")

    id = values[0] + "_" + values[1] + "_" + values[3] + "_" + values[4]

    if type == "ALL":
        if id in total:
            outfile.write("\t".join(values[:6]) + "\t" + values[7] + ";GERMLINE_PATHO=YES" + "\t" + values[8] + "\t" + values[9] + "\n" )
        else:
            outfile.write("\t".join(values[:6]) + "\t" + values[7] + ";GERMLINE_PATHO=-" + "\t" + values[8] + "\t" + values[9] + "\n" )


    else:

        if id in patho:
            outfile.write("\t".join(values[:6]) + "\t" + values[7] + ";GERMLINE_PATHO=YES" + "\t" + values[8] + "\t" + values[9] + "\n" )
        else:
            outfile.write("\t".join(values[:6]) + "\t" + values[7] + ";GERMLINE_PATHO=-" + "\t" + values[8] + "\t" + values[9] + "\n" )

outfile.close()
