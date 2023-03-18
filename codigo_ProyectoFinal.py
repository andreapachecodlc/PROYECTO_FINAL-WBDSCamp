#Importación de bibliotecas necesarias
#el archivo credentials.py contiene el correo electrónico y la API key que proporciona NCBI a sus usuarios registrados):

import credentials
import matplotlib.pyplot as plt
import numpy  as np
import pandas as pd
import pyrodigal
import requests
import seaborn as sns
import subprocess
import sys
from Bio import SeqIO
from Bio import Entrez
from io                 import StringIO
from matplotlib.patches import Patch
from pycirclize         import Circos
from pycirclize.parser  import Gff
from requests.adapters  import HTTPAdapter, Retry

#Obtención de una secuencia genómica 
#Se usó el genoma de Bacillus paralicheniformis cepa CPL618 con número de accesión CP090833.1

accession = "CP090833.1"
genome = Entrez.efetch(db="nucleotide",
                       id=accession,
                       format="gb",
                       rettype="text")
record = SeqIO.read(genome, "genbank")


#longitud en bp del genoma descargado

genome_length = len(record.seq)
genome_length

#Predicción de genes usando pyrodigal

#búsqueda de marcos abiertos de lectura en el genoma descargado

orf_finder = pyrodigal.OrfFinder()
orf_finder.train(bytes(record.seq))
orf_genes  = orf_finder.find_genes(bytes(record.seq))

#Con el siguiente código se almacenarán las secuencias aminoacídicas de los genes predichos en un nuevo archivo CP090833.1.faa
#Se usó "BacParCP" como prefijo del identificador de cada secuencia obtenida

aa_file = accession + ".faa"
prefix  = "BacParCP"
with open(aa_file, "w") as orf_gene:
    orf_genes.write_translations(orf_gene,sequence_id=prefix)

#También almacenaremos las características y coordenadas de los genes predichos en un nuevo archivo CP090833.1.gff
gff_file = accession + ".gff"
prefix  = "BacParCP"
with open(gff_file, "w") as orf_gene:
    orf_genes.write_gff(orf_gene,sequence_id=prefix)

#Obtención de un set de secuencias de referencia de la la API de UniProt mediante requests. 
#La consulta (query) será la palabra "biosurfactant" y también se solicitará que las secuencias tengan el status "reviewed"

uniprot_api_url  = "https://rest.uniprot.org/uniprotkb/stream"
uniprot_api_args = {"compressed" : "false",
                    "format"     : "fasta",
                    "query"      : "(biosurfactant) AND (reviewed:true)"}
uniprot_ref_seqs = requests.get(uniprot_api_url,params=uniprot_api_args).text

#El siguiente código pasará el objeto uniprot_ref_seqs al archivo uniprot_sequences.fasta:

uniprot_seqs_file = open("uniprot_sequences.fasta", "wt")
uniprot_seqs_file.write(uniprot_ref_seqs)
uniprot_seqs_file.close()

#Uso de BLAST para la comparación de secuencias
#Indicar adecuadamente el path en el que se encuentra instalados los programas makeblastdb y blastp

makeblastdb_path = "/bin/makeblastdb" #Path de programa makeblastdb instalado previamente en PC
makeblastdb_command = [makeblastdb_path,'-in',aa_file,'-dbtype','prot']
subprocess.call(makeblastdb_command)

blastp_path       = "/bin/blastp" #Path de programa blastp instalado previamente en PC
blastp_out_format = "6 qseqid sseqid qlen slen qstart sstart qend send score evalue length positive"
blastp_out_file   = accession + ".blast.tsv"
blastp_command    = [blastp_path,
                     "-db",          aa_file,
                     "-query",       "uniprot_sequences.fasta",
                     "-evalue",      "1e-6",
                     "-out",         blastp_out_file,
                     "-outfmt",      blastp_out_format,
                     "-num_threads", "12"]
subprocess.call(blastp_command)

#resultados de BLAST 

blastp_column_names = blastp_out_format.split(" ")[1:]
blastp_df = pd.read_csv(blastp_out_file,sep="\t",names=blastp_column_names)

#Segundo análisis BLAST

makeblastdb_path = "/bin/makeblastdb" #Path de programa makeblastdb instalado previamente en PC
makeblastdb_command = [makeblastdb_path,'-in',"uniprot_sequences.fasta",'-dbtype','prot']
subprocess.call(makeblastdb_command)

blastp_path      = "/bin/blastp" #Path de programa blastp instalado previamente en PC
blast_out_format = "6 qseqid sseqid qlen slen qstart sstart qend send score evalue length positive"
blast_out_file   = "uniprot_sequences.blast.tsv"
blastp_command   = [blastp_path,
                    "-db",          "uniprot_sequences.fasta",
                    "-query",       aa_file,
                    "-evalue",      "1e-6",
                    "-out",         blast_out_file,
                    "-outfmt",      blast_out_format,
                    "-num_threads", "12"]
subprocess.call(blastp_command)

blastp_column_names = blast_out_format.split(" ")[1:]
blastp_rev_df = pd.read_csv(blast_out_file,sep="\t",names=blastp_column_names)


#dataframe consenso

blastp_df['both'] = blastp_rev_df.qseqid.isin(blastp_df.sseqid)
blastp_cons_df = blastp_df.drop(blastp_df[blastp_df['both']==False].index)
candidate_genes=blastp_cons_df["sseqid"].unique().tolist()
print('Genes candidatos resultantes de BLAST:'+' '+str(len(candidate_genes)))

#Visualización de datos obtenidos:
#Visualización de los genes de interés en el genoma bacteriano usando pycirclize

gff_columns     = ["chr","source","feature_type","start","end","score","strand","phase","info"]
gff_df          = pd.read_csv(gff_file,sep="\t",comment="#",header=None,names=gff_columns)
gff_df["start"] = gff_df["start"].astype(int)
gff_df["end"]   = gff_df["end"].astype(int)

def get_gff_info(info_str):
    out_dict = {}
    info_arr = info_str.split(";")
    for line in info_arr:
        if "=" in line:
            line_arr    = line.split("=")
            field_name  = line_arr[0]
            field_value = line_arr[1]
            out_dict[field_name] = field_value
    return out_dict

gff_df["annotation"] = gff_df["info"].apply(lambda x: get_gff_info(x))

gff_df["candidate"] = gff_df["annotation"].apply(lambda x: "include" if x["ID"] in candidate_genes else "exclude")

candidate_df = gff_df.copy()
candidate_df = candidate_df[candidate_df["candidate"]=="include"][gff_columns]
candidate_df.to_csv("candidates.gff",sep="\t",header=False,index=False)

circos = Circos(sectors={accession: genome_length})
circos.text("Bacillus paralicheniformis")
circos_gff = Gff(gff_file="candidates.gff")
sector = circos.get_sector(accession)
sector = circos.sectors[0]
cds_track = sector.add_track((80, 100))
cds_track.axis(fc="#EEEEEE", ec="none")
cds_track.genomic_features(circos_gff.extract_features("CDS", target_strand =  1), r_lim=(90, 100),fc="red" )
cds_track.genomic_features(circos_gff.extract_features("CDS", target_strand = -1), r_lim=(80,  90),fc="blue")
pos_list, labels = [], []
cds_track.xticks_by_interval(
    interval=500000,
    label_formatter=lambda label_value: f"{label_value/ 1000000:.1f} Mb",
    label_orientation="vertical")
fig = circos.plotfig().set_figwidth(5)
circos.savefig('Cromosoma-GenesInteres.png')


#Visualización de los datos con seaborn

num_bins = 25
counter_1 = 0
counter_2 = 0
fig, axes = plt.subplots(5,5,figsize=(30,30))
bin_len  = (genome_length - (genome_length % (num_bins - 1))) / (num_bins)
for bin_num in range(num_bins):
    start_pos = bin_num * bin_len
    end_pos   = (bin_num + 1) * bin_len
    mb_df = gff_df.copy()
    mb_df = mb_df[(mb_df["start"]>start_pos) & (mb_df["end"]<=end_pos)]
    sns.swarmplot(ax = axes[counter_1,counter_2],data = mb_df,y="candidate",x="start",hue="strand",dodge=True,palette = 'Set1',order=["exclude","include"],hue_order=["+","-"])
    axes[counter_1,counter_2].set(ylabel=None)
    counter_2 += 1
    if (counter_2%5 == 0):
        counter_2 = 0
        counter_1 += 1
plt.show()
plt.savefig('Swarmplots_operones.png')

#Examinación a detalle del operón seleccionado y sus dominios funcionales

operon_df = gff_df.copy()
operon_df = operon_df[(operon_df["start"]     >= 650000) &
                      (operon_df["end"]       <= 725000) &
                      (operon_df["strand"]    == "-")     &
                      (operon_df["candidate"] == "include")]
operon_df.reset_index(drop=True, inplace=True)
print('Cantidad de genes candidatos en operón seleccionado:'+' '+str(len(operon_df)))

#Archivo fasta con genes candidatos de operón
operon_gene_list = []
for index in operon_df.index.tolist():
    gene_id = operon_df["annotation"][index]["ID"]
    operon_gene_list.append(gene_id)

query_str = ""
for record in SeqIO.parse(aa_file, "fasta"):
    seq_id  = record.id
    if(seq_id in operon_gene_list):
        seq_str = str(record.seq)
        query_str+=">"+seq_id+"\n"+seq_str+"\n"
query_str = query_str.replace("*","")

#búsqueda de dominios se hará a través de la API de InterProScan
#submit_url: Envío de las secuencias
#progress_url: Consulta del status del envío
#results_url: Descarga de resultados

submit_url   = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
progress_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status"
results_url  = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result"

submit_headers   = {"Accept":"text/plain"}
progress_headers = {"Accept":"text/plain"}
results_headers  = {"Accept":"text/tab-separated-values"}

submit_data = {"email":"andrea.pacheco@gmail.com", #Cambiar por correo personal para hacer seguimiento de resultados
               "title":"operon_065_0725",
               "goterms":"false",
               "pathways":"false",
               "stype":"p",
               "sequence":query_str}

submit_request = requests.post(submit_url,data=submit_data,headers=submit_headers)

submit_status_code = submit_request.status_code
submit_job_id      = submit_request.text

#status de envio a InterProScan

progress_request     = requests.get(progress_url+"/"+submit_job_id,headers=progress_headers)
progress_status_code = progress_request.status_code
progress_status      = progress_request.text

#resultado del análisis de dominios funcionales en InterProScan
#Se agregó directamente el job_id generado en la primera corrida y así evitar esperar resultados de InterProScan para continuar con análisis 
results_log_request = requests.get(results_url+"/"+'iprscan5-R20230317-004939-0038-41840473-p1m'+"/log",headers=results_headers)
results_tsv_request = requests.get(results_url+"/"+'iprscan5-R20230317-004939-0038-41840473-p1m'+"/tsv",headers=results_headers)

results_tsv_str = StringIO(results_tsv_request.text)

results_column_names = ["sequence","md5","length","database","accession","description","start","end","evalue","post_processed","date","entry","name"]
results_df = pd.read_csv(results_tsv_str,sep="\t",names=results_column_names)

#Dominios más frecuentes en las columnas 'description' y 'name' del dataframe

results_df_description = results_df[results_df['description'] != '-']
results_df_description.description.value_counts().sort_values().plot(kind='barh', fontsize=12,figsize=(13,13),colormap='Paired')
dominio_description = (results_df_description['description'].mode()).tolist()

results_df_name = results_df[results_df['name'] != '-']
results_df_name.name.value_counts().sort_values().plot(kind='barh',colormap='Paired')
dominio_name = (results_df_name['name'].mode()).tolist()
print('Los dominios funcionales más frecuentes en las secuencias seleccionadas son:')
for x in range (len(dominio_description)):
    print(dominio_description[x])
for j in range (len(dominio_name)):
    print(dominio_name[j])




    







