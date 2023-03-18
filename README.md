# ProyectoFinal-WBDSCamp
**A. Instalación:**

Las bibliotecas necesarias para que el código funcione adecuadamente son:

- El archivo credentials.py con claves de acceso para la API de NCBI
- pandas 
- pyCirclize 
- pyrodigal 
- requests 
- seaborn 
- subprocess 
- BioPython 
- io 
- Adicionalmente se debe tener instalado el sowftware de NCBI Blast+ en la máquina local

Para instalar las bibliotecas mencionadas, puedes ejecutar los siguientes comandos
pip3 install matplotlib
pip3 install pandas
pip3 install pycirclize
pip3 install pyrodigal
pip3 install requests
pip3 install seaborn
pip3 install biopython

**B. Corrida:**

1. Al correr el código que emplea la biblioteca **pyrodigal** otendremos 2 archivos:

`CP090833.1.faa` que contiene las secuencias aminoacídicas de los genes predichos en el genoma estudiado

`CP090833.1.gff` que contiene las características y coordenadas de los genes predichos en el genoma estudiado

2. Luego de correr el código que extrae información de la API de UniProt se otendrá un archivo de salida:

`uniprot_sequences.fasta`

3. **En los análisis BLAST es necesario cambiar el path los programas makeblastdb y blastp según estén ubicados en la máquina en la que se correrá el codigo**

Luego de hacer los 2 análisis BLAST complementarios se obtendrán 8 archivos:

- Los correspondientes al primer análisis BLAST, es decir en el que se empleó como base de datos a `CP090833.1.faa` y como query a `uniprot_sequences.fasta`:

 `CP090833.1.faa.phr`
 `CP090833.1.faa.pin`
 `CP090833.1.faa.psq`
 
 Siendo el archivo son los resultados finales de este análisis: `CP090833.1.blast.tsv`
 
- Los correspondientes al segundo  análisis BLAST, es decir en el que se empleó como base de datos a    `uniprot_sequences.fasta` y como query a `CP090833.1.faa`:

 `uniprot_sequences.fasta.phr`
 `uniprot_sequences.fasta.pin`
 `uniprot_sequences.fasta.psq`
 
 Siendo el archivo son los resultados finales de este análisis: `uniprot_sequences.blast.tsv`

4. Como resultado de correr el código correspondiente a pycirclize se obtienen los archivos:

`candidates.gff` con las coordenadas de los genes candidatos obtenidos por BLAST

`Cromosoma-GenesInteres.png` con el gráfico del cromosoma bacteriano y los genes de interés

5. Al correr el código correspondiente a la visualización de swarmplots con seaborn se obtiene el archivo:

`Swarmplots_operones.png` con la representación gráfica de la acumulación de genes a lo largo del genoma, para la elección de coordenadas correspondientes a un potencial operón (en este caso se eligieron las coordenadas  0.65 Mbp y 0.725 Mbp en ña variable `operon_df`)

6. Al subir el archivo fasta a InterProScan se requiere cambiar el correo personal correspondiente en la variable `submit_data`

7. En las variables `results_log_request` y `results_tsv_request` se agregó el job_id generado en la primera corrida del código (iprscan5-R20230317-004939-0038-41840473-p1m) para no repetir el análisis cada que se corra el código y así ahorrar tiempo, pero se puede reemplazar por la variable `submit_job_id` para trabajar con los resultados que se acaban de enviar a InterProScan

8. Finalmente a lo largo del código se irá imprimiendo información reelevante respecto a nuestro análisis: 
**El número de genes candidatos resultantes de BLAST, la cantidad de genes candidatos en operón seleccionado y los dominios funcionales más frecuentes en las secuencias seleccionadas**
