import logging
import os


logger = logging.getLogger("== analyze_marker_genes.py ==")

def get_file_name_without_extension(file_path):
    file_name_with_extension = os.path.basename(file_path)
    file_name, _ = os.path.splitext(file_name_with_extension)
    return file_name

def generate_gene_faa_file(fasta_file, out_dir, exp_dir):
    gene_file = f"{out_dir}/fasta_reads_genes.faa"
    fragCmd = (f"prodigal -i {fasta_file} -a {gene_file} 1> {exp_dir}/prodigalResult.out 2> {exp_dir}/prodigalResult.err")
    print("Prodigal started generating gene file...")
    logger.debug(f"exec cmd: {fragCmd}")
    os.system(fragCmd)
    return gene_file
    
def generate_hmmout_file(gene_file, out_dir, marker_files, kingdoms, exp_dir):
    
    hmmout_file_list = []
    for marker, kingdom in zip(marker_files, kingdoms):
        hmmout_file = f"{out_dir}/{kingdom}.hmmout"
        hmmCmd = (f"hmmsearch --domtblout {hmmout_file} {marker} {gene_file} 1> {exp_dir}/hmmsearchResult.out 2> {exp_dir}/hmmsearchResult.err")
        print(f"HMMER started hmmsearch...{kingdom}")
        logger.debug(f"exec cmd: {hmmCmd}")
        os.system(hmmCmd)
        hmmout_file_list.append(hmmout_file)
    return hmmout_file_list

# ---------------------------------------------MMSEQS2----------------------------------------------

def generate_phrog_db(out_dir,exp_dir):
    file_path = "phrog_db/phrogdb"
    mkdir_cmd = (f"mkdir phrog_db")
    logger.debug(f"exec cmd: {mkdir_cmd}")
    os.system(mkdir_cmd)
    dw_hmm_cmd = (f"wget https://phrogs.lmge.uca.fr/downloads_from_website/HMM_phrog.tar.gz")
    logger.debug(f"exec cmd: {dw_hmm_cmd}")
    os.system(dw_hmm_cmd)
    open_cmd = (f"tar -xvf HMM_phrog.tar.gz 1> {exp_dir}/unzip.out 2> {exp_dir}/unzip.err")
    logger.debug(f"exec cmd: {open_cmd}")
    os.system(open_cmd)
    rm_unzip_logs = (f"rm {exp_dir}/unzip.out {exp_dir}/unzip.err")
    logger.debug(f"exec cmd: {rm_unzip_logs}")
    os.system(rm_unzip_logs)
    print("Generating PHROG DB...")
    create_db_cmd = (f"ffindex_build -a phrog_hhm_db phrog_hhm_db.index HMM_phrog/ 1> {exp_dir}/ffindexResult.out 2> {exp_dir}/ffindexResult.err")
    logger.debug(f"exec cmd: {create_db_cmd}")
    os.system(create_db_cmd)
    con_db_cmd = (f"mmseqs convertprofiledb phrog_hhm_db phrog_db/phrogdb 1> {exp_dir}/convertprofiledb.out 2> {exp_dir}/convertprofiledb.err")
    logger.debug(f"exec cmd: {con_db_cmd}")
    os.system(con_db_cmd)
    mv_dir = (f"mv phrog_db {out_dir}")
    logger.debug(f"exec cmd: {mv_dir}")
    os.system(mv_dir)
    rm_files_cmd = (f"rm HMM_phrog.tar.gz phrog_hhm_db phrog_hhm_db.index")
    logger.debug(f"exec cmd: {rm_files_cmd}")
    os.system(rm_files_cmd)
    rm_dir_cmd = (f"rm -r HMM_phrog")
    logger.debug(f"exec cmd: {rm_dir_cmd}")
    os.system(rm_dir_cmd)
    
    return file_path

def generate_vog_db(out_dir,exp_dir):
    file_path = f"{out_dir}/vog_db/vogdb"
    mkdir_cmd = (f"mkdir {out_dir}/vog_db")
    logger.debug(f"exec cmd: {mkdir_cmd}")
    os.system(mkdir_cmd)
    
    query_db_cmd = (f"mmseqs databases VOGDB {file_path} tmp 1> {exp_dir}/vogdbsetup.out 2> {exp_dir}/vogdbsetup.err")
    logger.debug(f"exec cmd: {query_db_cmd}")
    print("Generating VOG DB...")
    os.system(query_db_cmd)
    
    rm_dir_cmd = (f"rm -r tmp")
    logger.debug(f"exec cmd: {rm_dir_cmd}")
    os.system(rm_dir_cmd)
    
    return file_path

def generate_fungi_db(out_dir,exp_dir):
    file_path = f"{out_dir}/fungi_db/fungidb"
    
    mkdir_cmd = (f"mkdir {out_dir}/fungi_db")
    logger.debug(f"exec cmd: {mkdir_cmd}")
    os.system(mkdir_cmd)
    
    create_tdb_cmd = (f" mmseqs createdb graphk/fungi_seq_markers.fa {file_path} 1> {exp_dir}/fungidbsetup.out 2> {exp_dir}/fungidbsetup.err")
    logger.debug(f"exec cmd: {create_tdb_cmd}")
    print("Generating Fungi DB...")
    os.system(create_tdb_cmd)
    
    return file_path

def generate_query_db(gene_file, out_dir,exp_dir):
    
    query_db = f"{out_dir}/query_db/querydb"
    mkdir_cmd = (f"mkdir {out_dir}/query_db")
    logger.debug(f"exec cmd: {mkdir_cmd}")
    os.system(mkdir_cmd)
    
    query_db_cmd = (f"mmseqs createdb {gene_file} {query_db} 1> {exp_dir}/querydbsetup.out 2> {exp_dir}/querydbsetup.err")
    print("Generating Query DB...")
    logger.debug(f"exec cmd: {query_db_cmd}")
    os.system(query_db_cmd)
    
    return query_db

def generate_mmseqs_result(query_db, out_dir, profile_databases, mmseqs2_kingdoms,exp_dir):
    
    create_output_dir = (f"mkdir {out_dir}/mmseqs2_files")
    logger.debug(f"exec cmd: {create_output_dir}")
    os.system(create_output_dir)
    
    mmseqs_tab_files = []

    for profile_db, kingdom in zip(profile_databases, mmseqs2_kingdoms):    
        if kingdom == "phrog":
            exhaustive_search = "--exhaustive-search"
        else:
            exhaustive_search = ""
            
        final_file_name = f"{kingdom}"
        
        mmseq_search_cmd = (f"mmseqs search {query_db} {profile_db} {out_dir}/mmseqs2_files/{final_file_name} {out_dir}/tmp_{final_file_name} --threads 16 {exhaustive_search} 1> {exp_dir}/mmseqs2search.out 2> {exp_dir}/mmseqs2search.err")
        print(f"MMseqs2 started searching... {kingdom}")
        logger.debug(f"exec cmd: {mmseq_search_cmd}")
        os.system(mmseq_search_cmd)
        
        final_db_file = (f"ls {out_dir}/mmseqs2_files/{final_file_name}* | grep -v \'{final_file_name}.dbtype \\|{final_file_name}.index\' > {out_dir}/mmseqs2_files/{final_file_name}") #renaming with _
        logger.debug(f"exec cmd: {final_db_file}")
        os.system(final_db_file)    
        
        convertalist = (f"mmseqs convertalis {query_db} {profile_db} {out_dir}/mmseqs2_files/{final_file_name} {out_dir}/mmseqs2_files/{final_file_name}.tab --format-output query,target,evalue,gapopen,pident,fident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen 1> {exp_dir}/converttab.out 2> {exp_dir}/converttab.err")
        logger.debug(f"exec cmd: {convertalist}")
        os.system(convertalist)
        
        mmseqs_tab_files.append(f"{out_dir}/mmseqs2_files/{final_file_name}.tab")
      
    remove_tmp = (f"rm -r {out_dir}/tmp*")
    logger.debug(f"exec cmd: {remove_tmp}")
    os.system(remove_tmp)
    
    return mmseqs_tab_files


def process_hmmout_file(filenames, score_data, threshold):
    try:
        for filename in filenames:
            with open(filename, 'r') as file:
                hmmout_kingdom = get_file_name_without_extension(filename)
                for line in file:
                    if line.startswith('#'):
                        continue
                    data = line.split()
                    read_id = "_".join(data[0].split("_")[:-1])
                    temp_score = (int(data[16]) - int(data[15]) + 1 )/int(data[5])*100
                    temp_marker_gene = data[3]
                    if read_id in score_data:
                        if score_data[read_id]['score'] < temp_score:
                            score_data[read_id]['score'] = temp_score
                            score_data[read_id]['marker_gene'] = temp_marker_gene
                            score_data[read_id]['kingdom'] = hmmout_kingdom
                    else:
                        if temp_score >= threshold:
                            score_data[read_id] = {'marker_gene': temp_marker_gene, 'score': temp_score, 'kingdom': hmmout_kingdom}
                        
    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    return f"Dictionary updated with hmmout file results"

def process_mmseqs_file(filenames, score_data, threshold):
    try:
        for filename in filenames:
            with open(filename, 'r') as f:
                mmseqs_kingdom = get_file_name_without_extension(filename)
                for line in f:
                    row = line.strip().split('\t')
                    
                    read_id = "_".join(row[0].split("_")[:-1])
                    
                    temp_score = 100 * (float(row[11]) - float(row[10])) / float(row[12])
                    
                    temp_marker_gene = row[1]
                    
                    if read_id in score_data:
                        if score_data[read_id]['score'] < temp_score:
                            score_data[read_id]['score'] = temp_score
                            score_data[read_id]['marker_gene'] = temp_marker_gene
                            score_data[read_id]['kingdom'] = mmseqs_kingdom
                            
                    else:
                        if temp_score >= threshold:
                            score_data[read_id] = {'marker_gene': temp_marker_gene, 'score': temp_score, 'kingdom': mmseqs_kingdom}

    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

    print(f"Dictionary updated with mmseqs .tab file results")

def remove_temp_files(out_dir):
    rm_dir_cmd = (f"rm -r {out_dir}/mmseqs2_files")
    logger.debug(f"exec cmd: {rm_dir_cmd}")
    os.system(rm_dir_cmd)

def generate_marker_scores(exp_dir, out_dir):
    
    generate_phrog_db(out_dir,exp_dir)
    generate_vog_db(out_dir,exp_dir)
    generate_fungi_db(out_dir,exp_dir)
    
    profile_databases = [f"{out_dir}/phrog_db/phrogdb", f"{out_dir}/vog_db/vogdb", f"{out_dir}/fungi_db/fungidb"]
    # profile_databases = [f"{out_dir}/vog_db/vogdb", f"{out_dir}/fungi_db/fungidb"]
    mmseqs2_kingdoms = ["phrog", "vog", "fungi"]
    # mmseqs2_kingdoms = ["vog", "fungi"]
    
    marker_files = ["graphk/markers/bacteria_archaea.hmm", "graphk/markers/protist.hmm"]
    marker_kingdoms = ["bacteria", "protist"]
    
    read_file = f"{exp_dir}/reads.fasta"
    
    gene_file = generate_gene_faa_file(read_file, out_dir,exp_dir)


    score_data = {}
    
    query_db = generate_query_db(gene_file, out_dir,exp_dir)
    hmmout_files = generate_hmmout_file( gene_file,out_dir, marker_files, marker_kingdoms,exp_dir)
    # threshold = 50
    process_hmmout_file(hmmout_files, score_data, 50)
    
    resultdb_files = generate_mmseqs_result(query_db,out_dir, profile_databases, mmseqs2_kingdoms,exp_dir)
    process_mmseqs_file(resultdb_files, score_data, 50)
    
    output_file_path = f"{out_dir}/marker_scores.txt"
    with open(output_file_path, 'w') as output_file:
        for key,value in score_data.items():
            output_file.write(f"{key.ljust(20)}\t{value['marker_gene'].ljust(20)}\t{value['score']}\t{value['kingdom']}\n")
    
    remove_temp_files(out_dir)
    
    
def run(exp_dir, out_dir):
    try:
        print("Starting step 2 ...")
        generate_marker_scores(exp_dir, out_dir)
        print("Step 2 completed!")
    except Exception as e:
        print(f"Error: {e}")
