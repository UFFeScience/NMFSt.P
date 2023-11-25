# %% [markdown]
# - O formato das Sequências de proteínas deve ser especificado em `DATA_FORMATT`

# %%
# Módulos de bibliotecas externas
import os, shutil, sys, time, hashlib, json

# Bibliotecas específicas
import pandas as pd
from Bio import Phylo
from Bio import *
from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpmax, fpgrowth
from dendropy import Tree


INPUT_PATH = '../data/out/Trees'
DATA_OUTPUT_PATH = '../data/out'
DATA_OUTPUT_PATH_SUBTREE = os.path.join(DATA_OUTPUT_PATH,'Subtrees')
OUTPUT_JSON = os.path.join(DATA_OUTPUT_PATH,'output.json')
OUTPUT_LOG = os.path.join(DATA_OUTPUT_PATH,'output_subtrees_log.txt')

COUNT_TREES = 0
COUNT_SUBTREES = 0

DATA_FORMAT = 'nexus' # newick ou nexus
match DATA_FORMAT:
    case 'nexus':
        EXTENTION_FORMAT = 'nexus' # Nexus: 'nexus'
   
    case 'nwk':
        EXTENTION_FORMAT = 'nwk' # Newick: 'nwk'

# Redirecione a saída padrão para o arquivo "output_log.txt"
sys.stdout = open(OUTPUT_LOG, "w")

start = time.time()

# %%
def wellcome():
    print('_'*100)
    print('')
    print('                                               NMFSt.P                                                   ')
    print('                                               -------                                                   ')
    print('')
    print('                           Notebook para Identificação em Paralelo de Subárvores                           ')
    print('                              Frequentes em Conjuntos de Árvores Filogenéticas                              ')
    print('')
    print('                                             - OUTPUT LOG -                                                   ')
    print('_'*100)
    
def resume(start, num_trees, num_subtrees):
    print('')
    print('_'*100)
    print('')
    print('                                               NMFSt.P                                                   ')
    print('                                               -------                                                   ')
    print('')
    print('                           Notebook para Identificação em Paralelo de Subárvores                           ')
    print('                              Frequentes em Conjuntos de Árvores Filogenéticas                              ')
    print('                           ------------------------------------------------------')
    print('')
    print(f'Relatório:\n')
    print(f'Tempo de execução: {"{:.2f}".format(time.time() - start)}s')
    print(f'Árvores analisadas com sucesso: {num_trees}')
    print(f'Número de subárvores produzidas e analisadas: {num_subtrees}\n')
    print('_'*100)
    print('')


wellcome()

# %% [markdown]
# <h3 align="center"> <b>Limpeza</b> </h3>
# 

# %%
def clean_dir(path_clean):
    if(os.path.exists(path_clean)):
        for name in os.listdir(path_clean):
            if(name != "file.gitkeep"):
                file_path = os.path.join(path_clean, name)
                os.remove(file_path)
        shutil.rmtree(path_clean)

def clean_files():
    dir_tmp = os.path.join(DATA_OUTPUT_PATH_SUBTREE) 
    arquivos_tmp = os.listdir(dir_tmp)
    for name_file in arquivos_tmp:
        if name_file != "file.gitkeep":
            os.remove(os.path.join(dir_tmp,name_file))

clean_files()



# %% [markdown]
# <h3 align="center"> <b>Exibe subárvore</b> </h3>
# 

# %%
def print_trees_in_directory(directory):
    for filename in os.listdir(directory):
        if filename.endswith(f'.{EXTENTION_FORMAT}'):
            filepath = os.path.join(directory, filename)
            
            print(filename.upper() + "\n")
            tree = Tree.get_from_path(filepath, DATA_FORMAT)
            tree.print_plot()


# %% [markdown]
# 
# <h4 align="center">Codificadores e Decodificadores </h4>
# 

# %%
def calculate_tree_hash(data):
    newick = data.format("newick")
    hash_object = hashlib.md5(newick.encode())
    return {
        'newick': newick,
        'terrminal_hash': int(hash_object.hexdigest()[:2], 16)
    }

def decode_tree_hash(encoded_data):
    newick = encoded_data['newick']
    original_hash = encoded_data['terrminal_hash']
    hash_object = hashlib.md5(newick.encode())
    if int(hash_object.hexdigest()[:2], 16) == original_hash:
        return newick
    else:
        return None

# %%

def encode_list_to_int(lst):
    lst_str = ",".join(map(str, lst))
    hash_object = hashlib.md5(lst_str.encode())
    return int(hash_object.hexdigest()[:2], 16)

def decode_int_to_list(encoded_int, original_list):
    lst_str = ",".join(map(str, original_list))
    hash_object = hashlib.md5(lst_str.encode())
    if int(hash_object.hexdigest()[:2], 16) == encoded_int:
        return original_list
    else:
        return None

# %%
def tree_to_dict(clade):
    node_dict = {
        "name": clade.name,
        "branch_length": clade.branch_length if clade.branch_length is not None else 0,
        "children": [tree_to_dict(child) for child in clade.clades],
    }
    return node_dict

def dict_to_tree(node_dict):
    clade = Phylo.BaseTree.Clade()
    clade.name = node_dict["name"]
    clade.branch_length = node_dict["branch_length"]

    for child_dict in node_dict["children"]:
        child_clade = dict_to_tree(child_dict)
        clade.clades.append(child_clade)

    return clade

# %% [markdown]
# 
# <h3 align="center"> <b>Construtor de subárvores</b> </h3>
# 

# %%
def sub_tree(path, name_tree):
    tree = Phylo.read(path, DATA_FORMAT)
    name_tree = name_tree.rsplit(".", 1)[0]
    
    print('='*100)
    
    tree_hash = calculate_tree_hash(tree)['terrminal_hash']
    print(f'\nTREE: HASH VALUE\n{name_tree}: {tree_hash}\n')

    dict_tree_terminals_hash = {}
    dict_aux = {}
    result_dict = {}

    for clade in tree.find_clades():
        hash_list = []
        subtree_list_termials = []

        subtree = Phylo.BaseTree.Tree(clade)
        name_subtree = f'{name_tree}_{clade.name}'

        if(clade.is_terminal()):
            clade_hash = calculate_tree_hash(clade.name)['terrminal_hash']
            dict_tree_terminals_hash[clade.name] = clade_hash
            
        if subtree.count_terminals() > 1:
            global COUNT_SUBTREES
            COUNT_SUBTREES += 1 
            print(f'\nSUBTREE:\n{name_subtree}\n')
            Phylo.draw_ascii(subtree)
            print('\n')

            filepath_out = os.path.join(DATA_OUTPUT_PATH_SUBTREE, f'{name_tree}_{clade.name}.{EXTENTION_FORMAT}')
            Phylo.write(subtree, filepath_out, DATA_FORMAT)
                           
            for subtree_clade in subtree.find_clades():
                name_terminal = subtree_clade.name
                if(subtree_clade.is_terminal()):

                    hash_result = calculate_tree_hash(name_terminal)
                    hash_list.append(hash_result)
                
                    subtree_list_termials.append(hash_result['terrminal_hash'])
                    
                    print(f'subtree clade: {subtree_clade}: {hash_result["terrminal_hash"]}')
                    print(f'Decode: {decode_tree_hash(hash_result)}')

            dict_aux[name_subtree] = {'Terminals':subtree_list_termials,
                                      'List_terminals_hash':encode_list_to_int(subtree_list_termials),
                                      'data_terminals':hash_list,
                                      'metadata': tree_to_dict(clade)
                                      }
            print(f'\nTerminais: {subtree_list_termials}:{encode_list_to_int(subtree_list_termials)}\n')

    result_dict[name_tree] = dict_aux
    print(f'CLADES DA ARVORE {name_tree}:\n')       
    for key, value in dict_tree_terminals_hash.items():
        print(f'{key}: {value}') 
    print('-'*100)
 
    return result_dict

# %% [markdown]
# <h3 align="center"> <b>Chamada de funções</b> </h3>
# 

# %% [markdown]
# <p align="center"> <b> ATENÇÃO: SE A ÁRVORE DE ENTRADA NÃO TIVER FORMATO NEXUS, MUDAR O SCHEMA EM Phylo.write( )
# </b> </p>
# 

# %%
# Diretório de entrada de árvores 
arquivos = os.listdir(INPUT_PATH)

# matriz com todas as subárvores
dados_brutos = []
matriz_subtree = []

print(f'\n{"-"*38} SUBÁRVORES PRODUZIDAS {"-"*38}\n')

for name_file in arquivos:
    if(name_file != "file.gitkeep"):
        # global COUNT_TREES
        COUNT_TREES += 1
        dir_path_out_epsecif = os.path.join(DATA_OUTPUT_PATH_SUBTREE,name_file.split(".")[0])
        file_path = os.path.join(INPUT_PATH,name_file)
        dados_brutos.append(sub_tree(file_path, name_file))


# %%
json_result = json.dumps(dados_brutos, indent=2)

with open(OUTPUT_JSON, 'w') as output_file:
    output_file.write(json_result)
    


# %%
def extract_terminais_list(data):
    all_terminais = []
    
    for _, tree_info in data.items():
        for _, subtree_info in tree_info.items():
            all_terminais.append(subtree_info['List_terminals_hash'])
                

    return all_terminais

# %%
for metadata in dados_brutos:
    matriz_subtree.append(extract_terminais_list(metadata))

# %%
print(f'\nMATRIZ DE SUBÁRVORES\n')

for linha in matriz_subtree:
    print(linha)

# %%
max_columns = max(len(row) for row in matriz_subtree)
# máximo de linhas
max_rows = len(matriz_subtree)

def preencher_matriz(matriz, valor_preenchimento):
    # Preencher as células vazias com o valor de preenchimento
    for row in matriz:
        while len(row) < max_columns:
            row.append(valor_preenchimento)

    return matriz

print(f'\nDIM: "{max_rows}x{max_columns}" ','-'*20)

# %% [markdown]
# 
# - [FPmax - python](https://rasbt.github.io/mlxtend/user_guide/frequent_patterns/fpmax/)

# %%
te = TransactionEncoder()
te_ary = te.fit(matriz_subtree).transform(matriz_subtree)
df = pd.DataFrame(te_ary, columns=te.columns_)
df

# %%
result_fpmax = fpmax(df, min_support=0.03, use_colnames=True)
result_fpmax

# %%
fpgrowth(df, min_support=0.2, use_colnames=True)


# %%
def find_exact_subsets(data, result_fpmax, rows):
    print('')
    print('- -'*33)
    print(f'                                       Subárvores frequentes:                                      \n')
    data_dict = result_fpmax.to_dict(orient='records')
    
    print(f'Número de Árvores que contém Subárvores frequêntes no dataset: {len(data_dict)} árvores\n')
    
    for key in data_dict:
        itemset_list = list(key['itemsets'])
        suport = key['support']
        for tree in itemset_list:
            tree_name = extract_name_tree(data, tree)
            print(f'\nÁrvore: { tree_name }\n')
            
            for subtree in itemset_list:
                subtree_name, terminals_list, metadata = extract_subtree_info(data, subtree)

                print(f'     - Subárvore: { subtree_name }')
                print(f'     - Número de aparições: {int(rows * suport)} de {COUNT_TREES} ')
                print(f'     - Suporte: {"{:.3f}".format(suport*100)}%')
                print(f'     - Terminais: { terminals_list }')
                print(f'     - Subárvore ASCII: \n')
                subtree_metadata = dict_to_tree(metadata)
                decode_subtree = Phylo.BaseTree.Tree(subtree_metadata)
                Phylo.draw_ascii(decode_subtree)
            break
    print('- -'*33)

def extract_name_tree(data, hash):

    tree_name = str()
    
    for info in data:
        for tree_name, subtree_info in info.items():
            for _, data in subtree_info.items():
                if "List_terminals_hash" in data and data["List_terminals_hash"] == hash:
                    return tree_name
                
def extract_subtree_info(data, hash):

    subtree_name = str()
    terminals_list = list()
    
    for info in data:
        for _, subtree_info in info.items():
            for subtree_name, data in subtree_info.items():
                if "List_terminals_hash" in data and data["List_terminals_hash"] == hash:
                    terminals_list = [item["newick"] for item in data["data_terminals"]]
                    metadata = data['metadata']
                    return subtree_name, terminals_list, metadata
   

find_exact_subsets(dados_brutos, result_fpmax, max_rows)

# %%
resume(start, COUNT_TREES, COUNT_SUBTREES)

# %% [markdown]
# <h3 align="center"> <b> Comparação das subárvores
# </b> </h3>
# 

# %% [markdown]
# 


