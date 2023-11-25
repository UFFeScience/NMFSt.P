# %%
import sqlite3, json, os
import pandas as pd

# PATH_JSON = os.path.join(DATA_OUTPUT_PATH,'outputs','data.json')

# with open(PATH_JSON, 'r') as file:
#         data_json = json.load(file)
        
def init_database(data_json):
    
    DATA_OUTPUT_PATH = '../data/out'
    OUTPUT_DB = os.path.join(DATA_OUTPUT_PATH,'outputs','database.db')
    OUTPUT_CSV = os.path.join(DATA_OUTPUT_PATH,'outputs','data.csv')

    conn = sqlite3.connect(OUTPUT_DB)
    cursor = conn.cursor()

    # Criar tabela
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS tree_data (
            tree_name TEXT,
            subtree_name TEXT,
            terminal_list TEXT,
            list_terminals_hash INTEGER,
            terminal_data TEXT
        )
    ''')

    # Inserir dados do JSON no banco de dados
    for data in data_json:
        for tree_name, inner_data in data.items():
            for inner_name, inner_details in inner_data.items():
                terminal_list = json.dumps(inner_details['Terminals'])
                terminal_data = json.dumps(inner_details['data_terminals'])

                cursor.execute('''
                    INSERT INTO tree_data (tree_name, subtree_name, terminal_list, list_terminals_hash, terminal_data)
                    VALUES (?, ?, ?, ?, ?)
                ''', (tree_name, inner_name, terminal_list, inner_details['List_terminals_hash'], terminal_data))

    conn.commit()

    print("Dados inseridos no banco de dados.")


    # %%

    # Executar uma consulta SQL para selecionar todos os dados da tabela
    consulta_sql = 'SELECT * FROM tree_data'
    cursor.execute(consulta_sql)

    # Obter os resultados da consulta e carregar em um DataFrame do pandas
    dados_df = pd.read_sql_query(consulta_sql, conn)

    # Fechar a conexão com o banco de dados
    conn.close()


    dados_df.to_csv(OUTPUT_CSV)
    dados_df


