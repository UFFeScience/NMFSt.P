from collections import *


def get_frequent_patterns(dataset, min_support):
    # Cria o dicionário para armazenar a contagem de cada item no dataset
    item_counts = defaultdict(int)


    # Conta a frequência de cada item no dataset
    for transaction in dataset:
        for item in transaction:
            item_counts[item] += 1


    # Filtra os itens que têm suporte mínimo
    frequent_items = {item for item, count in item_counts.items() if count >= min_support}


    # Inicializa o conjunto de itemsets frequentes
    frequent_itemsets = set()


    # Gera itemsets frequentes de tamanho 1
    for item in frequent_items:
        frequent_itemsets.add(frozenset([item]))


    # Gera itemsets frequentes maiores que 1
    k = 2
    while True:
        # Gera candidatos para o próximo conjunto de itemsets frequentes
        candidate_itemsets = set()
        for itemset1 in frequent_itemsets:
            for itemset2 in frequent_itemsets:
                # Verifica se os primeiros k-2 itens são iguais
                if len(itemset1.intersection(itemset2)) == k - 2:
                    # Combina os dois itemsets
                    candidate = itemset1.union(itemset2)
                    candidate_itemsets.add(candidate)


        # Conta a frequência dos candidatos no dataset
        itemset_counts = defaultdict(int)
        for transaction in dataset:
            for candidate in candidate_itemsets:
                if candidate.issubset(transaction):
                    itemset_counts[candidate] += 1


        # Filtra os candidatos que têm suporte mínimo
        frequent_itemsets = {itemset for itemset, count in itemset_counts.items() if count >= min_support}


        # Se não há itemsets frequentes, encerra o algoritmo
        if not frequent_itemsets:
            break


        # Adiciona os itemsets frequentes encontrados nessa iteração
        frequent_itemsets.update(frequent_itemsets)


        k += 1


    return frequent_itemsets


# Exemplo de uso do algoritmo FPmax
dataset = [
    [1, 2, 3],
    [1, 2, 4],
    [1, 3, 4],
    [2, 3, 4],
    [2, 3, 5],
    [3, 4, 5]
]
min_support = 2


frequent_itemsets = get_frequent_patterns(dataset, min_support)


# Exibe os itemsets frequentes encontrados
for itemset in frequent_itemsets:
    print(itemset)