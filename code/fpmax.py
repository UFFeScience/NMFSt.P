from collections import defaultdict


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
    [198, 85, 226, 33, 104, 99, 165, 225, 198, 85, 226, 33, 104, 85, 226, 33, 104, 85, 226, 33, 104, 99, 165, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [123, 117, 44, 175, 103, 16, 162, 10, 171, 123, 117, 44, 175, 103, 16, 117, 44, 175, 103, 16, 117, 44, 175, 44, 175, 103, 16, 162, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [59, 35, 191, 47, 63, 59, 35, 47, 63, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [43, 73, 186, 3, 218, 59, 131, 21, 181, 108, 43, 73, 186, 3, 218, 59, 131, 21, 73, 186, 3, 218, 59, 131, 21, 186, 3, 218, 59, 131, 21, 186, 3, 218, 59, 3, 218, 59, 218, 59, 131, 21],
    [181, 164, 17, 176, 120, 122, 33, 27, 42, 181, 164, 17, 176, 120, 122, 181, 164, 17, 176, 164, 17, 176, 17, 176, 120, 122, 33, 27, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [65, 6, 97, 29, 253, 142, 52, 64, 66, 65, 6, 97, 29, 253, 65, 6, 97, 29, 253, 97, 29, 142, 52, 64, 52, 64, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [3, 181, 8, 58, 165, 79, 21, 254, 212, 3, 181, 8, 58, 165, 79, 21, 254, 212, 165, 79, 21, 254, 79, 21, 254, 21, 254, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    [43, 73, 186, 141, 218, 59, 178, 21, 181, 108, 43, 73, 186, 141, 218, 59, 178, 21, 73, 186, 141, 218, 59, 178, 21, 186, 141, 218, 59, 178, 21, 186, 141, 218, 59, 141, 218, 59, 218, 59, 178, 21]

]
min_support = 2


frequent_itemsets = get_frequent_patterns(dataset, min_support)


# Exibe os itemsets frequentes encontrados
for itemset in frequent_itemsets:
    print(itemset)