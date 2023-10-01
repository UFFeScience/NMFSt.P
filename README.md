# NMFSt.P

<a href="https://sol.sbc.org.br/index.php/bresci/article/view/25492">NMFSt.P: um Notebook para Identificação em Paralelo de Subárvores Frequentes em Conjuntos de Árvores Filogenéticas</a>

## Pré-requisitos

- Python 3.10.12
- Clustalw (versão 2.1) 
- Arquivo FASTA com sequências de proteínas



## Instalação de Dependências

Antes de executar o projeto, é necessário instalar as dependências do Python especificadas no arquivo "requirements.txt". Para isso, execute o seguinte comando no terminal:

```
pip install -r requirements.txt
```

Para instalação do Clustalw no linux ( Ubuntu ):

```
sudo apt update
sudo apt-get install clustalw
```
## Organização dos Arquivos
Certifique-se de que todos os arquivos necessários, incluindo as sequências de proteínas no formato FASTA, estejam no diretório especificado em 'input_path'.

## Etapas do Workflow

## 1. Construção de Árvores Filogenéticas

O primeiro passo do workflow é a construção de árvores filogenéticas a partir das sequências de proteínas fornecidas. Para isso, execute o script "Constructor.ipynb" no terminal:

```
python Constructor.ipynb
```
Esse script realiza o alinhamento múltiplo das sequências usando o ClustalW e, em seguida, constrói a árvore filogenética utilizando o método Neighbor-Joining (NJ).

## 2. Construção de Subárvores e Análise de MAF

Após a construção das árvores filogenéticas, o próximo passo é gerar subárvores a partir das árvores principais e realizar a análise de MAF (matriz de frequência de pares de subárvores).

Execute o script "sub_find.ipynb" no terminal:

```
python sub_find.ipynb
```
Esse script irá gerar todas as subárvores a partir das árvores filogenéticas e, em seguida, calcular a matriz de frequência de pares de subárvores (MAF). O resultado será exibido no terminal.

## Saída

As subárvores geradas serão salvas no diretório "out/Subtrees". Além disso, a matriz de frequência de pares de subárvores (MAF) será exibida no terminal durante a execução do script "sub_find.ipynb".

## Limpeza de Arquivos Temporários

O script "Constructor.ipynb" e "sub_find.ipynb" realizarão automaticamente a limpeza dos arquivos temporários gerados durante o processo. Os arquivos temporários serão excluídos do diretório "out/tmp/".

## Considerações Finais

Este guia fornece uma visão geral do workflow para construção de árvores filogenéticas e análise de subárvores. Certifique-se de que os arquivos de entrada estejam corretamente organizados nos diretórios indicados e execute os scripts conforme as etapas descritas.

## Licença

Este projeto está licenciado sob a [MIT License](https://opensource.org/licenses/MIT).
