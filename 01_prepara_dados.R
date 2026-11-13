## Script para preparar os dados completos ##
## Carrega pacotes necessários
library(tidyverse)
library(lubridate)

################################################################################
## Leitura dos dados
################################################################################
## Baixa os arquivos do GitHub do monitora.
monitora <- read_csv("data/dados_relatorio_florestal_mastoaves_2026.csv") 
################################################################################
## Preparar um data.frame com os registros de mamiferos apenas
################################################################################
## Retira Cachorro domestico dos dados, corrige grafia de Alouatta sp
## e transfere duas categorias spp para sp.
## e elimina os poucos registros de 2023
## Seleciona bioma Amazônia e mamiferos)
monitora.mam <-
    monitora |>
    filter(taxon_validado!="Canis familiaris",
           classe=="Mammalia", bioma=="Amazônia",
           nivel_taxon == "Espécie",
           distancia < 100) |>
    mutate(cnuc = as.character(cnuc),
           ea = paste(cnuc,ea, sep="_"),
           data = mdy(data)) |>
    select(cnuc,uc, ano, ea, data:distancia)|>
    rename(dist.m = distancia) |>
    na.omit()
################################################################################
## Data frame com esforco em cada  UC, e data
################################################################################
## Data frame com o esforco por cnuc e data
esforco <-
    monitora |>
    filter(taxon_validado!="Canis familiaris",
           classe=="Mammalia", bioma=="Amazônia",
           nivel_taxon == "Espécie",
           distancia < 100) |>
    mutate(cnuc = as.character(cnuc),
           ea = paste(cnuc,ea, sep="_"),
           data = mdy(data)) |>
    ## seleciona variaveis necessarias
    select(cnuc, ea, data, esforco) |>
    ## remove linhas com NAs
    na.omit() |>
    group_by(cnuc, ea, data) |>
    ## Retem o valor maximo de esforço de cada UC, EA e data, convertido em metros
    summarise(tlength = max(esforco)*1000, .groups = "drop")
    

