## Script para preparar os dados completos ##
## Carrega pacotes necessários
library(tidyverse)
library(lubridate)
source(Sys.getenv("CONFIG_ANALISE", "00_config.R"))  ## configuracoes (00_config.R por padrao)

################################################################################
## Leitura dos dados
################################################################################
monitora <- read_csv("data/dados_relatorio_florestal_mastoaves_2026.csv",
                     show_col_types = FALSE)

## Tabela de tamanho máximo plausível de grupo por espécie (conferir à mão).
## Registros com n_animais acima deste valor são marcados como anômalos.
tam_max <- read_csv("data/tamanho_max_grupo.csv", show_col_types = FALSE)

## Todos os registros do bioma Amazônia, de qualquer classe (aves, mamíferos
## e linhas de "visita sem registro"), com identificadores padronizados.
monitora.amz <-
    monitora |>
    filter(bioma == "Amazônia") |>
    mutate(cnuc = as.character(cnuc),
           ea = paste(cnuc, ea, sep = "_"),
           data = mdy(data),
           cont_tot_par = str_trim(cont_tot_par))

## Coerência entre 'ano' (ano da campanha) e a data.
## O ano da campanha pode diferir em 1 do ano da data (campanhas que atravessam
## a virada do ano, ex. Parna do Viruá). Duas correções:
## (a) diferença > 1 ano: considera erro de digitação no ano da data e o
##     substitui pelo ano da campanha;
## (b) mesma UC e data com mais de um 'ano': usa o ano mais frequente nas
##     linhas daquela UC e data (empate: o ano da data).
correcoes_ano <- list(
    a = monitora.amz |>
        filter(abs(year(data) - ano) > 1) |>
        distinct(cnuc, uc, ea, data, ano))
monitora.amz <-
    monitora.amz |>
    mutate(data = if_else(abs(year(data) - ano) > 1,
                          make_date(ano, month(data), day(data)), data))
ano_por_data <-
    monitora.amz |>
    count(cnuc, data, ano) |>
    arrange(cnuc, data, desc(n), desc(ano == year(data))) |>
    group_by(cnuc, data) |>
    mutate(n_anos = n()) |>
    slice(1) |>
    ungroup() |>
    select(cnuc, data, ano_corrigido = ano, n_anos)
correcoes_ano$b <-
    monitora.amz |>
    left_join(ano_por_data, by = c("cnuc", "data")) |>
    filter(ano != ano_corrigido) |>
    distinct(cnuc, uc, ea, data, ano, ano_corrigido)
monitora.amz <-
    monitora.amz |>
    left_join(select(ano_por_data, -n_anos), by = c("cnuc", "data")) |>
    mutate(ano = ano_corrigido) |>
    select(-ano_corrigido)

## Período do estudo (usado para padronizar o ano em 0-1)
ano_min <- min(monitora.amz$ano)
ano_max <- max(monitora.amz$ano)

################################################################################
## Esforço por transecto (ea) e data
################################################################################
## O esforço (km) é anotado em uma única linha por transecto e data, que pode
## ser de qualquer registro (ave, mamífero ou linha sem registro). Por isso o
## esforço é calculado com a tabela completa, ANTES de qualquer filtro
## taxonômico. Cada linha desta tabela é uma visita a um transecto.
esforco <-
    monitora.amz |>
    group_by(cnuc, ea, data, ano) |>
    summarise(n_valores = sum(!is.na(esforco)),
              tlength = if (all(is.na(esforco))) NA_real_
                        else max(esforco, na.rm = TRUE) * 1000, ## em metros
              .groups = "drop")
## Visitas sem esforço anotado não podem ser usadas (nem seus registros)
visitas_sem_esforco <- filter(esforco, is.na(tlength))
esforco <- filter(esforco, !is.na(tlength)) |> select(-n_valores)

################################################################################
## Registros de mamíferos
################################################################################
## - Exclui cachorro doméstico
## - Usa o nível taxonômico VALIDADO (class_taxon_val), não o de campo
## - Exige distância anotada e < 100 m
## - Mantém registros com n_animais ausente: cada linha é um avistamento
##   (grupo), válido para a contagem de grupos mesmo sem o tamanho
## - Descarta registros de visitas sem esforço anotado
monitora.mam <-
    monitora.amz |>
    filter(classe == "Mammalia",
           class_taxon_val == "Espécie",
           taxon_validado != "Canis familiaris",
           !is.na(distancia), distancia < 100) |>
    semi_join(esforco, by = c("cnuc", "ea", "data")) |>
    select(cnuc, uc, ano, ea, data, estacao, taxon_validado,
           n_animais, cont_tot_par, distancia, plaqueta) |>
    rename(dist.m = distancia) |>
    left_join(tam_max, by = "taxon_validado") |>
    ## Marcação dos registros que podem entrar na estimativa de tamanho de grupo
    mutate(
        anomalo = !is.na(n_animais) & !is.na(tam_max_grupo) &
                  n_animais > tam_max_grupo,
        tam_valido = !is.na(n_animais) & n_animais > 0 &
                     cont_tot_par %in% "total" & !anomalo,
        motivo_tam_invalido = case_when(
            is.na(n_animais) | n_animais == 0 ~ "n_animais ausente ou zero",
            anomalo                           ~ "acima do máximo plausível",
            !(cont_tot_par %in% "total")      ~ "contagem parcial ou sem informação",
            TRUE                              ~ NA_character_)
    )

## Registros idênticos em todos os campos (possíveis duplicatas de digitação).
## Apenas relatados; não são removidos.
duplicatas <-
    monitora.mam |>
    filter(duplicated(pick(everything())) | duplicated(pick(everything()), fromLast = TRUE))

message(sprintf(
    "01_prepara_dados: correções de ano/data: %d visitas (data) e %d visitas (ano); %d visitas transecto-data com esforço (%.0f km); %d sem esforço descartadas; %d registros de mamíferos; %d com tamanho válido; %d linhas em possíveis duplicatas.",
    nrow(correcoes_ano$a), nrow(correcoes_ano$b),
    nrow(esforco), sum(esforco$tlength) / 1000, nrow(visitas_sem_esforco),
    nrow(monitora.mam), sum(monitora.mam$tam_valido), nrow(duplicatas)))
