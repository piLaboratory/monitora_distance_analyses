## Carrega o arquivo com os dados
source("01_prepara_dados.R")
## Um dataframe que soma o numero de encontros e n de anos em que foi encontrada para cada especie
n.cnuc.ano <-
    monitora.mam |>
    mutate(cnuc = as.character(cnuc)) |>
    summarise(n_encontros=n(), .by=c(taxon_validado, cnuc, ano))
## Data frame com o esforco total por ano
esforco.cnuc.ano <-
    monitora.mam |>
    ## seleciona variaveis necessarias
    select(cnuc, data, esforco, ano) |> 
    mutate(cnuc = as.character(cnuc)) |>
    ## remove linhas com NAs
    na.omit() |>
    ## Retem o valor maximo de esforço de cada UC, EA e data
    summarise(esf.max = max(esforco), .by = c(cnuc, ano, data)) |>
    ## Calcula numero de dias e total de distancia percorrida nos transectos
    summarise (esforco_anual = sum(esf.max), .by=c(cnuc, ano))
## Anexa esforco anual ao dataframe de n de encontros
n.cnuc.ano.esf <-
    left_join(n.cnuc.ano, esforco.cnuc.ano, by=join_by(cnuc, ano))
## Total de registros, de esforco e de anos
tot.encontros <- 
    n.cnuc.ano.esf |>
    group_by(taxon_validado) |>
    summarise(anos=length(unique(ano)),
              sites = length(unique(cnuc)),
              tot_encontros = sum(n_encontros),
              tot_esforco=sum(esforco_anual)) |>
    ungroup() |>
    mutate(taxa_enc_10km=10*tot_encontros/tot_esforco) 

## Nosso criterio : pelo menos 100 registros em pelo menos 14 UCs
selecionadas <-
    filter(tot.encontros, tot_encontros>=150, sites>=14)|>
    arrange(-tot_encontros, -sites, -anos) |>
    as.data.frame() |>
    filter(taxon_validado!="Dasyprocta sp.") ## remove Dasyprocta sem identificacao especifica
## Cria uma tabela com todas as UCs amazonicas e as especies selecionadas
ocorrencia_nas_UCs <-
    monitora.mam |>
    filter(taxon_validado %in% selecionadas$taxon_validado) |>
    mutate(cnuc = as.character(cnuc)) |>
    summarise(presente=n()>0, .by=c(cnuc, uc, taxon_validado)) |>
    pivot_wider(names_from=taxon_validado, values_from=presente) |>
    select(cnuc, uc, sort(selecionadas$taxon_validado)) |>
    arrange(uc)
## exporta a tabela para um csv
## nao executar mais, já foi feito uma vez
## write.csv(ocorrencia_nas_UCs, file="data/conferencia_ocorrencia_especies_nas_ucs_novo_para_conferencia.csv", row.names=FALSE)

## Criterio usado no artigo do Monitora, mas com outro modelo e premissas
## Especies de mamiferos com pelo menos 0,5 encontros/10km de esforço
## filter(tot.encontros, taxa_enc_10km>0.5)|>
##     arrange(-anos, -taxa_enc_10km) |>
##     as.data.frame()
