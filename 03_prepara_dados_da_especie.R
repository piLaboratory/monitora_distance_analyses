library(unmarked)
## Prepara os dados da especie a ser analisada ##
## Requer, antes de executar:
##   especie  : nome valido da especie, ex. "Tayassu pecari"
##   contagem : opcional, "individuos" (padrao) ou "grupos" (n de avistamentos)
if (!exists("contagem")) contagem <- "individuos"
stopifnot(contagem %in% c("individuos", "grupos"))

## Executa o script que prepara os dados com todos os registros
source("01_prepara_dados.R")

################################################################################
## UCs onde a especie pode ocorrer
################################################################################
## Leitura do csv com as ocorrencias (conferido a mao)
ocorrencias <-
    read.csv("data/conferencia_ocorrencia_especies_nas_ucs.csv", na.strings = "NA")[, c("cnuc", gsub(" ", ".", especie))]
ucs_especie <- as.character(ocorrencias$cnuc[ocorrencias[, 2] %in% TRUE])

################################################################################
## Dados com avistamentos da especie
################################################################################
## Obs.: .env$especie e necessario porque monitora.mam pode ter uma coluna
## chamada "especie"; sem isso o filtro compararia com a coluna.
dados_todos <- filter(monitora.mam, taxon_validado == .env$especie)
## Registros fora das UCs de ocorrencia (devem ser zero; se nao, conferir o csv)
fora_ucs <- filter(dados_todos, !cnuc %in% ucs_especie)
if (nrow(fora_ucs) > 0)
    warning(sprintf("%d registros de %s em UCs fora da tabela de ocorrencia; descartados.",
                    nrow(fora_ucs), especie))
## Truncagem: descarta as 10% distancias mais longas (recomendado para este modelo)
dist_trunc <- unname(quantile(dados_todos$dist.m, 0.9))
dados <-
    dados_todos |>
    filter(cnuc %in% ucs_especie, dist.m <= dist_trunc)

## Classes de distancia (necessario para o modelo): quintis das distancias.
## Quando ha empates (ex.: muitas distancias 0), os limites repetidos sao
## fundidos e o numero de classes fica menor que 5.
dist_breaks <- unique(unname(quantile(dados$dist.m, probs = seq(0, 1, length.out = 6))))
if (length(dist_breaks) < 4) {
    ## Menos de 3 classes: usa limites fixos dentro da amplitude observada
    dist_breaks <- unique(c(0, c(1, 2, 5, 10)[c(1, 2, 5, 10) < max(dados$dist.m)],
                            max(dados$dist.m)))
}
## O ultimo limite e o maximo observado; cut(include.lowest = TRUE) inclui 0 e o maximo
## Converte em km
dist_breaks <- dist_breaks / 1000
nbins <- length(dist_breaks) - 1
bin_labels <- 1:nbins
## Adiciona a classe de distancia ao dataframe de avistamentos
dados$dist_bin <- cut(dados$dist.m / 1000, breaks = dist_breaks,
                      include.lowest = TRUE, labels = bin_labels)
stopifnot(!any(is.na(dados$dist_bin)))

################################################################################
## Visitas: UC x data, com o esforco somado de todos os transectos percorridos
################################################################################
## Esforco em km de transecto para cada UC e data (tabela esforco do script 01)
visits <-
    esforco |>
    filter(cnuc %in% ucs_especie) |>
    summarise(tlength = sum(tlength) / 1000, .by = c(cnuc, data, ano))
stopifnot(!anyDuplicated(visits[, c("cnuc", "data")]))

################################################################################
## Contagens por visita e classe de distancia
################################################################################
## n_ind: soma de n_animais (registros sem n_animais contam 0 individuos)
## n_grp: numero de avistamentos (grupos), independente do tamanho
y_long <-
    dados |>
    group_by(cnuc, data, dist_bin) |>
    summarise(n_ind = sum(n_animais, na.rm = TRUE),
              n_grp = n(), .groups = "drop") |>
    mutate(dist_bin = factor(dist_bin, levels = bin_labels))
## Formata para uma tabela com uma coluna por classe de distancia
## (todas as classes presentes, mesmo as sem nenhum registro)
para_matriz <- function(var) {
    visits |>
        left_join(select(y_long, cnuc, data, dist_bin, n = all_of(var)),
                  by = c("cnuc", "data")) |>
        complete(nesting(cnuc, data, ano, tlength), dist_bin = factor(bin_labels, levels = bin_labels),
                 fill = list(n = 0)) |>
        filter(!is.na(dist_bin)) |>
        pivot_wider(names_from = dist_bin, values_from = n,
                    names_prefix = "bin_", names_sort = TRUE) |>
        arrange(cnuc, data)
}
visit_ind <- para_matriz("n_ind")
visit_grp <- para_matriz("n_grp")
stopifnot(identical(visit_ind[, c("cnuc", "data")], visit_grp[, c("cnuc", "data")]))
bin_cols <- paste0("bin_", bin_labels)
y_ind <- as.matrix(visit_ind[, bin_cols]); storage.mode(y_ind) <- "integer"
y_grp <- as.matrix(visit_grp[, bin_cols]); storage.mode(y_grp) <- "integer"
## Matriz usada no modelo: cada linha e uma visita (UC em uma data)
y <- if (contagem == "grupos") y_grp else y_ind

################################################################################
## Monta a tabela com as covariaveis de site (visita)
################################################################################
site_covs <-
    visit_ind |>
    mutate(
        ## Cada ano em proporcao ao intervalo ano minimo-ano maximo do estudo
        ano_c = (ano - ano_min) / (ano_max - ano_min),
        cnuc_ano = factor(paste(cnuc, ano, sep = "_")),
        cnuc  = factor(cnuc),
        ano_f = factor(ano)) |>
    select(cnuc, data, ano, ano_f, ano_c, cnuc_ano, tlength) |>
    as.data.frame()

## Algumas verificacoes ##
stopifnot(nrow(site_covs) == nrow(y),
          !any(is.na(site_covs)),
          ## todos os avistamentos retidos entraram na matriz
          sum(y_grp) == nrow(dados),
          sum(y_ind) == sum(dados$n_animais, na.rm = TRUE))

################################################################################
## Monta o objeto necessário para o ajuste do modelo
## da classe unmarkedFrameDS
umf <- unmarkedFrameDS(
  y           = y,
  siteCovs    = site_covs,
  dist.breaks = dist_breaks,
  survey      = "line",
  tlength     = site_covs$tlength,
  unitsIn     = "km"   ## distancias de avistamento e comprimentos de transecto em km
)
#verifica
summary(umf)
