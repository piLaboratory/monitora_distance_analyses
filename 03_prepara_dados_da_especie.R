library(unmarked)
## Prepara os dados da especie a ser analisada ##

## Executa o script que baixa e prepara os dados com todos os registros
source("01_prepara_dados.R")

################################################################################
## Dados com avistamentos da especie
################################################################################
## Dataframe com os avistamentos apenas da especie
dados <-
    filter(monitora.mam,
           ## seleciona registros da especie
           taxon_validado==especie,
           ## descarta as 10% distancias mais longas (recomendado para este modelo)
           dist.m<=quantile(dist.m, 0.9)) |>
    na.omit() |>
    select(dist.m, data, ano, cnuc, n_animais)
## Define classes de distancia de avistamento (necessario para o modelo)
## Classes como quintis de distancia, para declarar no modelo
dist_breaks <- quantile(dados$dist.m, probs = seq(0, 1, length.out = 6), na.rm = TRUE)
dist_breaks[5] <- dist_breaks[5] + 0.1
if(any(diff(dist_breaks) <= 0))
    dist_breaks <- c(0,1,2,5,10,max(dados$dist.m) + 0.1)
## Converte em KM
dist_breaks <- unname(dist_breaks)/1000
nbins <- length(dist_breaks) - 1
bin_labels <- 1:nbins  
## Adiciona a classe de distancia ao dataframe de avistamentos
dados$dist_bin <- cut(dados$dist.m/1000, breaks = dist_breaks,
                      include.lowest=TRUE, labels=bin_labels)

## Dados de esforco (em Km de transecto) para cada UC, e data
visits <- 
  mutate(esforco,
         data     = as.Date(data)) |>
    summarise(tlength=sum(tlength)/1000, .by=c(cnuc, data)) |>
    mutate(ano= as.integer(format(data, "%Y")))


## Filtra apenas paras UCS onde a especie pode ocorrer
## Leitura do csv que tem os dados de ocorrencias (conferido a mão, deve estar no mesmo diretorio)
ocorrencias <-
    read.csv("data/conferencia_ocorrencia_especies_nas_ucs.csv", na.string="NA")[, c("cnuc",gsub(" ", ".", especie))]
## Filtra apenas as UCS com registro de ocorrencia
visits <- filter(visits, cnuc %in% ocorrencias[ocorrencias[,2],1])

## Dataframe com n de animais avistados
## Agrega avistamentos (individuos ao invés de eventos)
## Soma de  n_animais por visita e classe de distancia
y_long <-
    dados |>
    mutate(data = as.Date(data)) |>
    group_by(cnuc, data, dist_bin) |>
    summarise(n_ind = sum(n_animais, na.rm = TRUE), .groups = "drop")
## Formata para uma tabela com uma coluna por classe de distancia
y_wide <-
    y_long |>
    mutate(dist_bin = factor(dist_bin, levels = bin_labels)) |>
    pivot_wider(names_from = dist_bin, values_from = n_ind,
                names_prefix = "bin_", values_fill = 0,
                names_sort = TRUE)
## Agrega aos dados de eventos os dados de visitas
visit_y <-
    visits |>
  left_join(y_wide, by = c("cnuc", "data")) |>
  mutate(across(starts_with("bin_"), ~ replace_na(.x, 0)))
## visit_y <-
##     visits |>
##   left_join(y_wide, by = c("cnuc", "ano")) |>
##   mutate(across(starts_with("bin_"), ~ replace_na(.x, 0)))

## Verifica se todas as colunas de distancia existem, mesmo se alguns bins têm zero por todo o dataset
missing_bins <- setdiff(paste0("bin_", bin_labels), names(visit_y))
for (b in missing_bins) visit_y[[b]] <- 0
## Monta a matriz para o modelo: cada linha é um evento de inspeção
## (um transecto em uma data) e as colunas indicam quantos animais
## foram avistados a cada intervalo de distancia
bin_cols <- paste0("bin_", bin_labels)
y <- as.matrix(visit_y[, bin_cols])
storage.mode(y) <- "integer"

################################################################################
## Monta a tabela com as covariaveis de site (UC)
################################################################################
site_covs <-
    visit_y |>
    mutate(
        ## Cada ano em proporcao ao intervalo ano minimo-ano-maximo
        ano_c = (ano - min(monitora.mam$ano))/(max(monitora.mam$ano)-min(monitora.mam$ano)), 
        cnuc_ano = factor(paste(cnuc,ano,sep="_")),
        cnuc  = factor(cnuc),
        ano_f = factor(ano)) |>
    select(cnuc, data, ano, ano_f, ano_c, cnuc_ano, tlength) |>
    na.omit()

## Algumas verificacoes ##
## Numero de linhas das duas tabelas deve ser igual
stopifnot(nrow(site_covs) == nrow(y))
## Sem valores NA para comprimento do transecto na tabela de covariaveis do site
all(!is.na(site_covs$tlength)) ## deve ser TRUE


################################################################################
## Monta o objeto necessário para o ajuste do modelo
## da classe unmarkedFrameDS
umf <- unmarkedFrameDS(
  y           = y,
  siteCovs    = site_covs,
  dist.breaks = dist_breaks,
  survey      = "line",
  tlength     = site_covs$tlength,
  unitsIn     = "km"   ## distancia de avistamento e de comprimento do transecto em metros
)
#verifica 
summary(umf)

