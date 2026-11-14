## Modelo de TAMANHO DE GRUPO, com contagens parciais censuradas a direita.
##   s = n_animais - 1 ~ binomial negativa
##   s | cens(censura) ~ ano_c + (ano_c || cnuc)
## Tamanho medio do grupo num ano: E[n] = 1 + exp(eta)
##
## Uso:
##   especie <- "Tayassu pecari"
##   source("06_ajusta_tamanho.R")
## Salva em ajustes_grupos/<Especie>_tam.rds (e _tam_dist.rds)

library(brms)
if (!exists("dados_tam") || !identical(attr(dados_tam, "especie"), especie)) {
    source("04_prepara_dados_uc_ano.R")
    attr(dados_tam, "especie") <- especie
}

dir.create("ajustes_grupos", showWarnings = FALSE)
arq_tam <- function(m) file.path("ajustes_grupos",
                                 paste0(gsub(" ", "_", especie), "_", m))

## Priors fracamente informativos (escala log do tamanho medio - 1)
priors <- c(prior(normal(0, 2), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(student_t(3, 0, 1), class = "sd"),
            prior(gamma(2, 0.5), class = "shape"))

m_tam <- brm(s_obs | cens(censura) ~ ano_c + (ano_c || cnuc),
             data = dados_tam, family = negbinomial(), prior = priors,
             chains = 6, iter = 2000, warmup = 1000,
             cores = min(6, parallel::detectCores()),
             control = list(adapt_delta = 0.95),
             backend = "rstan", seed = 1,
             file = arq_tam("tam"))   ## carrega o ajuste se o arquivo ja existir

print(summary(m_tam))

################################################################################
## Checagens
################################################################################
## 1. Tendencia do tamanho medio do grupo numa UC tipica (efeitos aleatorios = 0)
R_anos <- ano_max - ano_min
dr <- as_draws_df(m_tam)
n_ini <- 1 + exp(dr$b_Intercept)                  ## tamanho medio em ano_min
n_fim <- 1 + exp(dr$b_Intercept + dr$b_ano_c)     ## tamanho medio em ano_max
pct_tam <- 100 * ((n_fim / n_ini)^(1 / R_anos) - 1)
print(c(tam_medio_ini = mean(n_ini), tam_medio_fim = mean(n_fim),
        pct_ano = mean(pct_tam), quantile(pct_tam, c(0.025, 0.975))))

## 2. Tamanho registrado depende da distancia? (grupos grandes detectados
##    mais longe indicam vies de tamanho na deteccao). Coeficiente de dist_c
##    claramente positivo = sinal de vies.
m_tam_dist <- update(m_tam, formula. = ~ . + dist_c, newdata = dados_tam,
                     file = arq_tam("tam_dist"))
print(fixef(m_tam_dist)["dist_c", ])

## 3. Fracao de contagens parciais ao longo dos anos (se muda, a censura
##    e importante para a tendencia do tamanho)
print(dados_tam |> summarise(n = n(), frac_parcial = mean(parcial), .by = ano) |>
      arrange(ano))
