# Como reproduzir a análise no Windows

Este guia leva você do zero até o relatório compilado (`relatorio_v2.html`), usando os modelos já ajustados. Com os ajustes prontos, nada precisa ser reajustado: os scripts reconhecem os arquivos salvos e os reaproveitam.

Tempo estimado: 1 hora para instalar tudo; de 15 a 30 minutos para gerar as tabelas e o relatório.

## 1. O que instalar

Instale nesta ordem, aceitando as opções padrão dos instaladores, exceto onde indicado.

1. **R**, versão 4.4 ou mais recente: <https://cran.r-project.org/bin/windows/base/>.
   - Versões anteriores à 4.2 não servem: a partir da 4.2 o R no Windows lê textos em UTF-8, e os scripts filtram nomes com acento (por exemplo, "Amazônia").
2. **Rtools**, da **mesma versão** do R (por exemplo, Rtools45 para o R 4.5): <https://cran.r-project.org/bin/windows/Rtools/>.
   - É o compilador usado pelo Stan. Só é indispensável se algum modelo precisar ser reajustado, mas evita erros na instalação dos pacotes.
3. **RStudio Desktop**: <https://posit.co/download/rstudio-desktop/>. Já inclui o pandoc, necessário para compilar o relatório.
4. **Git for Windows**: <https://git-scm.com/download/win>. Opcional: dá para baixar o repositório como `.zip` (passo 3).

## 2. Escolher a pasta do projeto

Use um caminho **curto, sem espaços, sem acentos e fora do OneDrive**, por exemplo:

```
C:\projetos\monitora_distance_analyses
```

O Stan e o Rtools costumam falhar com caminhos que têm espaços ou acentos (por exemplo, `C:\Users\João Silva\OneDrive\Área de Trabalho\`). Pastas sincronizadas pelo OneDrive também podem travar arquivos grandes durante a leitura.

## 3. Baixar o código e os dados

**Com o Git** (no "Git Bash", instalado com o Git for Windows):

```bash
mkdir -p /c/projetos
cd /c/projetos
git clone https://github.com/piLaboratory/monitora_distance_analyses.git
```

**Sem o Git**: no GitHub, botão verde "Code" → "Download ZIP"; descompacte em `C:\projetos\` e renomeie a pasta para `monitora_distance_analyses`, se necessário.

O repositório já contém os dados (`data/`), os scripts e as tabelas do relatório (`relatorio_v2/tabelas/`).

## 4. Baixar os modelos ajustados (Google Drive)

Os modelos ajustados (cerca de 6 GB) e os resultados combinados não estão no GitHub, por serem grandes. Estão no Google Drive do laboratório: <!-- incluir o link da pasta do drive -->

Baixe as duas pastas e coloque-as **dentro da pasta do projeto**, com estes nomes exatos:

```
C:\projetos\monitora_distance_analyses\ajustes_v2\       (arquivos .rds e dois .csv)
C:\projetos\monitora_distance_analyses\resultados_v2\    (arquivos .csv e _concluido.txt)
```

Cuidados:

- O Google Drive compacta pastas grandes em vários `.zip` ao baixar. Descompacte todos e junte o conteúdo numa única pasta `ajustes_v2`, sem subpastas a mais (os `.rds` devem ficar diretamente dentro de `ajustes_v2`).
- Confira se o tamanho total de `ajustes_v2` ficou perto de 6 GB.

Se você só quer **ler e adaptar o relatório** (tabelas e figuras), os passos 4 e 7 podem ser pulados: as tabelas já estão no repositório. Os modelos só são necessários para refazer as tabelas (passo 7) ou para explorar os ajustes.

## 5. Abrir o projeto no RStudio

1. Abra o RStudio.
2. Menu **File → New Project → Existing Directory**, escolha `C:\projetos\monitora_distance_analyses` e clique em "Create Project".
3. Confira no console que a pasta de trabalho é a do projeto:

```r
getwd()
## [1] "C:/projetos/monitora_distance_analyses"
```

Todos os scripts usam caminhos relativos a essa pasta. Nas próximas vezes, abra o projeto pelo arquivo `monitora_distance_analyses.Rproj`.

## 6. Instalar os pacotes do R

No console do RStudio:

```r
install.packages(c("tidyverse", "unmarked", "ubms", "brms", "rstan",
                   "loo", "posterior", "rmarkdown", "knitr", "pkgbuild"))
```

Depois, confira se o compilador foi encontrado (deve responder `TRUE`):

```r
pkgbuild::has_build_tools(debug = TRUE)
```

E se os pacotes principais carregam sem erro:

```r
library(ubms); library(brms); library(tidyverse)
```

Mensagens como "The following objects are masked from ..." são normais. Se a instalação do `rstan` falhar, siga <https://github.com/stan-dev/rstan/wiki/Configuring-C-Toolchain-for-Windows>.

## 7. Gerar as tabelas do relatório

No console do RStudio:

```r
source("09_tabelas_relatorio.R")
```

- Lê os modelos de `ajustes_v2/` e os resultados de `resultados_v2/` e grava as tabelas em `relatorio_v2/tabelas/`.
- Leva de 10 a 30 minutos, quase tudo no teste de qualidade de ajuste (`gof`). O progresso aparece no console ("09: Alouatta macconnelli ...").
- É retomável: se for interrompido, rode de novo; as espécies já calculadas ficam em `relatorio_v2/cache/`.
- Para continuar usando o RStudio enquanto roda: no painel do script aberto, botão **Source → Source as Background Job**.

## 8. Compilar o relatório

Abra `relatorio_v2.Rmd` no RStudio e clique em **Knit** (ou, no console, `rmarkdown::render("relatorio_v2.Rmd")`).

- Resultado: `relatorio_v2.html` na pasta do projeto, e as figuras em `relatorio_v2/figuras/`.
- Para gerar um Word: setinha ao lado de "Knit" → "Knit to Word", ou `rmarkdown::render("relatorio_v2.Rmd", "word_document")`.
- Para mudar uma figura ou tabela, edite o bloco de código correspondente no `.Rmd` e compile de novo. Cada bloco lê um `.csv` de `relatorio_v2/tabelas/`, então não é preciso rodar o passo 7 de novo.

## 9. Se precisar reajustar modelos (opcional, demorado)

Só é necessário se mudar algo nos dados ou nos modelos. Leva horas para todas as espécies.

- **Uma espécie**: no console (roda num processo R separado, como o `08`):

  ```r
  system2(file.path(R.home("bin"), "Rscript"),
          c("roda_especie.R", shQuote("Tayassu pecari")),
          stdout = "logs/Tayassu_pecari.log", stderr = "logs/Tayassu_pecari.log")
  ```

  Crie antes a pasta `logs` (`dir.create("logs")`). O andamento fica no arquivo de log.
- **Todas as espécies**: no console, `source("08_roda_todas.R")`, de preferência como **Background Job**. As espécies já concluídas (com `resultados_v2/<Especie>_concluido.txt`) são puladas; para refazer uma, apague esse arquivo e os `.rds` dela em `ajustes_v2/`.
- O comando `nohup ... &` citado em alguns comentários dos scripts é do Linux; no Windows, use o Background Job do RStudio.
- Configurações (critérios, iterações, pastas): `00_config.R`. Mudar `versao` cria pastas novas (`ajustes_v3/`, `resultados_v3/`) sem apagar as anteriores.

## 10. Problemas comuns

| Sintoma | Causa provável | Solução |
|---|---|---|
| `cannot open file 'data/...': No such file or directory` | pasta de trabalho errada | abrir pelo projeto do RStudio; conferir `getwd()` |
| `cannot open compressed file 'ajustes_v2/...rds'` | modelos não baixados, ou em subpasta a mais | conferir o passo 4 |
| `todas_tend_media.csv` não encontrado | falta `resultados_v2/` | conferir o passo 4 |
| filtro de bioma retorna 0 linhas, ou acentos aparecem como `Ã´` | R anterior à 4.2 | atualizar o R |
| erro ao compilar modelos do Stan / `g++` não encontrado | Rtools ausente ou de outra versão | instalar o Rtools da mesma versão do R; `pkgbuild::has_build_tools(debug = TRUE)` |
| R trava ou "cannot allocate vector" | pouca memória (alguns modelos têm 400 MB) | fechar outros programas; recomenda-se 16 GB de RAM |
| o R encerra com erro **depois** de terminar uma espécie | liberação dos objetos do Stan | inofensivo, se o arquivo `_concluido.txt` foi gravado |
| `Error in pandoc` ao compilar | compilando fora do RStudio sem pandoc | usar o botão Knit do RStudio |

## 11. Onde está cada coisa

| Pasta ou arquivo | Conteúdo |
|---|---|
| `00_config.R` a `10_...R` | scripts, na ordem de execução (ver comentários no início de cada um) |
| `roda_especie.R` | roda 05 → 06 → 07 para uma espécie |
| `data/` | dados do Monitora, tabela de ocorrência por UC, tamanhos máximos de grupo, resultados da versão anterior |
| `ajustes_v2/` | modelos ajustados (`.rds`) e escolha dos modelos (Google Drive) |
| `resultados_v2/` | tendências e densidades combinadas por espécie (Google Drive) |
| `relatorio_v2.Rmd` | relatório de métodos e resultados |
| `relatorio_v2/tabelas/`, `relatorio_v2/figuras/` | tabelas (csv) e figuras (png) do relatório |
| `NOTAS_relatorio.md` | notas das decisões de análise, em ordem cronológica |
