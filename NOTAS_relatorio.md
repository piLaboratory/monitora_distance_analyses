# Notas para o relatório

Registro das decisões de análise e dos achados que precisam entrar no texto.
Atualizar ao longo das análises.

## Preparação dos dados

- **Esforço**: calculado a partir da tabela completa (todas as classes), uma
  linha por transecto e data. Na versão anterior o esforço era tirado só de
  linhas de mamíferos identificados até espécie, o que perdia ~40% do esforço
  e as visitas sem detecção de mamíferos.
- **Filtro taxonômico**: pelo nível validado (`class_taxon_val`), não pelo de campo.
- **Correções de ano/data**: 1 data com ano digitado errado (2014 em vez de
  2024, Parna da Amazônia); 2 casos de mesma UC e data com `ano` diferente
  (Ouro Preto, Tapajós-Arapiuns), corrigidos para o ano mais frequente no dia.
  `ano` é o ano da campanha (pode diferir do ano da data em campanhas que
  atravessam a virada do ano, ex. Viruá).
- **Período**: 2014–2025 (R = 11 anos).
- **Tamanho de grupo**: valores acima do máximo plausível por espécie
  (`data/tamanho_max_grupo.csv`, com referências) marcados como anômalos e
  excluídos da estimativa de tamanho, mas mantidos como 1 avistamento na
  contagem de grupos. Anômalos concentrados em poucas UCs (ex. Hadrosciurus em
  Cazumbá-Iracema; Mazama americana em Tapajós-Arapiuns) → possível diferença de
  registro entre equipes; vale consultar a equipe do Monitora.

## Modelo

- **Unidade**: grupo (avistamento), não indivíduo. Com indivíduos havia forte
  sobredispersão (φ até ~32 em Tayassu), explicada pelo agrupamento. Com grupos,
  φ ajustado ≈ 1–2.
- **Agregação**: contagens e esforço por UC-ano.
- **Abundância de grupos** (`ubms::stan_distsamp`, half-normal):
  `~ ano_c ~ ano_c + (ano_c || cnuc)` (g_sem) e o mesmo com `+ (1 | cnuc_ano)`
  (g_com). Principal: g_sem se φ ajustado ≤ 1,5, senão g_com; o outro é
  sensibilidade.
- **Detecção**: `~ ano_c`, sem efeito de UC (efeito de UC na detecção não é
  identificável com poucas detecções por UC: divergências concentradas em
  σ_det alto).
- **Inclinação por UC**: mantida se o g_sem com inclinação convergir
  (Rhat ≤ 1,01 e BFMI ≥ 0,2); senão, modelos só com tendência média.
- **Faixas de distância**: 0–2 m, depois tercis acima de 2 m. ~20% das
  distâncias são 0 (arredondamento); com faixas por quantis a 1ª faixa ficava
  estreita (0–0,4 a 0–1 m) e o modelo half-normal não ajustava.
- **Tamanho de grupo** (`brms`, binomial negativa em n − 1): com distância
  (tamanho avaliado à distância zero, corrige o viés de detecção de grupos
  grandes); contagens parciais censuradas à direita; sensibilidade sem censura
  (só contagens totais).
- **UC classificável**: ≥ 10 grupos em ≥ 3 anos com detecção.

## Achado importante: contagens parciais e tendência de tamanho

- **A tendência de GRUPOS é o resultado principal.** Não depende de
  `n_animais` nem do campo de contagem parcial.
- A fração de contagens marcadas como **parciais aumentou em todas as 15
  espécies** entre 2014–17 e 2022–25, inclusive nas solitárias (Tapirus
  0 → 13%, Eira 5 → 17%, Mazama americana 7 → 18%), e em 15 de 17 combinações
  espécie-UC com dados suficientes nos primatas.
- Nas espécies gregárias (Ateles, Cebus, Sapajus apella, S. macrocephalus,
  Dicotyles, Nasua, Tayassu) a fração subiu 16–32 pontos, enquanto a média das
  contagens totais **caiu** e o mínimo contado nas parciais **não subiu** (caiu
  na maioria).
- Interpretação: mudança na forma de registrar (protocolo/treinamento), não
  aumento real dos grupos. Com censura e distribuição de tamanho de cauda
  pesada, cada parcial "puxa" o tamanho estimado para cima; mais parciais ao
  longo do tempo → tendência de tamanho positiva espúria (v1: Ateles +6 %/ano,
  Cebus +6, S. apella +8; sem censura, Cebus −3).
- Alternativa não descartada: grupos maiores serem mais difíceis de contar
  (aí o modelo sem censura subestimaria). O padrão do mínimo contado pesa
  contra, mas não elimina.
- **No relatório**: indivíduos nas espécies gregárias reportados com as duas
  estimativas (com e sem censura) como faixa de sensibilidade; nas solitárias
  as duas devem coincidir.
- **A fazer**: perguntar à equipe do Monitora se houve mudança de orientação
  sobre marcar contagens como parciais (por volta de 2017–2018).

## Outros pontos para a discussão

- Tendências por UC: com poucas detecções a estimativa é encolhida para a
  média geral; só UCs classificáveis devem ser interpretadas.
- φ do `ubms` (`posterior_predict`, `gof`) é condicional aos dados (abundância
  latente sorteada dado y) → φ e `gof` pouco sensíveis; φ calculado com o
  valor esperado não condicional.
- LOO com efeito por observação (cnuc_ano) tem muitos k de Pareto altos; não
  usado para escolher modelos.
