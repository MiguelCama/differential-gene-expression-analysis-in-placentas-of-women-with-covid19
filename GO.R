# ==============================================================================
# R SCRIPT PARA ANÁLISE DE ENRIQUECIMENTO DE GENE ONTOLOGY (GO)
# E GERAÇÃO DE NUVEM DE PALAVRAS (WORD CLOUD)
# Requer: Resultados anotados do DESeq2 (CSV) e um banco de dados de organismo.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP E CARREGAMENTO DE PACOTES
# ------------------------------------------------------------------------------
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(wordcloud)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(ragg) # NOVO PACOTE: Para renderização robusta de gráficos PNG

# !!! ATENÇÃO: SUBSTITUA PELO PACOTE DE ORGANISMO CORRETO !!!
library(org.Hs.eg.db) # Placeholder: Usando Homo sapiens
# ------------------------------------------------------------------------------

# --- File Configurations ---
# Assumimos que os arquivos CSV gerados pelo script DESeq2 estão no mesmo diretório base
BASE_DIR <- "featurecounts/star_results/"
ANNOTATION_FILE <- "kallisto/output/description.tsv" # Annotation file path
OUTPUT_DIR <- file.path(BASE_DIR, "GO_ANALYSIS_OUTPUT")

# Cria o diretório de saída
if (!dir.exists(OUTPUT_DIR)) { dir.create(OUTPUT_DIR) }

# ------------------------------------------------------------------------------
# 2. CONFIGURAÇÕES DE FILTRO
# ------------------------------------------------------------------------------
ALPHA_GO <- 0.05       # Limite de padj para DEGs (padj < 0.05)
LOG2FC_GO <- 1         # Limite de log2FC para DEGs (|log2FC| > 1)
P_ADJ_GO_FILTER <- 0.05 # Limite de padj para termos GO enriquecidos

# ------------------------------------------------------------------------------
# 3. FUNÇÃO PARA RODAR ENRIQUECIMENTO DE GO
# ------------------------------------------------------------------------------

run_go_enrichment <- function(results_file, comparison_name, organism_db, alpha, log2fc_thresh) {
    cat(paste("\n--- Processando GO para:", comparison_name, "---\n"))

    # Carrega os resultados anotados do DESeq2
    res_df <- read.csv(file.path(BASE_DIR, results_file), stringsAsFactors = FALSE)

    # 1. Filtra por DEGs (genes significativamente regulados)
    # A significância estatística é definida por padj < alpha E abs(log2FoldChange) > log2fc_thresh
    DEGs <- res_df %>%
        filter(padj < alpha & abs(log2FoldChange) > log2fc_thresh) %>%
        arrange(padj)

    if(nrow(DEGs) == 0) {
        cat("Nenhum DEG encontrado para esta comparação. Pulando GO.\n")
        return(NULL)
    }
    
    # 2. Mapeamento de IDs: De Gene Symbol para Entrez ID (OBRIGATÓRIO para clusterProfiler)
    
    gene_map <- tryCatch({
        # Usa SYMBOL (o nome do gene, e.g., LACTB) para mapear para Entrez ID
        bitr(DEGs$gene_id, 
             fromType = "SYMBOL", 
             toType = "ENTREZID",
             OrgDb = organism_db)
    }, error = function(e) {
        cat("\nERRO DE MAPEAMENTO: O tipo de ID 'SYMBOL' falhou. Verifique se os IDs no arquivo de DEG são Símbolos de Gene.")
        return(NULL)
    })

    if(is.null(gene_map)) {
        return(NULL)
    }
    
    # 3. Filtra os DEGs que foram mapeados com sucesso
    mapped_genes <- DEGs %>% filter(gene_id %in% gene_map$SYMBOL) 
    entrez_ids <- unique(gene_map$ENTREZID)

    if (length(entrez_ids) < 10) {
        cat("Poucos genes mapeados para Entrez ID. Pulando GO.\n")
        return(NULL)
    }
    cat(paste(length(entrez_ids), "Genes mapeados para Entrez ID para enriquecimento.\n"))


    # 4. Executa o enriquecimento de GO
    ego <- enrichGO(gene          = entrez_ids,
                    OrgDb         = organism_db,
                    keyType       = "ENTREZID",
                    ont           = "BP",       
                    pAdjustMethod = "BH",       
                    pvalueCutoff  = P_ADJ_GO_FILTER,
                    readable      = TRUE)       
    
    # 5. Salva e plota resultados
    
    # CORREÇÃO: as.as.data.frame -> as.data.frame
    go_results_df <- as.data.frame(ego) 
    output_csv <- file.path(OUTPUT_DIR, paste0("GO_results_", comparison_name, "_BP.csv"))
    write.csv(go_results_df, output_csv, row.names = FALSE)
    cat(paste("Resultados GO salvos em:", output_csv, "\n"))
    
    if (nrow(go_results_df) > 0) {
      
        # Plot 1: Dotplot (Melhores 10 termos)
        dot_plot <- dotplot(ego, showCategory = 10, title = paste("Top 10 GO Biological Process:", comparison_name))
        ggsave(file.path(OUTPUT_DIR, paste0("Dotplot_GO_", comparison_name, "_BP.png")), plot = dot_plot, width = 10, height = 6)
      
        # Plot 2: Cnetplot (Rede de Termos e Genes) - Amostra
        if (nrow(go_results_df) >= 3) {
            # É crucial que mapped_genes$log2FoldChange seja um vetor nomeado.
            fc_vector <- mapped_genes$log2FoldChange
            names(fc_vector) <- mapped_genes$gene_id
            
            cnet_plot <- cnetplot(ego, showCategory = 3, categorySize="pvalue", foldChange=fc_vector, 
                                  layout = "circle", title = paste("Cnetplot GO (Top 3 Termos):", comparison_name))
            ggsave(file.path(OUTPUT_DIR, paste0("Cnetplot_GO_", comparison_name, "_BP.png")), plot = cnet_plot, width = 12, height = 10)
        }
        
    } else {
        cat("Nenhum termo GO significativo encontrado para p.ajustado <", P_ADJ_GO_FILTER, "\n")
    }

    # Adiciona Entrez ID aos DEGs antes de retornar
    DEGs_with_Entrez <- DEGs %>%
        left_join(gene_map, by = c("gene_id" = "SYMBOL")) %>%
        # Remove duplicatas que podem ter ocorrido no join (quando um símbolo mapeia para vários Entrez IDs, o que é raro, mas acontece)
        distinct(gene_id, .keep_all = TRUE) 
    
    return(list(ego_object = ego, deg_data = DEGs_with_Entrez))
}

# ------------------------------------------------------------------------------
# 4. FUNÇÃO PARA GERAÇÃO DE NUVEM DE PALAVRAS (WORD CLOUD) E EXPORTAÇÃO
# ------------------------------------------------------------------------------

# 4.1. Nuvem de Palavras de Termos GO Enriquecidos
create_go_word_cloud <- function(ego_object, comparison_name) {
    # CORREÇÃO: as.as.data.frame -> as.data.frame
    if (is.null(ego_object) || nrow(as.data.frame(ego_object)) == 0) {
        cat("Nenhum termo GO para gerar a nuvem de palavras.\n")
        return()
    }
    
    go_df <- as.data.frame(ego_object)
    
    # Filtra os termos GO pela significância estatística (p.adjust < 0.05)
    go_df_top <- go_df %>% 
        filter(p.adjust < P_ADJ_GO_FILTER) %>%
        arrange(p.adjust) %>%
        head(100)
    
    if (nrow(go_df_top) == 0) {
        cat("Nenhum termo significativo após o filtro de P.ajustado para a nuvem de palavras.\n")
        return()
    }

    # O peso das palavras é baseado na significância estatística (-log10 do p.adjust)
    go_df_top$weight <- -log10(go_df_top$p.adjust)
    
    # Limpeza de descrição
    go_df_top$Description <- gsub("process", "", go_df_top$Description, ignore.case = TRUE)
    go_df_top$Description <- gsub("biological", "", go_df_top$Description, ignore.case = TRUE)
    go_df_top$Description <- gsub("regulation", "regul", go_df_top$Description, ignore.case = TRUE)
    
    # Cria o Word Cloud
    png_file <- file.path(OUTPUT_DIR, paste0("WordCloud_GO_Terms_", comparison_name, ".png"))
    
    agg_png(png_file, width = 800, height = 800, res = 150, units = "px")
    
    wordcloud(words = go_df_top$Description, 
              freq = go_df_top$weight, 
              min.freq = quantile(go_df_top$weight, 0.5), 
              max.words = 75,
              random.order = FALSE, 
              rot.per = 0.35, 
              colors = brewer.pal(8, "Dark2"))
    
    # === AJUSTE PARA DESLOCAR O TÍTULO PARA CIMA (line=1) ===
    title(paste("Nuvem de Palavras: Termos GO Enriquecidos (", comparison_name, ")", sep=""), 
          cex.main = 1.2, 
          line = 1) 
    # =======================================================
    dev.off()
    cat(paste("Nuvem de Palavras (GO Terms) salva em:", png_file, "\n"))
}

# 4.2. Nuvem de Palavras de Descrições de Produtos de Genes
create_description_word_cloud <- function(deg_data, comparison_name) {
    
    # deg_data JÁ CONTÉM SOMENTE GENES SIGNIFICATIVOS (DEGs)
    if (nrow(deg_data) == 0) {
        cat("Nenhum gene para gerar a nuvem de palavras de descrição.\n")
        return()
    }
    
    # Carrega o arquivo de anotação - USANDO read.table com sep="\t" para TSV
    # quote="" garante que aspas não causem problemas na leitura de colunas
    annotation_data <- read.table(ANNOTATION_FILE, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    
    # --- BLOCO DE IDENTIFICAÇÃO E CRIAÇÃO DA COLUNA PADRÃO ---
    
    # 1. Limpa nomes de colunas para remover espaços em branco e caracteres invisíveis
    colnames(annotation_data) <- trimws(colnames(annotation_data)) 
    
    cat("\n--- DEBUG: Colunas lidas no ANNOTATION_FILE (limpas) ---\n")
    print(colnames(annotation_data)) # Print 1: Colunas lidas
    
    # 2. Tenta identificar a coluna de descrição de forma flexível (ex: Description, Product, etc.)
    desc_col_index <- which(grepl("description|product", colnames(annotation_data), ignore.case = TRUE))
    
    if (length(desc_col_index) == 0) {
        cat("\n!!! AVISO DE ERRO DE COLUNA !!! As colunas lidas no arquivo de anotação são:\n")
        print(colnames(annotation_data))
        stop("Não foi possível encontrar uma coluna de descrição ('description' ou 'product'). Verifique o nome da coluna no arquivo TSV.")
    }
    
    # Assume a primeira coluna encontrada como a coluna de descrição
    desc_col_name_found <- colnames(annotation_data)[desc_col_index[1]]
    
    # 3. Cria a coluna 'product_description' padronizada usando a sintaxe base do R ($)
    annotation_data$product_description <- gsub("%2C", ",", annotation_data[[desc_col_name_found]])
    
    # 4. Seleciona apenas as colunas essenciais para a junção, garantindo que 'product_description' seja mantida.
    annotation_data <- annotation_data %>%
        dplyr::select(gene_id, product_description)

    cat("\n--- DEBUG: HEAD de ANNOTATION_DATA após criação de 'product_description' ---\n")
    print(head(annotation_data)) # Print 2: Verifica a estrutura antes do join
    
    # --- FIM DO BLOCO DE CORREÇÃO ---
    
    # Mescla as descrições nos DEGs usando 'gene_id'
    deg_descriptions <- deg_data %>%
        left_join(annotation_data, by = "gene_id")

    # >>> CORREÇÃO CRÍTICA PARA O ERRO DE COLUNA NÃO ENCONTRADA (CONFLITO DE NOMES) <<<
    # O join criou product_description.x (do deg_data) e product_description.y (do annotation_data).
    # Precisamos renomear o .y para o nome padrão e remover o .x
    deg_descriptions <- deg_descriptions %>%
        # Garante que a coluna de anotação limpa seja a principal
        rename(product_description = product_description.y) %>%
        # Remove a coluna conflitante que veio do deg_data
        dplyr::select(-any_of("product_description.x"))
    
    cat("\n--- DEBUG: HEAD de DEG_DESCRIPTIONS após o left_join e RENAME (antes do filtro) ---\n")
    # Agora este print deve mostrar a coluna 'product_description' sem sufixos.
    print(head(deg_descriptions)) 
        
    # FIX FINAL: Filtra por descrições não nulas
    deg_descriptions <- deg_descriptions %>%
        filter(!is.na(.data$product_description))
        
    if (nrow(deg_descriptions) == 0) {
        cat("Descrições não disponíveis para os DEGs após a junção. Verifique se os Gene Symbols correspondem.\n")
        return()
    }
    
    # Combina todas as descrições em um único texto
    text_corpus <- paste(deg_descriptions$product_description, collapse = " ")
    
    # Tokenização e contagem de frequência
    words <- unlist(str_split(tolower(text_corpus), "\\W+"))
    words <- words[nchar(words) > 2 & !grepl("^[0-9]+$", words)]
    stop_words <- c("protein", "like", "putative", "family", "member", "containing", "similar", "chain", "domain", "associated", "predicted", "gene", "product", "description", "uncharacterized", "transcript", "variant")
    words <- words[!words %in% stop_words]
    
    word_counts <- table(words)
    word_counts_df <- as.data.frame(word_counts, stringsAsFactors = FALSE)
    colnames(word_counts_df) <- c("word", "freq")
    word_counts_df <- arrange(word_counts_df, desc(freq))
    
    # Cria o Word Cloud
    png_file <- file.path(OUTPUT_DIR, paste0("WordCloud_Gene_Descriptions_", comparison_name, ".png"))
    
    agg_png(png_file, width = 800, height = 800, res = 150, units = "px")
    
    wordcloud(words = word_counts_df$word, 
              freq = word_counts_df$freq, 
              min.freq = quantile(word_counts_df$freq, 0.75), 
              max.words = 100,
              random.order = FALSE, 
              rot.per = 0.35, 
              colors = brewer.pal(8, "Spectral"))
    
    # === AJUSTE PARA DESLOCAR O TÍTULO PARA CIMA (line=1) ===
    title(paste("Nuvem de Palavras: Descrição de Genes (", comparison_name, ")", sep=""), 
          cex.main = 1.2,
          line = 1)
    # =======================================================
    dev.off()
    cat(paste("Nuvem de Palavras (Descriptions) salva em:", png_file, "\n"))
}

# 4.3. FUNÇÃO PARA EXPORTAR DADOS COMPLETOS DE DEGS
# ------------------------------------------------------------------------------
export_deg_data <- function(deg_data, comparison_name) {
    if (nrow(deg_data) == 0) {
        cat("Nenhum DEG para exportar dados.\n")
        return()
    }
    
    # Carrega o arquivo de anotação novamente para obter o GenBank ID
    annotation_data <- read.table(ANNOTATION_FILE, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    colnames(annotation_data) <- trimws(colnames(annotation_data)) 
    
    # 1. Tenta identificar a coluna do GenBank ID de forma flexível (ex: GenBank_ID, Accession, etc.)
    genbank_col_index <- which(grepl("genbank|accession|protein_id", colnames(annotation_data), ignore.case = TRUE))
    
    if (length(genbank_col_index) == 0) {
        # Se não for encontrado, emite um aviso e usa o gene_id
        warning("Não foi possível encontrar a coluna 'GenBank_ID' ou similar no arquivo de anotação. Exportando Gene Symbol (gene_id) em vez disso.")
        genbank_col_name <- "gene_id"
        deg_data_to_export <- deg_data
        output_name <- "Gene Symbol"
        
    } else {
        # Coluna encontrada. Faz a junção.
        genbank_col_name <- colnames(annotation_data)[genbank_col_index[1]]
        output_name <- genbank_col_name
        
        # Seleciona apenas as colunas de anotação relevantes (gene_id e o GenBank ID encontrado)
        annotation_subset <- annotation_data %>%
            dplyr::select(gene_id, !!genbank_col_name) %>%
            distinct(gene_id, .keep_all = TRUE) # Remove duplicatas se houver
        
        # Faz a junção com os dados de DEG
        deg_data_to_export <- deg_data %>%
            left_join(annotation_subset, by = "gene_id")
    }

    # Seleciona as colunas a serem exportadas
    deg_output <- deg_data_to_export %>%
        dplyr::select(!!genbank_col_name, gene_id, log2FoldChange, padj, ENTREZID) %>%
        # Renomeia a coluna GenBank para um nome padrão para a exportação
        rename(Exported_ID = !!genbank_col_name) %>% 
        # Garante que não exporta linhas onde o Entrez ID ou o GenBank ID falharam no mapeamento/junção
        filter(!is.na(Exported_ID) & !is.na(ENTREZID))
        
    
    # Cria o nome do arquivo de saída (agora é CSV)
    output_csv <- file.path(OUTPUT_DIR, paste0("DEG_Data_", comparison_name, ".csv"))
    
    # Escreve os dados no arquivo CSV
    write.csv(deg_output, output_csv, row.names = FALSE)
    
    cat(paste("Dados de", nrow(deg_output), "DEGs (", output_name, ", Gene Symbol, log2FC, padj, Entrez ID) salvos em:", output_csv, "\n"))
    
    # --- NOVO BLOCO: EXPORTAÇÃO SEPARADA APENAS DOS IDs SOLICITADOS (GenBank ID) ---
    output_txt <- file.path(OUTPUT_DIR, paste0("DEG_IDs_", comparison_name, ".txt"))
    writeLines(unique(deg_output$Exported_ID), output_txt)
    cat(paste("Lista de IDs solicitados (", output_name, ") salva em:", output_txt, "\n"))
    # ---------------------------------------------------------------------------------
}


# ------------------------------------------------------------------------------
# 5. EXECUÇÃO PRINCIPAL
# ------------------------------------------------------------------------------

# Lista de comparações a serem processadas
comparisons_to_run <- list(
    list(file = "DESeq2_results_Low_vs_Control_annotated.csv", name = "Low_vs_Control"),
    list(file = "DESeq2_results_High_vs_Control_annotated.csv", name = "High_vs_Control"),
    list(file = "DESeq2_results_High_vs_Low_annotated.csv", name = "High_vs_Low"),
    list(file = "DESeq2_results_High_vs_LowControlCombined_annotated.csv", name = "High_vs_Baseline")
)

all_results <- list()

for (comp in comparisons_to_run) {
    
    # Roda a análise de GO para a comparação (inclui mapeamento de Entrez ID)
    go_result <- run_go_enrichment(comp$file, comp$name, org.Hs.eg.db, ALPHA_GO, LOG2FC_GO)
    
    if (!is.null(go_result)) {
        all_results[[comp$name]] <- go_result
        
        # EXPORTAÇÃO: Dados de DEGs (agora DEG_Data_*.csv e DEG_IDs_*.txt com GenBank ID)
        export_deg_data(go_result$deg_data, comp$name)

        # Geração de Nuvem de Palavras
        # Nuvem de Palavras 1: Termos GO Enriquecidos
        create_go_word_cloud(go_result$ego_object, comp$name)
        
        # Nuvem de Palavras 2: Descrições de Genes (Se houver DEGs)
        create_description_word_cloud(go_result$deg_data, comp$name)
    }
}

cat("\n==============================================================================")
cat("\nProcesso de Análise de GO e Nuvem de Palavras Concluído.")
cat("\nVerifique a pasta 'GO_ANALYSIS_OUTPUT' para resultados CSV, Dotplots, Cnetplots, Word Clouds E os dados completos dos DEGs ('DEG_Data_*.csv').")
cat("\n==============================================================================\n")