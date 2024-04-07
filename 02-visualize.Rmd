# (PART) VISUALIZATION TERMS {-}


# Visualization Terminologies
This compilation provides concise definitions of essential bioinformatics concepts, complemented by simple R code illustrations demonstrating various visualization techniques. These illustrations cover a wide range of visualization types, offering a comprehensive overview of bioinformatics visualization practices.


## ...IN PROGRESS...

<!-- ## Term: Phylogenetic Tree -->
<!-- ### Definition -->
<!-- A phylogenetic tree is a branching diagram that represents the evolutionary relationships among a group of organisms based on their genetic similarities and differences. -->

<!-- ### R Code Illustration -->

<!-- ```{r phytree} -->
<!-- # Load required packages -->
<!-- library(ape) -->

<!-- # Create a simple phylogenetic tree -->
<!-- tree <- read.tree(text = "(((Human, Chimp), Gorilla), Orangutan);") -->

<!-- # Plot the phylogenetic tree -->
<!-- plot(tree) -->


<!-- ``` -->

<!-- ### Interpretation: -->
<!-- - Phylogenetic trees visually represent evolutionary relationships among organisms. -->
<!-- - Constructed from genetic data, such as DNA or protein sequences. -->
<!-- - Branching patterns illustrate common ancestry and divergence over time. -->
<!-- - Each node represents a common ancestor, with branch lengths indicating genetic divergence. -->
<!-- - Used in evolutionary biology, biodiversity studies, and comparative genomics. -->
<!-- - Enables inference of evolutionary relationships, tracing species origins, and understanding diversification patterns. -->


<!-- ## Term Nucleotide sequence logo -->

<!-- ### Definition -->
<!-- A nucleotide sequence logo is a graphical representation of the consensus sequence derived from multiple sequence alignment, displaying the relative frequency of each nucleotide (A, T, G, C) at each position. It provides a visual summary of sequence conservation, highlighting important regions and motifs within DNA or RNA sequences. -->

<!-- ### R Code Illustration -->

<!-- ```{r nuclogo} -->
<!-- # Load required package -->
<!-- library(ggseqlogo) -->

<!-- # Example aligned sequences -->
<!-- aligned_sequences <- c("ATGC", "ATGC", "ATGC", "ATGC") -->

<!-- # Create sequence logo -->
<!-- logo <- ggseqlogo(aligned_sequences) -->

<!-- print(logo) -->


<!-- ``` -->

<!-- ### Interpretation: -->
<!-- - Nucleotide sequence logos depict the frequency of each nucleotide at each position in a sequence alignment. -->
<!-- - Conserved regions are represented by tall letters, indicating a high frequency of a particular nucleotide at that position. -->
<!-- - Variable positions have shorter letters or gaps, indicating lower sequence conservation. -->
<!-- - The height of the letters at each position reflects the relative frequency of the corresponding nucleotide. -->
<!-- - Sequence logos help identify conserved motifs, regulatory elements, and functional sites within DNA or RNA sequences. -->
<!-- - They are widely used in bioinformatics to visualize sequence conservation and infer functional significance from sequence data. -->


<!-- ## Amino Acid Sequence Logo -->

<!-- ### Definition -->
<!-- An amino acid sequence logo is a graphical representation of the consensus sequence derived from multiple sequence alignment of protein sequences, displaying the relative frequency of each amino acid at each position. It provides a visual summary of sequence conservation, highlighting important regions and motifs within protein sequences. -->

<!-- ### R Code Illustration -->

<!-- ```{r aalogo} -->
<!-- # Load required library -->
<!-- library(ggseqlogo) -->

<!-- # Example protein sequences with more similarities -->
<!-- protein_sequences <- c( -->
<!--   "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKK", -->
<!--   "MTDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRV", -->
<!--   "MVKVIVTDCTHTKHPGASRSGMVRHGNSGQQRQISLLKDLPEVKNRDIMVVFPNEITFEYGI", -->
<!--   "MVKVIVTDCTHTKHPGASRSGMVRHGNSGQQRQISLLKDLPEVKNRDIMVVFPNEITFEYGI", -->
<!--   "MVKVIVTDCTHTKHPGASRSGMVRHGNSGQQRQISLLKDLPEVKNRDIMVVFPNEITFEYGI" -->
<!-- ) -->

<!-- # Create protein sequence logo -->
<!-- ggseqlogo(protein_sequences) -->


<!-- ``` -->

<!-- ### Interpretation -->
<!-- - Amino acid sequence logos depict the frequency of each amino acid at each position in a sequence alignment. -->
<!-- - Conserved regions are represented by tall letters, indicating a high frequency of a particular amino acid at that position. -->
<!-- - Variable positions have shorter letters or gaps, indicating lower sequence conservation. -->
<!-- - The height of the letters at each position reflects the relative frequency of the corresponding amino acid. -->
<!-- - Sequence logos help identify conserved motifs, functional domains, and structurally important regions within protein sequences. -->
<!-- - They are widely used in bioinformatics to visualize sequence conservation and infer functional significance from protein sequence data. -->


<!-- ## Term: Genome Browser -->

<!-- ### Definition -->
<!-- A genome browser is a bioinformatics tool used to visualize and explore genomic data, allowing researchers to navigate, analyze, and interpret various genomic features within a genome assembly. Genome browsers provide an interactive interface for displaying genomic sequences, gene annotations, regulatory elements, genetic variation, and other relevant genomic information. -->

<!-- ### R Code Illustration -->

<!-- ```R -->
<!-- # Open a genome browser in a web browser -->
<!-- browseURL("https://genome.ucsc.edu/cgi-bin/hgGateway") -->

<!-- ``` -->

<!-- ```R -->
<!-- # Example code for using the Gviz package to visualize genomic data -->
<!-- library(Gviz) -->
<!-- library(BSgenome.Hsapiens.UCSC.hg19) -->
<!-- library(TxDb.Hsapiens.UCSC.hg19.knownGene) -->
<!-- txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene -->
<!-- gene <- "ENSG00000157764" -->
<!-- chr <- "chr6" -->
<!-- start <- 26025651 -->
<!-- end <- 26078124 -->
<!-- plotTracks( -->
<!--   GenomeAxisTrack(), -->
<!--   AnnotationTrack( -->
<!--     txdb, -->
<!--     name = "Genes", -->
<!--     chromosome = chr, -->
<!--     start = start, -->
<!--     end = end, -->
<!--     isLabelInside = TRUE, -->
<!--     fontsize = 8 -->
<!--   ) -->
<!-- ) -->
<!-- ``` -->

<!-- ## Term: Heatmap -->

<!-- ### Definition -->
<!-- A heatmap is a graphical representation of data where values in a matrix are represented as colors. It is commonly used in bioinformatics to visualize high-dimensional data, such as gene expression levels, protein abundance, or genomic features. In a heatmap, rows and columns of the matrix represent individual samples or features, and the color intensity indicates the magnitude of the values. Heatmaps are useful for identifying patterns, clusters, or correlations within large datasets. -->


<!-- ```{r heatmap, message=FALSE, warning=FALSE} -->
<!-- # Load required library -->
<!-- library(tidyverse) -->

<!-- # Example data -->
<!-- data_matrix <- matrix(rnorm(100, mean = 0, sd = 1), nrow = 10, ncol = 10) -->

<!-- # Convert matrix to dataframe -->
<!-- data_df <- as.data.frame(data_matrix) -->

<!-- # Reshape data for ggplot -->
<!-- data_df <- data_df %>% -->
<!--   mutate(row = row_number()) %>% -->
<!--   pivot_longer(-row, names_to = "column", values_to = "value") -->

<!-- # Create heatmap -->
<!-- ggplot(data = data_df, aes(x = column, y = row, fill = value)) + -->
<!--   geom_tile() + -->
<!--   scale_fill_gradient(low = "green", high = "red") + -->
<!--   theme_minimal() -->


<!-- ``` -->

<!-- ### Interpretation -->
<!-- - Heatmaps provide a visual representation of complex datasets, allowing researchers to quickly identify patterns and trends. -->
<!-- - Each cell in the heatmap corresponds to a value in the data matrix, and the color intensity represents the magnitude of the value. -->
<!-- - Heatmaps are commonly used in gene expression analysis to visualize gene expression profiles across samples or conditions. -->
<!-- - Clustering algorithms can be applied to heatmap data to group samples or features with similar expression patterns. -->
<!-- - Heatmaps are versatile and can be customized with various color scales, annotations, and clustering methods to enhance data interpretation. -->


<!-- ## Term: Volcano Plots -->

<!-- ### Definition -->
<!-- A volcano plot is a type of scatter plot commonly used in bioinformatics to visualize the results of statistical tests, such as differential gene expression analysis. The plot displays the relationship between statistical significance (usually represented as the negative logarithm of the p-value) on the y-axis and the magnitude of change (e.g., fold change) on the x-axis. Each point on the plot represents a gene or feature, with significantly differentially expressed genes typically highlighted based on predefined thresholds for statistical significance and fold change. -->

<!-- ### R Code Illustration -->
<!-- ```{r volcanoplot} -->
<!-- # Load required library -->
<!-- library(ggplot2) -->

<!-- # Generate example data1 -->
<!-- log_fc <- rnorm(100, mean = 0, sd = 1) -->
<!-- p_values <- runif(100, min = 0, max = 0.05) -->

<!-- # Create volcano plot -->
<!-- df <- data.frame(log_fc = log_fc, p_values = p_values) -->
<!-- ggplot(df, aes(x = log_fc, y = -log10(p_values))) + -->
<!--   geom_point() + -->
<!--   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + -->
<!--   geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") + -->
<!--   labs(x = "Log2 Fold Change", y = "-log10(P-value)") -->


<!-- # Generate example data2 -->
<!-- set.seed(42)  # for reproducibility -->
<!-- n <- 100 -->
<!-- log_fc <- c(runif(n/3, min = 1, max = 2), runif(n/3, min = -2, max = -1), runif(n/3, min = -0.5, max = 0.5))  # UP, DOWN, and non-regulated genes -->
<!-- p_values <- c(runif(n/3, min = 0.0001, max = 0.01), runif(n/3, min = 0.0001, max = 0.01), runif(n/3, min = 0.5, max = 1))  # Significant and non-significant genes -->

<!-- # Create dataframe -->
<!-- df <- data.frame(log_fc = log_fc, p_values = p_values) -->

<!-- # Define significance level -->
<!-- significance_level <- 0.05 -->

<!-- # Create a new column for significance -->
<!-- df$significance <- ifelse(df$p_values <= significance_level, "Significant", "Not Significant") -->
<!-- df$direction <- ifelse(df$log_fc > 0.5, "Up", ifelse(df$log_fc < -0.5, "Down", "Non-significant")) -->

<!-- # Create volcano plot with color differentiation -->
<!-- ggplot(df, aes(x = log_fc, y = -log10(p_values), color = significance, shape = direction)) + -->
<!--   geom_point(size = 2) + -->
<!--   geom_hline(yintercept = -log10(significance_level), linetype = "dashed", color = "red") + -->
<!--   labs(x = "Log2 Fold Change", y = "-log10(P-value)", color = "Significance", shape = "Direction") + -->
<!--   scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) + -->
<!--   scale_shape_manual(values = c("Up" = 17, "Down" = 15, "Non-significant" = 16)) + -->
<!--   theme_minimal() -->


<!-- ``` -->

<!-- ### Interpretation -->
<!-- - Volcano plots provide a visual representation of differential expression analysis results, highlighting genes or features that exhibit significant changes in expression levels. -->
<!-- - The x-axis represents the magnitude of change (e.g., fold change), while the y-axis represents statistical significance (usually represented as the negative logarithm of the p-value). -->
<!-- - Genes with high statistical significance and large fold changes are often considered as potential candidates for further investigation. -->
<!-- - Thresholds for statistical significance and fold change can be adjusted based on study-specific criteria or biological context. -->


<!-- ## Term: Principal Component Analysis (PCA) Plot -->

<!-- ### Definition -->
<!-- A Principal Component Analysis (PCA) plot is a visualization technique used to reduce the dimensionality of multivariate data while preserving its variance. It is commonly employed in bioinformatics to explore patterns and relationships within high-dimensional datasets, such as gene expression or metabolomics data. PCA identifies orthogonal axes, called principal components, that capture the maximum variance in the data. By projecting data points onto these components, PCA plots provide insights into the underlying structure of the dataset and facilitate the identification of clusters or trends. -->

<!-- ### R Code Illustration -->

<!-- ```{r pca} -->
<!-- # Load required libraries -->
<!-- library(ggplot2) -->

<!-- # Set seed for reproducibility -->
<!-- set.seed(123) -->

<!-- # Define parameters -->
<!-- num_samples <- 150    # Number of samples -->
<!-- num_genes <- 500      # Number of genes -->
<!-- num_clusters <- 3     # Number of clusters -->

<!-- # Generate synthetic gene expression data with distinct clusters -->
<!-- gene_expr <- matrix(0, nrow = num_samples, ncol = num_genes) -->

<!-- # Assign clusters to samples -->
<!-- cluster_assignments <- sample(1:num_clusters, num_samples, replace = TRUE) -->

<!-- # Generate gene expression data for each cluster -->
<!-- for (i in 1:num_clusters) { -->
<!--   cluster_size <- sum(cluster_assignments == i) -->
<!--   cluster_mean <- rnorm(num_genes, mean = runif(1, min = -1, max = 1), sd = 0.5) -->
<!--   gene_expr[cluster_assignments == i, ] <- matrix(rnorm(cluster_size * num_genes, mean = cluster_mean, sd = 0.5), nrow = cluster_size, ncol = num_genes) -->
<!-- } -->

<!-- # Perform PCA -->
<!-- pca_result <- prcomp(gene_expr, scale. = TRUE)  # Scale the data for PCA -->

<!-- # Create PCA plot -->
<!-- pca_data <- as.data.frame(pca_result$x) -->

<!-- # Add cluster assignments to the PCA data -->
<!-- pca_data$Cluster <- factor(cluster_assignments) -->

<!-- # Create PCA plot with color grouping by cluster -->
<!-- ggplot(data = pca_data, aes(x = PC1, y = PC2, color = Cluster)) + -->
<!--   geom_point() + -->
<!--   labs(x = "Principal Component 1", y = "Principal Component 2", color = "Cluster") + -->
<!--   theme_bw() -->


<!-- ``` -->

<!-- ### Interpretation -->
<!-- - PCA plots visualize high-dimensional data in a lower-dimensional space, capturing the most significant variations. -->
<!-- - Each data point represents a sample or observation in the dataset, projected onto the principal components. -->
<!-- - Clustering or patterns in the PCA plot indicate similarities or differences among samples, revealing underlying structure or trends. -->
<!-- - The first few principal components often explain the majority of the variance in the data, making them useful for data exploration and visualization. -->


<!-- ## Term: Network Graphs -->

<!-- ### Definition: -->
<!-- Network graphs, also known as graphs or networks, are mathematical representations of relationships between entities, often visualized as nodes (vertices) connected by edges (links or relationships). In bioinformatics, network graphs are widely used to model complex biological systems, such as protein-protein interactions, gene regulatory networks, metabolic pathways, and phylogenetic relationships. They provide a powerful framework for analyzing and visualizing biological data in the context of interconnected networks, enabling researchers to uncover functional modules, identify key nodes or hubs, and infer biological processes. -->

<!-- ### R Code Illustration: -->
<!-- ```{r netgraph} -->
<!-- # Load required library -->
<!-- library(igraph) -->

<!-- # Example data1 -->
<!-- nodes <- letters[1:20] -->
<!-- edges <- data.frame(from = sample(nodes, size = 20, replace = TRUE), -->
<!--                     to = sample(nodes, size = 20, replace = TRUE)) -->

<!-- # Create network graph -->
<!-- graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes) -->

<!-- plot(graph) -->

<!-- # Example data2 -->
<!-- # Load required libraries -->
<!-- library(igraph) -->
<!-- library(ggraph) -->
<!-- library(ggrepel) -->

<!-- # Set seed for reproducibility -->
<!-- set.seed(123) -->

<!-- # Define parameters -->
<!-- num_genes <- 20  # Number of genes -->
<!-- num_edges <- 30  # Number of edges -->

<!-- # Generate random gene names -->
<!-- genes <- paste0("Gene", 1:num_genes) -->

<!-- # Generate random edges between genes -->
<!-- edges <- data.frame( -->
<!--   from = sample(genes, size = num_edges, replace = TRUE), -->
<!--   to = sample(genes, size = num_edges, replace = TRUE) -->
<!-- ) -->

<!-- # Create igraph graph object -->
<!-- graph <- graph_from_data_frame(edges, directed = FALSE, vertices = genes) -->

<!-- # Plot network graph using ggraph -->
<!-- g <- ggraph(graph, layout = "fr") + -->
<!--   geom_edge_link() + -->
<!--   geom_node_point(color = "steelblue", size = 5) + -->
<!--   geom_node_text(aes(label = name), repel = TRUE) +  # Add node labels with repulsion -->
<!--   theme_graph() -->
<!-- g -->

<!-- ``` -->

<!-- ### Interpretation -->
<!-- - Network graphs represent relationships between entities as interconnected nodes and edges. -->
<!-- - Nodes typically represent biological entities such as genes, proteins, metabolites, or individuals, while edges denote interactions or relationships between them. -->
<!-- - The topology of a network graph, including node degree, centrality, and clustering coefficient, provides insights into the structure and organization of biological systems. -->
<!-- -nAnalysis of network graphs can reveal important biological properties, such as modular organization, network motifs, and key regulatory elements. -->


<!-- ## Term: Manhattan Plots -->

<!-- ### Definition -->
<!-- Manhattan plots are graphical representations used in genome-wide association studies to visualize the p-values resulting from statistical tests for association between genetic variants and a trait of interest. In these plots, each data point represents a single nucleotide polymorphism (SNP) or genetic variant, plotted according to its genomic position on the x-axis and its corresponding -log10(p-value) on the y-axis. -->

<!-- ### R Code Illustration -->
<!-- ```{r manhattanplot, message=FALSE, warning=FALSE} -->

<!-- # Generate example data for a Manhattan plot -->
<!-- chromosomes <- rep(1:5, each = 10) -->
<!-- positions <- rep(1:10, 5) -->
<!-- p_values <- c(rnorm(20, mean = 0, sd = 1), rnorm(30, mean = 0, sd = 2), rnorm(20, mean = 0, sd = 3), rnorm(20, mean = 0, sd = 4), rnorm(10, mean = 0, sd = 5)) -->
<!-- manhattan_data <- data.frame(Chromosome = chromosomes, Position = positions, P_Value = p_values) -->

<!-- # Plot Manhattan plot -->
<!-- library(ggplot2) -->
<!-- ggplot(manhattan_data, aes(x = Position, y = -log10(P_Value), color = factor(Chromosome))) + -->
<!--   geom_point(size = 2) + -->
<!--   scale_color_manual(values = c("blue", "red", "green", "orange", "purple")) + -->
<!--   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") + -->
<!--   labs(x = "Genomic Position", y = "-log10(P-value)", title = "Manhattan Plot", color="Chromosome") -->


<!-- # Load required library -->
<!-- library(ggplot2) -->

<!-- # Set seed for reproducibility -->
<!-- set.seed(123) -->

<!-- # Number of markers -->
<!-- num_markers <- 200 -->

<!-- # Generate simulated genomic data -->
<!-- genomic_position <- 1:num_markers  # Assuming each marker is located at consecutive genomic positions -->
<!-- p_values <- c(runif(num_markers - 10, min = 0, max = 0.1), runif(10, min = 0, max = 1e-6))  # Simulate p-values -->

<!-- # Create a data frame -->
<!-- manhattan_data <- data.frame( -->
<!--   Genomic_Position = genomic_position, -->
<!--   P_Value = p_values -->
<!-- ) -->

<!-- # Plot Manhattan plot with skyscrapers -->
<!-- ggplot(manhattan_data, aes(x = Genomic_Position, y = -log10(P_Value))) + -->
<!--   geom_point(size = 2, color="black") +  # Points representing markers -->
<!--   geom_segment(aes(xend = Genomic_Position, yend = 0), color = "steelblue") +  # Skyscrapers -->
<!--   geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red", size=1.1) +  # Significance threshold line -->
<!--   labs(x = "Genomic Position", y = "-log10(P-value)", title = "Manhattan Plot") + -->
<!--   theme_minimal() -->


<!-- ``` -->

<!-- ### Interpretation -->
<!-- - Manhattan plots are commonly used in genome-wide association studies (GWAS) to identify genetic variants associated with complex traits or diseases. -->
<!-- - SNPs or genetic variants that surpass a significance threshold (-log10(p-value)) are considered potentially associated with the trait under investigation. -->
<!-- - The distinctive peaks in a Manhattan plot indicate regions of the genome where significant genetic associations may be located. -->
<!-- - Researchers use Manhattan plots to prioritize regions of interest for further investigation, such as candidate gene identification or functional annotation. -->

