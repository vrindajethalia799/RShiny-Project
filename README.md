# RShiny Application for Bioinformatics Analysis

## About
This application is a user-friendly RShiny application designed to simplify bioinformatic analyses of mRNA sequencing data. Tailored for researchers, biologists, and bioinformaticians, this Bioinformatics Analysis webapp aims to streamline the tasks of processing genomic data and extracting meaningful insights.

## App Features and Layout
1. **Effortless Data Upload and Visualisation:** Simply upload the files from your computer and visualise the data in a searchable table in the application. <br>
2. **Sample Exploration:** Load and examine the metadata. View distribution of all continous variables in the sample information matrix as density plots. <br>
3. **Counts Exploration:** Input normalized counts matrix by any method as a CSV and be able to choose different gene filtering thresholds and assess their effects using diagnostic plots of the counts matrix. Available analyses include: <br>
* Diagnostic scatter plots of median count vs variance and median count vs number of zeros.
* Clustered heatmap of counts remaining after filtering.
* Scatter plot of principal component analysis projections where user can visualize the top N principal components.
4. **Differential Expression:** Simply upload differential expression dataset to visualise it as a volcano plot with customizable X axis and Y axis. View the filtered DEGs in a sortable table.
5. **Gene Set Enrichment Analysis:** Upload the gene set enrichment analysis result file to view the top 'n' enriched pathways, view the data based on different p-adjusted values and visualize a scatter plot of padj vs NES (normalized enrichment score).

