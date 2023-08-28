# RShiny Application for Some Bioinformatics Analysis

## About
This application is a user-friendly RShiny application designed to simplify bioinformatic analyses of mRNA sequencing data. Tailored for researchers, biologists, and bioinformaticians, this Bioinformatics Analysis webapp aims to streamline the tasks of processing genomic data and extracting meaningful insights.

## App Layout
Effortless Data Upload and Visualisation: Simply upload the files from your computer and visualise the data in a searchable table in the application.
Sample Information Exploration: Load and examine sample information matrix. View distribution of all continous variables in the sample information matrix as density plots.
Counts Matrix Exploration: Input normalized counts matrix by any method as a CSV and be able to choose different gene filtering thresholds and assess their effects using diagnostic plots of the counts matrix. Available analyses include:
Diagnostic scatter plots of median count vs variance and median count vs number of zeros.
Clustered heatmap of counts remaining after filtering.
Scatter plot of principal component analysis projections where user can visualize the top N principal components in a beeswarm plot.
Explore Differential Expression: Simply upload differential expression dataset to visualise it as a volcano plot with customizable X axis and Y axis. View the filtered DEGs in a sortable table.
Visualization of Individual Gene Expression(s): Selects counts from an arbitrary gene and visualize them broken out by a desired sample information variable as a bar plot, boxplot, violin plot, or beeswarm plot.
Future Directions
Develop a Test Suite to test all back-end functions.
Allow user to upload raw data and perform normalization and differential expression analysis in the back-end for a more user friendly experience.
