## Email from Hye Seon

I have 6 samples from endothelial cells (EC) from a mouse aorta.

Samples are 
1. Y_A_P4
2. Y_A_P5

3. O_A_P4
4. O_A_P5

5. O_T_P4
6. O_T_P5

Y, young vs. O, old 
A, EC from aortic arch vs. T, EC from thoracic aorta; 
P4 and P5 means passage 4 & 5 from each cell line.

I thought P4 and P5 are biologic replicates at first, because they are from the same cell line (only different passage). 
But, PCA analysis showed that both are quite different.

So, when I group P4 and P5 together, and compare young vs. old , or Arch vs. Thoracic aorta, there were a small number of genes that showed significant differences, resulting in no significant pathway. 
(for data analysis, I used DEseq2 as Kelsey taught in BI lesson)

Volcano plot and MA plot are as I attached.

I & Dr. Issa thought P4 and P5 should not be biological replicates, so he recommended me to do an analysis just 1 vs. 1 for each comparison (age or region).

For example, 
Y_A_P4 vs. O_A_P4 (for age difference)
O_A_P4 vs. O_T_P4 (for region difference)

No combination of P4 and P5

But, DESeq2 was not working for 1 vs. 1 comparison. 

Could you let me know the method (or code) for 1 vs.1 comparison?

In case you need the R code that I did for 6 samples,
I uploaded Rmd files in the server, my data folder as follow,
/mnt/data/data_hsj/RNAseq_data_process_Hyeseon/2020-06-24_rnaseq_mouse_hyeseon

you also can find count.txt files in the processed data folder.
/mnt/data/data_hsj/RNAseq_data_process_Hyeseon/2020-06-24_rnaseq_mouse_hyeseon/processed_data/03_count