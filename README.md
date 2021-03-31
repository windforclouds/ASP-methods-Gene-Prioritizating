## Gene Prioritizating by PS-V2N methods(*by PS-V2N.R*)

### This is a shortest path-based algorithm named PS-V2N (Proximity Score of Vertex to Network) which was proposed for the target identification. Here is our work pipeline:

![](https://github.com/windforclouds/PS-V2N/blob/master/pictures/pipeline.png)

### **1**. work before 
Before calculation，ready for one txt/csv file(such as **MSPPIN_vertex.csv**), first column is gene name, and the second colname is the type identification code 0, 1, and 2 (0 means groupB, 1 maens groupA, 2 means share genes); get your edges-network file ready (**such as MSPPIN_edge.csv**)

### **2**. Our method could be summarized into four steps：

①  Divided genes in network into subsets **A** and **B**;

②  Choose gene *vi* from **A** and all genes from **B**;

③  Keep the edges among the selected genes in step 2) to construct the sub-network ***Gsub***;

④  Calculate the proximity score for gene *vi* in ***Gsub***.

![](https://github.com/windforclouds/PS-V2N/blob/master/pictures/PS-V2N.jpg)

### **3**. Run PS-V2N.R and then you'll get a table list ordered by PS-V2N value,results are as follows:

![](https://github.com/windforclouds/PS-V2N/blob/master/pictures/results.png)

