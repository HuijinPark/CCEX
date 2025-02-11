
## Variable of Filenames in "data" directory

You can find the following naming in filename in data directory

When you calculate the coherence function with CCE method, 
we have to divde the coherence function of super cluster to the coherence functions of its sub-clusters. 

During the process, when you divde the coherence into **zero** due to its modulation effect, then you can encounter the numerical error which is larger than 1.0.

To correct the data, we have to compare the final result to a result that didn't divde the coherence into its subcluster' coherence.

- `wD`: The final data is called **wD**
- `nD`: And not divided data is called **nD**
