## Data

We have access to two publicly available studies on Metabolights:

- **[MTBLS2542: The trans-omics landscape of COVID-19](https://www.ebi.ac.uk/metabolights/editor/MTBLS2542/descriptors)**
- **[MTBLS1866: Large-Scale Plasma Analysis Revealed New Mechanisms and Molecules Associated with the Host Response to SARS-CoV-2](https://www.ebi.ac.uk/metabolights/editor/MTBLS1866/descriptors)**

### Data Overview

- Under the **`Metabolites`** tab, you'll find the annotated metabolites commonly used in pathway analysis.
- Under the **`Files`** tab, you can access the raw `.raw` files. Untransparently, the preprocessed data is typically not uploaded. Thankfully we have the preprocessed data located in the **`./data`** directory. This dataset is an abundance matrix where each row corresponds to a measured metabolite (referred to as a "feature" due to the ambiguity in identifying it), and each column represents the abundance of that feature across different samples.

While preprocessing is not a focus for this project, if you're interested in understanding how the feature table was generated, this is a good [tutorial](https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.html).

Please note that I need to give the preprocessed data a sanity check and tidy it up with some imputation and removing features with too many missing values first. 

---

## Reading list: 

* Mummichog: [the original paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003123). There have been several improvemtns which you can find on their [website](https://shuzhao-li.github.io/mummichog.org/publications.html) along with notebook tutorials to implement it.
* Pathway analysis with a focus on over representation analysis in metabolomics: [Pathway analysis in metabolomics: Recommendations for the use of over-representation analysis](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009105)
* Pathway analysis using Single Sample Pathway Analysis in metabolomics: [Single sample pathway analysis in metabolomics: performance evaluation and application
](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05005-1)
* A useful [colab tutorial](https://colab.research.google.com/drive/1rUVW7tYKRdVBikpAO2CUk2RkSwrFubLi?usp=sharing) covering pathway analysis in Python. 
* This is more a review for approaches to help in biomarker identification in metabolomics but could shed some insights on the nature of metabolomic datasets: [Statistical analysis in metabolic phenotyping](https://www.nature.com/articles/s41596-021-00579-1#Abs1)


