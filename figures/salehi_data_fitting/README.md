Fitting ALFA-K model to passaging data
================

Experimental data is available here
<https://osf.io/tvzw9/?view_only=0c216deb154b4625a3bbb7db280882ca>

Code chunks are not shown in this version of the README. See README.Rmd
for more detail.

Code in this folder covers

1)  Processing experimental data
2)  Training ALFA-K on experimental data
3)  Predicting karyotype evolution with experimental data
4)  Generating the following figures:

![Figure 3. Fits to cell passaging data A) Maps showing the relationship
cell passages that were sequenced and used as inputs to ALFA-K. Each
filled circle represents a sequenced sample. Lines connect parental
passages to children, with line type indicating the presence or absence
of cisplatin treatment in the intervening period. Circles are colored
according to the
![R^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R%5E2 "R^2")
metric value for fitted landscapes including that sample and all its
ancestors. B)
![R^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R%5E2 "R^2")
metric scores for all sublineages in the dataset, grouped according to
the number of samples used in the fitting. C) Angle metric scores for
all lineages with at least one descendent, grouped according to the
number of samples used in the fitting and the
![R^2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R%5E2 "R^2")
metric score. Red horizontal line indicates the expected value of the
angle under the null hypothesis (i.e. that the fitted and true landscape
are different).](figures/figure.png)
