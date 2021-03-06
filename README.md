# Spatial-Transcriptomics
## Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution

Spatial Transcriptomics is an overarching term for a range of methods designed
for assigning cell types(identified by mRNA readouts) to their locations in the
histological sections. This method can also be used to determine subcellular
localisation of mRNA molecules.


The Stahl method implies positioning individual tissue samples on the arrays
of spatially barcoded reverse transcription primers able to capture mRNA with
oligo(dT) tails. Besides oligo(dT) tails and spatial barcode, which indicates
the x and y position on the arrayed slide, the probe contains a cleavage site,
amplification and sequencing handle, and unique molecular identifier.


In the broader meaning of this term, spatial transcriptomics include methods
that can be divided into five principal approaches to resolving spatial distribution
of transcripts. They are microdissection techniques, Fluorescent techniques
in situ hybridisation methods, in situ sequencing, in situ capture protocols and
in silico approaches.

What we do at our labs is in silico analysis. We are provided with the
Single-cell data and the spatial transcriptomic data. Single cell data contain
the information regarding the distribution of mRNA counts in specific cells.
These data are mainly generated by 10X genomics. The spatial transcriptomic
data contain RNA distribution at all the spots in the tissue section.
We run those datasets on the published methods like DestVI, Stereoscope,
Seurat, CARD, DSTG, Autogenes and many more. DestVI and Autogenes, and
almost all other methods, have two model sc â model and st â model. DestVI
posits that for each gene g and each cell n , the number of observed transcripts
follows a negative binomial distribution.The distribution is parametrized
as (rng, pg) with mean (rngpg)
1âpg
and where pg is the gene specific parameter determining
the mean-variance relationship at each spot.Parameter rng = lnÏng
of the negative binomial depends on the type assigned to cn, and its overall
number of detected molecules ln and a low dimensional latent vector Î³n which
captures the variability within its respective cell types. A neural framework
maps Î³n and cn to Ïng. The other st-model also relies on Negative Binomial
distribution. Finally, we deconvolute the sc-model with the st-model getting
the final disttribution at each spot.
