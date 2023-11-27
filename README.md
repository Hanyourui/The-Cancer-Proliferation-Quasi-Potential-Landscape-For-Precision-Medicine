# The-Cancer-Proliferation-Quasi-Potential-Landscape-For-Precision-Medicine

Cancer is an evolutionary disease, despite breakthroughs in the understanding of cancer genetics, the mechanisms of cancer evolution still remain unknown. 
Cancer is accompanied by a large number of gene mutations and dysfunction in the process of its initiation and progression, which means the pathological state and evolutionary direction of cancer are very complex. 
That makes it difficult to predict a given patient as a cancer carrier stand where in the evolution of cancer, and in what direction? 
In this study, the focus was on depicting the evolution process of cancer, which shows cancer evolution as a ball falling from the peak of high proliferation along the trend of surface into the valley of low proliferation. 
The cancer evolves in its favor in a way that enhances its capacity to proliferate to characterize this capacity of cancer in the form of quasi-potential landscape, the layer-wise narrowing down process is first used to obtain differentially expressed genes strongly related to survival(DESGs). 
Then proliferation entropy is defined to describe the functional relationship between DESGs and proliferation capacity z of sample. After that, the self-organizing constraint map method is designed to characterize the functional relationship between DESGs and geographic coordinate (x,y) of sample. 
These two function relationships are eventually compounded together in the form of quasi-potential landscape to depict the evolution process of cancer.
This study deepens the understanding of cancer evolution. 
The geographic coordinate of landscape represents the pathological state of patient during the whole cancer progression, while the terrain trend represents the evolution direction of the pathological state over time. 
The landscape takes into account the proliferation of cancer, which is main characteristic of cancer. 
It achieves a balance between inter-sample variability and inter-sample heterogeneity by gradually narrowing down the process through layers and utilizing a constraint map to facilitate small-scale clustering.
And it offers potential strategies for diagnosis and treatment, providing valuable insights into the underlying mechanisms.


Description

"TCGA_LUAD_all_deg_sur_expr_sub_smoke.xlsx" is the data of 560 samples' DESGs gene expression data from TCGA-LUAD.
"The Cancer Proliferation Quasi-Potential Landscape For Precision Medicine.py" is the main code to obtain the cancer proliferation quasi-potential landscape of LUAD.
The 3D figure of the landscape is drew by matlab"landscape3D.m".
The Accuracy, Precision,Recall,f1-measure and AUC is calculated by "Subtype_classification.m".
