saturation_nmds_plot.R:
Contains R script used to produce saturation plot and NMDS plot. Requires intermediate files votus_cov75thres.txt and metadata.csv. 
This script requires the following packages:
        ggplot2 3.4.2
        reshape2 1.4.4
        vegan 2.5-6
        grid 4.0.2
        gridExtra 2.3

abundance.R:
Contains R script used to produce viral and prokaryotic abundance plots and vOTU count and rank order plots. Requires intermediate files votus_cov75thres.txt, metadata.csv, and qPCR_data.csv. 
This script requires the following packages:
        ggplot2 3.4.2
        reshape2 1.4.4
        readr 2.1.4
        grid 4.0.2
        gridExtra 2.3
        ggpubr 0.6.0

host_plot.R:
Contains R script used to produce host phyla plot and biogeochemical properties. Requires intermediate files votus_cov75thres.txt, metadata.csv, HostPrediction.tsv, and hostproperties.tsv.
This script requires the following packages:
        ggplot2 3.4.2
        readr 2.1.4
        gridExtra 2.3
        ggpubr 0.6.0
        ggvenn 0.1.10
